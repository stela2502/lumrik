//! REF/ALT single-cell SNP evidence and optional H5AD annotation.
use anyhow::{bail, Context, Result};
use clap::Args;
use flate2::read::MultiGzDecoder;
use scdata::{FeatureIndex, MexFeatureIndex, load_mtx_feature_matrix, read_mtx_barcodes};
use std::{collections::{BTreeMap, HashMap, HashSet}, fs::{self, File}, io::{BufRead, BufReader, BufWriter, Write}, path::{Path, PathBuf}};

#[derive(Debug, Clone, Args)]
pub struct SnpCli {
    /// Matrix Market folder containing reference allele counts.
    #[arg(long)] pub r#ref: PathBuf,
    /// Matrix Market folder containing alternate allele counts.
    #[arg(long)] pub alt: PathBuf,
    /// Optional GEX H5AD to annotate; never overwritten.
    #[arg(long)] pub data: Option<PathBuf>,
    /// Optional text/TSV list of variant IDs to analyze (first column).
    #[arg(long)] pub analyzed_snps: Option<PathBuf>,
    /// Minimum REF counts to call an individual SNP WT (when ALT=0).
    #[arg(long, default_value_t=4)] pub min_ref_for_wt: u64,
    /// Minimum WT SNP positions needed to call a cell WT overall.
    #[arg(long, default_value_t=2)] pub min_wt_snps: usize,
    #[arg(long, default_value="snp_analysis")] pub out: PathBuf,
}

#[derive(Default, Clone, Copy)]
struct Counts { r: u64, a: u64 }

fn normalize(raw: &str) -> String {
    // Match the notebook: the first DNA-only barcode before a numeric 10x suffix;
    // otherwise use the first 16 characters for padded ONT identifiers.
    let b=raw.as_bytes();
    for start in 0..b.len() {
        if !matches!(b[start], b'A'|b'C'|b'G'|b'T') { continue; }
        let mut end=start;
        while end<b.len() && matches!(b[end],b'A'|b'C'|b'G'|b'T') {end+=1;}
        if end+1<b.len() && b[end]==b'-' && b[end+1].is_ascii_digit() {
            return raw[start..end].to_owned();
        }
    }
    raw.chars().take(16).collect()
}
fn text_lines(path: &Path) -> Result<Box<dyn BufRead>> {
    let f=File::open(path).with_context(|| format!("opening {}",path.display()))?;
    if path.extension().is_some_and(|e| e=="gz") { Ok(Box::new(BufReader::new(MultiGzDecoder::new(f)))) }
    else { Ok(Box::new(BufReader::new(f))) }
}
fn feature_type(dir: &Path) -> Result<String> {
    let path=["features.tsv.gz","features.tsv","genes.tsv.gz","genes.tsv"].iter().map(|p|dir.join(p)).find(|p|p.is_file())
        .with_context(|| format!("no features file in {}",dir.display()))?;
    for line in text_lines(&path)?.lines() {
        let line=line?; if line.trim().is_empty(){continue;}
        return Ok(line.split('\t').nth(2).unwrap_or("Gene Expression").to_owned());
    }
    bail!("empty features file: {}",path.display())
}
fn read_counts(dir: &Path) -> Result<(Vec<String>, HashMap<(String,String),u64>,HashSet<String>)> {
    let kind=feature_type(dir)?;
    let idx=MexFeatureIndex::from_dir(dir,&kind)?;
    let features=idx.ordered_feature_ids().iter().map(|&id|idx.feature_name(id).to_owned()).collect::<Vec<_>>();
    let barcodes=read_mtx_barcodes(dir)?;
    let (data,_,_)=load_mtx_feature_matrix(dir,&kind,1)?;
    let mat=data.as_sprs().map_err(anyhow::Error::msg)?;
    if mat.rows()!=features.len() || mat.cols()!=barcodes.len() { bail!("unexpected matrix dimensions in {}",dir.display()); }
    let mut counts=HashMap::new(); let mut cells=HashSet::new();
    for (i,(barcode,_)) in barcodes.iter().enumerate() {
        let norm=normalize(barcode); cells.insert(norm.clone());
        for (j,v) in mat.outer_iterator().enumerate().filter_map(|(j,row)|row.get(i).map(|v|(j,*v))) {
            if !v.is_finite() || v<0.0 || v.fract()!=0.0 {bail!("invalid SNP count {v}");}
            *counts.entry((norm.clone(),features[j].clone())).or_insert(0)+=v as u64;
        }
    }
    Ok((features,counts,cells))
}
fn requested(path: Option<&Path>) -> Result<Option<HashSet<String>>> {
    let Some(path)=path else {return Ok(None)};
    let mut out=HashSet::new();
    for line in text_lines(path)?.lines() {
        let line=line?; let s=line.split('\t').next().unwrap_or("").trim();
        if s.is_empty() || s.starts_with('#') || matches!(s,"variant"|"variant_id"|"snp") {continue;}
        out.insert(s.to_owned());
    }
    if out.is_empty(){bail!("--analyzed-snps contains no variants")}
    Ok(Some(out))
}
fn call(c: Counts, min_ref: u64) -> &'static str {
    if c.a>0 {"MUT"} else if c.r>=min_ref {"WT"} else {"NO_INFO"}
}
fn safe_name(s: &str) -> String {
    s.chars().map(|c|if c.is_ascii_alphanumeric() || c=='_' {c}else{'_'}).collect()
}

pub fn run(cli: SnpCli) -> Result<()> {
    if cli.min_ref_for_wt==0 || cli.min_wt_snps==0 {bail!("WT thresholds must be positive")}
    fs::create_dir_all(&cli.out)?;
    let (rf,rc,mut cells)=read_counts(&cli.r#ref)?;
    let (af,ac,alt_cells)=read_counts(&cli.alt)?;
    cells.extend(alt_cells);
    let rset:HashSet<_>=rf.iter().cloned().collect(); let aset:HashSet<_>=af.iter().cloned().collect();
    if rset!=aset {bail!("REF and ALT feature sets differ; refusing positional matching")}
    if rset.len()!=rf.len() || aset.len()!=af.len(){bail!("duplicate SNP feature names")}
    let filter=requested(cli.analyzed_snps.as_deref())?;
    let mut snps=rf.into_iter().filter(|s|filter.as_ref().is_none_or(|f|f.contains(s))).collect::<Vec<_>>();
    if let Some(f)=&filter { for s in f {if !rset.contains(s){bail!("requested SNP {s:?} not found")}} }
    if snps.is_empty(){bail!("no SNPs selected")}
    snps.sort();
    let mut cell_ids=cells.into_iter().collect::<Vec<_>>(); cell_ids.sort();
    let mut per=BufWriter::new(File::create(cli.out.join("snp_variant_evidence.tsv"))?);
    writeln!(per,"barcode\tvariant\tref_count\talt_count\talt_fraction\tgenotype")?;
    let mut summary=BufWriter::new(File::create(cli.out.join("snp_cell_genotypes.tsv"))?);
    writeln!(summary,"barcode\tgenotype\tref_total\talt_total\tmut_snps\twt_snps\tno_info_snps")?;
    let mut observed=BTreeMap::new();
    for cell in &cell_ids {
        let mut total=Counts::default(); let mut mut_n=0; let mut wt_n=0;
        let mut evidence=Vec::new();
        for snp in &snps {
            let c=Counts {r:*rc.get(&(cell.clone(),snp.clone())).unwrap_or(&0),a:*ac.get(&(cell.clone(),snp.clone())).unwrap_or(&0)};
            total.r+=c.r;total.a+=c.a;
            let state=call(c,cli.min_ref_for_wt);
            if state=="MUT" {mut_n+=1} else if state=="WT" {wt_n+=1}
            let fraction=if c.r+c.a>0 {c.a as f64/(c.r+c.a) as f64} else {0.0};
            writeln!(per,"{cell}\t{snp}\t{}\t{}\t{fraction}\t{state}",c.r,c.a)?;
            evidence.push(c);
        }
        let genotype=if mut_n>0 {"MUT"} else if wt_n>=cli.min_wt_snps {"WT"} else {"NO_INFO"};
        writeln!(summary,"{cell}\t{genotype}\t{}\t{}\t{mut_n}\t{wt_n}\t{}",total.r,total.a,snps.len()-mut_n-wt_n)?;
        observed.insert(cell.clone(),(genotype.to_owned(),total,mut_n,wt_n,evidence));
    }
    per.flush()?;summary.flush()?;
    if let Some(input)=&cli.data {
        enrich_h5ad(input,&cli.out.join("snp_enriched.h5ad"),&snps,&observed)?;
    }
    eprintln!("[sc-analyse snp] {} SNPs, {} cells; output {}",snps.len(),cell_ids.len(),cli.out.display());
    Ok(())
}

fn enrich_h5ad(input:&Path, output:&Path, snps:&[String], observed:&BTreeMap<String,(String,Counts,usize,usize,Vec<Counts>)>) -> Result<()> {
    use hdf5::{File as H5File, Group, types::{VarLenUnicode, VarLenAscii}};
    use std::str::FromStr;

    fn str_attr(group:&Group, name:&str) -> Result<String> {
        let a=group.attr(name)?;
        if let Ok(v)=a.read_scalar::<VarLenUnicode>() { return Ok(v.as_str().to_owned()); }
        Ok(a.read_scalar::<VarLenAscii>()?.as_str().to_owned())
    }
    fn index_values(group:&Group, name:&str) -> Result<Vec<String>> {
        let ds=group.dataset(name)?;
        if let Ok(v)=ds.read_raw::<VarLenUnicode>() {
            return Ok(v.iter().map(|x|x.as_str().to_owned()).collect());
        }
        Ok(ds.read_raw::<VarLenAscii>()?.iter().map(|x|x.as_str().to_owned()).collect())
    }
    fn add_str(group:&Group, name:&str, values:Vec<String>) -> Result<()> {
        let v=values.iter().map(|s|VarLenUnicode::from_str(s)).collect::<std::result::Result<Vec<_>,_>>()?;
        let ds=group.new_dataset::<VarLenUnicode>().shape(v.len()).create(name)?;
        ds.write(&v)?;
        ds.new_attr::<VarLenUnicode>().shape(()).create("encoding-type")?
            .write_scalar(&VarLenUnicode::from_str("string-array")?)?;
        ds.new_attr::<VarLenUnicode>().shape(()).create("encoding-version")?
            .write_scalar(&VarLenUnicode::from_str("0.2.0")?)?;
        Ok(())
    }
    fn add_int(group:&Group, name:&str, values:Vec<i64>) -> Result<()> {
        let ds=group.new_dataset::<i64>().shape(values.len()).create(name)?;
        ds.write(&values)?;
        ds.new_attr::<VarLenUnicode>().shape(()).create("encoding-type")?
            .write_scalar(&VarLenUnicode::from_str("array")?)?;
        ds.new_attr::<VarLenUnicode>().shape(()).create("encoding-version")?
            .write_scalar(&VarLenUnicode::from_str("0.2.0")?)?;
        Ok(())
    }

    if input==output { bail!("refusing to overwrite input H5AD"); }
    if output.exists() { bail!("output H5AD already exists: {}",output.display()); }
    // Preserve all unknown AnnData extensions, layers, matrices, metadata, and
    // compression byte-for-byte at the HDF5-object level: modify only /obs.
    let staging=output.with_extension("h5ad.partial");
    if staging.exists() { bail!("staging file already exists: {}",staging.display()); }
    fs::copy(input,&staging).with_context(||format!("copying H5AD {}",input.display()))?;
    let result=(|| -> Result<()> {
        let f=H5File::open_rw(&staging)?;
        let obs=f.group("obs")?;
        let index_name=str_attr(&obs,"_index").unwrap_or_else(|_|"_index".into());
        let index=index_values(&obs,&index_name)?;
        let keys=index.iter().map(|s|normalize(s)).collect::<Vec<_>>();
        let names=std::cell::RefCell::new(Vec::<String>::new());
        let add_s=|name:String, vals:Vec<String>| -> Result<()> {
            if obs.link_exists(&name) { bail!("existing obs column {name:?}; refusing overwrite"); }
            add_str(&obs,&name,vals)?; names.borrow_mut().push(name); Ok(())
        };
        let add_i=|name:String, vals:Vec<i64>| -> Result<()> {
            if obs.link_exists(&name) { bail!("existing obs column {name:?}; refusing overwrite"); }
            add_int(&obs,&name,vals)?; names.borrow_mut().push(name); Ok(())
        };
        add_s("snp_genotype".into(),keys.iter().map(|k|observed.get(k).map(|x|x.0.clone()).unwrap_or_else(||"NO_INFO".into())).collect())?;
        add_i("snp_ref_total".into(),keys.iter().map(|k|observed.get(k).map_or(0,|x|x.1.r as i64)).collect())?;
        add_i("snp_alt_total".into(),keys.iter().map(|k|observed.get(k).map_or(0,|x|x.1.a as i64)).collect())?;
        add_i("snp_mut_snps".into(),keys.iter().map(|k|observed.get(k).map_or(0,|x|x.2 as i64)).collect())?;
        add_i("snp_wt_snps".into(),keys.iter().map(|k|observed.get(k).map_or(0,|x|x.3 as i64)).collect())?;
        for (i,snp) in snps.iter().enumerate() {
            let name=safe_name(snp);
            add_i(format!("snp_{name}_ref_count"),keys.iter().map(|k|observed.get(k).map_or(0,|x|x.4[i].r as i64)).collect())?;
            add_i(format!("snp_{name}_alt_count"),keys.iter().map(|k|observed.get(k).map_or(0,|x|x.4[i].a as i64)).collect())?;
        }
        // AnnData dataframe readers use column-order, not just HDF5 member names.
        let attr=obs.attr("column-order")?;
        let mut order=if let Ok(v)=attr.read_raw::<VarLenUnicode>() {
            v.into_iter().map(|x|x.as_str().to_owned()).collect::<Vec<_>>()
        } else {
            attr.read_raw::<VarLenAscii>()?.into_iter().map(|x|x.as_str().to_owned()).collect()
        };
        order.extend(names.into_inner());
        // HDF5 attributes are fixed-size: remove and recreate with the new length.
        // hdf5-rust does not expose attr deletion directly; H5Adelete is required.
        let cname=std::ffi::CString::new("column-order")?;
        let rc=unsafe { hdf5_sys::h5a::H5Adelete(obs.id(),cname.as_ptr()) };
        if rc<0 { bail!("failed to replace /obs column-order"); }
        let values=order.iter().map(|s|VarLenUnicode::from_str(s)).collect::<std::result::Result<Vec<_>,_>>()?;
        obs.new_attr::<VarLenUnicode>().shape(values.len()).create("column-order")?.write(&values)?;
        f.flush()?;
        Ok(())
    })();
    if let Err(err)=result { let _=fs::remove_file(&staging); return Err(err); }
    fs::rename(&staging,output)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test] fn calls() {
        assert_eq!(call(Counts{r:0,a:1},4),"MUT");
        assert_eq!(call(Counts{r:4,a:0},4),"WT");
        assert_eq!(call(Counts{r:0,a:0},4),"NO_INFO");
    }
    #[test] fn barcodes() {
        assert_eq!(normalize("AAACCCGGGTTTAAAA-1"),"AAACCCGGGTTTAAAA");
    }
}
