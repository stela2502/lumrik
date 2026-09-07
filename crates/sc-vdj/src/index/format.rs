use super::{Chain, SegmentKind, Strand, VdjIndex, VdjSegment};
use anyhow::{bail, Context, Result};
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;

const MAGIC: &[u8; 8] = b"LVDJIDX5";
const VERSION: u32 = 5;

pub(super) fn save(index: &VdjIndex, path: &Path) -> Result<()> {
    let mut w =
        BufWriter::new(File::create(path).with_context(|| format!("creating {}", path.display()))?);
    w.write_all(MAGIC)?;
    write_u32(&mut w, VERSION)?;
    write_u32(&mut w, index.segments.len() as u32)?;
    for s in &index.segments {
        write_str(&mut w, &s.name)?;
        write_str(&mut w, &s.transcript_id)?;
        write_str(&mut w, &s.gene_id)?;
        write_u8(&mut w, s.chain.code())?;
        write_u8(
            &mut w,
            match s.kind {
                SegmentKind::V => 0,
                SegmentKind::D => 1,
                SegmentKind::J => 2,
                SegmentKind::C => 3,
            },
        )?;
        write_str(&mut w, &s.chromosome)?;
        write_u32(&mut w, s.start)?;
        write_u32(&mut w, s.end)?;
        write_u8(&mut w, if s.strand == Strand::Minus { 1 } else { 0 })?;
        write_bytes(&mut w, &s.sequence)?;
    }
    w.flush()?;
    Ok(())
}

pub(super) fn load(path: &Path) -> Result<VdjIndex> {
    let mut r =
        BufReader::new(File::open(path).with_context(|| format!("opening {}", path.display()))?);
    let mut magic = [0u8; 8];
    r.read_exact(&mut magic)?;
    if &magic != MAGIC {
        bail!(
            "{} is not a clean sc-vdj v5 index; regenerate it with vdj-index",
            path.display()
        )
    }
    let version = read_u32(&mut r)?;
    if version != VERSION {
        bail!("unsupported VDJ index version {version}")
    }
    let n = read_u32(&mut r)? as usize;
    let mut segments = Vec::with_capacity(n);
    for i in 0..n {
        let name = read_str(&mut r)?;
        let transcript_id = read_str(&mut r)?;
        let gene_id = read_str(&mut r)?;
        let chain = Chain::from_code(read_u8(&mut r)?)?;
        let kind = match read_u8(&mut r)? {
            0 => SegmentKind::V,
            1 => SegmentKind::D,
            2 => SegmentKind::J,
            3 => SegmentKind::C,
            x => bail!("invalid segment kind {x}"),
        };
        let chromosome = read_str(&mut r)?;
        let start = read_u32(&mut r)?;
        let end = read_u32(&mut r)?;
        let strand = if read_u8(&mut r)? == 1 {
            Strand::Minus
        } else {
            Strand::Plus
        };
        let sequence = read_bytes(&mut r)?;
        segments.push(VdjSegment {
            id: i as u16,
            name,
            transcript_id,
            gene_id,
            chain,
            kind,
            chromosome,
            start,
            end,
            strand,
            sequence,
        });
    }
    VdjIndex::from_segments(segments)
}
fn write_u8<W: Write>(w: &mut W, x: u8) -> Result<()> {
    w.write_all(&[x])?;
    Ok(())
}
fn read_u8<R: Read>(r: &mut R) -> Result<u8> {
    let mut b = [0];
    r.read_exact(&mut b)?;
    Ok(b[0])
}
fn write_u32<W: Write>(w: &mut W, x: u32) -> Result<()> {
    w.write_all(&x.to_le_bytes())?;
    Ok(())
}
fn read_u32<R: Read>(r: &mut R) -> Result<u32> {
    let mut b = [0; 4];
    r.read_exact(&mut b)?;
    Ok(u32::from_le_bytes(b))
}
fn write_bytes<W: Write>(w: &mut W, x: &[u8]) -> Result<()> {
    write_u32(w, x.len() as u32)?;
    w.write_all(x)?;
    Ok(())
}
fn read_bytes<R: Read>(r: &mut R) -> Result<Vec<u8>> {
    let n = read_u32(r)? as usize;
    let mut x = vec![0; n];
    r.read_exact(&mut x)?;
    Ok(x)
}
fn write_str<W: Write>(w: &mut W, x: &str) -> Result<()> {
    write_bytes(w, x.as_bytes())
}
fn read_str<R: Read>(r: &mut R) -> Result<String> {
    String::from_utf8(read_bytes(r)?).context("invalid UTF-8 in VDJ index")
}
