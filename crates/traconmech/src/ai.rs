//! Ollama investigation loop: the model requests evidence; Rust computes it.
use anyhow::{bail, Context, Result};
use sc_analysis::{SingleCellData, write_group_stats};
use serde_json::{json, Value};
use std::fs::{self, OpenOptions};
use std::io::{self, BufRead, Write};
use std::path::Path;
use std::collections::BTreeMap;

const MAX_TOOL_ROUNDS: usize = 8;

fn append(path: &Path, text: &str) -> Result<()> {
    let mut file = OpenOptions::new().create(true).append(true).open(path)?;
    file.write_all(text.as_bytes())?;
    Ok(())
}
fn json_line(path: &Path, value: &Value) -> Result<()> {
    append(path, &format!("{}\n", serde_json::to_string(value)?))
}
fn failure(out: &Path, stage: &str, raw: &str, error: &str) {
    let _ = append(&out.join("error.log"),
        &format!("\n[{stage}] {error}\nRaw response/request:\n{raw}\n"));
}

fn read_cluster_summary(out: &Path, group: &str) -> Result<Value> {
    let dir = out.join("stats");
    let labels = fs::read_to_string(dir.join("cluster_labels.tsv"))?;
    let mut names = BTreeMap::new();
    for line in labels.lines().skip(1) {
        if let Some((id, label)) = line.split_once('\t') { names.insert(id.to_string(), label.to_string()); }
    }
    let summary = fs::read_to_string(dir.join("marker_summary.tsv"))?;
    let mut clusters = Vec::new();
    for line in summary.lines().skip(1) {
        let fields: Vec<_> = line.split('\t').collect();
        anyhow::ensure!(fields.len() == 5, "invalid marker_summary.tsv row: {line}");
        clusters.push(json!({"cluster":names.get(fields[0]).context("unmapped cluster label")?,
            "n_cells":fields[1].parse::<usize>()?,
            "fraction_cells":fields[2].parse::<f64>()?,
            "top20_up":fields[3].split(',').filter(|x| !x.is_empty()).collect::<Vec<_>>(),
            "top20_down":fields[4].split(',').filter(|x| !x.is_empty()).collect::<Vec<_>>() }));
    }
    Ok(json!({"grouping":group,"clusters":clusters,
        "source":"sc-analysis standard marker_summary.tsv; ranked by adjusted p-value and absolute log2 fold change",
        "statistics":"pseudo-sample Wilcoxon; not biological replicate-level differential expression"}))
}
fn evidence(data: &SingleCellData, group: &str, out: &Path, name: &str, args: &Value) -> Result<Value> {
    match name {
        "get_cluster_signatures" => read_cluster_summary(out, group),
        "list_populations" => {
            let groups = data.annotations.groups(group)?;
            Ok(json!({"grouping": group, "populations": groups.iter()
                .map(|(label, cells)| json!({"population":label,"cells":cells.len()}))
                .collect::<Vec<_>>(), "available_annotations":data.annotations.headers()}))
        }
        "get_marker_expression" => {
            let genes = args.get("genes").and_then(Value::as_array)
                .context("get_marker_expression requires genes array")?;
            anyhow::ensure!(!genes.is_empty() && genes.len() <= 100,
                "genes must contain 1 to 100 names");
            let groups = data.annotations.groups(group)?;
            let mut results = Vec::new();
            for gene in genes {
                let name = gene.as_str().context("gene must be a string")?;
                let Some(idx) = data.feature_index(name) else {
                    results.push(json!({"gene":name,"status":"not_found"}));
                    continue;
                };
                let row = data.matrix.outer_view(idx).context("missing matrix feature row")?;
                let mut per_group = Vec::new();
                for (label, cells) in &groups {
                    let mut detected = 0usize;
                    let mut sum = 0.0f64;
                    let mut sum_sq = 0.0f64;
                    for &cell in cells {
                        let value = row.get(cell).copied().unwrap_or(0.0) as f64;
                        if value > 0.0 { detected += 1; }
                        sum += value;
                        sum_sq += value * value;
                    }
                    let n = cells.len() as f64;
                    let mean = sum / n;
                    let sd = if cells.len() > 1 {
                        ((sum_sq - sum * sum / n).max(0.0) / (n - 1.0)).sqrt()
                    } else { 0.0 };
                    per_group.push(json!({"population":label,"cells":cells.len(),
                        "mean":mean,"sd":sd,"detection_fraction":detected as f64 / n}));
                }
                results.push(json!({"gene":name,"status":"found","populations":per_group}));
            }
            Ok(json!({"grouping":group, "measurements":results,
                "note":"mean and sample SD include zero-valued cells and use input H5AD X values, which may be normalized/log-transformed"}))
        }
        _ => bail!("unsupported tool {name:?}"),
    }
}
fn tool_specs() -> Value {
    json!([
      {"type":"function","function":{"name":"get_cluster_signatures",
        "description":"Read the standard sc-analysis marker_summary.tsv: top 20 enriched and depleted genes for each existing cluster. Use this FIRST to identify cell types and developmental relationships.",
        "parameters":{"type":"object","properties":{}}}},
      {"type":"function","function":{"name":"list_populations",
        "description":"List measured populations, cell counts and available cell annotations.",
        "parameters":{"type":"object","properties":{}}}},
      {"type":"function","function":{"name":"get_marker_expression",
        "description":"Follow-up measurements: mean and sample standard deviation (including zeros), plus detection fraction, per cluster. Start with get_cluster_signatures instead.",
        "parameters":{"type":"object","properties":{"genes":{"type":"array","items":{"type":"string"}}},
          "required":["genes"]}}}
    ])
}
pub fn investigate(data: &SingleCellData, group: &str, out: &Path,
                   model: &str, url: &str, initial: Option<&str>) -> Result<()> {
    fs::create_dir_all(out).with_context(|| format!("creating {}", out.display()))?;
    let conversation = out.join("conversation.md");
    let evidence_log = out.join("evidence.jsonl");
    let errors = out.join("error.log");
    // Create error.log even for a successful session, so scripts can check it.
    OpenOptions::new().create(true).append(true).open(&errors)?;
    let groups = data.annotations.groups(group)?;
    let stats_dir = out.join("stats");
    if !stats_dir.join("marker_summary.tsv").exists() || !stats_dir.join("cluster_labels.tsv").exists() {
        println!("Computing standard sc-analysis marker statistics for {group} ({} populations)...", groups.len());
        write_group_stats(&stats_dir, data, group, 10)?;
    }
    let signatures = read_cluster_summary(out, group)?;
    json_line(&evidence_log, &json!({"tool":"get_cluster_signatures","result":signatures}))?;
    append(&conversation, "\n## Dataset cluster signatures\n\nStandard sc-analysis marker_summary.tsv was computed for the supplied grouping. See stats/marker_summary.tsv and evidence.jsonl.\n")?;
    let mut messages = vec![json!({"role":"system","content":format!(
        "You are TraConMech, an evidence-driven biological investigator. \
         Dataset: {} cells, {} genes, {} populations grouped by {}. \
         Never invent measurements, population identities, or sample diagnoses. \
         The supplied cluster signatures are your primary evidence. First assign tentative cell identities \
         using top enriched/depleted genes, then propose a directed, possibly branching developmental \
         order using the ACTUAL cluster identifiers. Explain evidence for every proposed edge, and \
         distinguish inferred order from measured lineage tracing. Request compact follow-up mean/SD \
         marker expression only where needed. Clearly distinguish observations from hypotheses. \
         If a tool is unavailable, say so rather than guessing.",
        data.n_cells(), data.n_features(), groups.len(), group)})];
    println!("TraConMech AI — {model} | {} populations | output {}", groups.len(), out.display());
    let stdin = io::stdin();
    let mut lines = stdin.lock().lines();
    let mut pending = initial.map(str::to_owned);
    loop {
        let question = match pending.take() {
            Some(question) => question,
            None => {
                print!("\nTraConMech> ");
                io::stdout().flush()?;
                match lines.next() {
                    Some(Ok(line)) => line,
                    Some(Err(e)) => return Err(e.into()),
                    None => break,
                }
            }
        };
        let question = question.trim();
        if question.eq_ignore_ascii_case("exit") || question.eq_ignore_ascii_case("quit") { break; }
        if question.is_empty() { continue; }
        append(&conversation, &format!("\n## Question\n\n{question}\n"))?;
        messages.push(json!({"role":"user","content":format!(
            "Question: {question}\n\nMeasured cluster signatures from sc-analysis (use these before other tools):\n{}",
            serde_json::to_string(&signatures)?)}));
        let mut answered = false;
        for round in 0..MAX_TOOL_ROUNDS {
            let payload = json!({"model":model,"stream":false,"messages":messages,
                                 "tools":tool_specs()});
            let (response, raw) = match ai_llm_con::chat(out, url, &payload) {
                Ok(value) => value,
                Err(e) => {
                    failure(out,"ollama_http_or_json",&payload.to_string(),&e.to_string());
                    return Err(e);
                }
            };
            let Some(message) = response.get("message").filter(|v| v.is_object()) else {
                failure(out,"ollama_schema",&raw,"missing message object");
                bail!("Ollama response missing message; details in {}", errors.display());
            };
            let Some(calls) = message.get("tool_calls").and_then(Value::as_array) else {
                if message.get("tool_calls").is_some() {
                    failure(out,"ollama_schema",&raw,"tool_calls is not an array");
                    bail!("Malformed tool_calls; details in {}", errors.display());
                }
                let Some(answer) = message.get("content").and_then(Value::as_str) else {
                    failure(out,"ollama_schema",&raw,"missing text content");
                    bail!("Missing Ollama answer; details in {}", errors.display());
                };
                if answer.trim().is_empty() {
                    failure(out,"ollama_schema",&raw,"empty answer without tool calls");
                    bail!("Empty Ollama answer; details in {}", errors.display());
                }
                println!("\n{answer}");
                if question.to_lowercase().contains("develop") || question.to_lowercase().contains("order") {
                    append(&out.join("developmental_order.md"), &format!("\n## {question}\n\n{answer}\n"))?;
                }
                append(&conversation, &format!("\n## Answer\n\n{answer}\n"))?;
                messages.push(message.clone());
                answered = true;
                break;
            };
            if calls.is_empty() {
                failure(out,"ollama_schema",&raw,"empty tool_calls array");
                bail!("Empty tool calls; details in {}", errors.display());
            }
            messages.push(message.clone());
            for call in calls {
                let Some(function) = call.get("function") else {
                    failure(out,"tool_schema",&raw,"tool call missing function");
                    bail!("Malformed tool call; details in {}", errors.display());
                };
                let Some(name) = function.get("name").and_then(Value::as_str) else {
                    failure(out,"tool_schema",&raw,"tool call missing function name");
                    bail!("Malformed tool call; details in {}", errors.display());
                };
                let args = function.get("arguments").cloned().unwrap_or(json!({}));
                let args = if let Some(s) = args.as_str() {
                    match serde_json::from_str::<Value>(s) {
                        Ok(value) => value,
                        Err(e) => {
                            failure(out,"tool_arguments",s,&e.to_string());
                            bail!("Malformed tool arguments; details in {}", errors.display());
                        }
                    }
                } else { args };
                println!("  → {name}({args})");
                let result = match evidence(data, group, out, name, &args) {
                    Ok(value) => value,
                    Err(e) => {
                        failure(out,"tool_execution",&call.to_string(),&e.to_string());
                        json!({"error":e.to_string()})
                    }
                };
                json_line(&evidence_log,&json!({"round":round,"tool":name,"arguments":args,"result":result}))?;
                append(&conversation,&format!("\n### Evidence: {name}\n\nArguments: `{args}`\n\nSee `evidence.jsonl` for measured results.\n"))?;
                if url == "openai" {
                    let id = call.get("id").and_then(Value::as_str)
                        .context("OpenAI tool call missing id")?;
                    messages.push(json!({"role":"tool","tool_call_id":id,"content":result.to_string()}));
                } else {
                    messages.push(json!({"role":"tool","tool_name":name,"content":result.to_string()}));
                }
            }
        }
        if !answered {
            failure(out,"tool_limit",&question.to_string(),"model exceeded tool-call round limit");
            eprintln!("Investigation reached the {MAX_TOOL_ROUNDS}-round limit; see {}",errors.display());
        }
        if initial.is_some() { break; }
    }
    Ok(())
}

/// Reuse complete validated annotations; otherwise call the existing Ollama client.
/// Do not silently accept an incomplete or stale developmental order.
pub fn developmental_order_cached(out: &Path, model: &str, url: &str,
                                  allowed: &[String], primer: Option<&str>) -> Result<Vec<String>> {
    let annotations = out.join("cluster_ai_annotations.tsv");
    let primer_path = out.join("developmental_ai_primer.txt");
    let requested_primer = primer.unwrap_or("");
    let same_primer = fs::read_to_string(&primer_path).ok().as_deref() == Some(requested_primer);
    let backend_identity = format!("{url}\n{model}");
    let same_backend = fs::read_to_string(out.join("developmental_ai_backend.txt"))
        .ok().as_deref() == Some(backend_identity.as_str());
    let order_file = out.join("developmental_display_order.tsv");
    if annotations.is_file() && order_file.is_file() && same_primer && same_backend {
        let lines = fs::read_to_string(&order_file)?;
        let mut order = Vec::new();
        for line in lines.lines().skip(1) {
            let parts: Vec<_> = line.split('\t').collect();
            anyhow::ensure!(parts.len() == 2, "invalid developmental order row");
            order.push(parts[1].to_owned());
        }
        let expected: std::collections::BTreeSet<_> = allowed.iter().collect();
        let unique: std::collections::BTreeSet<_> = order.iter().collect();
        let text = fs::read_to_string(&annotations)?;
        let annotated: std::collections::BTreeSet<_> = text.lines().skip(1)
            .filter_map(|line| line.split('\t').next()).collect();
        let expected_names: std::collections::BTreeSet<_> = allowed.iter().map(String::as_str).collect();
        if order.len() == allowed.len() && unique == expected && annotated == expected_names {
            eprintln!("[traconmech] developmental AI: reusing {}", annotations.display());
            return Ok(order);
        }
        eprintln!("[traconmech] developmental AI: cached clusters differ; querying Ollama");
    }
    eprintln!("[traconmech] developmental AI: requesting cluster identities and order from Ollama");
    developmental_order(out, model, url, allowed, primer)
}

/// Reuse the existing Ollama transport and sc-analysis marker evidence.
/// AI proposes annotations and a *linear display order*; it never changes
/// measured cluster identities, module membership or reference matching.
pub fn developmental_order(out: &Path, model: &str, url: &str,
                           allowed: &[String], primer: Option<&str>) -> Result<Vec<String>> {
    let signatures = read_cluster_summary_from_markers(&out.join("cluster_markers"))?;
    let context = primer.filter(|p| !p.trim().is_empty()).unwrap_or("No experiment description supplied.");
    eprintln!("[traconmech] developmental AI: sending request to Ollama (timeout: 1200 seconds)");
    let prompt = json!({"role":"user","content":format!(
        "Experiment background (context, not marker evidence):\n{}\n\nUsing ONLY these sc-analysis cluster markers, propose tentative cell identities and a developmental display order.\n\n{}\n\nAllowed cluster IDs: {}\n\nReturn ONLY a JSON object with keys 'annotations' (array of objects each with cluster, cell_type, lineage, confidence, rationale) and 'order' (array of cluster ID strings). Include each allowed cluster exactly once in order. Group related lineages; do NOT imply that distinct lineages form a single differentiation path. Label stress/cycling states as states, not stages. Do not fabricate certainty.",
        context, serde_json::to_string(&signatures)?, serde_json::to_string(allowed)?) });
    let payload = json!({"model":model,"stream":false,"format":"json","messages":[
        {"role":"system","content":"You are an evidence-driven hematopoiesis cell annotation assistant. Return strictly valid JSON. Developmental order is a tentative display arrangement, not measured pseudotime or lineage tracing."},prompt]});
    let (response, raw) = ai_llm_con::chat(out, url, &payload)?;
    fs::write(out.join("developmental_ai_response.json"), &raw)?;
    let content = response.get("message").and_then(|m| m.get("content"))
        .and_then(Value::as_str).context("Ollama missing developmental annotation content")?;
    let parsed: Value = serde_json::from_str(content).context("invalid AI developmental JSON")?;
    let order = parsed.get("order").and_then(Value::as_array).context("missing order array")?
        .iter().map(|v| v.as_str().map(str::to_owned).context("order must contain strings"))
        .collect::<Result<Vec<_>>>()?;
    let unique: std::collections::BTreeSet<_> = order.iter().collect();
    let expected: std::collections::BTreeSet<_> = allowed.iter().collect();
    anyhow::ensure!(order.len() == allowed.len() && unique == expected,
        "AI order must contain every plotted cluster exactly once; raw response saved");
    let annotations = parsed.get("annotations").and_then(Value::as_array)
        .context("missing annotations array")?;
    let mut rows = BTreeMap::new();
    for entry in annotations {
        let field = |name: &str| -> Result<&str> {
            entry.get(name).and_then(Value::as_str)
                .with_context(|| format!("missing annotation field {name}"))
        };
        let cluster = field("cluster")?;
        anyhow::ensure!(expected.contains(&cluster.to_owned()), "unknown AI cluster {cluster}");
        anyhow::ensure!(!rows.contains_key(cluster), "duplicate annotation for {cluster}");
        let values = [field("cell_type")?, field("lineage")?, field("confidence")?, field("rationale")?];
        anyhow::ensure!(values.iter().all(|v| !v.chars().any(|c| c == '\t' || c == '\n' || c == '\r')), "invalid TSV annotation");
        rows.insert(cluster.to_owned(), values.map(str::to_owned));
    }
    anyhow::ensure!(rows.len() == allowed.len(), "AI did not annotate every cluster");
    let mut tsv = String::from("cluster\tcell_type\tlineage\tconfidence\trationale\n");
    for (cluster, values) in rows {
        tsv.push_str(&format!("{}\t{}\n", cluster, values.join("\t")));
    }
    fs::write(out.join("cluster_ai_annotations.tsv"), tsv)?;
    fs::write(out.join("developmental_ai_primer.txt"), primer.unwrap_or(""))?;
    fs::write(out.join("developmental_ai_backend.txt"), format!("{url}\n{model}"))?;
    fs::write(out.join("developmental_display_order.tsv"),
        format!("display_rank\tcluster\n{}", order.iter().enumerate()
            .map(|(i,c)| format!("{}\t{}\n", i+1,c)).collect::<String>()))?;
    Ok(order)
}

fn read_cluster_summary_from_markers(dir: &Path) -> Result<Value> {
    let labels = fs::read_to_string(dir.join("cluster_labels.tsv"))?;
    let names: BTreeMap<_,_> = labels.lines().skip(1)
        .filter_map(|line| line.split_once('\t'))
        .map(|(id,name)| (id.to_owned(),name.to_owned())).collect();
    let summary = fs::read_to_string(dir.join("marker_summary.tsv"))?;
    let mut clusters = Vec::new();
    for line in summary.lines().skip(1) {
        let fields: Vec<_> = line.split('\t').collect();
        anyhow::ensure!(fields.len()==5, "invalid marker summary row");
        clusters.push(json!({"cluster":names.get(fields[0]).context("unmapped cluster")?,
            "n_cells":fields[1],"top20_up":fields[3],"top20_down":fields[4]}));
    }
    Ok(json!({"clusters":clusters,"source":"sc-analysis marker_summary.tsv"}))
}
