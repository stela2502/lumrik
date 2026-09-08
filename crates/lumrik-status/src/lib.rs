use std::io::{Read, Write};
use std::net::{SocketAddr, TcpListener, TcpStream};
use std::process::Command;
use std::sync::{Arc, RwLock};
use std::thread::{self, JoinHandle};

use anyhow::{Context, Result};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StatusMetric {
    pub label: String,
    pub value: String,
}

impl StatusMetric {
    pub fn new(label: impl Into<String>, value: impl Into<String>) -> Self {
        Self {
            label: label.into(),
            value: value.into(),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StatusSection {
    pub title: String,
    pub metrics: Vec<StatusMetric>,
}

impl StatusSection {
    pub fn new(title: impl Into<String>, metrics: Vec<StatusMetric>) -> Self {
        Self {
            title: title.into(),
            metrics,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ServerSnapshot {
    pub title: String,
    pub subtitle: String,
    pub started_unix_ms: u128,
    pub finished_unix_ms: Option<u128>,
    pub stage: String,
    pub public_url: Option<String>,
    pub sections: Vec<StatusSection>,
}

pub trait ServerContent {
    fn server_snapshot(&self) -> ServerSnapshot;
}

pub struct StatusServer {
    handle: JoinHandle<()>,
    addr: SocketAddr,
}

impl StatusServer {
    pub fn addr(&self) -> SocketAddr {
        self.addr
    }

    pub fn join(self) -> thread::Result<()> {
        self.handle.join()
    }
}

pub fn spawn_status_server<C>(state: Arc<RwLock<C>>, addr: SocketAddr) -> Result<StatusServer>
where
    C: ServerContent + Send + Sync + 'static,
{
    let listener = TcpListener::bind(addr)
        .with_context(|| format!("binding Lumrik status server to {addr}"))?;
    let addr = listener
        .local_addr()
        .context("reading Lumrik status server address")?;

    let handle = thread::spawn(move || {
        for connection in listener.incoming() {
            match connection {
                Ok(stream) => {
                    if let Err(error) = handle_connection(stream, &state) {
                        eprintln!("[lumrik-status] request failed: {error}");
                    }
                }
                Err(error) => {
                    eprintln!("[lumrik-status] connection failed: {error}");
                }
            }
        }
    });

    Ok(StatusServer { handle, addr })
}

fn handle_connection<C>(mut stream: TcpStream, state: &Arc<RwLock<C>>) -> Result<()>
where
    C: ServerContent,
{
    let mut buffer = [0u8; 4096];
    let n = stream
        .read(&mut buffer)
        .context("reading status-server request")?;
    if n == 0 {
        return Ok(());
    }

    let request = String::from_utf8_lossy(&buffer[..n]);
    let request_line = request.lines().next().unwrap_or_default();
    let mut parts = request_line.split_whitespace();
    let method = parts.next().unwrap_or_default();
    let path = parts.next().unwrap_or_default();

    match (method, path) {
        ("GET", "/") => write_response(
            &mut stream,
            "200 OK",
            "text/html; charset=utf-8",
            DASHBOARD_HTML,
        )?,
        ("GET", "/health") => {
            write_response(&mut stream, "200 OK", "text/plain; charset=utf-8", "OK\n")?
        }
        ("GET", "/status") => {
            let body = status_json(state)?;
            write_response(
                &mut stream,
                "200 OK",
                "application/json; charset=utf-8",
                &body,
            )?;
        }
        _ => write_response(
            &mut stream,
            "404 Not Found",
            "text/plain; charset=utf-8",
            "Not Found\n",
        )?,
    }

    Ok(())
}

fn status_json<C>(state: &Arc<RwLock<C>>) -> Result<String>
where
    C: ServerContent,
{
    let state = state
        .read()
        .map_err(|_| anyhow::anyhow!("status-server state lock poisoned"))?;
    Ok(snapshot_json(&state.server_snapshot()))
}

fn snapshot_json(snapshot: &ServerSnapshot) -> String {
    let finished = snapshot
        .finished_unix_ms
        .map(|value| value.to_string())
        .unwrap_or_else(|| "null".to_string());
    let public_url = json_optional_string(snapshot.public_url.as_deref());
    let sections = snapshot
        .sections
        .iter()
        .map(section_json)
        .collect::<Vec<_>>()
        .join(",");

    format!(
        concat!(
            "{{",
            "\"title\":\"{}\",",
            "\"subtitle\":\"{}\",",
            "\"started_unix_ms\":{},",
            "\"finished_unix_ms\":{},",
            "\"stage\":\"{}\",",
            "\"public_url\":{},",
            "\"sections\":[{}]",
            "}}\n"
        ),
        escape_json(&snapshot.title),
        escape_json(&snapshot.subtitle),
        snapshot.started_unix_ms,
        finished,
        escape_json(&snapshot.stage),
        public_url,
        sections,
    )
}

fn section_json(section: &StatusSection) -> String {
    let metrics = section
        .metrics
        .iter()
        .map(|metric| {
            format!(
                "{{\"label\":\"{}\",\"value\":\"{}\"}}",
                escape_json(&metric.label),
                escape_json(&metric.value),
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    format!(
        "{{\"title\":\"{}\",\"metrics\":[{}]}}",
        escape_json(&section.title),
        metrics,
    )
}

fn json_optional_string(value: Option<&str>) -> String {
    match value {
        Some(value) => format!("\"{}\"", escape_json(value)),
        None => "null".to_string(),
    }
}

fn write_response(
    stream: &mut TcpStream,
    status: &str,
    content_type: &str,
    body: &str,
) -> Result<()> {
    let response = format!(
        concat!(
            "HTTP/1.1 {}\r\n",
            "Content-Type: {}\r\n",
            "Content-Length: {}\r\n",
            "Connection: close\r\n",
            "Cache-Control: no-store\r\n",
            "\r\n",
            "{}"
        ),
        status,
        content_type,
        body.len(),
        body,
    );

    stream
        .write_all(response.as_bytes())
        .context("writing status-server response")?;
    stream.flush().context("flushing status-server response")?;
    Ok(())
}

fn escape_json(value: &str) -> String {
    let mut out = String::with_capacity(value.len());
    for ch in value.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            ch if ch.is_control() => {
                use std::fmt::Write as _;
                let _ = write!(out, "\\u{:04x}", ch as u32);
            }
            ch => out.push(ch),
        }
    }
    out
}


pub fn snapshot_html(snapshot: &ServerSnapshot) -> String {
    let elapsed_ms = snapshot
        .finished_unix_ms
        .unwrap_or(snapshot.started_unix_ms)
        .saturating_sub(snapshot.started_unix_ms);
    let elapsed = format_elapsed_ms(elapsed_ms);
    let mut sections = String::new();
    for section in &snapshot.sections {
        sections.push_str("<section class=\"panel\"><h2>");
        sections.push_str(&escape_html(&section.title));
        sections.push_str("</h2><div class=\"rows\">");
        for metric in &section.metrics {
            sections.push_str("<div class=\"row\"><div class=\"label\">");
            sections.push_str(&escape_html(&metric.label));
            sections.push_str("</div><div class=\"metric\">");
            sections.push_str(&escape_html(&metric.value));
            sections.push_str("</div></div>");
        }
        sections.push_str("</div></section>");
    }

    format!(
        r#"<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{}</title>
<style>{}</style>
</head>
<body>
<main>
<header><div><h1>{}</h1><p class="subtitle">{}</p></div><div class="badge">final report</div></header>
<div class="summary"><div><span>Stage</span><strong>{}</strong></div><div><span>Elapsed</span><strong>{}</strong></div></div>
<div class="sections">{}</div>
<footer>Static Lumrik run report</footer>
</main>
</body>
</html>
"#,
        escape_html(&snapshot.title),
        DASHBOARD_CSS,
        escape_html(&snapshot.title),
        escape_html(&snapshot.subtitle),
        escape_html(&snapshot.stage),
        elapsed,
        sections,
    )
}

fn format_elapsed_ms(ms: u128) -> String {
    let total = ms / 1000;
    let seconds = total % 60;
    let minutes = (total / 60) % 60;
    let hours = total / 3600;
    format!("{hours:02}:{minutes:02}:{seconds:02}")
}

fn escape_html(value: &str) -> String {
    let mut out = String::with_capacity(value.len());
    for ch in value.chars() {
        match ch {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '\"' => out.push_str("&quot;"),
            '\'' => out.push_str("&#39;"),
            _ => out.push(ch),
        }
    }
    out
}

#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct MemoryStatus {
    pub process_rss_mib: f64,
    pub process_peak_rss_mib: f64,
    pub system_available_mib: f64,
}

pub fn memory_status() -> MemoryStatus {
    let (process_rss_mib, process_peak_rss_mib) = process_memory_mib();
    MemoryStatus {
        process_rss_mib,
        process_peak_rss_mib,
        system_available_mib: system_available_memory_mib(),
    }
}

pub fn public_hostname(override_hostname: Option<&str>) -> String {
    if let Some(hostname) = override_hostname {
        if !hostname.trim().is_empty() {
            return hostname.trim().to_string();
        }
    }

    if let Ok(hostname) = std::env::var("SLURMD_NODENAME") {
        if !hostname.trim().is_empty() {
            return hostname;
        }
    }

    if let Ok(output) = Command::new("hostname").arg("-f").output() {
        if output.status.success() {
            let hostname = String::from_utf8_lossy(&output.stdout).trim().to_string();
            if !hostname.is_empty() {
                return hostname;
            }
        }
    }

    if let Ok(hostname) = std::env::var("HOSTNAME") {
        if !hostname.trim().is_empty() {
            return hostname;
        }
    }

    "localhost".to_string()
}

fn process_memory_mib() -> (f64, f64) {
    let Ok(status) = std::fs::read_to_string("/proc/self/status") else {
        return (0.0, 0.0);
    };
    let mut rss_kib = 0u64;
    let mut peak_kib = 0u64;
    for line in status.lines() {
        if let Some(value) = parse_proc_kib(line, "VmRSS:") {
            rss_kib = value;
        } else if let Some(value) = parse_proc_kib(line, "VmHWM:") {
            peak_kib = value;
        }
    }
    (rss_kib as f64 / 1024.0, peak_kib as f64 / 1024.0)
}

fn system_available_memory_mib() -> f64 {
    let Ok(meminfo) = std::fs::read_to_string("/proc/meminfo") else {
        return 0.0;
    };
    meminfo
        .lines()
        .find_map(|line| parse_proc_kib(line, "MemAvailable:"))
        .map(|kib| kib as f64 / 1024.0)
        .unwrap_or(0.0)
}

fn parse_proc_kib(line: &str, key: &str) -> Option<u64> {
    line.strip_prefix(key)?
        .split_whitespace()
        .next()?
        .parse()
        .ok()
}

const DASHBOARD_CSS: &str = r#"
:root{color-scheme:dark}*{box-sizing:border-box}body{margin:0;background:#0e1116;color:#e8edf3;font-family:Inter,ui-sans-serif,system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}main{max-width:1040px;margin:0 auto;padding:32px 22px 48px}header{display:flex;justify-content:space-between;gap:20px;align-items:flex-start;margin-bottom:20px}h1{margin:0 0 5px;font-size:1.9rem;letter-spacing:-.025em}.subtitle{margin:0;color:#95a0ad}.badge{font-size:.78rem;color:#a9b4c0;border:1px solid #2a3440;border-radius:999px;padding:6px 10px;white-space:nowrap}.summary{display:grid;grid-template-columns:minmax(0,2fr) minmax(130px,1fr);gap:12px;margin:18px 0 20px}.summary>div,.panel{background:#151a21;border:1px solid #262e39;border-radius:10px}.summary>div{padding:13px 15px}.summary span{display:block;color:#8995a3;font-size:.76rem;text-transform:uppercase;letter-spacing:.06em;margin-bottom:5px}.summary strong{font-size:1rem;font-weight:650}.sections{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:14px}.panel{overflow:hidden}.panel h2{font-size:.9rem;margin:0;padding:12px 14px;border-bottom:1px solid #252d37;color:#c7d0da;background:#171d25}.rows{padding:3px 0}.row{display:grid;grid-template-columns:minmax(0,1fr) auto;gap:16px;padding:8px 14px;border-bottom:1px solid #202731;align-items:baseline}.row:last-child{border-bottom:0}.label{color:#909ba8;font-size:.84rem}.metric{font-variant-numeric:tabular-nums;text-align:right;font-weight:620;font-size:.88rem;overflow-wrap:anywhere}footer{margin-top:22px;color:#66717e;font-size:.78rem}@media(max-width:760px){main{padding:22px 14px 36px}.sections{grid-template-columns:1fr}.summary{grid-template-columns:1fr}header{display:block}.badge{display:inline-block;margin-top:10px}.row{grid-template-columns:1fr;gap:3px}.metric{text-align:left}}
"#;

const DASHBOARD_HTML: &str = r##"<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Lumrik</title>
<style>
:root{color-scheme:dark}*{box-sizing:border-box}body{margin:0;background:#0e1116;color:#e8edf3;font-family:Inter,ui-sans-serif,system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}main{max-width:1040px;margin:0 auto;padding:32px 22px 48px}header{display:flex;justify-content:space-between;gap:20px;align-items:flex-start;margin-bottom:20px}h1{margin:0 0 5px;font-size:1.9rem;letter-spacing:-.025em}.subtitle{margin:0;color:#95a0ad}.badge{font-size:.78rem;color:#a9b4c0;border:1px solid #2a3440;border-radius:999px;padding:6px 10px;white-space:nowrap}.summary{display:grid;grid-template-columns:minmax(0,2fr) minmax(130px,1fr);gap:12px;margin:18px 0 20px}.summary>div,.panel{background:#151a21;border:1px solid #262e39;border-radius:10px}.summary>div{padding:13px 15px}.summary span{display:block;color:#8995a3;font-size:.76rem;text-transform:uppercase;letter-spacing:.06em;margin-bottom:5px}.summary strong{font-size:1rem;font-weight:650}.sections{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:14px}.panel{overflow:hidden}.panel h2{font-size:.9rem;margin:0;padding:12px 14px;border-bottom:1px solid #252d37;color:#c7d0da;background:#171d25}.rows{padding:3px 0}.row{display:grid;grid-template-columns:minmax(0,1fr) auto;gap:16px;padding:8px 14px;border-bottom:1px solid #202731;align-items:baseline}.row:last-child{border-bottom:0}.label{color:#909ba8;font-size:.84rem}.metric{font-variant-numeric:tabular-nums;text-align:right;font-weight:620;font-size:.88rem;overflow-wrap:anywhere}footer{margin-top:22px;color:#66717e;font-size:.78rem}@media(max-width:760px){main{padding:22px 14px 36px}.sections{grid-template-columns:1fr}.summary{grid-template-columns:1fr}header{display:block}.badge{display:inline-block;margin-top:10px}.row{grid-template-columns:1fr;gap:3px}.metric{text-align:left}}
</style>
</head>
<body>
<main>
<header><div><h1 id="title">Lumrik</h1><p class="subtitle" id="subtitle">Live processing status</p></div><div class="badge">live</div></header>
<div class="summary"><div><span>Stage</span><strong id="stage">startup</strong></div><div><span>Elapsed</span><strong id="elapsed">00:00:00</strong></div></div>
<div class="sections" id="sections"></div>
<footer id="updated">Waiting for status...</footer>
</main>
<script>
let runStartedMs=null,runFinishedMs=null;
function updateElapsed(){if(runStartedMs===null)return;const end=runFinishedMs??Date.now();const s=Math.floor(Math.max(0,end-runStartedMs)/1000);const sec=s%60,min=Math.floor(s/60)%60,h=Math.floor(s/3600);document.getElementById("elapsed").textContent=String(h).padStart(2,"0")+":"+String(min).padStart(2,"0")+":"+String(sec).padStart(2,"0")}
function renderSections(sections){const root=document.getElementById("sections");root.replaceChildren();for(const section of sections){const panel=document.createElement("section");panel.className="panel";const heading=document.createElement("h2");heading.textContent=section.title;panel.appendChild(heading);const rows=document.createElement("div");rows.className="rows";for(const item of section.metrics){const row=document.createElement("div");row.className="row";const label=document.createElement("div");label.className="label";label.textContent=item.label;const metric=document.createElement("div");metric.className="metric";metric.textContent=item.value;row.append(label,metric);rows.appendChild(row)}panel.appendChild(rows);root.appendChild(panel)}}
async function updateStatus(){try{const response=await fetch("/status",{cache:"no-store"});if(!response.ok)throw new Error("HTTP "+response.status);const s=await response.json();document.title=s.title;document.getElementById("title").textContent=s.title;document.getElementById("subtitle").textContent=s.subtitle;document.getElementById("stage").textContent=s.stage;runStartedMs=Number(s.started_unix_ms);runFinishedMs=s.finished_unix_ms===null?null:Number(s.finished_unix_ms);renderSections(s.sections);document.getElementById("updated").textContent="Updated "+new Date().toLocaleTimeString();updateElapsed()}catch(error){document.getElementById("updated").textContent="Status unavailable: "+error}}
setInterval(updateElapsed,1000);setInterval(updateStatus,2000);updateStatus();
</script>
</body>
</html>"##;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn snapshot_json_preserves_sections_and_escapes_values() {
        let snapshot = ServerSnapshot {
            title: "Lumrik \"test\"".into(),
            subtitle: "status".into(),
            started_unix_ms: 42,
            finished_unix_ms: None,
            stage: "running".into(),
            public_url: Some("http://host:8787".into()),
            sections: vec![StatusSection::new(
                "Work",
                vec![StatusMetric::new("Cells", "209")],
            )],
        };
        let json = snapshot_json(&snapshot);
        assert!(json.contains("Lumrik \\\"test\\\""));
        assert!(json.contains("\"label\":\"Cells\""));
        assert!(json.contains("\"value\":\"209\""));
        assert!(json.contains("\"finished_unix_ms\":null"));
    }

    #[test]
    fn static_html_preserves_sections_and_escapes_values() {
        let snapshot = ServerSnapshot {
            title: "Lumrik <test>".into(),
            subtitle: "status & report".into(),
            started_unix_ms: 1_000,
            finished_unix_ms: Some(3_000),
            stage: "finished".into(),
            public_url: None,
            sections: vec![StatusSection::new(
                "Evidence",
                vec![StatusMetric::new("J→C", "42 < 100")],
            )],
        };
        let html = snapshot_html(&snapshot);
        assert!(html.contains("Lumrik &lt;test&gt;"));
        assert!(html.contains("status &amp; report"));
        assert!(html.contains("J→C"));
        assert!(html.contains("42 &lt; 100"));
        assert!(html.contains("00:00:02"));
    }

}
