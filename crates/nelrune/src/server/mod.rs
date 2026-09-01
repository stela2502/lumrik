use std::io::{Read, Write};
use std::net::{SocketAddr, TcpListener, TcpStream};
use std::sync::{Arc, RwLock};
use std::thread::{self, JoinHandle};

use anyhow::{Context, Result};

use crate::progress::RunStatus;

mod dashboard;

pub struct HealthServer {
    handle: JoinHandle<()>,
    addr: SocketAddr,
}

impl HealthServer {
    pub fn addr(&self) -> SocketAddr {
        self.addr
    }

    pub fn join(self) -> thread::Result<()> {
        self.handle.join()
    }
}

pub fn spawn_health_server(
    state: Arc<RwLock<RunStatus>>,
    addr: SocketAddr,
) -> Result<HealthServer> {
    let listener = TcpListener::bind(addr)
        .with_context(|| format!("binding Nelrune health server to {addr}"))?;

    let addr = listener
        .local_addr()
        .context("reading Nelrune health server address")?;

    let handle = thread::spawn(move || {
        for connection in listener.incoming() {
            match connection {
                Ok(stream) => {
                    if let Err(error) = handle_connection(stream, &state) {
                        eprintln!("[nelrune] health server request failed: {error}");
                    }
                }
                Err(error) => {
                    eprintln!("[nelrune] health server connection failed: {error}");
                }
            }
        }
    });

    Ok(HealthServer { handle, addr })
}

fn handle_connection(mut stream: TcpStream, state: &Arc<RwLock<RunStatus>>) -> Result<()> {
    let mut buffer = [0u8; 4096];
    let n = stream
        .read(&mut buffer)
        .context("reading health-server request")?;
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
            dashboard::HTML,
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

fn status_json(state: &Arc<RwLock<RunStatus>>) -> Result<String> {
    let state = state
        .read()
        .map_err(|_| anyhow::anyhow!("RunStatus lock poisoned"))?;

    let finished_unix_ms = state
        .finished_unix_ms
        .map(|value| value.to_string())
        .unwrap_or_else(|| "null".to_string());
    let input_file = json_optional_string(state.input_file.as_deref());
    let public_url = json_optional_string(state.public_url.as_deref());

    Ok(format!(
        concat!(
            "{{",
            "\"started_unix_ms\":{},",
            "\"finished_unix_ms\":{},",
            "\"stage\":\"{}\",",
            "\"reads_processed\":{},",
            "\"reads_per_second\":{:.3},",
            "\"no_cell_umi\":{},",
            "\"duplicates\":{},",
            "\"unique_genomic\":{},",
            "\"unique_feature\":{},",
            "\"duplicate_pct\":{:.3},",
            "\"unique_yield_pct\":{:.3},",
            "\"process_rss_mib\":{:.3},",
            "\"process_peak_rss_mib\":{:.3},",
            "\"system_available_mib\":{:.3},",
            "\"input_file\":{},",
            "\"public_url\":{}",
            "}}\n"
        ),
        state.started_unix_ms,
        finished_unix_ms,
        escape_json(&state.stage),
        state.reads_processed,
        state.reads_per_second,
        state.no_cell_umi,
        state.duplicates,
        state.unique_genomic,
        state.unique_feature,
        state.duplicate_pct,
        state.unique_yield_pct,
        state.process_rss_mib,
        state.process_peak_rss_mib,
        state.system_available_mib,
        input_file,
        public_url,
    ))
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
        .context("writing health-server response")?;
    stream.flush().context("flushing health-server response")?;
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
                use std::fmt::Write;
                let _ = write!(out, "\\u{:04x}", ch as u32);
            }
            ch => out.push(ch),
        }
    }
    out
}
