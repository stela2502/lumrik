//! Shared AI transport and audit logging for Lumrik.
//! The caller owns prompts, tool semantics and interpretation.
use anyhow::{Context, Result, bail};
use serde_json::{Value, json};
use std::fs::{self, OpenOptions};
use std::io::Write;
use std::path::Path;
use std::process::{Command, Stdio};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

fn json_line(path: &Path, value: &Value) -> Result<()> {
    let mut file = OpenOptions::new().create(true).append(true).open(path)?;
    writeln!(file, "{}", serde_json::to_string(value)?)?;
    Ok(())
}

pub fn chat(out: &Path, url: &str, payload: &Value) -> Result<(Value, String)> {
    // Audit the complete request BEFORE contacting Ollama, including failed/time-out calls.
    let audit = out.join("ai_audit");
    fs::create_dir_all(&audit).context("creating AI audit directory")?;
    let openai = url == "openai";
    let endpoint = if openai {
        "https://api.openai.com/v1/chat/completions".to_string()
    } else {
        format!("{}/api/chat", url.trim_end_matches('/'))
    };
    let api_payload = if openai {
        let mut p = payload.clone();
        let object = p.as_object_mut().context("AI payload must be an object")?;
        object.remove("stream");
        if object.remove("format").is_some() {
            object.insert("response_format".into(), json!({"type":"json_object"}));
        }
        // Ollama uses tool_name for tool responses; OpenAI requires a call ID.
        // The interactive loop provides that ID when using the OpenAI backend.
        p
    } else {
        payload.clone()
    };
    let timestamp_ms = SystemTime::now().duration_since(UNIX_EPOCH)?.as_millis();
    let request = json!({"timestamp_unix_ms":timestamp_ms,"endpoint":endpoint,
        "payload":api_payload});
    fs::write(
        audit.join("ai_last_request.json"),
        serde_json::to_vec_pretty(&request)?,
    )?;
    json_line(
        &audit.join("ai_history.jsonl"),
        &json!({"event":"request", "request":request}),
    )?;
    let started = Instant::now();
    let result = (|| -> Result<(Value, String)> {
        let mut command = Command::new("curl");
        command.args([
            "--fail-with-body",
            "--silent",
            "--show-error",
            "--max-time",
            "1200",
            "-H",
            "Content-Type: application/json",
            "-X",
            "POST",
        ]);
        if openai {
            let key = std::env::var("OPENAI_API_KEY")
                .context("OPENAI_API_KEY must be set for --ai-provider openai")?;
            anyhow::ensure!(!key.trim().is_empty(), "OPENAI_API_KEY is empty");
            command
                .arg("-H")
                .arg(format!("Authorization: Bearer {key}"));
        }
        let mut child = command
            .args(["--data-binary", "@-", &endpoint])
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .context("starting AI HTTP request")?;
        child
            .stdin
            .take()
            .context("curl stdin")?
            .write_all(serde_json::to_string(&api_payload)?.as_bytes())?;
        let output = child.wait_with_output().context("waiting for Ollama")?;
        let raw = String::from_utf8_lossy(&output.stdout).into_owned();
        if !output.status.success() {
            bail!(
                "Ollama HTTP request failed: {} (body: {})",
                String::from_utf8_lossy(&output.stderr),
                raw
            );
        }
        let parsed: Value = serde_json::from_str(&raw).context("invalid AI response JSON")?;
        if openai {
            let choice = parsed
                .pointer("/choices/0/message")
                .context("OpenAI response missing choices[0].message")?;
            let mut message = choice.clone();
            // Normalize OpenAI tool calls to the existing Ollama interface.
            if let Some(calls) = message.get_mut("tool_calls").and_then(Value::as_array_mut) {
                for call in calls.iter_mut() {
                    if let Some(arguments) =
                        call.pointer("/function/arguments").and_then(Value::as_str)
                    {
                        let decoded: Value = serde_json::from_str(arguments)
                            .context("invalid OpenAI tool arguments")?;
                        call["function"]["arguments"] = decoded;
                    }
                }
            }
            Ok((json!({"message":message}), raw))
        } else {
            Ok((parsed, raw))
        }
    })();
    let elapsed_ms = started.elapsed().as_millis();
    match &result {
        Ok((response, raw)) => {
            let record = json!({"timestamp_unix_ms":timestamp_ms,"elapsed_ms":elapsed_ms,
                "status":"ok","response":response,"raw_response":raw});
            fs::write(
                audit.join("ai_last_response.json"),
                serde_json::to_vec_pretty(&record)?,
            )?;
            json_line(
                &audit.join("ai_history.jsonl"),
                &json!({"event":"response","record":record}),
            )?;
        }
        Err(error) => {
            let record = json!({"timestamp_unix_ms":timestamp_ms,"elapsed_ms":elapsed_ms,
                "status":"error","error":format!("{error:#}")});
            fs::write(
                audit.join("ai_last_error.json"),
                serde_json::to_vec_pretty(&record)?,
            )?;
            json_line(
                &audit.join("ai_history.jsonl"),
                &json!({"event":"error","record":record}),
            )?;
        }
    }
    result
}
