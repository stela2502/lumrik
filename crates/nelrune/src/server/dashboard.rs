pub const HTML: &str = r##"<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">

<title>Nelrune</title>

<style>
body {
    font-family:
        system-ui,
        sans-serif;

    max-width: 900px;
    margin: 40px auto;
    padding: 0 24px;

    background: #111;
    color: #eee;
}

h1 {
    margin-bottom: 4px;
}

.subtitle {
    color: #999;
    margin-bottom: 32px;
}

.grid {
    display: grid;

    grid-template-columns:
        repeat(
            auto-fit,
            minmax(220px, 1fr)
        );

    gap: 16px;
}

.card {
    background: #1c1c1c;

    border:
        1px solid #333;

    border-radius: 10px;

    padding: 20px;
}

.label {
    color: #999;

    font-size: 0.9rem;

    margin-bottom: 8px;
}

.value {
    font-size: 1.5rem;
    font-weight: 600;

    overflow-wrap: anywhere;
}

.file {
    font-size: 1rem;
    font-family: monospace;
}

.footer {
    margin-top: 32px;
    color: #666;
    font-size: 0.85rem;
}
</style>
</head>

<body>

<h1>Nelrune</h1>

<div class="subtitle">
Live processing status
</div>

<div class="card" style="margin-bottom:16px">
    <div class="label">External URL</div>
    <div class="value file" id="public_url">-</div>
</div>

<div class="grid">

<div class="card">
    <div class="label">
        Stage
    </div>

    <div
        class="value"
        id="stage"
    >
        startup
    </div>
</div>

<div class="card">
    <div class="label">Elapsed</div>
    <div class="value" id="elapsed">00:00:00</div>
</div>

<div class="card">
    <div class="label">
        Reads processed
    </div>

    <div
        class="value"
        id="reads"
    >
        0
    </div>
</div>

<div class="card">
    <div class="label">
        Reads / second
    </div>

    <div
        class="value"
        id="rate"
    >
        0
    </div>
</div>


<div class="card">
    <div class="label">Cell / UMI not detected</div>
    <div class="value" id="no_cell_umi">0</div>
</div>

<div class="card">
    <div class="label">Duplicates</div>
    <div class="value" id="duplicates">0</div>
</div>

<div class="card">
    <div class="label">Unique genomic</div>
    <div class="value" id="unique_genomic">0</div>
</div>

<div class="card">
    <div class="label">Unique feature</div>
    <div class="value" id="unique_feature">0</div>
</div>

<div class="card">
    <div class="label">Duplicate fraction</div>
    <div class="value" id="duplicate_pct">0%</div>
</div>

<div class="card">
    <div class="label">Unique molecule yield</div>
    <div class="value" id="unique_yield_pct">0%</div>
</div>

<div class="card">
    <div class="label">Nelrune RSS / peak</div>
    <div class="value" id="process_memory">0 / 0 MiB</div>
</div>

<div class="card">
    <div class="label">System memory available</div>
    <div class="value" id="system_available">0 MiB</div>
</div>

<div class="card">
    <div class="label">
        Current FASTQ
    </div>

    <div
        class="value file"
        id="file"
    >
        -
    </div>
</div>

</div>

<div
    class="footer"
    id="updated"
>
    Waiting for status...
</div>

<script>
let runStartedMs = null;
let runFinishedMs = null;

function updateElapsed() {
    if (runStartedMs === null) return;
    const endMs = runFinishedMs ?? Date.now();
    const totalSeconds = Math.floor(Math.max(0, endMs - runStartedMs) / 1000);
    const seconds = totalSeconds % 60;
    const minutes = Math.floor(totalSeconds / 60) % 60;
    const hours = Math.floor(totalSeconds / 3600);
    document.getElementById("elapsed").textContent =
        String(hours).padStart(2, "0") + ":" +
        String(minutes).padStart(2, "0") + ":" +
        String(seconds).padStart(2, "0");
}

async function updateStatus() {
    try {
        const response =
            await fetch(
                "/status",
                {
                    cache: "no-store"
                }
            );

        if (!response.ok) {
            throw new Error(
                "HTTP " + response.status
            );
        }

        const status =
            await response.json();

        if (runStartedMs === null) {
            runStartedMs = Number(status.started_unix_ms);
        }
        runFinishedMs = status.finished_unix_ms === null
            ? null
            : Number(status.finished_unix_ms);
        updateElapsed();
        document.querySelector("#elapsed").previousElementSibling.textContent =
            runFinishedMs === null ? "Elapsed" : "Total elapsed";

        document
            .getElementById("stage")
            .textContent =
            status.stage;

        document
            .getElementById("reads")
            .textContent =
            status
                .reads_processed
                .toLocaleString();

        document
            .getElementById("rate")
            .textContent =
            Math.round(
                status.reads_per_second
            ).toLocaleString();


        document
            .getElementById("no_cell_umi")
            .textContent =
            status.no_cell_umi.toLocaleString();

        document
            .getElementById("duplicates")
            .textContent =
            status.duplicates.toLocaleString();

        document
            .getElementById("unique_genomic")
            .textContent =
            status.unique_genomic.toLocaleString();

        document
            .getElementById("unique_feature")
            .textContent =
            status.unique_feature.toLocaleString();

        document.getElementById("duplicate_pct").textContent =
            status.duplicate_pct.toFixed(2) + "%";

        document.getElementById("unique_yield_pct").textContent =
            status.unique_yield_pct.toFixed(2) + "%";

        document.getElementById("process_memory").textContent =
            status.process_rss_mib.toFixed(0) + " / " +
            status.process_peak_rss_mib.toFixed(0) + " MiB";

        document.getElementById("system_available").textContent =
            status.system_available_mib.toFixed(0) + " MiB";

        document
            .getElementById("file")
            .textContent =
            status.input_file ?? "-";

        document
            .getElementById("public_url")
            .textContent =
            status.public_url ?? "-";

        document
            .getElementById("updated")
            .textContent =
            "Updated " +
            new Date()
                .toLocaleTimeString();
    }
    catch (error) {
        document
            .getElementById("updated")
            .textContent =
            "Connection lost";
    }
}

updateElapsed();
updateStatus();
setInterval(updateElapsed, 1000);
setInterval(updateStatus, 1000);
</script>

</body>
</html>
"##;
