pub const HTML: &str = r#"<!doctype html>
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
            .getElementById("file")
            .textContent =
            status.input_file ?? "-";

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

updateStatus();

setInterval(
    updateStatus,
    1000
);
</script>

</body>
</html>
"#;