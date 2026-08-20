use std::time::{Duration, Instant};

#[derive(Debug)]
pub struct RunProgress {
    started: Instant,

    last_report: Instant,
    last_reads: usize,

    reads_seen: usize,

    report_every: usize,
}

impl Default for RunProgress {
    fn default() -> Self {
        Self::new()
    }
}

impl RunProgress {
    pub fn new() -> Self {
        let now = Instant::now();

        Self {
            started: now,
            last_report: now,
            last_reads: 0,
            reads_seen: 0,

            report_every: 100_000,
        }
    }

    pub fn with_report_every(
        mut self,
        report_every: usize,
    ) -> Self {
        self.report_every = report_every.max(1);
        self
    }

    /// Print a major pipeline transition.
    pub fn stage(
        &self,
        message: impl AsRef<str>,
    ) {
        eprintln!(
            "[nelrune] {}",
            message.as_ref(),
        );
    }

    /// Record one successfully accepted input molecule/read pair.
    pub fn read_processed(
        &mut self,
    ) {
        self.reads_seen += 1;

        if self.reads_seen % self.report_every != 0 {
            return;
        }

        self.print_read_progress();
    }

    pub fn reads_seen(&self) -> usize {
        self.reads_seen
    }

    pub fn elapsed(&self) -> Duration {
        self.started.elapsed()
    }

    pub fn elapsed_seconds(&self) -> f64 {
        self.elapsed().as_secs_f64()
    }

    pub fn average_reads_per_second(
        &self,
    ) -> f64 {
        let elapsed =
            self.elapsed_seconds();

        if elapsed <= 0.0 {
            return 0.0;
        }

        self.reads_seen as f64 / elapsed
    }

    fn print_read_progress(
        &mut self,
    ) {
        let now = Instant::now();

        let elapsed =
            now
                .duration_since(
                    self.last_report,
                )
                .as_secs_f64();

        let reads =
            self.reads_seen
                - self.last_reads;

        let rate =
            if elapsed > 0.0 {
                reads as f64 / elapsed
            } else {
                0.0
            };

        eprintln!(
            "[nelrune] {:>12} reads | {:>10.0} reads/s",
            self.reads_seen,
            rate,
        );

        self.last_report = now;
        self.last_reads = self.reads_seen;
    }
}