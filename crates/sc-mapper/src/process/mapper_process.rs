use anyhow::{bail, Context, Result};
use std::fs::{File, OpenOptions};
use std::io::Write;
use std::os::fd::AsRawFd;
use std::path::{Path, PathBuf};
use std::process::{Child, ChildStdin, Command, Stdio};
use std::sync::mpsc::{Receiver, TryRecvError};
use std::thread::{self, JoinHandle};
use std::io::stdin;

use bam_tide::fastq::FastqRecord;
use rust_htslib::bam::{self, Header, Read};
use tempfile::TempDir;

use crate::core::MapperProcessLike;
use crate::process::{
    sam_cluster_channel, MapperRecord, SamClusterBuffer, SamClusterReceiver, SamReadCluster,
};

const DEFAULT_CLUSTER_CHANNEL_SIZE: usize = 10_000;
const DEFAULT_CLUSTER_MAX_GAP: u64 = 100_000;
const DEFAULT_CLUSTER_FLUSH_EVERY: u64 = 1_000;

pub enum FastqInput {
    SingleStdin(ChildStdin),
    SingleFifo {
        r1: File,
    },
    PairedFifo {
        r1: File,
        r2: File,
    },
}

pub struct MapperProcess {
    child: Child,
    input: Option<FastqInput>,

    // The BAM/SAM header is owned by MapperProcess once the
    // stdout reader thread has published it.
    header: Option<Header>,
    header_rx: Receiver<Header>,

    // Completed read clusters produced by the stdout reader thread.
    clusters: SamClusterReceiver,
    stdout_thread: Option<JoinHandle<Result<()>>>,
}

impl MapperProcess {
    pub fn spawn_single_stdin(binary: &Path, args: &[String]) -> Result<Self> {

        let mut child = Command::new(binary)
            .args(args)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::inherit())
            .spawn()
            .with_context(|| {
                format!("failed to spawn mapper `{}`", binary.display())
            })?;

        let stdout = child
            .stdout
            .take()
            .context("failed to capture mapper stdout")?;

        let stdin = child
            .stdin
            .take()
            .context("failed to open mapper stdin")?;

        let (clusters, header_rx, stdout_thread) =
            spawn_stdout_cluster_thread(stdout);

        Ok(Self {
            child,
            input: Some(FastqInput::SingleStdin(stdin)),
            header: None,
            header_rx,
            clusters,
            stdout_thread: Some(stdout_thread),
        })
    }

    pub fn spawn_single_fifo(
        binary: &Path,
        base_args: &[String],
    ) -> Result<Self> {
        let input_tmpdir = tempfile::tempdir()
            .context("failed to create temporary FASTQ FIFO directory")?;

        let r1_path = input_tmpdir.path().join("r1.fq");

        make_fifo(&r1_path)?;

        // The caller has already added the mapper-specific argument
        // preceding the input path, e.g. STAR's "--readFilesIn".
        let mut args = base_args.to_vec();
        args.push(r1_path.to_string_lossy().to_string());

        let mut child = Command::new(binary)
            .args(&args)
            .stdin(Stdio::null())
            .stdout(Stdio::piped())
            .stderr(Stdio::inherit())
            .spawn()
            .with_context(|| {
                format!(
                    "failed to spawn mapper `{}`",
                    binary.display()
                )
            })?;

        let stdout = child
            .stdout
            .take()
            .context("failed to capture mapper stdout")?;

        let (clusters, header_rx, stdout_thread) =
            spawn_stdout_cluster_thread(stdout);

        // This blocks until the mapper opens the FIFO for reading.
        // That is fine for a single FIFO: there is no R1/R2 ordering
        // problem like there is for paired input.
        let r1 = OpenOptions::new()
            .write(true)
            .open(&r1_path)
            .with_context(|| {
                format!(
                    "failed to open FASTQ FIFO `{}`",
                    r1_path.display()
                )
            })?;

        Ok(Self {
            child,

            input: Some(FastqInput::SingleFifo {
                r1,
            }),

            header: None,
            header_rx,
            clusters,
            stdout_thread: Some(stdout_thread),
        })
    }

    pub fn spawn_paired_fifo(
        binary: &Path,
        base_args: &[String],
    ) -> Result<Self> {
        let input_tmpdir = tempfile::tempdir()
            .context("failed to create temporary paired FASTQ FIFO directory")?;

        let r1_path = input_tmpdir.path().join("r1.fq");
        let r2_path = input_tmpdir.path().join("r2.fq");

        make_fifo(&r1_path)?;
        make_fifo(&r2_path)?;

        let mut args = base_args.to_vec();
        args.push(r1_path.to_string_lossy().to_string());
        args.push(r2_path.to_string_lossy().to_string());

        let mut child = Command::new(binary)
            .args(&args)
            .stdin(Stdio::null())
            .stdout(Stdio::piped())
            .stderr(Stdio::inherit())
            .spawn()
            .with_context(|| {
                format!(
                    "failed to spawn paired mapper `{}`",
                    binary.display()
                )
            })?;

        let stdout = child
            .stdout
            .take()
            .context("failed to capture mapper stdout")?;

        let (clusters, header_rx, stdout_thread) =
            spawn_stdout_cluster_thread(stdout);

        let r1_open_path = r1_path.clone();
        let r2_open_path = r2_path.clone();

        let r1_open = thread::spawn(move || -> Result<File> {
            OpenOptions::new()
                .write(true)
                .open(&r1_open_path)
                .with_context(|| {
                    format!(
                        "failed to open R1 FIFO `{}`",
                        r1_open_path.display()
                    )
                })
        });

        let r2_open = thread::spawn(move || -> Result<File> {
            OpenOptions::new()
                .write(true)
                .open(&r2_open_path)
                .with_context(|| {
                    format!(
                        "failed to open R2 FIFO `{}`",
                        r2_open_path.display()
                    )
                })
        });

        let r1 = r1_open
            .join()
            .map_err(|_| anyhow::anyhow!("R1 FIFO opener thread panicked"))??;

        let r2 = r2_open
            .join()
            .map_err(|_| anyhow::anyhow!("R2 FIFO opener thread panicked"))??;

        Ok(Self {
            child,

            input: Some(FastqInput::PairedFifo {
                r1,
                r2,
            }),

            header: None,
            header_rx,

            clusters,
            stdout_thread: Some(stdout_thread),
        })
    }
    pub fn spawn_paired_fifo_with_arg_insertion(
        binary: &Path,
        base_args: &[String],
    ) -> Result<Self> {
        Self::spawn_paired_fifo(binary, base_args)
    }

    fn join_stdout_thread(&mut self) -> Result<()> {
        let Some(handle) = self.stdout_thread.take() else {
            return Ok(());
        };

        handle
            .join()
            .map_err(|_| anyhow::anyhow!("mapper stdout reader thread panicked"))?
            .context("mapper stdout reader thread failed")
    }

    fn close_input(&mut self) -> Result<()> {
        let Some(mut input) = self.input.take() else {
            return Ok(());
        };

        match &mut input {
            FastqInput::SingleStdin(stdin) => {
                stdin.flush()?;
            }

            FastqInput::SingleFifo { r1 } => {
                r1.flush()?
            }

            FastqInput::PairedFifo { r1, r2, .. } => {
                r1.flush()?;
                r2.flush()?;
            }
        }

        // Dropping stdin/FIFO writers signals EOF to the mapper.
        drop(input);

        Ok(())
    }
}

impl MapperProcessLike for MapperProcess {
    fn write_fastq(
        &mut self,
        r1: &FastqRecord,
        r2: Option<&FastqRecord>,
    ) -> Result<()> {
        let input = self
            .input
            .as_mut()
            .context("cannot write FASTQ: mapper input has already been closed")?;

        match input {
            FastqInput::SingleStdin(stdin) => {
                write!(stdin, "{r1}")?;
                stdin.flush()?;
            }

            FastqInput::SingleFifo { r1: w1 } => {
                write!(w1, "{r1}")?;
                w1.flush()?;
            }

            FastqInput::PairedFifo {
                r1: w1,
                r2: w2,
                ..
            } => {
                let Some(r2) = r2 else {
                    bail!("paired mapper requires R2 but got None");
                };

                write!(w1, "{r1}")?;
                write!(w2, "{r2}")?;

                w1.flush()?;
                w2.flush()?;
            }
        }

        Ok(())
    }

    /// Return any completed mapper result currently available.
    ///
    /// This intentionally does NOT try to match the result to the FASTQ read
    /// that was most recently submitted. External mappers may complete reads
    /// out of order.
    fn next_cluster(&mut self) -> Result<Option<SamReadCluster>> {
        match self.clusters.try_recv() {
            Ok(cluster) => Ok(Some(cluster)),

            Err(TryRecvError::Empty) => Ok(None),

            Err(TryRecvError::Disconnected) => {
                self.join_stdout_thread()?;
                Ok(None)
            }
        }
    }

    /// Close mapper input, drain all outstanding mapping results, then wait
    /// for the mapper and stdout reader to terminate.
    ///
    /// Draining happens before waiting on the child because the stdout reader
    /// publishes through a bounded channel. Waiting first could deadlock if
    /// that channel fills while the mapper is still producing output.
    fn finish(mut self: Box<Self>) -> Result<Vec<SamReadCluster>> {
        self.close_input()?;

        let mut remaining = Vec::new();

        // Once input is closed, the mapper will eventually finish and its
        // stdout reader will drop the final sender. recv() is appropriate here:
        // shutdown is deliberately blocking until all remaining results arrive.
        while let Ok(cluster) = self.clusters.recv() {
            remaining.push(cluster);
        }

        // Propagate any parsing/buffering error from the stdout thread.
        self.join_stdout_thread()?;

        let status = self
            .child
            .wait()
            .context("failed while waiting for mapper process")?;

        if !status.success() {
            bail!("mapper process failed with exit status {status}");
        }

        Ok(remaining)
    }

    fn header(&mut self) -> Result<&Header> {
        if self.header.is_none() {
            let header = self
                .header_rx
                .recv()
                .context("mapper stdout closed before providing a BAM header")?;

            self.header = Some(header);
        }

        Ok(self
            .header
            .as_ref()
            .expect("header was initialized above"))
    }
}

struct MapperStdoutPipe {
    tmpdir: TempDir,
    reader_path: PathBuf,
}

impl MapperStdoutPipe {
    fn new() -> Result<Self> {
        let tmpdir =
            tempfile::tempdir().context("failed to create temporary mapper stdout FIFO directory")?;

        let reader_path = tmpdir.path().join("mapper.stdout.sam_or_bam");

        make_fifo(&reader_path)?;

        Ok(Self {
            tmpdir,
            reader_path,
        })
    }

    fn open_writer(&self) -> Result<File> {
        OpenOptions::new()
            .write(true)
            .open(&self.reader_path)
            .with_context(|| {
                format!(
                    "failed to open mapper stdout FIFO writer `{}`",
                    self.reader_path.display()
                )
            })
    }
}

pub fn check_binary(path: &Path, display_name: &str) -> Result<()> {
    match std::process::Command::new(path).arg("--version").output() {
        Ok(_) => Ok(()),

        Err(err) if err.kind() == std::io::ErrorKind::NotFound => {
            bail!("{} binary not found: {}", display_name, path.display())
        }

        Err(err) => {
            bail!(
                "failed to execute {} binary {}: {}",
                display_name,
                path.display(),
                err
            )
        }
    }
}

fn spawn_stdout_cluster_thread(
    stdout: std::process::ChildStdout,
) -> (
    SamClusterReceiver,
    Receiver<Header>,
    JoinHandle<Result<()>>,
) {
    let (tx, rx) =
        sam_cluster_channel(DEFAULT_CLUSTER_CHANNEL_SIZE);

    let (header_tx, header_rx) =
        std::sync::mpsc::channel();

    let handle = thread::spawn(move || -> Result<()> {
        let mut reader =
            bam_reader_from_child_stdout(stdout)
                .context(
                    "failed to create BAM reader for mapper stdout",
                )?;

        // HeaderView belongs to the reader.
        // Convert it into an owned Header that we can retain and
        // later use for one or many BAM writers.
        let header =
            Header::from_template(reader.header());

        header_tx
            .send(header)
            .context(
                "failed to publish mapper BAM header",
            )?;

        let mut buffer = SamClusterBuffer::new(
            tx,
            DEFAULT_CLUSTER_MAX_GAP,
            DEFAULT_CLUSTER_FLUSH_EVERY,
        );

        for record in reader.records() {
            let record = record.context(
                "failed to read mapper SAM/BAM record",
            )?;

            buffer.push(
                MapperRecord::new(record),
            )?;
        }

        // Flush all clusters still buffered when mapper stdout
        // reaches EOF.
        buffer.finish()?;

        Ok(())
    });

    (rx, header_rx, handle)
}

#[cfg(unix)]
fn bam_reader_from_child_stdout(
    stdout: std::process::ChildStdout,
) -> Result<bam::Reader> {
    use std::os::fd::AsRawFd;

    let fd = stdout.as_raw_fd();
    let path = format!("/proc/self/fd/{fd}");

    // Important:
    // keep `stdout` alive until Reader has duplicated/opened the fd.
    let reader = bam::Reader::from_path(&path)
        .with_context(|| {
            format!(
                "failed to create BAM reader from mapper stdout fd {fd}"
            )
        })?;

    drop(stdout);

    Ok(reader)
}


#[cfg(unix)]
fn make_fifo(path: &Path) -> Result<()> {
    use nix::sys::stat::Mode;
    use nix::unistd::mkfifo;

    mkfifo(path, Mode::S_IRUSR | Mode::S_IWUSR)
        .with_context(|| format!("failed to create FIFO `{}`", path.display()))
}

#[cfg(not(unix))]
compile_error!("sc-mapper currently requires a Unix platform.");