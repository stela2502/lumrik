use anyhow::{bail, Context, Result};

use std::fs::{File, OpenOptions};
use std::io::Write;
use std::path::Path;
use std::process::{
    Child,
    ChildStdin,
    Command,
    Stdio,
};
use std::sync::mpsc::{
    Receiver,
    TryRecvError,
};
use std::thread::{
    self,
    JoinHandle,
};

use bam_tide::fastq::FastqRecord;

use rust_htslib::bam::{
    self,
    Header,
    Read,
};

use tempfile::TempDir;

use crate::core::MapperProcessLike;

use crate::process::{
    sam_cluster_channel,
    MapperRecord,
    SamClusterBuffer,
    SamClusterReceiver,
    SamReadCluster,
};


const DEFAULT_CLUSTER_CHANNEL_SIZE: usize = 10_000;
const DEFAULT_CLUSTER_MAX_GAP: u64 = 100_000;
const DEFAULT_CLUSTER_FLUSH_EVERY: u64 = 1_000;


pub enum FastqInput {
    SingleStdin(
        ChildStdin,
    ),

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

    /*
     * ------------------------------------------------------------
     * Temporary resources
     * ------------------------------------------------------------
     *
     * These MUST remain alive for the complete lifetime of the
     * mapper process.
     *
     * In particular STAR receives the FIFO path through
     * --readFilesIn and uses the private work directory for files
     * such as _STARtmp.
     *
     * Dropping these TempDirs in spawn_*() would remove the paths
     * while the external mapper is still running.
     */
    _input_tmpdir: Option<TempDir>,
    _work_tmpdir: Option<TempDir>,

    /*
     * ------------------------------------------------------------
     * Header
     * ------------------------------------------------------------
     *
     * The stdout reader thread publishes one owned BAM/SAM header.
     *
     * header_loaded() may receive it non-blockingly and cache it
     * here. header() returns the cached object or waits for it.
     */
    header: Option<Header>,
    header_rx: Receiver<Header>,

    /*
     * ------------------------------------------------------------
     * Mapper output
     * ------------------------------------------------------------
     */
    clusters: SamClusterReceiver,
    stdout_thread: Option<JoinHandle<Result<()>>>,
}


impl MapperProcess {
    /*
     * ============================================================
     * Single stdin
     * ============================================================
     */

    pub fn spawn_single_stdin(
        binary: &Path,
        args: &[String],
        header: Option<Header>,
    ) -> Result<Self> {
        let mut child =
            Command::new(binary)
                .args(args)
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::inherit())
                .spawn()
                .with_context(|| {
                    format!(
                        "failed to spawn mapper `{}`",
                        binary.display(),
                    )
                })?;

        let stdout =
            child
                .stdout
                .take()
                .context(
                    "failed to capture mapper stdout",
                )?;

        let stdin =
            child
                .stdin
                .take()
                .context(
                    "failed to open mapper stdin",
                )?;

        let (
            clusters,
            header_rx,
            stdout_thread,
        ) =
            spawn_stdout_cluster_thread(
                stdout,
            );

        Ok(Self {
            child,

            input:
                Some(
                    FastqInput::SingleStdin(
                        stdin,
                    ),
                ),

            _input_tmpdir: None,
            _work_tmpdir: None,

            header,
            header_rx,

            clusters,
            stdout_thread:
                Some(stdout_thread),
        })
    }


    /*
     * ============================================================
     * Single FIFO
     * ============================================================
     */

    pub fn spawn_single_fifo(
        binary: &Path,
        base_args: &[String],
        header: Option<Header>,
    ) -> Result<Self> {
        /*
         * Keep the FIFO directory alive inside MapperProcess.
         */
        let input_tmpdir =
            tempfile::tempdir()
                .context(
                    "failed to create temporary FASTQ FIFO directory",
                )?;

        /*
         * Keep the mapper's private working directory alive too.
         *
         * This is important for STAR, which creates _STARtmp and
         * other runtime files relative to its working directory.
         */
        let work_tmpdir =
            tempfile::tempdir()
                .context(
                    "failed to create temporary mapper working directory",
                )?;

        let r1_path =
            input_tmpdir
                .path()
                .join("r1.fq");

        make_fifo(
            &r1_path,
        )?;

        /*
         * Caller already added the mapper-specific option preceding
         * the read path, e.g.
         *
         *     STAR ... --readFilesIn
         *
         * We append the actual named FIFO.
         */
        let mut args =
            base_args.to_vec();

        args.push(
            r1_path
                .to_string_lossy()
                .to_string(),
        );

        let mut child =
            Command::new(binary)
                .args(&args)
                .current_dir(
                    work_tmpdir.path(),
                )
                .stdin(Stdio::null())
                .stdout(Stdio::piped())
                .stderr(Stdio::inherit())
                .spawn()
                .with_context(|| {
                    format!(
                        "failed to spawn mapper `{}`",
                        binary.display(),
                    )
                })?;

        let stdout =
            child
                .stdout
                .take()
                .context(
                    "failed to capture mapper stdout",
                )?;

        let (
            clusters,
            header_rx,
            stdout_thread,
        ) =
            spawn_stdout_cluster_thread(
                stdout,
            );

        /*
         * Opening a FIFO for writing blocks until the mapper has
         * opened the other side for reading.
         *
         * For one FIFO this is fine and also gives us a useful
         * synchronization point: STAR has at least opened its input.
         */
        let r1 =
            OpenOptions::new()
                .write(true)
                .open(&r1_path)
                .with_context(|| {
                    format!(
                        "failed to open FASTQ FIFO `{}`",
                        r1_path.display(),
                    )
                })?;

        Ok(Self {
            child,

            input:
                Some(
                    FastqInput::SingleFifo {
                        r1,
                    },
                ),

            /*
             * CRITICAL:
             * move both TempDirs into MapperProcess.
             */
            _input_tmpdir:
                Some(input_tmpdir),

            _work_tmpdir:
                Some(work_tmpdir),

            header,
            header_rx,

            clusters,
            stdout_thread:
                Some(stdout_thread),
        })
    }


    /*
     * ============================================================
     * Paired FIFO
     * ============================================================
     */

    pub fn spawn_paired_fifo(
        binary: &Path,
        base_args: &[String],
        header: Option<Header>,
    ) -> Result<Self> {
        let input_tmpdir =
            tempfile::tempdir()
                .context(
                    "failed to create temporary paired FASTQ FIFO directory",
                )?;

        let work_tmpdir =
            tempfile::tempdir()
                .context(
                    "failed to create temporary mapper working directory",
                )?;

        let r1_path =
            input_tmpdir
                .path()
                .join("r1.fq");

        let r2_path =
            input_tmpdir
                .path()
                .join("r2.fq");

        make_fifo(
            &r1_path,
        )?;

        make_fifo(
            &r2_path,
        )?;

        let mut args =
            base_args.to_vec();

        args.push(
            r1_path
                .to_string_lossy()
                .to_string(),
        );

        args.push(
            r2_path
                .to_string_lossy()
                .to_string(),
        );

        let mut child =
            Command::new(binary)
                .args(&args)
                .current_dir(
                    work_tmpdir.path(),
                )
                .stdin(Stdio::null())
                .stdout(Stdio::piped())
                .stderr(Stdio::inherit())
                .spawn()
                .with_context(|| {
                    format!(
                        "failed to spawn paired mapper `{}`",
                        binary.display(),
                    )
                })?;

        let stdout =
            child
                .stdout
                .take()
                .context(
                    "failed to capture mapper stdout",
                )?;

        let (
            clusters,
            header_rx,
            stdout_thread,
        ) =
            spawn_stdout_cluster_thread(
                stdout,
            );

        /*
         * Opening paired FIFOs sequentially can deadlock:
         *
         * STAR may open R1 and then R2.
         * Parent opening R1 might wait while STAR waits elsewhere.
         *
         * Therefore open both writers concurrently.
         */
        let r1_open_path =
            r1_path.clone();

        let r2_open_path =
            r2_path.clone();

        let r1_open =
            thread::spawn(
                move || -> Result<File> {
                    OpenOptions::new()
                        .write(true)
                        .open(
                            &r1_open_path,
                        )
                        .with_context(|| {
                            format!(
                                "failed to open R1 FIFO `{}`",
                                r1_open_path.display(),
                            )
                        })
                },
            );

        let r2_open =
            thread::spawn(
                move || -> Result<File> {
                    OpenOptions::new()
                        .write(true)
                        .open(
                            &r2_open_path,
                        )
                        .with_context(|| {
                            format!(
                                "failed to open R2 FIFO `{}`",
                                r2_open_path.display(),
                            )
                        })
                },
            );

        let r1 =
            r1_open
                .join()
                .map_err(|_| {
                    anyhow::anyhow!(
                        "R1 FIFO opener thread panicked"
                    )
                })??;

        let r2 =
            r2_open
                .join()
                .map_err(|_| {
                    anyhow::anyhow!(
                        "R2 FIFO opener thread panicked"
                    )
                })??;

        Ok(Self {
            child,

            input:
                Some(
                    FastqInput::PairedFifo {
                        r1,
                        r2,
                    },
                ),

            _input_tmpdir:
                Some(input_tmpdir),

            _work_tmpdir:
                Some(work_tmpdir),

            header,
            header_rx,

            clusters,
            stdout_thread:
                Some(stdout_thread),
        })
    }


    pub fn spawn_paired_fifo_with_arg_insertion(
        binary: &Path,
        base_args: &[String],
        header: Option<Header>,
    ) -> Result<Self> {
        Self::spawn_paired_fifo(
            binary,
            base_args,
            header,
        )
    }


    /*
     * ============================================================
     * Process/thread helpers
     * ============================================================
     */

    fn join_stdout_thread(
        &mut self,
    ) -> Result<()> {
        let Some(handle) =
            self.stdout_thread.take()
        else {
            return Ok(());
        };

        handle
            .join()
            .map_err(|_| {
                anyhow::anyhow!(
                    "mapper stdout reader thread panicked"
                )
            })?
            .context(
                "mapper stdout reader thread failed",
            )
    }


    fn close_input(
        &mut self,
    ) -> Result<()> {
        let Some(mut input) =
            self.input.take()
        else {
            return Ok(());
        };

        /*
         * Flush before dropping the writers.
         */
        match &mut input {
            FastqInput::SingleStdin(
                stdin,
            ) => {
                stdin.flush()?;
            }

            FastqInput::SingleFifo {
                r1,
            } => {
                r1.flush()?;
            }

            FastqInput::PairedFifo {
                r1,
                r2,
            } => {
                r1.flush()?;
                r2.flush()?;
            }
        }

        /*
         * Dropping stdin/FIFO writers signals EOF to the mapper.
         */
        drop(input);

        Ok(())
    }
}


impl MapperProcessLike for MapperProcess {
    /*
     * ============================================================
     * FASTQ input
     * ============================================================
     */

    fn write_fastq(
        &mut self,
        r1: &FastqRecord,
        r2: Option<&FastqRecord>,
    ) -> Result<()> {
        let input =
            self
                .input
                .as_mut()
                .context(
                    "cannot write FASTQ: mapper input has already been closed",
                )?;

        match input {
            FastqInput::SingleStdin(
                stdin,
            ) => {
                write!(
                    stdin,
                    "{r1}"
                )?;

                stdin.flush()?;
            }

            FastqInput::SingleFifo {
                r1: w1,
            } => {
                write!(
                    w1,
                    "{r1}"
                )?;

                w1.flush()?;
            }

            FastqInput::PairedFifo {
                r1: w1,
                r2: w2,
            } => {
                let Some(r2) =
                    r2
                else {
                    bail!(
                        "paired mapper requires R2 but got None"
                    );
                };

                write!(
                    w1,
                    "{r1}"
                )?;

                write!(
                    w2,
                    "{r2}"
                )?;

                w1.flush()?;
                w2.flush()?;
            }
        }

        Ok(())
    }


    /*
     * ============================================================
     * Child liveness
     * ============================================================
     */

    fn is_running(
        &mut self,
    ) -> Result<bool> {
        match self
            .child
            .try_wait()
            .context(
                "failed to query mapper process state",
            )?
        {
            None => {
                Ok(true)
            }

            Some(status) => {
                eprintln!(
                    "[sc-mapper] mapper exited with status {status}"
                );

                Ok(false)
            }
        }
    }


    /*
     * ============================================================
     * Streaming results
     * ============================================================
     */

    /// Return any completed mapper result currently available.
    ///
    /// This intentionally does NOT try to match the result to the
    /// FASTQ read that was most recently submitted. External
    /// mappers may complete reads out of order.
    fn next_cluster(
        &mut self,
    ) -> Result<Option<SamReadCluster>> {
        match self
            .clusters
            .try_recv()
        {
            Ok(cluster) => {
                Ok(
                    Some(cluster)
                )
            }

            Err(
                TryRecvError::Empty,
            ) => {
                Ok(None)
            }

            Err(
                TryRecvError::Disconnected,
            ) => {
                /*
                 * If stdout processing has terminated, propagate any
                 * parser/thread error rather than silently hiding it.
                 */
                self.join_stdout_thread()?;

                Ok(None)
            }
        }
    }


    /*
     * ============================================================
     * Shutdown
     * ============================================================
     */

    /// Close mapper input, drain all outstanding mapping results,
    /// then wait for the mapper and stdout reader to terminate.
    ///
    /// Draining happens before waiting on the child because the
    /// stdout reader publishes through a bounded channel. Waiting
    /// first could deadlock if that channel fills while the mapper
    /// is still producing output.
    fn finish(
        mut self: Box<Self>,
    ) -> Result<Vec<SamReadCluster>> {
        /*
         * Signal EOF first.
         */
        self.close_input()?;

        let mut remaining =
            Vec::new();

        /*
         * Once input is closed, the mapper will eventually finish
         * and the stdout reader will drop the final sender.
         *
         * Blocking recv() is intentional here.
         */
        while let Ok(cluster) =
            self.clusters.recv()
        {
            remaining.push(
                cluster,
            );
        }

        /*
         * Propagate any SAM/BAM parsing or buffering error.
         */
        self.join_stdout_thread()?;

        /*
         * Finally wait for the actual mapper process.
         */
        let status =
            self
                .child
                .wait()
                .context(
                    "failed while waiting for mapper process",
                )?;

        if !status.success() {
            bail!(
                "mapper process failed with exit status {status}"
            );
        }

        Ok(remaining)
    }


    /*
     * ============================================================
     * Header
     * ============================================================
     */

    /// Blocking header retrieval.
    ///
    /// If header_loaded() has already consumed the header from the
    /// channel, it will be returned immediately from self.header.
    fn header(
        &mut self,
    ) -> Result<&Header> {
        if self.header.is_none() {
            let header =
                self
                    .header_rx
                    .recv()
                    .context(
                        "mapper stdout closed before providing a BAM header",
                    )?;

            self.header =
                Some(header);
        }

        Ok(
            self
                .header
                .as_ref()
                .expect(
                    "header was initialized above"
                ),
        )
    }


    /// Non-blocking check for mapper header availability.
    ///
    /// If a header has arrived, this consumes it from the channel
    /// and caches it in self.header so a later header() call remains
    /// lossless and immediate.
    fn header_loaded(
        &mut self,
    ) -> bool {
        if self.header.is_some() {
            return true;
        }

        match self
            .header_rx
            .try_recv()
        {
            Ok(header) => {
                self.header =
                    Some(header);

                true
            }

            Err(
                TryRecvError::Empty,
            ) => {
                false
            }

            Err(
                TryRecvError::Disconnected,
            ) => {
                false
            }
        }
    }
}


/*
 * ================================================================
 * Binary availability
 * ================================================================
 */

pub fn check_binary(
    path: &Path,
    display_name: &str,
) -> Result<()> {
    match Command::new(path)
        .arg("--version")
        .output()
    {
        Ok(output) => {
            if output.status.success() {
                Ok(())
            } else {
                bail!(
                    "{} binary `{}` executed but returned non-zero status {}",
                    display_name,
                    path.display(),
                    output.status,
                );
            }
        }

        Err(err)
            if err.kind()
                == std::io::ErrorKind::NotFound =>
        {
            bail!(
                "{} binary not found: {}",
                display_name,
                path.display(),
            );
        }

        Err(err) => {
            bail!(
                "failed to execute {} binary {}: {}",
                display_name,
                path.display(),
                err,
            );
        }
    }
}


/*
 * ================================================================
 * Mapper stdout reader
 * ================================================================
 */

fn spawn_stdout_cluster_thread(
    stdout:
        std::process::ChildStdout,
) -> (
    SamClusterReceiver,
    Receiver<Header>,
    JoinHandle<Result<()>>,
) {
    let (
        tx,
        rx,
    ) =
        sam_cluster_channel(
            DEFAULT_CLUSTER_CHANNEL_SIZE,
        );

    let (
        header_tx,
        header_rx,
    ) =
        std::sync::mpsc::channel();

    let handle =
        thread::spawn(
            move || -> Result<()> {
                /*
                 * This call is intentionally allowed to block until
                 * the external mapper starts producing SAM/BAM.
                 *
                 * It happens in this dedicated stdout thread, never
                 * in the main Nelrune thread.
                 */
                let mut reader =
                    bam_reader_from_child_stdout(
                        stdout,
                    )
                    .context(
                        "failed to create BAM reader for mapper stdout",
                    )?;

                /*
                 * HeaderView belongs to Reader. Convert it into an
                 * owned Header that MapperProcess can retain.
                 */
                let header =
                    Header::from_template(
                        reader.header(),
                    );

                header_tx
                    .send(header)
                    .context(
                        "failed to publish mapper BAM header",
                    )?;

                let mut buffer =
                    SamClusterBuffer::new(
                        tx,
                        DEFAULT_CLUSTER_MAX_GAP,
                        DEFAULT_CLUSTER_FLUSH_EVERY,
                    );

                for record in
                    reader.records()
                {
                    let record =
                        record.context(
                            "failed to read mapper SAM/BAM record",
                        )?;

                    buffer.push(
                        MapperRecord::new(
                            record,
                        ),
                    )?;
                }

                /*
                 * Flush clusters remaining when mapper stdout reaches
                 * EOF.
                 */
                buffer.finish()?;

                Ok(())
            },
        );

    (
        rx,
        header_rx,
        handle,
    )
}


/*
 * ================================================================
 * rust-htslib reader from child stdout
 * ================================================================
 */

#[cfg(unix)]
fn bam_reader_from_child_stdout(
    stdout:
        std::process::ChildStdout,
) -> Result<bam::Reader> {
    use std::os::fd::AsRawFd;

    let fd =
        stdout.as_raw_fd();

    let path =
        format!(
            "/proc/self/fd/{fd}"
        );

    /*
     * Keep ChildStdout alive until rust-htslib has opened/duplicated
     * the descriptor.
     */
    let reader =
        bam::Reader::from_path(
            &path,
        )
        .with_context(|| {
            format!(
                "failed to create BAM reader from mapper stdout fd {fd}"
            )
        })?;

    /*
     * Reader now owns its own open descriptor.
     */
    drop(stdout);

    Ok(reader)
}


/*
 * ================================================================
 * FIFO helper
 * ================================================================
 */

#[cfg(unix)]
fn make_fifo(
    path: &Path,
) -> Result<()> {
    use nix::sys::stat::Mode;
    use nix::unistd::mkfifo;

    mkfifo(
        path,
        Mode::S_IRUSR
            | Mode::S_IWUSR,
    )
    .with_context(|| {
        format!(
            "failed to create FIFO `{}`",
            path.display(),
        )
    })
}


#[cfg(not(unix))]
compile_error!(
    "sc-mapper currently requires a Unix platform."
);