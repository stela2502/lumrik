//! Lumrik's CLI argument types, built on upstream `clap`.
//!
//! `ReadableFile` checks filesystem accessibility while clap parses arguments,
//! before the application starts loading reference indices or processing data.
//! This checks that the file can be opened, not that its contents are valid.

use std::{fmt, fs::{self, File}, ops::Deref, path::{Path, PathBuf}, str::FromStr};

pub use clap;
pub use clap::{Args, CommandFactory, FromArgMatches, Parser, Subcommand, ValueEnum};

/// An input file successfully opened for reading during CLI parsing.
///
/// Use this type for required and optional input files:
///
/// ```
/// use lumrik_clap::{Parser, ReadableFile};
///
/// #[derive(Parser)]
/// struct Cli {
///     #[arg(long)]
///     index: ReadableFile,
///     #[arg(long)]
///     vcf: Option<ReadableFile>,
/// }
/// ```
///
/// Clap includes the argument name in errors from `FromStr`.
/// The path is not held open; consumers should still handle subsequent I/O errors.
#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub struct ReadableFile(PathBuf);

impl ReadableFile {
    pub fn as_path(&self) -> &Path {
        &self.0
    }

    pub fn into_path_buf(self) -> PathBuf {
        self.0
    }
}

impl FromStr for ReadableFile {
    type Err = String;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        let path = PathBuf::from(value);
        let file = File::open(&path)
            .map_err(|error| format!("cannot read input file '{}': {error}", path.display()))?;
        let metadata = file.metadata()
            .map_err(|error| format!("cannot inspect input file '{}': {error}", path.display()))?;
        if !metadata.is_file() {
            return Err(format!("input path '{}' is not a regular file", path.display()));
        }
        Ok(Self(path))
    }
}

impl AsRef<Path> for ReadableFile {
    fn as_ref(&self) -> &Path {
        self.as_path()
    }
}

impl Deref for ReadableFile {
    type Target = Path;

    fn deref(&self) -> &Self::Target {
        self.as_path()
    }
}

impl fmt::Display for ReadableFile {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.display().fmt(f)
    }
}

/// A required CLI path that may identify an existing readable file or a new file.
///
/// Unlike `Option<ReadableFile>`, which permits *omitting the argument* but
/// rejects a supplied path that does not exist, `ExistingOrNewFile` requires
/// the argument while permitting its target to be absent.
///
/// - Existing regular file: must be readable.
/// - Missing target: accepted without creating it.
/// - Existing directory or other non-file: rejected.
/// - Any error other than `NotFound`: rejected.
///
/// Parsing never creates or modifies the target. The application remains
/// responsible for opening or creating it, including handling I/O races.
#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub struct ExistingOrNewFile(PathBuf);

impl FromStr for ExistingOrNewFile {
    type Err = String;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        if value.is_empty() {
            return Err("existing-or-new file path must not be empty".into());
        }
        let path = PathBuf::from(value);
        match File::open(&path) {
            Ok(file) => {
                let metadata = file.metadata().map_err(|error| {
                    format!("cannot inspect existing file '{}': {error}", path.display())
                })?;
                if !metadata.is_file() {
                    return Err(format!("path '{}' is not a regular file", path.display()));
                }
            }
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
                // A missing target is valid; never create it during parsing.
            }
            Err(error) => {
                return Err(format!("cannot read existing file '{}': {error}", path.display()));
            }
        }
        Ok(Self(path))
    }
}

/// An existing input directory that can be opened for listing.
#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub struct ReadableDir(PathBuf);

impl FromStr for ReadableDir {
    type Err = String;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        let path = PathBuf::from(value);
        fs::read_dir(&path)
            .map_err(|error| format!("cannot read input directory '{}': {error}", path.display()))?;
        Ok(Self(path))
    }
}

/// A destination filename. Parsing does not touch the filesystem.
/// Call `prepare()` after successful argument parsing to create its parent directories.
#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub struct OutputFile(PathBuf);

/// A filename prefix used to construct several output files.
///
/// Parsing validates the prefix without creating anything. `prepare()` creates
/// only its parent directory, never the prefix itself. A trailing path separator,
/// an empty filename, or a path already pointing to a directory is rejected.
/// `with_suffix(".results.tsv")` appends rather than replacing an extension.
#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub struct OutputPrefix(PathBuf);

impl OutputPrefix {
    pub fn prepare(&self) -> std::io::Result<()> {
        if let Some(parent) = self.0.parent().filter(|p| !p.as_os_str().is_empty()) {
            fs::create_dir_all(parent)?;
        }
        Ok(())
    }

    pub fn with_suffix(&self, suffix: &str) -> PathBuf {
        let mut name = self.0.as_os_str().to_os_string();
        name.push(suffix);
        PathBuf::from(name)
    }
}

impl FromStr for OutputPrefix {
    type Err = String;
    fn from_str(value: &str) -> Result<Self, Self::Err> {
        if value.is_empty() || value.ends_with('/') || value.ends_with('\\') {
            return Err("output prefix must be a nonempty filename without a trailing separator".into());
        }
        let path = PathBuf::from(value);
        if path.file_name().is_none() || path.is_dir() {
            return Err(format!("output prefix '{}' must be a filename, not a directory", path.display()));
        }
        Ok(Self(path))
    }
}

/// A destination directory. Parsing does not touch the filesystem.
/// Call `prepare()` after successful argument parsing to create it recursively.
#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub struct OutputDir(PathBuf);

impl OutputFile {
    pub fn prepare(&self) -> std::io::Result<()> {
        if let Some(parent) = self.0.parent().filter(|parent| !parent.as_os_str().is_empty()) {
            fs::create_dir_all(parent)?;
        }
        Ok(())
    }
}

impl OutputDir {
    pub fn prepare(&self) -> std::io::Result<()> {
        fs::create_dir_all(&self.0)
    }
}

/// Storage formats accepted by a particular single-cell CLI argument.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum SingleCellFormat { Mtx, H5, H5ad }

impl SingleCellFormat {
    fn label(self) -> &'static str {
        match self { Self::Mtx => "mtx", Self::H5 => "h5", Self::H5ad => "h5ad" }
    }
}

fn classify_single_cell(path: &Path, input: bool) -> Result<SingleCellFormat, String> {
    if input && path.is_dir() {
        fs::read_dir(path).map_err(|e| format!("cannot read single-cell directory '{}': {e}", path.display()))?;
        return Ok(SingleCellFormat::Mtx);
    }
    let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("").to_ascii_lowercase();
    let format = match ext.as_str() {
        "h5ad" => SingleCellFormat::H5ad,
        "h5" => SingleCellFormat::H5,
        _ if !input && !path.exists() || !input && path.is_dir() => SingleCellFormat::Mtx,
        _ => return Err(format!("cannot determine single-cell format of '{}': expected an MTX directory, .h5 or .h5ad file", path.display())),
    };
    if input {
        File::open(path).map_err(|e| format!("cannot read single-cell file '{}': {e}", path.display()))?
            .metadata().map_err(|e| format!("cannot inspect '{}': {e}", path.display()))?
            .is_file().then_some(()).ok_or_else(|| format!("'{}' is not a regular file", path.display()))?;
    } else if path.exists() && path.is_dir() && format != SingleCellFormat::Mtx {
        return Err(format!("'{}' is a directory, not a {} file", path.display(), format.label()));
    }
    Ok(format)
}

fn check_allowed(format: SingleCellFormat, allowed: &[SingleCellFormat]) -> Result<(), String> {
    if allowed.contains(&format) { return Ok(()); }
    Err(format!("detected {}, but this argument supports only {}", format.label(),
        allowed.iter().map(|f| f.label()).collect::<Vec<_>>().join(", ")))
}

/// An existing single-cell file or MTX/MEX directory. Parsing does not load matrices.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ReadableSingleCellData { path: PathBuf, format: SingleCellFormat }
impl ReadableSingleCellData {
    pub fn format(&self) -> SingleCellFormat { self.format }
    pub fn as_path(&self) -> &Path { &self.path }
    pub fn into_path_buf(self) -> PathBuf { self.path }
}
impl AsRef<Path> for ReadableSingleCellData { fn as_ref(&self) -> &Path { &self.path } }
impl Deref for ReadableSingleCellData { type Target = Path; fn deref(&self) -> &Path { &self.path } }
impl FromStr for ReadableSingleCellData {
    type Err = String;
    fn from_str(value: &str) -> Result<Self, Self::Err> {
        let path = PathBuf::from(value);
        let format = classify_single_cell(&path, true)?;
        Ok(Self { path, format })
    }
}

/// An output single-cell file or directory; parsing never creates it.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct OutputSingleCellData { path: PathBuf, format: SingleCellFormat }
impl OutputSingleCellData {
    pub fn format(&self) -> SingleCellFormat { self.format }
    pub fn as_path(&self) -> &Path { &self.path }
    pub fn into_path_buf(self) -> PathBuf { self.path }
    pub fn prepare(&self) -> std::io::Result<()> {
        match self.format {
            SingleCellFormat::Mtx => fs::create_dir_all(&self.path),
            _ => match self.path.parent().filter(|p| !p.as_os_str().is_empty()) {
                Some(parent) => fs::create_dir_all(parent), None => Ok(()),
            }
        }
    }
}
impl AsRef<Path> for OutputSingleCellData { fn as_ref(&self) -> &Path { &self.path } }
impl Deref for OutputSingleCellData { type Target = Path; fn deref(&self) -> &Path { &self.path } }
impl FromStr for OutputSingleCellData {
    type Err = String;
    fn from_str(value: &str) -> Result<Self, Self::Err> {
        if value.is_empty() { return Err("single-cell output path must not be empty".into()); }
        let path = PathBuf::from(value);
        let format = classify_single_cell(&path, false)?;
        Ok(Self { path, format })
    }
}

/// Clap parser restricting an input argument to the storage formats its reader supports.
pub fn single_cell_input<const N: usize>(allowed: [SingleCellFormat; N])
    -> impl Fn(&str) -> Result<ReadableSingleCellData, String> + Clone {
    move |value: &str| {
        let data: ReadableSingleCellData = value.parse()?;
        check_allowed(data.format(), &allowed)?;
        Ok(data)
    }
}

/// Clap parser restricting the output formats supported by a writer.
pub fn single_cell_output<const N: usize>(allowed: [SingleCellFormat; N])
    -> impl Fn(&str) -> Result<OutputSingleCellData, String> + Clone {
    move |value: &str| {
        let data: OutputSingleCellData = value.parse()?;
        check_allowed(data.format(), &allowed)?;
        Ok(data)
    }
}

// Keep the wrapper types ergonomic for libraries expecting paths, while ensuring
// callers explicitly choose whether to prepare an output destination.
macro_rules! path_accessors {
    ($type:ident) => {
        impl $type {
            pub fn as_path(&self) -> &Path { &self.0 }
            pub fn into_path_buf(self) -> PathBuf { self.0 }
        }
        impl AsRef<Path> for $type {
            fn as_ref(&self) -> &Path { self.as_path() }
        }
        impl Deref for $type {
            type Target = Path;
            fn deref(&self) -> &Path { self.as_path() }
        }
        impl fmt::Display for $type {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                self.0.display().fmt(f)
            }
        }
    };
}

path_accessors!(ExistingOrNewFile);
path_accessors!(ReadableDir);
path_accessors!(OutputFile);
path_accessors!(OutputDir);
path_accessors!(OutputPrefix);

impl FromStr for OutputFile {
    type Err = String;
    fn from_str(value: &str) -> Result<Self, Self::Err> {
        if value.is_empty() { return Err("output file path must not be empty".into()); }
        let path = PathBuf::from(value);
        if path.file_name().is_none() {
            return Err(format!("output file path '{}' must include a filename", path.display()));
        }
        Ok(Self(path))
    }
}

impl FromStr for OutputDir {
    type Err = String;
    fn from_str(value: &str) -> Result<Self, Self::Err> {
        if value.is_empty() { return Err("output directory path must not be empty".into()); }
        Ok(Self(PathBuf::from(value)))
    }
}


/// Clap value parsers for APIs that still store filesystem paths as `PathBuf`.
/// These enforce the same validation as the typed wrappers without changing
/// downstream library APIs or internal path ownership.
pub fn parse_readable_file(value: &str) -> Result<PathBuf, String> {
    value.parse::<ReadableFile>().map(ReadableFile::into_path_buf)
}
pub fn parse_readable_dir(value: &str) -> Result<PathBuf, String> {
    value.parse::<ReadableDir>().map(ReadableDir::into_path_buf)
}
pub fn parse_output_file(value: &str) -> Result<PathBuf, String> {
    value.parse::<OutputFile>().map(OutputFile::into_path_buf)
}
pub fn parse_output_dir(value: &str) -> Result<PathBuf, String> {
    value.parse::<OutputDir>().map(OutputDir::into_path_buf)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Debug, Parser)]
    struct Cli {
        #[arg(long)]
        index: ReadableFile,
        #[arg(long)]
        vcf: Option<ReadableFile>,
    }

    #[derive(Debug, Parser)]
    struct PathsCli {
        #[arg(long)]
        input_dir: ReadableDir,
        #[arg(long)]
        output_file: OutputFile,
        #[arg(long)]
        output_dir: OutputDir,
    }

    #[test]
    fn missing_directory_reports_option_and_path() {
        let missing = std::env::temp_dir().join(format!(
            "lumrik_clap_missing_dir_{}_{}", std::process::id(), line!()
        ));
        let err = ReadableDir::from_str(missing.to_str().unwrap()).unwrap_err();
        assert!(err.contains(&missing.display().to_string()), "{err}");
    }

    #[test]
    fn output_paths_do_not_create_directories_during_parsing() {
        let base = std::env::temp_dir().join(format!(
            "lumrik_clap_output_{}_{}", std::process::id(), line!()
        ));
        let output_file = base.join("nested").join("results.tsv");
        let output_dir = base.join("figures");
        assert!(!base.exists());
        let args = PathsCli::try_parse_from([
            "test", "--input-dir", std::env::temp_dir().to_str().unwrap(),
            "--output-file", output_file.to_str().unwrap(),
            "--output-dir", output_dir.to_str().unwrap(),
        ]).unwrap();
        assert!(!base.exists());
        args.output_file.prepare().unwrap();
        args.output_dir.prepare().unwrap();
        assert!(output_file.parent().unwrap().is_dir());
        assert!(!output_file.exists());
        assert!(output_dir.is_dir());
        fs::remove_dir_all(&base).unwrap();
    }

    #[test]
    fn output_file_never_truncates_existing_file() {
        let path = std::env::temp_dir().join(format!(
            "lumrik_clap_existing_{}_{}.txt", std::process::id(), line!()
        ));
        fs::write(&path, b"keep me").unwrap();
        let output = OutputFile::from_str(path.to_str().unwrap()).unwrap();
        output.prepare().unwrap();
        assert_eq!(fs::read(&path).unwrap(), b"keep me");
        fs::remove_file(&path).unwrap();
    }

    #[test]
    fn missing_file_reports_option_and_path() {
        let missing = std::env::temp_dir().join(format!(
            "lumrik_clap_missing_{}_{}", std::process::id(), line!()
        ));
        let err = Cli::try_parse_from(["test", "--index", missing.to_str().unwrap()])
            .unwrap_err().to_string();
        assert!(err.contains("--index"), "{err}");
        assert!(err.contains(&missing.display().to_string()), "{err}");
    }

    #[test]
    fn directory_is_not_an_input_file() {
        let directory = std::env::temp_dir();
        let err = Cli::try_parse_from(["test", "--index", directory.to_str().unwrap()])
            .unwrap_err().to_string();
        assert!(err.contains("--index"), "{err}");
        assert!(err.contains("not a regular file"), "{err}");
    }

    #[test]
    fn existing_file_is_accepted_and_exposes_path() {
        let path = PathBuf::from(concat!(env!("CARGO_MANIFEST_DIR"), "/src/lib.rs"));
        let parsed = Cli::try_parse_from(["test", "--index", path.to_str().unwrap()])
            .unwrap();
        assert_eq!(parsed.index.as_path(), path);
        assert!(parsed.vcf.is_none());
    }

    #[test]
    fn optional_input_is_validated_when_present() {
        let path = PathBuf::from(concat!(env!("CARGO_MANIFEST_DIR"), "/src/lib.rs"));
        let missing = std::env::temp_dir().join(format!(
            "lumrik_clap_optional_missing_{}_{}", std::process::id(), line!()
        ));
        let err = Cli::try_parse_from([
            "test", "--index", path.to_str().unwrap(), "--vcf", missing.to_str().unwrap()
        ]).unwrap_err().to_string();
        assert!(err.contains("--vcf"), "{err}");
        assert!(err.contains(&missing.display().to_string()), "{err}");
    }

    #[derive(Debug, Parser)]
    struct ExistingOrNewCli {
        #[arg(long)]
        store: ExistingOrNewFile,
    }

    #[test]
    fn existing_or_new_file_accepts_missing_path_without_creating_it() {
        let dir = std::env::temp_dir().join(format!(
            "lumrik_clap_new_store_{}_{}", std::process::id(), line!()
        ));
        let path = dir.join("store.bin");
        assert!(!path.exists());
        let args = ExistingOrNewCli::try_parse_from([
            "test", "--store", path.to_str().unwrap()
        ]).unwrap();
        assert_eq!(args.store.as_path(), path);
        assert!(!dir.exists());
    }

    #[test]
    fn existing_or_new_file_accepts_readable_file_but_rejects_directory() {
        let path = PathBuf::from(concat!(env!("CARGO_MANIFEST_DIR"), "/Cargo.toml"));
        ExistingOrNewCli::try_parse_from([
            "test", "--store", path.to_str().unwrap()
        ]).unwrap();
        let err = ExistingOrNewCli::try_parse_from([
            "test", "--store", std::env::temp_dir().to_str().unwrap()
        ]).unwrap_err().to_string();
        assert!(err.contains("--store"), "{err}");
        assert!(err.contains("not a regular file"), "{err}");
    }

    #[test]
    fn existing_or_new_file_is_not_an_optional_argument() {
        assert!(ExistingOrNewCli::try_parse_from(["test"]).is_err());
    }

    #[derive(Debug, Parser)]
    struct PrefixCli {
        #[arg(long)]
        out: OutputPrefix,
    }

    #[test]
    fn output_prefix_appends_suffix_and_prepares_only_parent() {
        let base = std::env::temp_dir().join(format!(
            "lumrik_clap_prefix_{}_{}", std::process::id(), line!()
        ));
        let prefix = base.join("nested").join("TP53_");
        let args = PrefixCli::try_parse_from([
            "test", "--out", prefix.to_str().unwrap()
        ]).unwrap();
        assert!(!base.exists());
        assert_eq!(args.out.with_suffix("genome.svg"), base.join("nested/TP53_genome.svg"));
        args.out.prepare().unwrap();
        assert!(base.join("nested").is_dir());
        assert!(!prefix.exists());
        fs::remove_dir_all(base).unwrap();
    }

    #[test]
    fn output_prefix_rejects_directories_and_trailing_separators() {
        assert!(OutputPrefix::from_str("").is_err());
        assert!(OutputPrefix::from_str("results/").is_err());
        assert!(OutputPrefix::from_str(std::env::temp_dir().to_str().unwrap()).is_err());
    }

}

#[cfg(test)]
mod single_cell_tests {
    use super::*;

    #[test]
    fn detects_input_formats_and_rejects_disallowed() {
        let base = std::env::temp_dir().join(format!("lumrik_sc_{}_{}", std::process::id(), line!()));
        fs::create_dir_all(&base).unwrap();
        let h5ad = base.join("data.h5ad");
        fs::write(&h5ad, b"fixture").unwrap();
        let data = single_cell_input([SingleCellFormat::H5ad])(h5ad.to_str().unwrap()).unwrap();
        assert_eq!(data.format(), SingleCellFormat::H5ad);
        assert!(single_cell_input([SingleCellFormat::Mtx])(h5ad.to_str().unwrap()).unwrap_err().contains("supports only mtx"));
        assert_eq!(single_cell_input([SingleCellFormat::Mtx])(base.to_str().unwrap()).unwrap().format(), SingleCellFormat::Mtx);
        fs::remove_dir_all(base).unwrap();
    }

    #[test]
    fn output_parser_does_not_create_files_and_prepares_correctly() {
        let base = std::env::temp_dir().join(format!("lumrik_sc_out_{}_{}", std::process::id(), line!()));
        let file = base.join("nested").join("data.h5ad");
        let output = single_cell_output([SingleCellFormat::H5ad])(file.to_str().unwrap()).unwrap();
        assert!(!base.exists());
        output.prepare().unwrap();
        assert!(file.parent().unwrap().is_dir());
        assert!(!file.exists());
        assert!(single_cell_output([SingleCellFormat::Mtx])(file.to_str().unwrap()).is_err());
        fs::remove_dir_all(base).unwrap();
    }
}


// Automatic command history for Lumrik executables. Linux startup hook:
// no explicit call is required from individual binaries.
#[cfg(target_os = "linux")]
#[used]
#[unsafe(link_section = ".init_array")]
static LUMRIK_COMMAND_HISTORY_INIT: extern "C" fn() = lumrik_command_history_init;

#[cfg(target_os = "linux")]
extern "C" fn lumrik_command_history_init() {
    // A logging failure must never prevent an analysis from starting.
    let _ = std::panic::catch_unwind(|| command_history::record());
}

#[cfg(target_os = "linux")]
mod command_history {
    use std::{
        env,
        ffi::OsStr,
        fs::{self, OpenOptions},
        io::{self, Write},
        os::{fd::AsRawFd, unix::ffi::OsStrExt},
        path::{Path, PathBuf},
        process::Command,
    };

    unsafe extern "C" {
        fn flock(fd: i32, operation: i32) -> i32;
    }
    const LOCK_EX: i32 = 2;
    const LOCK_UN: i32 = 8;
    const MAX_BYTES: u64 = 10 * 1024 * 1024;
    const ROTATIONS: usize = 5;

    fn shell_quote(value: &OsStr) -> String {
        let bytes = value.as_bytes();
        if !bytes.is_empty()
            && bytes.iter().all(|b| b.is_ascii_alphanumeric() || b"@%_+=:,./-".contains(b))
        {
            return String::from_utf8_lossy(bytes).into_owned();
        }
        // Shell-quote every argument without interpreting its contents.
        // Most command arguments are UTF-8; non-UTF8 bytes are represented
        // using Bash ANSI-C octal quoting to preserve them exactly.
        if let Ok(value) = std::str::from_utf8(bytes) {
            return format!("'{}'", value.replace('\'', "'\\''"));
        }
        let mut escaped = String::from("$'");
        for &byte in bytes {
            match byte {
                b'\\' => escaped.push_str("\\\\"),
                b'\'' => escaped.push_str("\\'"),
                0x20..=0x7e => escaped.push(byte as char),
                _ => escaped.push_str(&format!("\\{:03o}", byte)),
            }
        }
        escaped.push('\'');
        escaped
    }

    fn history_path() -> Option<PathBuf> {
        if env::var_os("LUMRIK_COMMAND_LOG_DISABLE").as_deref() == Some(OsStr::new("1")) {
            return None;
        }
        if let Some(path) = env::var_os("LUMRIK_COMMAND_LOG") {
            return Some(PathBuf::from(path));
        }
        let home = env::var_os("HOME")?;
        Some(PathBuf::from(home).join(".lumrik/commands.db"))
    }

    fn rotate(path: &Path) -> io::Result<()> {
        for i in (1..=ROTATIONS).rev() {
            let source = if i == 1 { path.to_path_buf() }
                else { PathBuf::from(format!("{}.{}", path.display(), i - 1)) };
            let dest = PathBuf::from(format!("{}.{}", path.display(), i));
            if source.exists() { fs::rename(source, dest)?; }
        }
        Ok(())
    }

    pub(super) fn record() -> io::Result<()> {
        let Some(path) = history_path() else { return Ok(()); };
        let cwd = env::current_dir()?;
        let args: Vec<_> = env::args_os().collect();
        if args.is_empty() { return Ok(()); }
        // No history for help/version invocations.
        if args.iter().any(|a| a == "--help" || a == "-h" || a == "--version" || a == "-V") {
            return Ok(());
        }
        let date = Command::new("date").arg("+%Y-%m-%d %H:%M:%S %z").output()
            .ok().filter(|output| output.status.success())
            .map(|output| String::from_utf8_lossy(&output.stdout).trim().to_owned())
            .unwrap_or_else(|| format!("{:?}", std::time::SystemTime::now()));
        let entry = format!(
            "# {date}\ncd {}\n{}\n\n",
            shell_quote(cwd.as_os_str()),
            args.iter().map(|a| shell_quote(a)).collect::<Vec<_>>().join(" ")
        );
        if let Some(parent) = path.parent() { fs::create_dir_all(parent)?; }
        // Separate stable lock file prevents concurrent writers from interleaving
        // entries even while the history is being rotated.
        let lock_path = PathBuf::from(format!("{}.lock", path.display()));
        let lock = OpenOptions::new().create(true).write(true).open(lock_path)?;
        if unsafe { flock(lock.as_raw_fd(), LOCK_EX) } != 0 {
            return Err(io::Error::last_os_error());
        }
        let result = (|| {
            if fs::metadata(&path).is_ok_and(|m| m.len() + entry.len() as u64 > MAX_BYTES) {
                rotate(&path)?;
            }
            let mut file = OpenOptions::new().create(true).append(true).open(&path)?;
            file.write_all(entry.as_bytes())
        })();
        unsafe { flock(lock.as_raw_fd(), LOCK_UN); }
        result
    }
}
