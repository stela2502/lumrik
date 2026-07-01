use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};
use directories::ProjectDirs;
use std::env;

pub struct ATACtoRNAMapper {
    mapping: HashMap<u32, u32>,
}

fn get_default_path() -> PathBuf {
    if let Ok(exe_path) = env::current_exe() {
        if let Some(parent) = exe_path.parent() {
            if parent.starts_with("/usr/local") {
                return PathBuf::from("/usr/local/share/bam_tide/whitelist_map.bin");
            }
        }
    }
    ProjectDirs::from("com", "stela2502", "bam_tide")
        .expect("Could not determine default data directory")
        .data_local_dir()
        .join("whitelist_map.bin")
}

impl ATACtoRNAMapper {
    /// Load the binary map from file
    pub fn from_binary_file<P: AsRef<Path>>(path: P) -> Result<Self, String> {
        let user_path = path.as_ref();



        // Try user-provided path
        let buffer = match File::open(user_path) {
            Ok(file) => {
                let mut reader = BufReader::new(file);
                let mut buffer = Vec::new();
                reader
                    .read_to_end(&mut buffer)
                    .map_err(|e| format!("Failed to read binary file: {}", e))?;
                buffer
            }
            Err(_) => {
                // Try default path
                let default_path = get_default_path();
                match File::open(&default_path) {
                    Ok(file) => {
                        let mut reader = BufReader::new(file);
                        let mut buffer = Vec::new();
                        reader
                            .read_to_end(&mut buffer)
                            .map_err(|e| format!("Failed to read binary file from default location: {}", e))?;
                        buffer
                    }
                    Err(_) => {
                        return Err(format!(
                            "Could not open file at {:?} and fallback default path {:?} does not exist.",
                            user_path, default_path
                        ));
                    }
                }
            }
        };

        if buffer.len() % 8 != 0 {
            return Err("Invalid binary format: not a multiple of 8 bytes".into());
        }

        let mut mapping = HashMap::with_capacity(buffer.len() / 8);
        for chunk in buffer.chunks_exact(8) {
            let key = u32::from_le_bytes(chunk[0..4].try_into().unwrap());
            let val = u32::from_le_bytes(chunk[4..8].try_into().unwrap());
            mapping.insert(key, val);
        }

        Ok(Self { mapping })
    }

    /// Translate an ATAC u32 ID to RNA u32 ID
    pub fn translate(&self, atac_id: u32) -> Option<u32> {
        self.mapping.get(&atac_id).copied()
    }
}