// crates/sc-mapper/src/process/sam_cluster_buffer.rs

use anyhow::{Context, Result};
use std::collections::HashMap;
use std::sync::mpsc::{sync_channel, Receiver, SyncSender};

use crate::process::{MapperRecord, SamReadCluster};

pub type SamClusterSender = SyncSender<SamReadCluster>;
pub type SamClusterReceiver = Receiver<SamReadCluster>;

pub fn sam_cluster_channel(buffer_size: usize) -> (SamClusterSender, SamClusterReceiver) {
    sync_channel(buffer_size)
}

#[derive(Debug)]
pub struct SamClusterBuffer {
    active: HashMap<String, ActiveCluster>,
    tick: u64,
    max_gap: u64,
    flush_every: u64,
    tx: SamClusterSender,
}

#[derive(Debug)]
struct ActiveCluster {
    records: Vec<MapperRecord>,
    last_seen_tick: u64,
}

impl SamClusterBuffer {
    pub fn new(tx: SamClusterSender, max_gap: u64, flush_every: u64) -> Self {
        Self {
            active: HashMap::new(),
            tick: 0,
            max_gap,
            flush_every,
            tx,
        }
    }

    pub fn push(&mut self, rec: MapperRecord) -> Result<()> {
        self.tick += 1;

        let read_id = rec.clean_id();

        let cluster = self
            .active
            .entry(read_id)
            .or_insert_with(|| ActiveCluster {
                records: Vec::new(),
                last_seen_tick: self.tick,
            });

        cluster.records.push(rec);
        cluster.last_seen_tick = self.tick;

        if self.flush_every > 0 && self.tick % self.flush_every == 0 {
            self.flush_old()?;
        }

        Ok(())
    }

    fn flush_old(&mut self) -> Result<()> {
        let ready_ids: Vec<String> = self
            .active
            .iter()
            .filter_map(|(read_id, cluster)| {
                if self.tick.saturating_sub(cluster.last_seen_tick) > self.max_gap {
                    Some(read_id.clone())
                } else {
                    None
                }
            })
            .collect();

        for read_id in ready_ids {
            self.emit(read_id)
                .context("failed to emit old mapper cluster")?;
        }

        Ok(())
    }

    fn emit(&mut self, read_id: String) -> Result<()> {
        let Some(cluster) = self.active.remove(&read_id) else {
            return Ok(());
        };

        self.tx
            .send(SamReadCluster {
                read_id,
                records: cluster.records,
            })
            .context("failed to send mapper cluster to consumer")?;

        Ok(())
    }

    pub fn finish(mut self) -> Result<()> {
        let ready_ids: Vec<String> = self.active.keys().cloned().collect();

        for read_id in ready_ids {
            self.emit(read_id)
                .context("failed to send final mapper cluster to consumer")?;
        }

        Ok(())
    }
}