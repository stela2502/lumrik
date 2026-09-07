# Nelrune integration

`sc-vdj` consumes the retained mapper BAM after normal Nelrune mapping. The BAM is authoritative for cell identity and exact genomic V/D/J/C overlap; the VDJ index is built from the same GTF and genome as the mapping run.

```text
Nelrune mapper BAM
      + matching .vdjidx
              |
              v
          VdjRunner
              |
              v
        CellEvidenceVdj
              |
              v
        Recombination
              |
              v
      AIRR-compatible TSV
```

Expression quantification is intentionally not part of receptor reconstruction. It can be joined downstream by cell ID without changing the V(D)J call.
