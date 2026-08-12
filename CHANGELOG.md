# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2026-08-12

### Added

- [GRD-1172](https://jira.oicr.on.ca/browse/GRD-1172)

- WDL wrapper around the biomodal duet Nextflow pipeline v1.7.0a1 in Ultima
  single-read CRAM-input mode, for Cromwell + UGE with Apptainer.
- `splitCram` task: partitions a single large Ultima CRAM into parts that the
  pipeline consumes as separate lanes, so `PRELUDE` and `BWA_MEM2` run per part in
  parallel and merge before deduplication. Used instead of the pipeline's own
  `split_reads_pre_resolution`, which round-trips through FASTQ and was measured at
  10.4 MB/s against samtools' 77.3 MB/s on the same data.
- Runtime patches appended to `nextflow_override.config` for OICR's UGE cluster:
  scheduling on `h_vmem` rather than `h_rss` to avoid cgroup memory kills, `penv`,
  project name, per-process cpu and memory clamps, and apptainer image caching on
  shared storage.
- Idempotent patches to the vendored Nextflow modules for two upstream defects: a
  malformed `@PG` header from `bwa-mem2` that crashed QUALIMAP, and a missing `$` in
  `samtools_cram_to_fastq.nf` that forced single-threaded CRAM to FASTQ conversion.
- Reads lacking the Ultima `tp`/`t0` flow tags are filtered out during splitting.
  `prelude` skips such reads, and its skip path leaves a truncated record in the
  resolved FASTQ that fails alignment.
- `bwaMem2Memory` is set independently of `maxMemory`: `bwa-mem2` loads a ~16GB index
  and is killed by the cgroup at the value that makes `PRELUDE` schedulable.
- Nextflow trace, report and timeline enabled through the config rather than the
  command line, which the vendor CLI does not forward.
