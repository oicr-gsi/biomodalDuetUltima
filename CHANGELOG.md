# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0] - 2026-09-03

### Added

- `scheduler`, `slurmPartition`, `slurmAccount`, `processBeforeScript` and
  `imagesStagingDir`, so the workflow runs on slurm as well as sge. `scheduler`
  accepts `sge` (the default, unchanged behaviour), `slurm`, `slurm-gcp`, or `auto`
  to resolve from the submit command the cluster provides. It is decided before any
  work is done, and the submit command is checked for, so a setting the cluster
  cannot satisfy fails at once rather than at the first process submission.
- `slurm-gcp` differs from `slurm` in its defaults only: it takes `compute` as the
  partition and sends no accounting group, which is how a cloud slurm deployment is
  normally configured. `slurm` requires `slurmPartition`, because the pipeline's own
  slurm profile names a queue that exists only at the vendor's site. `auto`
  distinguishes the two from the machine's firmware identity, which is a local file
  read; resolving the cloud metadata hostname is deliberately not used, as a
  resolver that answers wildcards reports a cloud instance where there is none.
- The task writes the executor settings itself and appends them last, so no config
  file has to be placed on the cluster to change scheduler. For slurm it also
  answers `data_path`, `work_path` and `reference_path`, which the pipeline's slurm
  profile points at the vendor's own installation, and sets a `beforeScript` that
  gives `TMPDIR` a value: Grid Engine always sets it and slurm does not, and the
  container run options bind it. `processBeforeScript` is appended to that rather
  than replacing it, so a site that has to put its container runtime on PATH keeps
  the guard.
- `clusterOptions` and `time` are replaced on every `withName` selector that sets
  them, not just on the generic scope, since a generic assignment does not reach a
  selector. The selectors are found rather than listed, so a pipeline upgrade that
  adds one is covered.

### Fixed

- `containerAutoMounts`, defaulting to true. Nextflow derives its container bind
  points from the paths a task touches and collapses them into their common parents.
  Where the execution tree and the installed software sit under a single top-level
  directory, that directory is what gets bound -- and if the image has its own copy,
  the host's hides it and every tool inside disappears. Set it false there. The task
  now lists its bind points individually in `runOptions` regardless, so turning the
  derivation off loses nothing: the working directory, the reference bundle, and both
  the staged and resolved locations of every input are bound by name.


- The JVM's cgroup probe throws a `NullPointerException` where a cgroup v2 hierarchy
  is present but exposes no controller the JDK recognises, which is how a
  scheduler-created job cgroup can look. Nextflow reaches it while printing system
  information at startup, so the run dies before submitting any process.
  `-XX:-UseContainerSupport` is now passed to the head JVM, by way of both `NXF_OPTS`
  and `JAVA_TOOL_OPTIONS`, since a launcher that does not build its command line from
  the former still honours the latter. The flag governs only cgroup-derived defaults
  for heap size and processor count, and the heap is set explicitly, so it is inert
  where the probe works. The container run options clear `JAVA_TOOL_OPTIONS` so the
  pipeline's own containerised JVMs keep their cgroup-derived defaults.

### Changed

- The `penv`/`h_vmem` cluster options and the `qsub` shim are now written only for
  sge. Both exist to work around `h_vmem` being a per-slot limit and Nextflow's SGE
  executor emitting RSS directives alongside it; slurm's `--mem` is per job and has
  no equivalent, so neither is generated there.
- The container image staging directory is an input rather than a fixed path, so a
  site other than the one the module was built for can supply its own.

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
