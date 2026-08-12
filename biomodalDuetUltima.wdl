version 1.0

workflow biomodalDuetUltima {
    input {
        Array[File] crams
        String sampleId
        String runName
        String outputFileNamePrefix
        String mode = "6bp"
        String additionalProfile = "deep_seq"
        String modules = "biomodal-duet-ultima/1.7.0a1 samtools/1.16.1"
        Int splitCramReads = -1
        Array[File]? craiIndexes
    }

    parameter_meta {
        crams: "Array of unaligned Ultima single-read CRAM files (one per lane) for a single sample"
        sampleId: "Sample identifier (used for naming input CRAMs and output files)"
        runName: "Sequencing run name / identifier (used in report file names)"
        outputFileNamePrefix: "Prefix for all output file names"
        mode: "Biomodal DUET mode: 6bp (duet evoC) or 5bp (duet +modC). Default: 6bp"
        additionalProfile: "Nextflow resource profile: deep_seq (<=500M reads), super_seq (>500M reads), or empty (<=50M reads). Default: deep_seq"
        modules: "Environment module providing the biomodal instance dir and its env vars (apptainer loads as a dependency)"
        craiIndexes: "Optional pre-built .crai, one per entry of crams, matched by basename. Only used when splitCram does not run -- when it does, its own indexes take precedence. Supply this when re-running against parts a previous splitCram already produced, otherwise step 6b re-indexes every part from scratch."
        splitCramReads: "If >0, split the single input CRAM into parts of this many reads and present them to the pipeline as separate lanes, so PRELUDE runs per part in parallel. This is the fast alternative to the pipeline's own split_reads_pre_resolution, which first converts CRAM->FASTQ and then splits with seqkit -- measured at 10.4 MB/s against samtools' 77.3 MB/s on the same data. Requires exactly one input CRAM. Default -1 (off)."
    }

    # Split mode: partition the one input CRAM into pseudo-lanes up front. The
    # pipeline parallelises PRELUDE and BWA_MEM2 per lane and merges before dedup
    # (resolve_align.nf), so this is the same parallelism its own splitter provides
    # without the CRAM->FASTQ round trip.
    if (splitCramReads > 0 && length(crams) == 1) {
        call splitCram {
            input:
                cram = crams[0],
                outputFileNamePrefix = outputFileNamePrefix,
                readsPerPart = splitCramReads,
                modules = modules
        }
    }

    call runDuet {
        input:
            crams = select_first([splitCram.parts, crams]),
            craiIndexes = if defined(splitCram.indexes) then splitCram.indexes else craiIndexes,
            sampleId = sampleId,
            runName  = runName,
            outputFileNamePrefix = outputFileNamePrefix,
            mode = mode,
            additionalProfile = additionalProfile,
            modules = modules
    }

    meta {
        author: "Gavin Peng"
        email: "gpeng@oicr.on.ca"
        description: "WDL wrapper for the Biomodal DUET evoC methylation-sequencing pipeline v1.7.0a1, running the Ultima single-read CRAM-input early-access mode on OICR's UGE/SGE cluster via Apptainer."
        dependencies: [
            {
                name: "biomodal-duet-ultima/1.7.0a1",
                url: "https://biomodal.com"
            }
        ]
        output_meta: {
            outputBam: {
                description: "Deduplicated, coordinate-sorted genome BAM file of aligned reads",
                vidarr_label: "outputBam"
            },
            outputBai:  {
                description: "BAM index (.bai) for random-access retrieval of the deduplicated BAM",
                vidarr_label: "outputBai"
            },
            hmc_cxreport:  {
                description: "Cytosine Report for 5-hydroxymethylcytosine (5hmC) at CpG sites. Tab-separated, one row per stranded CpG position; columns report chromosome, position, strand, methylated-read count, unmethylated-read count, and context (CG). Suitable for downstream epigenetic analysis tools (e.g. methylKit, DSS). Gzip-compressed.",
                vidarr_label: "hmc_cxreport"
            },
            hmc_cxreportIndex: {
                description: "Tabix index (.tbi) for the 5hmC Cytosine Report, enabling fast random-access queries by genomic region",
                vidarr_label: "hmc_cxreportIndex"
            },
            mc_cxreport:  {
                description: "Cytosine Report for 5-methylcytosine (5mC) at CpG sites. Same tab-separated, per-stranded-CpG format as the 5hmC report; columns give chromosome, position, strand, methylated-read count, unmethylated-read count, and context (CG). Suitable for downstream epigenetic analysis tools (e.g. methylKit, DSS). Gzip-compressed.",
                vidarr_label: "mc_cxreport"
            },
            mc_cxreportIndex: {
                description: "Tabix index (.tbi) for the 5mC Cytosine Report, enabling fast random-access queries by genomic region",
                vidarr_label: "mc_cxreportIndex"
            },
            modc_cxreport: {
                description: "Cytosine Report for total modified cytosine (5mC + 5hmC combined, modC) at CpG sites. Same tab-separated, per-stranded-CpG format; provides an aggregate modification signal across both marks. Gzip-compressed.",
                vidarr_label: "modc_cxreport"
            },
            modc_cxreportIndex: {
                description: "Tabix index (.tbi) for the modC Cytosine Report, enabling fast random-access queries by genomic region",
                vidarr_label: "modc_cxreportIndex"
            },
            vcf: {
                description: "Germline variant calls VCF (DeepVariant by default; optional, absent when no variants are called)",
                vidarr_label: "vcf"
            },
            vcfIndex: {
                description: "Tabix index (.tbi) for the germline VCF (optional)",
                vidarr_label: "vcfIndex"
            },
            summaryCsv: {
                description: "Run-level DUET summary metrics in CSV format",
                vidarr_label: "summaryCsv"
            },
            summaryHtml: {
                description: "Run-level DUET summary metrics as an interactive HTML report",
                vidarr_label: "summaryHtml"
            },
            summaryXlsx: {
                description: "Run-level DUET summary metrics in Excel format",
                vidarr_label: "summaryXlsx"
            },
            multiqcReport: {
                description: "MultiQC HTML report aggregating QC metrics across all pipeline steps",
                vidarr_label: "multiqcReport"
            },
            metricsDefinitions: {
                description: "CSV file defining and describing each metric reported in the summary outputs",
                vidarr_label: "metricsDefinitions"
            }
        }
    }

    output {
        File  outputBam = runDuet.outputBam
        File  outputBai = runDuet.outputBai
        File  hmc_cxreport = runDuet.hmc_cxreport
        File  hmc_cxreportIndex = runDuet.hmc_cxreportIndex
        File  mc_cxreport = runDuet.mc_cxreport
        File  mc_cxreportIndex = runDuet.mc_cxreportIndex
        File  modc_cxreport = runDuet.modc_cxreport
        File  modc_cxreportIndex = runDuet.modc_cxreportIndex
        File? vcf = runDuet.vcf
        File? vcfIndex = runDuet.vcfIndex
        File  summaryCsv = runDuet.summaryCsv
        File  summaryHtml = runDuet.summaryHtml
        File  summaryXlsx = runDuet.summaryXlsx
        File  multiqcReport = runDuet.multiqcReport
        File  metricsDefinitions = runDuet.metricsDefinitions
    }
}

task splitCram {
    input {
        File cram
        String outputFileNamePrefix
        Int readsPerPart
        String modules
        Int cores = 8
        Int ioSlots = 2
        Int jobMemory = 16
        Int timeout = 24
    }
    parameter_meta {
        cram: "The single unaligned Ultima CRAM to partition"
        outputFileNamePrefix: "Prefix for the part file names"
        readsPerPart: "Reads per part. 2.27B reads / 150000000 gives 16 parts."
        modules: "Environment modules; only samtools is used here"
        cores: "samtools thread count, NOT an SGE slot request: cores-2 threads decode and 2 encode. The UGE backend declares a cpu runtime attribute but never references it in its submit command, so no -pe smp is emitted and the task always gets one slot. Measured 4.7 effective cores at cores=8."
        ioSlots: "gsi_io_slots to request; the one resource lever this backend actually wires through. all.q carries 10 per queue instance (per host) and an ordinary task takes 1, so 2 declares this pass as roughly twice the I/O load of a normal task without crowding a node. Raising it reduces how many other io-slot-consuming jobs share the host, but it cannot help with contention from jobs on other hosts hitting the same filesystem."
        jobMemory: "Memory in GB. The pass is streaming, so this is generous."
        timeout: "Wall-time limit in hours"
    }

    command <<<
        set -euo pipefail

        mkdir -p parts

        # Every part needs the input's header. --no-PG keeps samtools from
        # appending a @PG line per part, which would make the parts differ
        # from each other for no reason.
        samtools view -H --no-PG "~{cram}" > header.sam

        # One streaming pass:
        #   samtools view   decodes CRAM -> SAM on stdout
        #   split -l        counts records and hands each block to its own filter
        #   the filter      prepends the header, re-encodes to CRAM, and indexes
        #
        # Only pipes sit between the three stages, so the disk cost is one read
        # of the input plus one write of the output -- no intermediate FASTQ.
        # Indexing happens inside the filter, while the part it just wrote is
        # still in page cache, which is far cheaper than a second pass over
        # every part afterwards.
        #
        # split runs one filter at a time, so the writer only ever needs a
        # couple of threads; the rest go to the decoder.
        #
        # -e drops reads lacking the Ultima flow tags. prelude cannot resolve those
        # reads and skips them, but its skip path also truncates a NEIGHBOURING
        # record in the resolved FASTQ, which then fails BWA_MEM2 with
        # "SEQ and QUAL are of different length" hours later. Removing them here
        # avoids the trigger entirely. Rejected reads are kept and counted rather
        # than silently discarded. See ref/biomodal_prelude_skip_bug.md; remove this
        # filter once biomodal fixes the skip path.
        #
        # The test is on t0 only. prelude reports the pair together ("missing tag(s)
        # tp t0") and the two are absent together in this data, but tp cannot be
        # tested here regardless: it is a B (array) tag and samtools filter
        # expressions reject those with "Aux type 'B' not yet supported by filters".
        # The assumption is self-checking -- if a later run still logs prelude skips
        # after this filter has run, then reads exist that carry t0 but not tp, and
        # this needs a `grep -F 'tp:B:c'` stage on the SAM stream as well.
        samtools view -@ ~{cores - 2} --no-PG \
                      -e 'exists([t0])' \
                      -U dropped_no_flow_tags.sam \
                      "~{cram}" \
            | split -l ~{readsPerPart} -d -a 4 --numeric-suffixes=1 \
                    --filter='{ cat header.sam; cat; } \
                              | samtools view -@ 2 --no-PG -C -o "${FILE}.cram" - \
                              && samtools index -@ 2 "${FILE}.cram" "${FILE}.cram.crai"' \
                    - "parts/~{outputFileNamePrefix}_part_"

        n_dropped=$(wc -l < dropped_no_flow_tags.sam)
        echo "splitCram: dropped ${n_dropped} read(s) lacking tp/t0; IDs in dropped_no_flow_tags.sam"

        # Fail loudly rather than silently handing the pipeline one lane.
        n=$(find parts -maxdepth 1 -name '*.cram' | wc -l)
        if [ "${n}" -lt 2 ]; then
            echo "ERROR: split produced ${n} part(s). Expected at least 2." >&2
            echo "       Check that readsPerPart (~{readsPerPart}) is smaller than the read count." >&2
            exit 1
        fi
        echo "Split into ${n} parts of ~{readsPerPart} reads:"
        ls -l parts/
    >>>

    output {
        Array[File] parts   = glob("parts/*.cram")
        Array[File] indexes = glob("parts/*.cram.crai")
    }

    runtime {
        # cpu is declared by the UGE backend's runtime-attributes but never used in
        # its submit command, so it emits no -pe smp and this task gets one slot
        # whatever it asks for. Kept so the request is correct if the backend is
        # fixed; until then `cores` only sets samtools -@ threading. io_slots IS
        # wired through (-l gsi_io_slots) and is the one resource lever available.
        cpu:      "~{cores}"
        io_slots: ioSlots
        memory:   "~{jobMemory} GB"
        timeout:  "~{timeout}"
        modules:  "~{modules}"
    }

    meta {
        output_meta: {
            parts: "CRAM parts, each presented to the pipeline as one lane",
            indexes: "CRAM index (.crai) for each part, required by PRELUDE_ULTIMA"
        }
    }
}

task runDuet {
    input {
        Array[File] crams
        Array[File]? craiIndexes
        String sampleId
        String runName
        String outputFileNamePrefix
        String mode
        String additionalProfile
        String modules
        Int hpMinOverlap = 10
        Float hpMaxErrorRate = 0.2
        Int frontQualityTrim = 15
        Int backQualityTrim = 15
        Int meanQualityR1 = 15
        Int meanQualityR2 = 15
        Int maskEndCs  = 5
        Boolean callGermlineVariants = true
        String variantCaller = "deepvariant"
        Int maxCpus = 16
        Int maxMemory = 64
        Int bwaMem2Memory = 40
        Int maxTime = 48
        Int splitReadsPreResolution = -1
        Int gvcfScatterCount = 20
        Int jobMemory = 16
        Int timeout = 96
    }

    Array[File] craiList = select_first([craiIndexes, []])

    parameter_meta {
        crams: "Array of unaligned Ultima single-read CRAM files (one per lane) for a single sample"
        craiIndexes: "Optional pre-built .crai for each entry of crams. Supplied by splitCram, which indexes each part while it is still in page cache. When absent, step 6b builds them here instead."
        sampleId: "Sample identifier (used for naming input CRAMs and output files)"
        runName: "Sequencing run name / identifier (used in report file names)"
        outputFileNamePrefix: "Prefix for all output file names"
        mode: "Biomodal DUET mode: 6bp (duet evoC) or 5bp (duet +modC)"
        additionalProfile: "Nextflow resource profile (deep_seq, super_seq, or empty)"
        modules: "Environment module providing the biomodal instance dir and its env vars"
        hpMinOverlap: "prelude.hp_min_overlap: overlap of the hairpin required to identify and remove it"
        hpMaxErrorRate: "prelude.hp_max_error_rate: error rate tolerated when identifying hairpin sequences for removal"
        frontQualityTrim: "prelude.front_quality_trim: minimum quality below which bases are trimmed from the start of reads"
        backQualityTrim: "prelude.back_quality_trim: minimum quality below which bases are trimmed from the end of reads"
        meanQualityR1: "prelude.mean_quality_R1: minimum mean quality below which an R1 read is discarded"
        meanQualityR2: "prelude.mean_quality_R2: minimum mean quality below which an R2 read is discarded"
        maskEndCs: "prelude.mask_end_cs: mask Cs in the last n bases at the tail of reads to improve methylation-calling sensitivity"
        callGermlineVariants: "Whether to run germline variant calling at all. Set false for a methylation-only run (e.g. when the DeepVariant model is unavailable)"
        variantCaller: "Germline variant caller: deepvariant (Ultima-trained model), gatk, or both"
        maxCpus: "Cap on per-process cpus/slots for heavy steps (PRELUDE, PRELUDE_ULTIMA, BIOMODAL_COLLAPSE, DEEPVARIANT_CALLER always; BWA_MEM2, MUTECT2 when a profile is set). OICR all.q offers at most 39 slots/node (31 on default nodes) but these steps hardcode/request 32-96, so they must be capped to schedule. Lower it (e.g. 8-16) for small test runs or to fit smaller/busier nodes; default 30 fits the 31-slot default nodes."
        maxMemory: "Cap (GB) on per-process memory for heavy steps (MUTECT2, HAPLOTYPE_CALLER, GENOMICS_DB_IMPORT, DEEPVARIANT_CALLER, PRELUDE, PRELUDE_ULTIMA, BIOMODAL_COLLAPSE). These hardcode 32-64GB, so a 64GB h_vmem request only fits the scarce big all.q nodes and can sit in 'qw'. Only reduces (min with the base value), so default 64 is a no-op that preserves production memory; lower it (e.g. 24) so PRELUDE fits the plentiful ~62GB nodes. Does NOT apply to BWA_MEM2 -- see bwaMem2Memory."
        bwaMem2Memory: "Memory in GB for BWA_MEM2, set independently of maxMemory. bwa-mem2 loads its entire index into RAM (~16GB for GRCh38Decoy+controls: .bwt.2bit.64 9.5GB + .0123 5.8GB + .pac 0.7GB) and peaks around 20GB resident, plus page cache from unpacking the index archive in the same cgroup. Sharing PRELUDE's lower cap kills it mid-alignment; biomodal's own 64GB does not fit the common ~62GB nodes. Default 40 clears the measured peak with headroom and still schedules."
        splitReadsPreResolution: "If >0, the pipeline splits the input into chunks of this many reads (seqkit) and runs PRELUDE per chunk in PARALLEL (converting CRAM->FASTQ first), merging before dedup. Essential for very large Ultima samples (~2.3B reads) where a single serial PRELUDE would exceed the wall-time limit. E.g. 285000000 -> ~8 chunks for a 2.3B-read sample. ONLY works with a single input CRAM (one lane) -- do not combine with multiple crams. Default -1 (off)."
        gvcfScatterCount: "Number of genomic intervals to scatter GATK variant calling across (SPLIT_INTERVALS -> per-interval HAPLOTYPE_CALLER/GENOMICS_DB_IMPORT, run in parallel). Increase for more variant-calling parallelism on large genomes/samples. Default 20 (pipeline default)."
        jobMemory: "Memory in GB for the head (Nextflow driver) task"
        maxTime: "Per-process wall-time limit in hours, applied to every Nextflow process as `time`. Distinct from timeout, which bounds the head task."
        timeout: "Wall-time limit in hours for the head (Nextflow driver) task. Must exceed the total runtime of the whole pipeline, not just any one process -- the driver stays alive until the last Nextflow task finishes."
    }

    command <<<
        set -euo pipefail

        # ---------------------------------------------------------------------------
        # 1. Build a writable instance directory. The biomodal CLI needs writable
        #    copies of the two config files it rewrites, plus a real (symlink-free)
        #    pipelines/ tree -- Nextflow cannot follow symlinks for includeConfig
        #    directives, which nextflow.config uses heavily. 
        # ---------------------------------------------------------------------------
        mkdir -p biomodal_instance

        cp -L --remove-destination "$BIOMODAL_INSTANCE_DIR/cli_config.yaml" ./biomodal_instance/cli_config.yaml
        cp -L --remove-destination "$BIOMODAL_INSTANCE_DIR/nextflow_override.config" ./biomodal_instance/nextflow_override.config
        chmod 770 ./biomodal_instance/cli_config.yaml ./biomodal_instance/nextflow_override.config

        cp -rL "$BIOMODAL_INSTANCE_DIR/pipelines" ./biomodal_instance/pipelines
        chmod -R u+w ./biomodal_instance/pipelines

        INSTANCE_DIR="$(pwd)/biomodal_instance"

        # ---------------------------------------------------------------------------
        # 1b. Patch bwa_mem2.nf for the biomodal @PG-header bug. bwa-mem2 writes its
        #     tab-delimited -R read group into the @PG CL: field, producing a
        #     malformed SAM header (duplicate ID tags) that crashes QUALIMAP_BAMQC
        #     under fail_fast. 
        # ---------------------------------------------------------------------------
        BWA_NF="${INSTANCE_DIR}/pipelines/duet/1.7.0a1/modules/bwa_mem2.nf" python3 <<'PYEOF'
import os, sys, pathlib
p = pathlib.Path(os.environ["BWA_NF"])
text = p.read_text()
marker = "OICR hotfix (biomodal @PG bug)"
anchor = "      samtools index -@"
if marker in text:
    print("bwa_mem2.nf: @PG hotfix already present, skipping")
elif anchor in text:
    insert = r'''      # --- OICR hotfix (biomodal @PG bug): bwa-mem2 writes its tab-delimited
      # -R read group into the @PG CL: field, corrupting the SAM header (duplicate
      # ID tags) and crashing qualimap. Collapse the tabs inside @PG CL: to spaces.
      samtools view -H ${bam_file_tag}.bam | sed '/^@PG/{:a; s/\\(CL:.*\\)\\t/\\1 /; ta}' > fixed_header.sam
      samtools reheader fixed_header.sam ${bam_file_tag}.bam > reheadered.bam && mv reheadered.bam ${bam_file_tag}.bam

'''
    text = text.replace(anchor, insert + anchor, 1)
    p.write_text(text)
    print("bwa_mem2.nf: applied @PG hotfix")
else:
    sys.stderr.write("ERROR: bwa_mem2.nf anchor '%s' not found; @PG hotfix NOT applied.\n" % anchor)
    sys.stderr.write("       biomodal may have restructured the module -- re-verify the fix.\n")
    sys.exit(1)
PYEOF

        # ---------------------------------------------------------------------------
        # 1c. Fix a biomodal typo in samtools_cram_to_fastq.nf: SAMTOOLS_CRAM_TO_FASTQ_SINGLE
        #     runs `samtools fastq -@ {task.cpus} ...` -- missing the '$', so Nextflow
        #     passes the literal string {task.cpus}, samtools parses it as 0 threads, and
        #     the CRAM->FASTQ conversion (used in split mode) runs SINGLE-THREADED. Restore ${task.cpus}.
        #     Idempotent (only edits if the buggy pattern is present).
        # ---------------------------------------------------------------------------
        CRAM2FQ_NF="${INSTANCE_DIR}/pipelines/duet/1.7.0a1/modules/samtools_cram_to_fastq.nf"
        if grep -q -- '-@ {task.cpus}' "${CRAM2FQ_NF}"; then
            sed -i 's/-@ {task.cpus}/-@ ${task.cpus}/g' "${CRAM2FQ_NF}"
            echo "Patched samtools_cram_to_fastq.nf: -@ {task.cpus} -> -@ \${task.cpus} (biomodal typo forced single-threaded CRAM->FASTQ)"
        fi

        # NOTE: do not bother patching seqkit flags in split_fastqs.nf. Benchmarked
        # against the module's own seqkit and the pipeline container (v2.9.0) on real
        # chunk data: with GZIPPED input -- which is what SPLIT_FASTQS always gets --
        # -j 8 gives 1.05x and --compress-level 1 gives nothing while writing 10-15%
        # more bytes. (-j 8 does give 2.34x on uncompressed input, which is why it
        # looks promising on paper.) See ref/devlog.txt for the measurements.

        # ---------------------------------------------------------------------------
        # 2. Rewrite cli_config.yaml with runtime paths from the module env vars.
        #    container_engine is apptainer for this release; work dir is task-local.
        # ---------------------------------------------------------------------------
        cat > "${INSTANCE_DIR}/cli_config.yaml" << CLIEOF
        cli:
            max_concurrent_transfers: 6
            max_retries: 3
        computing_platform:
            container_engine: apptainer
            error_strategy: fail_fast
            images_registry_location: ${BIOMODAL_IMAGES_DIR}
            nextflow_work_directory_location: $(pwd)/work
            reference_files_location: ${BIOMODAL_REF_DATA_DIR}
            type: sge
        pipelines:
            duet:
                version: 1.7.0a1
        telemetry:
            share_events: false
            share_metrics: false
CLIEOF

        # ---------------------------------------------------------------------------
        # 3. Append OICR runtime patches to nextflow_override.config.
        # ---------------------------------------------------------------------------

        # 3a. Apptainer image lookup, so containers are reused rather than re-pulled:
        #       libraryDir -> the module's image set, read-only
        #       cacheDir   -> GSI staging dir; writable, and holds images not yet built
        #                     into the module (currently seqkit, needed by split mode)
        #     Nextflow resolves libraryDir first and only falls back to cacheDir, so
        #     module images always win and the staging dir just fills the gaps. Point
        #     cacheDir back at ${BIOMODAL_IMAGES_DIR} once the module ships every image.
        #     runOptions replaces the biomodal-shipped one: keep the $TMPDIR->/tmp
        #     bind, and additionally force TMPDIR=/tmp INSIDE the container.
        IMAGES_STAGING_DIR="/.mounts/labs/gsi/src/biomodal/duet_ultima/images"

        cat >> "${INSTANCE_DIR}/nextflow_override.config" << NFEOF

// ---- OICR WDL runtime patches (env-var expanded) ----
apptainer {
    libraryDir = "${BIOMODAL_IMAGES_DIR}"
    cacheDir   = "${IMAGES_STAGING_DIR}"
    runOptions = '--bind "\$TMPDIR:/tmp" --env TMPDIR=/tmp'
}
NFEOF

        # 3b. Literal Groovy block (no shell expansion). Appended last so it is the
        #     final word on penv / clusterOptions, overriding any biomodal-shipped
        #     process{} defaults.
        #       - NUMBA_CACHE_DIR container env (required by this release).
        #       - penv = 'smp'                       (UGE parallel-environment policy, 1.2)
        #       - clusterOptions schedules on h_vmem  (cgroup memory-kill fix)
        #         and retains -S /bin/bash + -P gsi   (UGE project-name policy)
        cat >> "${INSTANCE_DIR}/nextflow_override.config" << 'NFEOF'

process {
    containerOptions = '--env NUMBA_CACHE_DIR=/tmp/numba_cache'
}

process {
    penv           = 'smp'
    time           = '~{maxTime}h'
    clusterOptions = { "-S /bin/bash -P gsi -l h_vmem=${task.memory.toMega().intdiv(task.cpus)}M" }
}
NFEOF

        # 3c. Clamp per-process cpus that exceed OICR's max smp slots.
        cat >> "${INSTANCE_DIR}/nextflow_override.config" << NFEOF

process {
    withName: 'PRELUDE'            { cpus = ~{maxCpus} }
    withName: 'PRELUDE_ULTIMA'     { cpus = ~{maxCpus} }
    withName: 'BIOMODAL_COLLAPSE'  { cpus = ~{maxCpus} }
    withName: 'DEEPVARIANT_CALLER' { cpus = ~{maxCpus} }
}
NFEOF

        if [ -n "~{additionalProfile}" ]; then
            cat >> "${INSTANCE_DIR}/nextflow_override.config" << NFEOF

process {
    withName: 'BWA_MEM2' { cpus = ~{maxCpus} }
    withName: 'MUTECT2'  { cpus = ~{maxCpus} }
}
NFEOF
        fi

        # 3d. Clamp per-process memory to maxMemory.
        #     BWA_MEM2 is deliberately NOT in this list -- it gets bwaMem2Memory.
        #     bwa-mem2 loads its whole index into RAM and needs far more than the
        #     value that makes PRELUDE schedulable, so one shared cap cannot serve
        #     both: clamped to the same number as PRELUDE it is killed by the cgroup
        #     ("failed 52 : cgroups enforced memory limit" in qacct) partway through
        #     alignment, which surfaces confusingly as "samtools sort: truncated
        #     file" because bwa dies mid-write and the pipe closes.
        {
            echo ""
            echo "process {"
            for spec in "MUTECT2:64" "HAPLOTYPE_CALLER:64" \
                        "GENOMICS_DB_IMPORT:64" "DEEPVARIANT_CALLER:64" \
                        "PRELUDE:32" "PRELUDE_ULTIMA:32" "BIOMODAL_COLLAPSE:32"; do
                pname="${spec%%:*}"; base="${spec##*:}"
                if [ "${base}" -gt "~{maxMemory}" ]; then
                    echo "    withName: '${pname}' { memory = '~{maxMemory}GB' }"
                fi
            done
            echo "    withName: 'BWA_MEM2' { memory = '~{bwaMem2Memory}GB' }"
            echo "}"
        } >> "${INSTANCE_DIR}/nextflow_override.config"

        # 3e. Enable execution tracing here rather than on the command line.
        #     -with-trace / -with-report are Nextflow CLI options (single dash); the
        #     biomodal CLI forwards --additional-params as pipeline params (double
        #     dash), so passing them that way silently does nothing -- .command.trace
        #     stays empty and no aggregate trace is written. Setting them in the
        #     config works regardless of how the CLI forwards arguments.
        #     nf_trace.tsv is the per-process record of runtime, cpu and peak memory,
        #     which is the only way to size maxMemory/bwaMem2Memory from evidence
        #     instead of guessing.
        cat >> "${INSTANCE_DIR}/nextflow_override.config" << NFEOF

trace {
    enabled   = true
    overwrite = true
    file      = "$(pwd)/nf_trace.tsv"
    fields    = 'task_id,name,status,exit,attempt,realtime,%cpu,peak_rss,peak_vmem,rchar,wchar'
}

report {
    enabled   = true
    overwrite = true
    file      = "$(pwd)/nf_report.html"
}

timeline {
    enabled   = true
    overwrite = true
    file      = "$(pwd)/nf_timeline.html"
}
NFEOF

        # ---------------------------------------------------------------------------
        # 4. qsub shim (cgroup memory-kill fix).
        #    clusterOptions above is CONCATENATED onto the s_rss/h_rss/mem_free
        #    directives Nextflow's SGE executor generates from the memory directive;
        #    those re-impose the cgroup RSS limit and must be stripped before qsub.
        #    Install a wrapper named qsub earlier on PATH than the real one; it edits
        #    each .command.run in place (next to the script, on NFS visible to exec
        #    nodes) and re-submits to the real qsub.
        # ---------------------------------------------------------------------------
        mkdir -p ./bin
        cat > ./bin/qsub << 'SHIMEOF'
#!/bin/bash
WRAPPER_DIR=$(dirname "$(readlink -f "$0")")
REAL_QSUB=$(PATH=$(echo "$PATH" | tr ':' '\n' | grep -v "^${WRAPPER_DIR}$" | tr '\n' ':') which qsub)
script_file=""
for arg in "$@"; do
    [[ -f "$arg" && "$arg" != -* ]] && script_file="$arg"
done
if [[ -n "$script_file" ]]; then
    script_dir=$(dirname "$(readlink -f "$script_file")")
    tmp="${script_dir}/.command.run.$$"
    # Strip the RSS/mem_free specs, then DELETE (not blank) any '#$ -l' line left
    # with no specs, so no empty line remains among the #$ directives -- some Grid
    # Engine builds stop scanning directives at the first non-# line.
    sed -e 's/h_rss=[^,]*,\?//g' \
        -e 's/s_rss=[^,]*,\?//g' \
        -e 's/mem_free=[^,]*,\?//g' \
        -e '/^#\$ -l[[:space:]]*,\?[[:space:]]*$/d' \
        "$script_file" > "$tmp"
    chmod +x "$tmp"
    exec "$REAL_QSUB" "${@/$script_file/$tmp}"
else
    exec "$REAL_QSUB" "$@"
fi
SHIMEOF
        chmod +x ./bin/qsub
        export PATH="$(pwd)/bin:$PATH"

        # ---------------------------------------------------------------------------
        # 5. Writable NXF_HOME, pre-seeded with the bundled Nextflow framework jar.
        # ---------------------------------------------------------------------------
        export NXF_HOME="$(pwd)/nxf_home"
        export NXF_OPTS="-Xms512m -Xmx8g"
        # Fully offline run: no network fetches at runtime.
        export NXF_OFFLINE=true
        # One line per event instead of an in-place ANSI progress block. The ANSI
        # renderer truncates process names to fit the terminal, which makes the
        # captured stdout nearly useless for telling PRELUDE from BWA_MEM2.
        # -ansi-log is a Nextflow CLI option, not a pipeline param, so it cannot be
        # passed through the biomodal CLI's --additional-params; the env var can.
        export NXF_ANSI_LOG=false
        JAR=$(find "${INSTANCE_DIR}/pipelines/duet/1.7.0a1" -name "nextflow-*-one.jar" 2>/dev/null | head -1 || true)
        if [ -n "${JAR}" ]; then
            JAR_VER=$(basename "${JAR}" | sed -E 's/^nextflow-(.*)-one\.jar$/\1/')
            mkdir -p "${NXF_HOME}/framework/${JAR_VER}"
            cp "${JAR}" "${NXF_HOME}/framework/${JAR_VER}/"
            echo "Staged Nextflow engine jar: ${JAR} -> ${NXF_HOME}/framework/${JAR_VER}/"
        else
            echo "ERROR: no nextflow-*-one.jar found under ${INSTANCE_DIR}/pipelines/duet/1.7.0a1." >&2
            echo "       The biomodal CLI would try to download the Nextflow engine, which" >&2
            echo "       fails on offline exec nodes. Bundle the jar (e.g." >&2
            echo "       nextflow-25.04.8-one.jar) into the module under pipelines/duet/1.7.0a1/." >&2
            exit 1
        fi

        # ---------------------------------------------------------------------------
        # 6. Stage Ultima CRAM inputs into nf-input using biomodal naming:
        #      {sample-id-no-underscores}_S1_L###_R1_001.cram
        #    Ultima reads are single-read, so there is no R2. The pipeline globs
        #    *.cram and parses the lane number from the _L### field.
        # ---------------------------------------------------------------------------
        SAMPLE_ID="~{sampleId}"
        RUN_NAME="~{runName}"
        SAMPLE_ID_DASH=$(echo "${SAMPLE_ID}" | tr '_' '-')
        mkdir -p nf-input

        sorted_crams=($(for f in ~{sep=' ' crams}; do echo "$f"; done | sort))
        crais=(~{sep=' ' craiList})

        for i in "${!sorted_crams[@]}"; do
            cram="${sorted_crams[$i]}"
            lane=$(printf 'L%03d' "$((i+1))")
            target="nf-input/${SAMPLE_ID_DASH}_S1_${lane}_R1_001.cram"
            ln -s "${cram}" "${target}"
            # Cromwell may localize a pre-built index into a different directory
            # than its CRAM, so match on basename rather than assuming adjacency.
            for crai in ${crais[@]+"${crais[@]}"}; do
                if [ "$(basename "${crai}")" = "$(basename "${cram}").crai" ]; then
                    ln -s "${crai}" "${target}.crai"
                    echo "Linked prebuilt index for lane ${lane}"
                    break
                fi
            done
            echo "Linked lane ${lane}: $(basename "${cram}")"
        done

        # ---------------------------------------------------------------------------
        # 6b. Index each CRAM (.crai). PRELUDE_ULTIMA runs `prelude -r1 <cram> -n N`,
        #     which reads the CRAM in parallel chunks and needs a CRAM index beside
        #     the file to seek to chunk boundaries; a large CRAM without one dies with
        #     "cram_index_load: Could not retrieve index file". SKIPPED ENTIRELY in 
        #     split mode (splitReadsPreResolution > 0): there the
        #     pipeline converts CRAM->FASTQ (sequential, no index) and PRELUDE reads
        #     FASTQ chunks, so the index is never used.
        # ---------------------------------------------------------------------------
        if [ "~{splitReadsPreResolution}" -le 0 ]; then
            for cram in nf-input/*.cram; do
                if [ ! -f "${cram}.crai" ]; then
                    echo "Indexing ${cram} ..."
                    samtools index -@ 4 "${cram}" "${cram}.crai"
                fi
            done
        else
            echo "Split mode (split_reads_pre_resolution=~{splitReadsPreResolution}): PRELUDE reads FASTQ chunks; skipping CRAM index."
        fi

        # ---------------------------------------------------------------------------
        # 7. Run biomodal DUET in Ultima mode.
        #    reference_path is <ref_data>/<ref_pipeline_version>_<ref_genome> =
        #    ${BIOMODAL_REF_DATA_DIR}/1.1.0_GRCh38Decoy 
        # ---------------------------------------------------------------------------
        mkdir -p nf-results
        REFERENCE_PATH="${BIOMODAL_REF_DATA_DIR}/1.1.0_GRCh38Decoy"

        # Fail fast if DeepVariant germline calling is requested but its Ultima model
        # is absent from the reference bundle.  Add the model under
        # ${REFERENCE_PATH}/duet/deepvariant/ultima_model/, set variantCaller=gatk, or
        # set callGermlineVariants=false.
        if [ "~{callGermlineVariants}" = "true" ] && \
           { [ "~{variantCaller}" = "deepvariant" ] || [ "~{variantCaller}" = "both" ]; }; then
            if ! ls "${REFERENCE_PATH}"/duet/deepvariant/ultima_model/checkpoint-* >/dev/null 2>&1; then
                echo "ERROR: variantCaller=~{variantCaller} needs the Ultima DeepVariant model at" >&2
                echo "       ${REFERENCE_PATH}/duet/deepvariant/ultima_model/checkpoint-*" >&2
                echo "       but it is not present in the reference bundle." >&2
                echo "       Add the model, or set variantCaller=gatk, or callGermlineVariants=false." >&2
                exit 1
            fi
        fi

        ADDITIONAL_PROFILE="~{additionalProfile}"
        PROFILE_ARGS=()
        if [ -n "${ADDITIONAL_PROFILE}" ]; then
            PROFILE_ARGS=(--additional-profile "${ADDITIONAL_PROFILE}")
        fi

        # Resolve the biomodal CLI. 
        if [ -x "$BIOMODAL_INSTANCE_DIR/biomodal" ]; then
            BIOMODAL="$BIOMODAL_INSTANCE_DIR/biomodal"
        elif command -v biomodal >/dev/null 2>&1; then
            BIOMODAL="$(command -v biomodal)"
        else
            echo "ERROR: biomodal CLI not found in \$BIOMODAL_INSTANCE_DIR or on PATH" >&2
            exit 1
        fi
        echo "Using biomodal CLI: ${BIOMODAL}"

        "${BIOMODAL}" run duet \
            --instance-directory "${INSTANCE_DIR}" \
            --work-dir "$(pwd)/work" \
            --input-path  "$(pwd)/nf-input" \
            --output-path "$(pwd)/nf-results" \
            --run-name    "${RUN_NAME}" \
            --tag         "${SAMPLE_ID_DASH}" \
            --additional-params "with-report=$(pwd)/nf_report.html" \
            --additional-params "with-trace=$(pwd)/nf_trace.tsv" \
            --additional-params "log=$(pwd)/nextflow.log" \
            "${PROFILE_ARGS[@]}" \
            --additional-params "reference_path=${REFERENCE_PATH}" \
            --additional-params ultima_cram_input=true \
            --additional-params ultima_single_end_input=true \
            --additional-params override_sequencer="ultima" \
            --additional-params input_file_pattern="*.cram" \
            --additional-params split_reads_pre_resolution=~{splitReadsPreResolution} \
            --additional-params gvcf_scatter_count=~{gvcfScatterCount} \
            --additional-params prelude.hp_min_overlap=~{hpMinOverlap} \
            --additional-params prelude.hp_max_error_rate=~{hpMaxErrorRate} \
            --additional-params prelude.front_quality_trim=~{frontQualityTrim} \
            --additional-params prelude.back_quality_trim=~{backQualityTrim} \
            --additional-params prelude.mean_quality_R1=~{meanQualityR1} \
            --additional-params prelude.mean_quality_R2=~{meanQualityR2} \
            --additional-params prelude.mask_end_cs=~{maskEndCs} \
            --additional-params call_germline_variants=~{callGermlineVariants} \
            --additional-params variant_caller="~{variantCaller}" \
            --mode ~{mode}

        # ---------------------------------------------------------------------------
        # 8. Locate the results subdirectory: nf-results/duet-1.7.0a1_<tag>_<mode>/
        # ---------------------------------------------------------------------------
        RESULTS_SUBDIR=$(find "$(pwd)/nf-results" -mindepth 1 -maxdepth 1 \
                           -type d -name "duet-*" | head -1 || true)
        echo "Results subdir: ${RESULTS_SUBDIR}"
        if [ -z "${RESULTS_SUBDIR}" ]; then
            echo "ERROR: could not find a duet-* results directory under nf-results" >&2
            exit 1
        fi

        OUTPUT_PREFIX="~{outputFileNamePrefix}"
        SAMPLE_OUT="${RESULTS_SUBDIR}/sample_outputs"

        # Helper: resolve exactly one file matching a find expression, or fail.
        find_one () {
            # $1 = search root ; remaining args = find predicates
            local root="$1"; shift
            local hit
            hit=$(find "${root}" "$@" 2>/dev/null | head -1 || true)
            if [ -z "${hit}" ]; then
                echo "ERROR: no file under ${root} matching: $*" >&2
                exit 1
            fi
            echo "${hit}"
        }

        # 8a. Deduplicated genome BAM (exclude decoy / non-primary-assembly BAMs).
        BAM=$(find_one "${SAMPLE_OUT}/bams" -maxdepth 1 -name "*.genome.*.dedup.bam")
        ln -s "${BAM}"        "${OUTPUT_PREFIX}.bam"
        ln -s "${BAM}.bai"    "${OUTPUT_PREFIX}.bam.bai"

        # 8b. modC quantification cytosine reports (genome, CpG context).
        MODC_DIR="${SAMPLE_OUT}/modc_quantification"
        HMC=$(find_one  "${MODC_DIR}" -name "*.hmc_cxreport.txt.gz")
        MC=$(find_one   "${MODC_DIR}" -name "*.mc_cxreport.txt.gz")
        MODC=$(find_one "${MODC_DIR}" -name "*.modc_cxreport.txt.gz")
        ln -s "${HMC}"       "${OUTPUT_PREFIX}.hmc_cxreport.txt.gz"
        ln -s "${HMC}.tbi"   "${OUTPUT_PREFIX}.hmc_cxreport.txt.gz.tbi"
        ln -s "${MC}"        "${OUTPUT_PREFIX}.mc_cxreport.txt.gz"
        ln -s "${MC}.tbi"    "${OUTPUT_PREFIX}.mc_cxreport.txt.gz.tbi"
        ln -s "${MODC}"      "${OUTPUT_PREFIX}.modc_cxreport.txt.gz"
        ln -s "${MODC}.tbi"  "${OUTPUT_PREFIX}.modc_cxreport.txt.gz.tbi"

        # 8c. Germline VCF (optional). DeepVariant publishes *.output.vcf.gz under
        #     variant_call_files/deepvariant/germline; GATK publishes
        #     *joint_genotyping.vcf.gz under variant_call_files/germline.
        VCF=$(find "${SAMPLE_OUT}/variant_call_files" \
                   \( -name "*.output.vcf.gz" -o -name "*joint_genotyping.vcf.gz" \) \
                   ! -name "*.g.vcf.gz" ! -name "*atomised*" 2>/dev/null | head -1 || true)
        if [ -n "${VCF}" ] && [ -f "${VCF}" ]; then
            ln -s "${VCF}"       "${OUTPUT_PREFIX}.vcf.gz"
            ln -s "${VCF}.tbi"   "${OUTPUT_PREFIX}.vcf.gz.tbi"
        else
            touch "${OUTPUT_PREFIX}.vcf.gz" "${OUTPUT_PREFIX}.vcf.gz.tbi"
        fi

        # 8d. Run-level reports.
        REPORTS="${RESULTS_SUBDIR}/reports"
        ln -s "$(find_one "${REPORTS}" -name "*Summary.csv")"              "${OUTPUT_PREFIX}.summary.csv"
        ln -s "$(find_one "${REPORTS}" -name "*Summary.html")"             "${OUTPUT_PREFIX}.summary.html"
        ln -s "$(find_one "${REPORTS}" -name "*Summary.xlsx")"             "${OUTPUT_PREFIX}.summary.xlsx"
        ln -s "$(find_one "${REPORTS}" -name "*multiqc_report.html")"      "${OUTPUT_PREFIX}.multiqc_report.html"
        ln -s "$(find_one "${REPORTS}" -name "*Metrics_Definitions.csv")"  "${OUTPUT_PREFIX}.metrics_definitions.csv"
    >>>

    runtime {
        memory:  "~{jobMemory} GB"
        timeout: "~{timeout}"
        modules: "~{modules}"
    }

    output {
        File  outputBam = "~{outputFileNamePrefix}.bam"
        File  outputBai = "~{outputFileNamePrefix}.bam.bai"
        File  hmc_cxreport = "~{outputFileNamePrefix}.hmc_cxreport.txt.gz"
        File  hmc_cxreportIndex  = "~{outputFileNamePrefix}.hmc_cxreport.txt.gz.tbi"
        File  mc_cxreport = "~{outputFileNamePrefix}.mc_cxreport.txt.gz"
        File  mc_cxreportIndex = "~{outputFileNamePrefix}.mc_cxreport.txt.gz.tbi"
        File  modc_cxreport  = "~{outputFileNamePrefix}.modc_cxreport.txt.gz"
        File  modc_cxreportIndex = "~{outputFileNamePrefix}.modc_cxreport.txt.gz.tbi"
        File? vcf = "~{outputFileNamePrefix}.vcf.gz"
        File? vcfIndex = "~{outputFileNamePrefix}.vcf.gz.tbi"
        File  summaryCsv = "~{outputFileNamePrefix}.summary.csv"
        File  summaryHtml = "~{outputFileNamePrefix}.summary.html"
        File  summaryXlsx = "~{outputFileNamePrefix}.summary.xlsx"
        File  multiqcReport = "~{outputFileNamePrefix}.multiqc_report.html"
        File  metricsDefinitions = "~{outputFileNamePrefix}.metrics_definitions.csv"
    }
}
