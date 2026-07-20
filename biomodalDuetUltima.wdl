version 1.0

workflow biomodalDuetUltima {
    input {
        Array[File] crams
        String sampleId
        String runName
        String outputFileNamePrefix
        String mode = "6bp"
        String additionalProfile = "deep_seq"
        String modules = "biomodal-duet-ultima/1.7.0a1"
    }

    parameter_meta {
        crams:                "Array of unaligned Ultima single-read CRAM files (one per lane) for a single sample"
        sampleId:             "Sample identifier (used for naming input CRAMs and output files)"
        runName:              "Sequencing run name / identifier (used in report file names)"
        outputFileNamePrefix: "Prefix for all output file names"
        mode:                 "Biomodal DUET mode: 6bp (duet evoC) or 5bp (duet +modC). Default: 6bp"
        additionalProfile:    "Nextflow resource profile: deep_seq (<=500M reads), super_seq (>500M reads), or empty (<=50M reads). Default: deep_seq"
        modules:              "Environment module providing the biomodal instance dir and its env vars (apptainer loads as a dependency)"
    }

    call runDuet {
        input:
            crams                = crams,
            sampleId             = sampleId,
            runName              = runName,
            outputFileNamePrefix = outputFileNamePrefix,
            mode                 = mode,
            additionalProfile    = additionalProfile,
            modules              = modules
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
        File    outputBam          = runDuet.outputBam
        File    outputBai          = runDuet.outputBai
        File    hmc_cxreport       = runDuet.hmc_cxreport
        File    hmc_cxreportIndex  = runDuet.hmc_cxreportIndex
        File    mc_cxreport        = runDuet.mc_cxreport
        File    mc_cxreportIndex   = runDuet.mc_cxreportIndex
        File    modc_cxreport      = runDuet.modc_cxreport
        File    modc_cxreportIndex = runDuet.modc_cxreportIndex
        File?   vcf                = runDuet.vcf
        File?   vcfIndex           = runDuet.vcfIndex
        File    summaryCsv         = runDuet.summaryCsv
        File    summaryHtml        = runDuet.summaryHtml
        File    summaryXlsx        = runDuet.summaryXlsx
        File    multiqcReport      = runDuet.multiqcReport
        File    metricsDefinitions = runDuet.metricsDefinitions
    }
}

task runDuet {
    input {
        Array[File] crams
        String sampleId
        String runName
        String outputFileNamePrefix
        String mode
        String additionalProfile
        String modules
        # Ultima recommended quality-filtering / trimming parameters (see
        # duet_on_Ultima_early_access_instructions_1-7-0a1). Defaults are the
        # biomodal-recommended values chosen to maximise recovered reads.
        Int    hpMinOverlap     = 10
        Float  hpMaxErrorRate   = 0.2
        Int    frontQualityTrim = 15
        Int    backQualityTrim  = 15
        Int    meanQualityR1    = 15
        Int    meanQualityR2    = 15
        Int    maskEndCs        = 5
        Boolean callGermlineVariants = true
        String variantCaller    = "deepvariant"
        Int    maxCpus          = 30
        Int    jobMemory        = 16
        Int    timeout          = 96
    }
    parameter_meta {
        crams:                "Array of unaligned Ultima single-read CRAM files (one per lane) for a single sample"
        sampleId:             "Sample identifier (used for naming input CRAMs and output files)"
        runName:              "Sequencing run name / identifier (used in report file names)"
        outputFileNamePrefix: "Prefix for all output file names"
        mode:                 "Biomodal DUET mode: 6bp (duet evoC) or 5bp (duet +modC)"
        additionalProfile:    "Nextflow resource profile (deep_seq, super_seq, or empty)"
        modules:              "Environment module providing the biomodal instance dir and its env vars"
        hpMinOverlap:         "prelude.hp_min_overlap: overlap of the hairpin required to identify and remove it"
        hpMaxErrorRate:       "prelude.hp_max_error_rate: error rate tolerated when identifying hairpin sequences for removal"
        frontQualityTrim:     "prelude.front_quality_trim: minimum quality below which bases are trimmed from the start of reads"
        backQualityTrim:      "prelude.back_quality_trim: minimum quality below which bases are trimmed from the end of reads"
        meanQualityR1:        "prelude.mean_quality_R1: minimum mean quality below which an R1 read is discarded"
        meanQualityR2:        "prelude.mean_quality_R2: minimum mean quality below which an R2 read is discarded"
        maskEndCs:            "prelude.mask_end_cs: mask Cs in the last n bases at the tail of reads to improve methylation-calling sensitivity"
        callGermlineVariants: "Whether to run germline variant calling at all. Set false for a methylation-only run (e.g. when the DeepVariant model is unavailable)"
        variantCaller:        "Germline variant caller: deepvariant (Ultima-trained model), gatk, or both"
        maxCpus:              "Cap on per-process cpus/slots for heavy steps (PRELUDE, PRELUDE_ULTIMA, BIOMODAL_COLLAPSE, DEEPVARIANT_CALLER always; BWA_MEM2, MUTECT2 when a profile is set). OICR all.q offers at most 39 slots/node (31 on default nodes) but these steps hardcode/request 32-96, so they must be capped to schedule. Lower it (e.g. 8-16) for small test runs or to fit smaller/busier nodes; default 30 fits the 31-slot default nodes."
        jobMemory:            "Memory in GB for the head (Nextflow driver) task"
        timeout:              "Timeout in hours"
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
        module use /.mounts/labs/gsiprojects/gsi/gsiusers/gpeng/modules/local/gsi/modulator/modulefiles/Ubuntu20.04
        module load biomodal-duet-ultima/1.7.0a1
        cp -L --remove-destination "$BIOMODAL_INSTANCE_DIR/cli_config.yaml" ./biomodal_instance/cli_config.yaml
        cp -L --remove-destination "$BIOMODAL_INSTANCE_DIR/nextflow_override.config" ./biomodal_instance/nextflow_override.config
        chmod 770 ./biomodal_instance/cli_config.yaml ./biomodal_instance/nextflow_override.config

        cp -rL "$BIOMODAL_INSTANCE_DIR/pipelines" ./biomodal_instance/pipelines

        INSTANCE_DIR="$(pwd)/biomodal_instance"

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

        # 3a. Point the apptainer image cache at the shared module images dir so
        #     containers are pulled once and reused (env-var expanded here).
        cat >> "${INSTANCE_DIR}/nextflow_override.config" << NFEOF

// ---- OICR WDL runtime patches (env-var expanded) ----
apptainer {
    libraryDir = "${BIOMODAL_IMAGES_DIR}"
    cacheDir   = "${BIOMODAL_IMAGES_DIR}"
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
    clusterOptions = { "-S /bin/bash -P gsi -l h_vmem=${task.memory.toGiga()}g" }
}
NFEOF

        # 3c. Clamp per-process cpus that exceed OICR's max smp slots. all.q offers
        #     at most 39 slots per node (31 on default nodes, 39 on 40-core
        #     hostgroups, 23 on 24-core); a request above a node's slot count sits
        #     in 'qw' forever. Two groups need capping to ~{maxCpus}:
        #      - PRELUDE, PRELUDE_ULTIMA, BIOMODAL_COLLAPSE, DEEPVARIANT_CALLER
        #        hardcode 'cpus 32' in their .nf, so they exceed the 31-slot default
        #        nodes in EVERY tier -> clamp unconditionally.
        #      - BWA_MEM2, MUTECT2 are small at base but the deep_seq/super_seq
        #        profiles push them to 32/64 -> clamp only when a profile is active
        #        (forcing their small base value up would be counterproductive).
        #     Memory needs no clamp: all.q has 256-768GB big-memory nodes that
        #     satisfy the 128GB dedup steps.
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
        #    The biomodal CLI bootstraps the Nextflow engine (nextflow-<ver>-one.jar)
        #    on run; on an offline exec node it cannot download it (curl 403). The
        #    module must therefore ship the jar under pipelines/duet/1.7.0a1/; we
        #    copy it into $NXF_HOME/framework/<ver>/ (where the launcher looks before
        #    downloading). If it is missing we fail fast here with a clear message
        #    rather than letting the CLI hit the network and emit an opaque 403.
        # ---------------------------------------------------------------------------
        export NXF_HOME="$(pwd)/nxf_home"
        export NXF_OPTS="-Xms512m -Xmx8g"
        # Fully offline run: no network fetches at runtime.
        export NXF_OFFLINE=true
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

        for i in "${!sorted_crams[@]}"; do
            cram="${sorted_crams[$i]}"
            lane=$(printf 'L%03d' "$((i+1))")
            ln -s "${cram}" "nf-input/${SAMPLE_ID_DASH}_S1_${lane}_R1_001.cram"
            echo "Linked lane ${lane}: $(basename "${cram}")"
        done

        # ---------------------------------------------------------------------------
        # 7. Run biomodal DUET in Ultima mode.
        #    reference_path is <ref_data>/<ref_pipeline_version>_<ref_genome> =
        #    ${BIOMODAL_REF_DATA_DIR}/1.1.0_GRCh38Decoy (same convention as the
        #    v1.5.0 WDL's 1.0.5_GRCh38Decoy). The pipeline reads TSS/prelude/
        #    deepvariant/blacklist from $reference_path/duet/... and the genome /
        #    control references from $reference_path/duet/duet-ref-1.1.0/... .
        #    ultima_cram_input / ultima_single_end_input / override_sequencer /
        #    input_file_pattern are the required Ultima flags.
        # ---------------------------------------------------------------------------
        mkdir -p nf-results
        REFERENCE_PATH="${BIOMODAL_REF_DATA_DIR}/1.1.0_GRCh38Decoy"

        # Fail fast if DeepVariant germline calling is requested but its Ultima model
        # is absent from the reference bundle. The pipeline passes the model path
        # straight to run_deepvariant (no startup file-existence check), so without
        # this guard a missing model would only surface after hours of alignment and
        # quantification. Add the model under
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

        # Resolve the biomodal CLI. It is NOT part of the downloaded pipeline files;
        # the module bundles it inside the instance dir ($BIOMODAL_INSTANCE_DIR/
        # biomodal). Prefer that, fall back to PATH (in case a future module puts it
        # there instead), and fail loudly if neither is found.
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
        File    outputBam          = "~{outputFileNamePrefix}.bam"
        File    outputBai          = "~{outputFileNamePrefix}.bam.bai"
        File    hmc_cxreport       = "~{outputFileNamePrefix}.hmc_cxreport.txt.gz"
        File    hmc_cxreportIndex  = "~{outputFileNamePrefix}.hmc_cxreport.txt.gz.tbi"
        File    mc_cxreport        = "~{outputFileNamePrefix}.mc_cxreport.txt.gz"
        File    mc_cxreportIndex   = "~{outputFileNamePrefix}.mc_cxreport.txt.gz.tbi"
        File    modc_cxreport      = "~{outputFileNamePrefix}.modc_cxreport.txt.gz"
        File    modc_cxreportIndex = "~{outputFileNamePrefix}.modc_cxreport.txt.gz.tbi"
        File?   vcf                = "~{outputFileNamePrefix}.vcf.gz"
        File?   vcfIndex           = "~{outputFileNamePrefix}.vcf.gz.tbi"
        File    summaryCsv         = "~{outputFileNamePrefix}.summary.csv"
        File    summaryHtml        = "~{outputFileNamePrefix}.summary.html"
        File    summaryXlsx        = "~{outputFileNamePrefix}.summary.xlsx"
        File    multiqcReport      = "~{outputFileNamePrefix}.multiqc_report.html"
        File    metricsDefinitions = "~{outputFileNamePrefix}.metrics_definitions.csv"
    }
}
