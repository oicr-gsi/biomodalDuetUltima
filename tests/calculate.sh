#!/bin/bash
set -euo pipefail

# Metrics for a biomodalDuetUltima run. Deliberately structural, not numeric: the
# pipeline's methylation percentages, coverage figures and control sensitivities all
# move with the biomodal release and with the sample, so pinning them would make every
# upgrade a test failure. What IS worth pinning is the SHAPE of the output -- which
# files were provisioned, and the schema of the tabular ones.
#
# Coreutils and zcat only. Nothing here loads an environment module, so the metric does
# not depend on what the Jenkins runner happens to have on PATH.

# Bytewise collation, so the baseline does not depend on the runner's locale.
export LC_ALL=C

cd "$1" || exit 1

echo "### provisioned files"
# `! -type d`, NOT `-type f`: Vidarr provisions outputs as SYMLINKS, so -type f matches
# nothing and this section would silently record an empty list -- which would then never
# detect a missing output.
find . -maxdepth 1 ! -type d -printf '%f\n' | sort

echo "### cytosine report schema"
for report in *cxreport.txt.gz; do
    [ -e "$report" ] || continue
    echo "-- $report"
    # Field count of the first record. The per-position calls underneath are expected to
    # drift; the column count is the schema and should not change silently.
    zcat "$report" | head -1 | awk -F'\t' '{print "fields: " NF}'
done

echo "### summary metrics present"
for summary in *.summary.csv; do
    [ -e "$summary" ] || continue
    echo "-- $summary"
    # Column names only. A new, renamed or dropped metric is a real change worth failing
    # on; the values underneath are expected to move with every sample and release.
    head -1 "$summary" | tr ',' '\n' | sort
done

echo "### metrics definitions"
for defs in *.metrics_definitions.csv; do
    [ -e "$defs" ] || continue
    echo "-- $defs"
    # The metric ids biomodal documents. Catches a metric being added or removed between
    # releases, which is exactly the sort of change that should be noticed deliberately.
    tail -n +2 "$defs" | cut -d, -f1 | sort
done
