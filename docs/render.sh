#!/bin/bash
set -euo pipefail

# Render the workflow diagram.
#
#   bash docs/render.sh          SVG only -- this is what is committed and what
#                                README.md references; GitHub renders it inline.
#   bash docs/render.sh --png    also render the PNG for a Confluence attachment.
#
# The PNG is gitignored on purpose: Confluence will not render an uploaded SVG inline, so
# the wiki needs a raster copy, but at half a megabyte of binary that would change on every
# edit it is rendered on demand rather than committed.

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
src="${here}/biomodalDuetUltima.dot"

if ! command -v dot >/dev/null 2>&1; then
    echo "error: graphviz 'dot' not found on PATH" >&2
    exit 1
fi

dot -Tsvg "${src}" -o "${here}/biomodalDuetUltima.svg"
echo "wrote ${here}/biomodalDuetUltima.svg"

if [ "${1:-}" = "--png" ]; then
    dot -Tpng -Gdpi=150 "${src}" -o "${here}/biomodalDuetUltima.png"
    echo "wrote ${here}/biomodalDuetUltima.png"
fi
