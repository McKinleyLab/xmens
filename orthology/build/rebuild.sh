#!/bin/bash
# Rebuilds ensembl_human_mouse_components.tsv from the frozen Ensembl 116 sources in build/sources/ and checks that the
# table and both validation records come out byte-identical to the ones in this folder.
#
# The builder writes its edge, summary and attrition tables next to its sources, so the rebuild runs on a copy in a
# scratch folder and leaves this folder untouched. It needs only Python 3 (standard library) and takes under a minute.
#
# Usage: bash build/rebuild.sh [scratch folder]      (default: a new folder in $TMPDIR, or in /tmp)

# Stop at the first failed command, unset variable or failed stage of a pipe.
set -euo pipefail

HERE=$(cd "$(dirname "$0")/.." && pwd)
PY=${PY:-/n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/rsc_uv/bin/python}
WORK=${1:-$(mktemp -d "${TMPDIR:-/tmp}/orthology_rebuild.XXXXXX")}

# The builder and validator expect the layout of the project they come from (human_ivd_repro):
# <root>/orthogroups/sources/ for the inputs, <root>/orthogroups/ensembl_human_mouse_components.tsv for the table.
mkdir -p "$WORK/orthogroups"
cp -r "$HERE/build/sources" "$WORK/orthogroups/sources"
(cd "$WORK/orthogroups/sources" && sha256sum --quiet -c SHA256SUMS)

"$PY" -I "$HERE/build/build_ensembl116_components.py" --repo-root "$WORK" > "$WORK/build.log"
"$PY" -I "$HERE/build/validate_ensembl116_components.py" --repo-root "$WORK" > "$WORK/validate.log"

cmp "$WORK/orthogroups/ensembl_human_mouse_components.tsv" "$HERE/ensembl_human_mouse_components.tsv"
for record in ensembl_116_component_validation.json ensembl_116_component_independent_validation.json; do
    cmp "$WORK/orthogroups/sources/$record" "$HERE/provenance/$record"
done
echo "rebuild identical to the table of record and its two validation records; scratch copy in $WORK"
