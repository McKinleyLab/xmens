#!/bin/bash
# Step 5 of 5: the SAMap figures, drawn by the figure notebooks of record (notebooks/5_figures/, kernel rsc), each
# run with papermill on the tables of step 4:
#   Supplementary Fig. S30, panels a-g and the assembled figure       -> $OUT/5_figures/sfig_samap/ (panels in panels/)
#   Fig. 6C (transferred pseudotime on the joint UMAP), 6E (in the two human sections) and the cell-identity twin
#   of 6C (transferred cell type), each with a .provenance.json      -> $OUT/5_figures/fig6/
#   Reads:  $OUT/4_statistics/ tables/, graph_block/, edge_share/, ptnull/ (step 4)
#   Logs:   $OUT/notebooks/5_figures/<notebook>.ipynb (one executed notebook per figure) and logs/5_figures_<job ID>.out
# A notebook whose executed copy exists is skipped; one that did not finish is rerun and overwrites what it wrote.
#
# Usage:
#   bash 5_figures.sbatch --check
#   sbatch 5_figures.sbatch

#SBATCH --job-name=5_figures
#SBATCH --partition=eddy
#SBATCH --account=eddy_lab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out  # job name and job ID; relative to the folder sbatch was run from

# Stop at the first failed command, unset variable or failed stage of a pipe.
set -euo pipefail

# The pipeline folder holds config.sh. Run with bash, this script sits in it; under sbatch, SLURM runs a copy from its
# spool folder, and the pipeline folder is then the one sbatch was run from.
HERE=$(cd "$(dirname "$0")" && pwd)
if [ ! -e "$HERE/config.sh" ]; then
    HERE=${SLURM_SUBMIT_DIR:-$PWD}
fi
if [ ! -e "$HERE/config.sh" ]; then
    echo "submit this from the samap folder (no config.sh in $HERE)" >&2
    exit 2
fi

# config.sh sets $OUT, run_notebook and exists.
source "$HERE/config.sh"

NB=$HERE/notebooks/5_figures
T=$OUT/4_statistics
F=$OUT/5_figures

# The tables each notebook reads, as papermill parameters (the notebooks' own names).
JOINT=(-p JOINT "$T/tables/joint.parquet")
CELLS=(-p MOUSE_CELLS "$T/tables/mouse_cells.csv" -p HUMAN_CELLS "$T/tables/human_cells.csv")
HUMAN_TABLE=(-p HUMAN_TABLE "$T/tables/human_plot_table.parquet")

if [ "${1:-}" = "--check" ]; then
    for f in "$T/tables/joint.parquet" "$T/tables/mouse_cells.csv" "$T/tables/human_cells.csv" \
             "$T/tables/human_plot_table.parquet" "$T/graph_block/C_hs_mm.npz" \
             "$T/edge_share/samap_compartment_edge_share.csv" "$T/ptnull/null_replicates.csv" \
             "$T/ptnull/example_transfers.csv"; do
        echo "reads:  $f ($(exists "$f"))"
    done
    for nb in "$NB"/*.ipynb; do
        name=$(basename "$nb" .ipynb)
        echo "$name: $(exists "$OUT/notebooks/5_figures/$name.ipynb")"
    done
    exit 0
fi

# figure <notebook> <output folder> [-p NAME VALUE ...]: draws one figure unless its executed copy exists. REPO is
# only the base of the notebooks' relative paths; every path given here is absolute.
figure() {
    local name=$1 out=$2
    shift 2
    local executed=$OUT/notebooks/5_figures/$name.ipynb
    if [ -e "$executed" ]; then
        echo "$name: $executed exists; skipped"
        return
    fi
    mkdir -p "$out"
    run_notebook "$NB/$name.ipynb" rsc "$executed" \
        -p REPO "$HERE" -p NOTEBOOK "$NB/$name.ipynb" -p OUT "$out" -p FORCE True "$@"
}

# Supplementary Fig. S30: the seven panels, then the assembled figure.
figure sfig_samap_a_umap_mouse_pseudotime       "$F/sfig_samap/panels" "${JOINT[@]}" "${CELLS[@]}"
figure sfig_samap_b_umap_transferred_pseudotime "$F/sfig_samap/panels" "${JOINT[@]}" "${CELLS[@]}"
figure sfig_samap_c_umap_human_pseudotime       "$F/sfig_samap/panels" "${JOINT[@]}" "${CELLS[@]}"
figure sfig_samap_d_edge_share                  "$F/sfig_samap/panels" -p EDGE_SHARE "$T/edge_share/samap_compartment_edge_share.csv"
figure sfig_samap_e_transferred_vs_human        "$F/sfig_samap/panels" -p HUMAN_CELLS "$T/tables/human_cells.csv"
figure sfig_samap_f_sections_shuffled           "$F/sfig_samap/panels" -p HUMAN_CELLS "$T/tables/human_cells.csv" \
    -p EXAMPLE_TRANSFERS "$T/ptnull/example_transfers.csv"
figure sfig_samap_g_permutation_statistics      "$F/sfig_samap/panels" -p NULL_REPLICATES "$T/ptnull/null_replicates.csv" \
    -p EXAMPLE_TRANSFERS "$T/ptnull/example_transfers.csv" -p HUMAN_CELLS "$T/tables/human_cells.csv"
figure sfig_samap_compose                       "$F/sfig_samap" -p PANELS "$F/sfig_samap/panels"

# Fig. 6C and 6E, and the cell-identity twin of 6C.
figure fig6c_umap_transferred_pseudotime        "$F/fig6" "${JOINT[@]}" "${HUMAN_TABLE[@]}"
figure fig6e_sections_transferred_pseudotime    "$F/fig6" "${JOINT[@]}" "${HUMAN_TABLE[@]}"
figure fig6_umap_transferred_cell_type          "$F/fig6" "${JOINT[@]}" "${CELLS[@]}" -p GRAPH_DIR "$T/graph_block"
