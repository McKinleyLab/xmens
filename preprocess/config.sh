# Paths shared by every step; each step script sources this file.
#
# OUT, the run folder, has no default, so a run never writes into an earlier one by accident. Set it once per run:
#   export OUT=/n/netscratch/eddy_lab/Everyone/nhilgert/$(date +%Y_%m_%d)_xmens_preprocess

OUT=${OUT:?export OUT, the run folder, before running a step (see README.md)}
NH=/n/lab_storage/eddy_lab/users/nhilgert

# Step 1: Cell Ranger 7.2.0, the FASTQs and the references (mouse: mm10-2020-A plus the GsD/GqD transgenes)
CELLRANGER=$NH/software/cellranger-7.2.0/cellranger
FASTQS=$NH/xmens_final/dataverse/fastqs
REF_MOUSE=$NH/reference/transgenic_mus
REF_HUMAN=$NH/reference/human/refdata-gex-GRCh38-2024-A

# Step 2: CellBender 0.2.2
CONDA=/n/sw/Miniforge3-25.3.1-0
CELLBENDER_ENV=$NH/.conda/envs/cellbender_022

# Steps 3-5: rsc_uv for ingest, concat and positions; the DoubletDetection 4.2 environment for QC
PY=$NH/.uv/envs/rsc_uv/bin/python
PY_QC=$NH/.conda/envs/cuda11.8_env/bin/python

# Slide-tags positions: Nick's run (attached in step 3) and Andy's original positioning (added in step 5)
COORDS_DIR=$NH/Slide-tags/coords
ANDY_DIR=/n/lab_storage/eddy_lab/Lab/mckinley_backup_2025-12-18/data/andy_raw_data

# Where steps 2 and 3 read their inputs; {lane} and {sample} are filled in per task. To start from earlier results,
# export either one before submitting, e.g. the July 2026 CellBender files behind the regen objects:
#   export CELLBENDER_H5=$NH/decidualization/xmens_inputs_regen/{sample}/{sample}_cellbender_filtered.h5
# (Plain ifs: a "}" inside ${VAR:-default} would end the default early.)
if [ -z "${CELLRANGER_RAW:-}" ]; then
    CELLRANGER_RAW="$OUT/cellranger/{lane}/outs/raw_feature_bc_matrix.h5"
fi
if [ -z "${CELLBENDER_H5:-}" ]; then
    CELLBENDER_H5="$OUT/cellbender/{sample}/{sample}_cellbender_filtered.h5"
fi

# Caches go to lab storage, never to the home directory.
CACHE=$NH/.cache_scratch
export HDF5_USE_FILE_LOCKING=FALSE
export MPLCONFIGDIR=$CACHE/mpl
export NUMBA_CACHE_DIR=$CACHE/numba
export XDG_CACHE_HOME=$CACHE/xdg

# field <row> <column>: one value from samples.tsv; row 1 is the first sample, as is array task 1.
field() {
    awk -F'\t' -v row="$1" -v col="$2" '
        NR == 1 { for (i = 1; i <= NF; i++) if ($i == col) c = i }
        NR == row + 1 { print $c }' "$HERE/samples.tsv"
}

N_SAMPLES=$(($(wc -l < "$HERE/samples.tsv") - 1))
