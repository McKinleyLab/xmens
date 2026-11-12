# Paths and settings shared by every step; each step script sources this file.
#
# OUT, the run folder, has no default, so a run never writes into an earlier one by accident. Set it once per run:
#   export OUT=/n/netscratch/eddy_lab/Everyone/nhilgert/$(date +%Y_%m_%d)_xmens_mi

OUT=${OUT:?export OUT, the run folder, before running a step (see README.md)}
NH=/n/lab_storage/eddy_lab/users/nhilgert
NS=/n/netscratch/eddy_lab/Everyone/nhilgert

# The Python environment of the run of record: Python 3.11, scikit-learn 1.8.0, numpy 2.2.0, pandas 2.3.3, scipy 1.17.1,
# statsmodels 0.14.6, anndata 0.12.11, openpyxl 3.1.5.
PY=$NH/.uv/envs/rsc_uv/bin/python

# The final objects (read only): fibroblast and decidual nuclei, raw counts in layers['counts'], and the rank-normalized
# pseudotime in obs['pseudotime'] (NaN for the 1,175 apoptotic mouse decidual nuclei). netscratch is purged after 90 days:
# point OBJECTS at a permanent copy of the same files (sha256 in record/OUTPUTS.tsv) once there is one.
# The fibroblast/decidual objects of 6-10 Nov 2026 with the pseudotime of 10 Nov. (The deposited objects of
# record: $NS/final_objects/Cevrim_XMens_Only_Fib_Dec_Cells.h5ad and Cevrim_Human_Only_Fib_Dec_Cells.h5ad,
# layer counts.) Step 1 checks the checksums before reading.
OBJECTS=/n/lab_storage/eddy_lab/Lab/mckinley/data
H5AD_MOUSE=$OBJECTS/20261110_mouse_only_eFib_dec_final_nick_pdt.h5ad
H5AD_HUMAN=$OBJECTS/20261107_Human_fibdec_final_nick_pdt.h5ad
H5AD_MOUSE_SHA256=e31372ba5fdd706ad620e8259e34216bc96b8f3b6394e4d2eedc0baeec8760dd
H5AD_HUMAN_SHA256=12429594a2e71ba39861e96edb0a19d23088fd1f026e11f7ade8856a10a4c425
COUNTS_LAYER=raw_counts  # the objects' raw-counts layer

# Ensembl Compara release 116 human-mouse homology components (read only).
ORTHOLOGS=$NH/human_ivd_repro/orthogroups/ensembl_human_mouse_components.tsv

# Gene chunks per species, one array task each in step 2. Each chunk draws its own 1,000 permutations of the pseudotime,
# seeded by the chunk number, so p and q depend on these numbers: 200 and 1,500 reproduce the run of record. If you
# change them, change the --array line of 2_mi.sbatch to 0-(N_CHUNKS_MOUSE + N_CHUNKS_HUMAN - 1).
N_CHUNKS_MOUSE=200
N_CHUNKS_HUMAN=1500

# Where step 3 reads the per-chunk results; {species} is filled in. To start at step 3 from the run of record's chunks,
# export this before submitting:
#   export CHUNKS=/n/lab_storage/eddy_lab/users/nhilgert/xmens_final/results/mi_dpt/pt_harmonized_{species}_nall/results
# (Plain if: a "}" inside ${VAR:-default} would end the default early.)
if [ -z "${CHUNKS:-}" ]; then
    CHUNKS="$OUT/chunks/{species}"
fi

# Caches go to lab storage, never to the home directory.
CACHE=$NH/.cache_scratch
export HDF5_USE_FILE_LOCKING=FALSE
export MPLCONFIGDIR=$CACHE/mpl
export NUMBA_CACHE_DIR=$CACHE/numba
export XDG_CACHE_HOME=$CACHE/xdg

# h5ad <species> and n_chunks <species> print the species' value from above.
h5ad() {
    case $1 in
        mouse) echo "$H5AD_MOUSE" ;;
        human) echo "$H5AD_HUMAN" ;;
    esac
}
n_chunks() {
    case $1 in
        mouse) echo "$N_CHUNKS_MOUSE" ;;
        human) echo "$N_CHUNKS_HUMAN" ;;
    esac
}

# Prints "found" or "missing", for --check.
exists() {
    if [ -e "$1" ]; then
        echo found
    else
        echo missing
    fi
}
