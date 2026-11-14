# Paths and settings shared by every step. Export OUT, the run folder, first; use a new one for each run:
#   export OUT=/n/netscratch/eddy_lab/Everyone/nhilgert/$(date +%Y_%m_%d)_xmens_samap_full

OUT=${OUT:?export OUT, the run folder, before running a step (see README.md)}
NH=/n/lab_storage/eddy_lab/users/nhilgert
NS=/n/netscratch/eddy_lab/Everyone/nhilgert

# The final fibroblast/decidual objects with pseudotime (steps 1, 3 and 4).
FINAL_OBJECTS=/n/lab_storage/eddy_lab/Lab/mckinley/data
MOUSE_H5AD=$FINAL_OBJECTS/20261110_mouse_only_eFib_dec_final_nick_pdt.h5ad
HUMAN_H5AD=$FINAL_OBJECTS/20261112_Human_fibdec_final_nick_pdt.h5ad

# The objects' raw-counts layer (step 1) and the human section coordinates (step 4).
COUNTS_LAYER=raw_counts
X_COLUMN=x_um_dbscan
Y_COLUMN=y_um_dbscan

# Step 0: the proteomes of Ensembl release 113, decompressed from
#   https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz
#   https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/pep/Mus_musculus.GRCm39.pep.all.fa.gz
# and phmmer of HMMER 3.4.
PROTEOMES=$NH/data/proteomes/fasta
HUMAN_PROTEOME=$PROTEOMES/Homo_sapiens.GRCh38.pep.all.fa
MOUSE_PROTEOME=$PROTEOMES/Mus_musculus.GRCm39.pep.all.fa
PHMMER=$NH/software/hmmer-3.4/src/phmmer

# Checksums of the four input files. A step that reads one checks it first and stops if it differs.
MOUSE_H5AD_SHA256=e31372ba5fdd706ad620e8259e34216bc96b8f3b6394e4d2eedc0baeec8760dd
HUMAN_H5AD_SHA256=c7d4bacb40f6472000c372cd005c6712a9c3f0e112ae5b23fb7df04b31dea82b
HUMAN_PROTEOME_SHA256=f38587d8e87564a08f99378c3416b47a866cd0aff5ed0bdeaeab0ca2d2a8ea3b
MOUSE_PROTEOME_SHA256=54cb3f64f6dd1d0a398713016e866c355cff8e21968f8d18e0d60f7b498712d6

# Where steps 0b-4 read their inputs: by default the outputs of the step before, in OUT. Export one of these to
# start from earlier results; steps whose output would then go unused skip themselves.
DEFAULT_PHMMER_HITS=$OUT/0_homology/phmmer  # hu_vs_mu/hu_<i>.tbl and mu_vs_hu/mu_<i>.tbl
DEFAULT_MAPS=$OUT/0_homology  # SAMap reads $MAPS/hsmm/hs_to_mm.txt and mm_to_hs.txt
DEFAULT_SAMAP_INPUTS=$OUT/1_inputs  # folder holding hs_fib_dec_all_raw.h5ad and mm_fib_dec_all_raw.h5ad
DEFAULT_SAMAP_JOINT=$OUT/2_samap/samap_joint.h5ad
DEFAULT_SAMAP_FINAL=$OUT/3_transfer/Samap_Final.h5ad
if [ -z "${PHMMER_HITS:-}" ]; then
    PHMMER_HITS=$DEFAULT_PHMMER_HITS
fi
if [ -z "${MAPS:-}" ]; then
    MAPS=$DEFAULT_MAPS
fi
if [ -z "${SAMAP_INPUTS:-}" ]; then
    SAMAP_INPUTS=$DEFAULT_SAMAP_INPUTS
fi
if [ -z "${SAMAP_JOINT:-}" ]; then
    SAMAP_JOINT=$DEFAULT_SAMAP_JOINT
fi
if [ -z "${SAMAP_FINAL:-}" ]; then
    SAMAP_FINAL=$DEFAULT_SAMAP_FINAL
fi

# Notebooks are run by papermill on the kernels in jupyter/kernels/: samap for steps 1 and 2, rsc for the rest.
PY=$NH/.uv/envs/rsc_uv/bin/python
PAPERMILL=$NH/.uv/envs/rsc_uv/bin/papermill
export JUPYTER_PATH=$HERE/jupyter

# hnswlib in samap_harmony_uv needs a newer libstdc++ (GLIBCXX_3.4.32) than the system's.
module load gcc/14.2.0-fasrc01

# Caches go to lab storage, not the home directory.
CACHE=$NH/.cache_scratch
export HDF5_USE_FILE_LOCKING=FALSE PYTHONDONTWRITEBYTECODE=1
export MPLCONFIGDIR=$CACHE/mpl NUMBA_CACHE_DIR=$CACHE/numba XDG_CACHE_HOME=$CACHE/xdg
export IPYTHONDIR=$CACHE/ipython JUPYTER_RUNTIME_DIR=$CACHE/jupyter_runtime
export JUPYTER_DATA_DIR=$CACHE/jupyter_data JUPYTER_CONFIG_DIR=$CACHE/jupyter_config
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR" "$XDG_CACHE_HOME" "$IPYTHONDIR" "$JUPYTER_RUNTIME_DIR"

# run_notebook <notebook> <kernel> <executed copy> [-p NAME VALUE ...]: runs a notebook with papermill. The executed
# copy is the step's log; it is named .partial.ipynb until the last cell has run.
run_notebook() {
    local notebook=$1 kernel=$2 executed=$3
    shift 3
    local partial=${executed%.ipynb}.partial.ipynb
    mkdir -p "$(dirname "$executed")"
    "$PAPERMILL" "$notebook" "$partial" -k "$kernel" --cwd "$(dirname "$notebook")" --log-output "$@"
    mv "$partial" "$executed"
}

# exists <path>: prints "found" or "missing", for --check.
exists() {
    if [ -e "$1" ]; then
        echo found
    else
        echo missing
    fi
}

# pinned <file> <checksum>: for --check, prints whether the file exists and has the pinned checksum.
pinned() {
    if [ ! -e "$1" ]; then
        echo missing
    elif [ "$(sha256sum "$1" | cut -c1-64)" = "$2" ]; then
        echo "found, sha256 as pinned"
    else
        echo "found, but its sha256 is not the one config.sh pins"
    fi
}

# verify_sha256 <file> <checksum>: stops the step when the file's sha256 is not the one config.sh pins.
verify_sha256() {
    local actual
    actual=$(sha256sum "$1" | cut -c1-64)
    if [ "$actual" != "$2" ]; then
        echo "$1 has sha256 $actual, not $2 as config.sh pins; stopping" >&2
        exit 1
    fi
    echo "$1: sha256 as pinned"
}

# skip_if_unused <step> <variable> ...: ends the step successfully when a later step's input was exported, because
# this step's output would then go unused.
skip_if_unused() {
    local step=$1 name default
    shift
    for name in "$@"; do
        default=DEFAULT_$name
        if [ "${!name}" != "${!default}" ]; then
            echo "$name is set to ${!name}, so the output of step $step would go unused; skipped"
            exit 0
        fi
    done
}

# made_from <json> <key> <path>: stops the step when an existing output was made from an input other than <path>.
made_from() {
    local noted
    noted=$("$PY" -c 'import json, sys
value = json.load(open(sys.argv[1]))
for key in sys.argv[2].split("/"):
    value = value[key]
print(value)' "$1" "$2")
    if [ "$(realpath -m "$noted")" != "$(realpath -m "$3")" ]; then
        echo "$(dirname "$1") was made from $noted, not $3; use a new OUT for a run with other inputs" >&2
        exit 1
    fi
}
