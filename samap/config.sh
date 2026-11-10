# Paths and settings shared by every step; each step script sources this file after setting HERE, the pipeline folder.
#
# OUT, the run folder, has no default, so a run never writes into an earlier one by accident. Use a new one for each
# run: a step whose output exists in OUT is skipped, or stops if that output was made from other inputs.
#   export OUT=/n/netscratch/eddy_lab/Everyone/nhilgert/$(date +%Y_%m_%d)_xmens_samap_full

OUT=${OUT:?export OUT, the run folder, before running a step (see README.md)}
NH=/n/lab_storage/eddy_lab/users/nhilgert
NS=/n/netscratch/eddy_lab/Everyone/nhilgert

# The deposited fibroblast/decidual objects (steps 1, 3 and 4): each species' counts, labels, positions and pseudotime
# of record.
FINAL_OBJECTS=$NS/final_objects
MOUSE_H5AD=$FINAL_OBJECTS/Cevrim_XMens_Only_Fib_Dec_Cells.h5ad
HUMAN_H5AD=$FINAL_OBJECTS/Cevrim_Human_Only_Fib_Dec_Cells.h5ad

# Step 0: the human and mouse proteomes of Ensembl release 113, byte-identical to the decompressed
#   https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz
#   https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/pep/Mus_musculus.GRCm39.pep.all.fa.gz
# (checked 7 Nov 2026); and phmmer of HMMER 3.4, as an absolute path (batch jobs do not read an interactive PATH).
PROTEOMES=$NH/data/proteomes/fasta
HUMAN_PROTEOME=$PROTEOMES/Homo_sapiens.GRCh38.pep.all.fa
MOUSE_PROTEOME=$PROTEOMES/Mus_musculus.GRCm39.pep.all.fa
PHMMER=$NH/software/hmmer-3.4/src/phmmer

# Checksums of the four input files above. Each step that reads one checks it first and stops if it differs, so a
# replaced file (netscratch is purgeable, and the final objects were rewritten once, on 2 Nov 2026) cannot flow
# through unnoticed. To run on other files, change their paths and checksums here.
MOUSE_H5AD_SHA256=e3e7ad3df4fd9ca58888a1e703ab70d27847be90d33fc46eae6e6860756344e7
HUMAN_H5AD_SHA256=e74d038e446049fe7bea8d18a5612379ebfbc22f1084eb104f75eb0695998ada
HUMAN_PROTEOME_SHA256=f38587d8e87564a08f99378c3416b47a866cd0aff5ed0bdeaeab0ca2d2a8ea3b
MOUSE_PROTEOME_SHA256=54cb3f64f6dd1d0a398713016e866c355cff8e21968f8d18e0d60f7b498712d6

# Where steps 0b-4 read their inputs: by default the outputs of the step before, in OUT. To start from earlier
# results, export one of these before submitting. Steps whose output would then go unused skip themselves, so the
# whole chain can still be queued. For example, the SAMap fit of record (steps 0-2 skip):
#   export SAMAP_JOINT=$NH/code/samap-scvi/full_fibdec/runs_raw/vanilla_raw/samap_transferred.h5ad
# or the homology maps of record (step 0 skips):
#   export MAPS=$NH/decidualization/revisions/samap_integrate/homology/hsmm_hmmer
# (Plain ifs, as in preprocess/.)
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

# Every step but the phmmer searches (0a) is a notebook, run by papermill (in rsc_uv) on one of the two kernels in
# jupyter/kernels/: samap (samap_harmony_uv, sc-samap 3.0.1) for steps 1 and 2, rsc (rsc_uv) for steps 0b and 3-5.
# PY runs compare_with_record.
PY=$NH/.uv/envs/rsc_uv/bin/python
PAPERMILL=$NH/.uv/envs/rsc_uv/bin/papermill
export JUPYTER_PATH=$HERE/jupyter

# hnswlib in samap_harmony_uv needs a newer libstdc++ (GLIBCXX_3.4.32) than the system's.
module load gcc/14.2.0-fasrc01

# Caches go to lab storage, never to the home directory, and no bytecode is written into this folder.
CACHE=$NH/.cache_scratch
export HDF5_USE_FILE_LOCKING=FALSE PYTHONDONTWRITEBYTECODE=1
export MPLCONFIGDIR=$CACHE/mpl NUMBA_CACHE_DIR=$CACHE/numba XDG_CACHE_HOME=$CACHE/xdg
export IPYTHONDIR=$CACHE/ipython JUPYTER_RUNTIME_DIR=$CACHE/jupyter_runtime
export JUPYTER_DATA_DIR=$CACHE/jupyter_data JUPYTER_CONFIG_DIR=$CACHE/jupyter_config
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR" "$XDG_CACHE_HOME" "$IPYTHONDIR" "$JUPYTER_RUNTIME_DIR"

# run_notebook <notebook> <kernel> <executed copy> [-p NAME VALUE ...]: runs a notebook with papermill from the
# notebook's folder, with the parameters given. The executed copy, holding everything the notebook printed, is the
# step's log; it is written as <name>.partial.ipynb and renamed when the last cell has run, so after a failure only
# the .partial copy exists, ending at the cell that failed.
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

# pinned <file> <checksum>: prints, for --check, whether the file exists and has the checksum config.sh pins.
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

# skip_if_unused <step> <variable> ...: the variables are inputs of later steps. When one of them was exported (it
# differs from its default), this step's output would go unused, so the step ends here, successfully, and the steps
# queued after it with afterok still start.
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

# made_from <json> <key> <path>: an existing output is reused only if it was made from the input in use. <json> is the
# file the step wrote beside its output, <key> the slash-separated key under which it noted the input's path. Stops
# the step when the noted path is another file.
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
