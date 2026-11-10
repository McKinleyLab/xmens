# Paths and settings shared by every step; each step script sources this file.
#
# OUT, the run folder, has no default, so a run never writes into an earlier one by accident. Set it once per run:
#   export OUT=/n/netscratch/eddy_lab/Everyone/nhilgert/$(date +%Y_%m_%d)_xmens_scenic

OUT=${OUT:?export OUT, the run folder, before running a step (see README.md)}
NH=/n/lab_storage/eddy_lab/users/nhilgert

# Step 1 input: the fibroblast/decidual objects SCENIC was run on, release 2026_09_15 (byte-identical to the
# netscratch 2026_08_28 objects the run of record read). Their nuclei, genes and raw counts (layers['counts']) are the
# same as the final objects'. The mouse final object carries a different scVI latent, and step 2 draws each run's
# sketch on obsm['X_scVI'], so these objects, not the final ones, reproduce the run of record.
SCENIC_INPUT_MOUSE=${SCENIC_INPUT_MOUSE:-$NH/xmens_final/data/deliverables/2026_09_15/Cevrim_XMens_Only_Fib_Dec_Cells.h5ad}
SCENIC_INPUT_HUMAN=${SCENIC_INPUT_HUMAN:-$NH/xmens_final/data/deliverables/2026_09_15/Cevrim_Human_Only_Fib_Dec_Cells.h5ad}
# The layer of the input objects that holds raw counts; the objects of record use "counts", Cagri's Nov 2026 objects
# "raw_counts". Export RAW_LAYER before submitting step 1 to change it.
RAW_LAYER=${RAW_LAYER:-counts}

# Step 5 input: the final objects. obs['pseudotime'] is the pseudotime of record (mouse: pt_harmonized_rank, NaN for
# the 907 apoptotic decidual nuclei; human: the 30 Oct refit on the fibroblast/decidual latent).
FINAL_MOUSE=/n/netscratch/eddy_lab/Everyone/nhilgert/final_objects/Cevrim_XMens_Only_Fib_Dec_Cells.h5ad
FINAL_HUMAN=/n/netscratch/eddy_lab/Everyone/nhilgert/final_objects/Cevrim_Human_Only_Fib_Dec_Cells.h5ad

# cisTarget databases (motif collection v10, clustered; mm10 and hg38 gene rankings 500 bp up/100 bp down and 10 kb
# around the TSS) and the TF lists allTFs_mm.txt and allTFs_hg38.txt
DB_DIR=$NH/decidualization/revisions/scenic/databases

# Ensembl release 116 one-to-one human-mouse orthologs: the pair table (step 5; protein-coding pairs with unique
# symbols, 16,275; a copy of xmens_final/results/mi_dpt/overlap_of_record/ortholog_pairs_one2one.csv) and the Compara
# component table it was made from (step 6)
ORTHOLOG_PAIRS=$HERE/resources/ortholog_pairs_one2one.csv
ORTHOLOG_COMPONENTS=$NH/xmens_final/resources/orthology/ensembl_human_mouse_components.tsv

# pySCENIC 0.12.1 (steps 1-4, and the regulon tables of step 6). Its pyscenic and arboreto scripts are run with this
# Python, because their own first lines point at the retired /n/eddy_lab path.
SCENIC_ENV=$NH/.conda/envs/pyscenic
PY_SCENIC=$SCENIC_ENV/bin/python
PYSCENIC_CLI=$SCENIC_ENV/bin/pyscenic
ARBORETO=$SCENIC_ENV/bin/arboreto_with_multiprocessing.py
# rsc_uv (scanpy 1.11): steps 5 and 6
PY=$NH/.uv/envs/rsc_uv/bin/python

# Settings of the run of record
SKETCH=15k                # every run trains GRNBoost2 on its own 15,000-nucleus geometric sketch
SKETCH_N=15000
TAG=sketch15k             # the folder tag build_consensus_aucell.py derives from --sketch 15k
MIN_OCCURRENCE=6          # consensus: a target is kept if its link to the TF recurs in at least 6 of the 100 runs
CONS6=consensus_${TAG}_weighted_soft_occ6   # consensus regulons of record (floor 6)
CONS2=consensus_${TAG}_weighted_soft        # the same at a floor of 2 runs; Fig. S28 panels C and D read its links
AUCELL=${CONS6}_pyscenic_cli                # regulon activity of record: pyscenic aucell --weights
# The 20 canonical decidualization TFs flagged in the Fig. 5E table, and the 13 named on the panel
LABELLED=CEBPB,CEBPD,CREB1,ETS1,FOS,FOSL2,FOXO1,FOXO3,GATA2,HOXA10,IRF1,JUN,JUNB,KLF4,KLF6,NFKB1,NR2F2,STAT1,TFAP2C,WT1
NAMED=JUN,WT1,FOSL2,GATA2,FOS,JUNB,KLF6,CEBPB,BPTF,NFIA,IRF1,ETS1,FOXO1
RHO_THRESHOLD=0.2         # |rho| needed in both species for a pair to be coloured in Fig. 5E

# The run of record (September-October 2026): compare_with_record.sbatch compares against it, and its runs can stand
# in for steps 2 and 3 (below).
RECORD=$NH/decidualization/revisions/scenic_v2_2026_08_28

# Where steps 3, 4 and 6 read the per-run results; each is a folder holding scenic_seeds/<species>/. Both default to
# the run folder. To start from the run of record, export either before submitting (see README.md):
#   GRN_FROM       step 3 reads looms_sketch15k/ and adjacencies_sketch15k/ from here
#   REGULONS_FROM  steps 4 and 6 read regulons_sketch15k/ (the regulons of the 100 runs) from here
GRN_FROM=${GRN_FROM:-$OUT}
REGULONS_FROM=${REGULONS_FROM:-$OUT}

# Caches and temporary files go to lab storage, never to the home directory.
CACHE=$NH/.cache_scratch
export HDF5_USE_FILE_LOCKING=FALSE
export MPLCONFIGDIR=$CACHE/mpl
export NUMBA_CACHE_DIR=$CACHE/numba
export XDG_CACHE_HOME=$CACHE/xdg
export TMPDIR=$NH/tmp
export PYTHONDONTWRITEBYTECODE=1

# species_paths <species>: sets SP, SP_DIR (the species' folder in the run folder), TF_LIST and SCENIC_INPUT.
species_paths() {
    SP=$1
    case "$SP" in
        mouse)
            TF_LIST=$DB_DIR/allTFs_mm.txt
            SCENIC_INPUT=$SCENIC_INPUT_MOUSE
            ;;
        human)
            TF_LIST=$DB_DIR/allTFs_hg38.txt
            SCENIC_INPUT=$SCENIC_INPUT_HUMAN
            ;;
        *)
            echo "species must be mouse or human, not '$SP'" >&2
            return 1
            ;;
    esac
    SP_DIR=$OUT/scenic_seeds/$SP
}

# task_species <task>: the species of array task 1 or 2 (steps 1 and 4).
task_species() {
    case "$1" in
        1) echo mouse ;;
        2) echo human ;;
        *) echo "array task $1: steps 1 and 4 have tasks 1 (mouse) and 2 (human)" >&2; return 1 ;;
    esac
}

# field <row> <column>: one value from runs.tsv; row 1 is the first run, as is array task 1.
field() {
    awk -F'\t' -v row="$1" -v col="$2" '
        NR == 1 { for (i = 1; i <= NF; i++) if ($i == col) c = i }
        NR == row + 1 { print $c }' "$HERE/runs.tsv"
}

# exists <path>: prints "found" or "missing", for --check.
exists() {
    if [ -e "$1" ]; then
        echo found
    else
        echo missing
    fi
}

mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR" "$XDG_CACHE_HOME" "$TMPDIR"
N_RUNS=$(($(wc -l < "$HERE/runs.tsv") - 1))
