#!/bin/bash

# Parse arguments
die() {
    printf '%s\n' "$1" >&2
    exit 1
}

_help() {
    printf "Usage: map_genes_blast.sh [--qname] [--sname] [--threads]\n\
    [--qname]: identifier for the query species (e.g., 'hs')\n\
    [--sname]: identifier for the subject species (e.g., 'mm')\n\
    [--threads]: number of CPU threads to use (default: 8)"
}

if [ $# -eq 0 ]; then
    _help
    exit 1
fi

# Parse command-line arguments
while :; do
    case $1 in
        -h|-?|--help)
            _help
            exit
            ;;
        --qname)
            if [ "$2" ]; then
                qname=$2
                shift
            else
                die 'ERROR: "--qname" requires a non-empty option argument.'
            fi
            ;;
        --sname)
            if [ "$2" ]; then
                sname=$2
                shift
            else
                die 'ERROR: "--sname" requires a non-empty option argument.'
            fi
            ;;
        --threads)
            if [ "$2" ]; then
                threads=$2
                shift
            else
                die 'ERROR: "--threads" requires a non-empty option argument.'
            fi
            ;;
        --)
            shift
            break
            ;;
        -?*)
            printf 'WARN: Unknown option (ignored): %s\n' "$1" >&2
            ;;
        *)
            break
    esac

    shift
done

# Set default number of threads if not provided
if [ -z ${threads+x} ]; then
    threads=8
fi

# Transcriptome paths
#HUMAN_TRANSCRIPTOME="/net/holy-nfsisilon/ifs/rc_labs/eddy_lab/users/nhilgert/reference/human/refdata-gex-GRCh38-2024-A/fasta/Homo_sapiens.GRCh38.cdna.all.fa"
#MOUSE_TRANSCRIPTOME="/net/holy-nfsisilon/ifs/rc_labs/eddy_lab/users/nhilgert/reference/mus/refdata-gex-GRCm39-2024-A/fasta/Mus_musculus.GRCm39.cdna.all.fa"

HUMAN_TRANSCRIPTOME="Homo_sapiens.GRCh38.cdna.all.fa"
MOUSE_TRANSCRIPTOME="Mus_musculus.GRCm39.cdna.all.fa"

# Set transcriptomes based on qname and sname
if [[ $qname == "hs" && $sname == "mm" ]]; then
    query=$HUMAN_TRANSCRIPTOME
    subject=$MOUSE_TRANSCRIPTOME
elif [[ $qname == "mm" && $sname == "hs" ]]; then
    query=$MOUSE_TRANSCRIPTOME
    subject=$HUMAN_TRANSCRIPTOME
else
    die "ERROR: Invalid species identifiers. Supported identifiers are 'hs' (human) and 'mm' (mouse)."
fi

# Output directory
OUTDIR="/n/eddy_lab/users/nhilgert/decidualization/code/python/samap_blast/${qname}_${sname}"
mkdir -p "${OUTDIR}"

# Run TBLASTX
echo "Running tblastx alignment from ${qname} to ${sname}..."
tblastx -query "${query}" -db "${subject}" -outfmt 6 \
    -out "${OUTDIR}/${qname}_to_${sname}.txt" \
    -num_threads ${threads} -max_hsps 1 -evalue 1e-6

echo "tblastx alignment completed. Results saved in ${OUTDIR}/${qname}_to_${sname}.txt"
