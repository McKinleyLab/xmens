#!/bin/bash
#SBATCH -p eddy
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -t 3-00:00
#SBATCH --mem=256gb
#SBATCH -o mm_to_hs_blast_%j.out
#SBATCH -e mm_to_hs_blast_%j.err

# Working directory
WORKDIR="/n/eddy_lab/users/nhilgert/decidualization/code/python/samap_blast/"
mkdir -p $WORKDIR
cd $WORKDIR

# Run tblastx from mouse to human
bash ./map_genes_blast.sh --qname mm --sname hs --threads $SLURM_CPUS_PER_TASK
