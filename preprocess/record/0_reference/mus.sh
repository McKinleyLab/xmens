#!/bin/bash
#SBATCH -c 1                # Number of cores (-c)
#SBATCH -t 5-00:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p eddy # Partition to submit to
#SBATCH --mem 100gb           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o mus_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e mus_%j.err  # File to which STDERR will be written, %j inserts jobid

# run code
/n/eddy_lab/users/nhilgert/cellranger-7.2.0/cellranger mkref --genome=transgenic_mus --fasta=/n/eddy_lab/users/nhilgert/reference/transgenic_mus_premkref/fasta/genome.fa --genes=/n/eddy_lab/users/nhilgert/reference/transgenic_mus_premkref/genes/genes.gtf
