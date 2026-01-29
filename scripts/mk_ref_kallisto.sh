#!/bin/bash
#SBATCH -A naiss2025-22-779
#SBATCH -p shared
#SBATCH -t 0-1:00:00
#SBATCH -n 1
#SBATCH -J mk_ref_kallisto
#SBATCH -e logs/mk_ref_kallisto.err
#SBATCH -o logs/mk_ref_kallisto.o
#SBATCH --mail-user emilio.skarwan@scilifelab.se
#SBATCH --mail-type=ALL


#cd /cfs/klemming/home/s/skarwan/EmilioTemp/rat/rat_bulk/local/reference
#wget https://ftp.ensembl.org/pub/release-109/fasta/rattus_norvegicus/cdna/Rattus_norvegicus.mRatBN7.2.cdna.all.fa.gz

module load PDC/24.11
module load kallisto/0.51.1-cpeCray-24.11

cd /cfs/klemming/home/s/skarwan/EmilioTemp/rat/rat_bulk/local/reference

# For 109
fasta_file=Rattus_norvegicus.mRatBN7.2.cdna.all.fa.gz
index_output_name=Rattus_norvegicus.mRatBN7.2.cdna.all.fa.index

kallisto index -i $index_output_name $fasta_file

# New assembly 115
#wget https://ftp.ensembl.org/pub/release-115/fasta/rattus_norvegicus/cdna/Rattus_norvegicus.GRCr8.cdna.all.fa.gz
fasta_file=Rattus_norvegicus.GRCr8.cdna.all.fa.gz
index_output_name=Rattus_norvegicus.GRCr8.cdna.all.fa.index

kallisto index -i $index_output_name $fasta_file