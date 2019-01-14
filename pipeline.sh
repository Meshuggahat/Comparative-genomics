#!/bin/bash

fasta_pep=SEQS/Solanum_tuberosum.SolTub_3.0.pep.all.fa
fasta_cds=SEQS/Solanum_tuberosum.SolTub_3.0.cds.all.fa
filtered=SEQS/filtered.fasta
blastp=RESULTS/blastp_result.txt
tmp_pep=tmp_pep.txt
tmp_cds=tmp_cds.txt


# ============================================= #
# Keep only the longest isoform for each gene
# ============================================= #

python3 ./SCRIPTS/filter_isoforms.py SEQS/Solanum_tuberosum.SolTub_3.0.pep.all.fa


# ============================================================== #
# Make database with fasta file of sequences of Vitis Vinifera
# ============================================================== #

makeblastdb -in SEQS/filtered.fasta -dbtype prot


# ========================================================================================= #
# Run blastp on previously created database and fasta file of sequences of Vitis Vinifera
# -evalue : threshold for saving hits
# -matrix : name of the scoring matrix to use (default BLOSUM62)
# -num_threads : nb of threads to use during the search
# -out : name of the file to write the output
# -outfmt : format of the output
# ========================================================================================= #

blastp -query SEQS/filtered.fasta \
	-db SEQS/filtered.fasta \
	-out RESULTS/blastp_result.txt \
	-outfmt "6 qseqid sseqid pident qlen length mismatch gapopen evalue bitscore qstart qend" \
	-num_threads 8


# ========================== #
# Ka/Ks ratio calculation
# ========================== #

python3 ./SCRIPTS/ka_ks.py \
	SEQS/Solanum_tuberosum.SolTub_3.0.pep.all.fa \
	SEQS/Solanum_tuberosum.SolTub_3.0.cds.all.fa  \
	RESULTS/MCL_random_families.tabular
	
python3 ./SCRIPTS/ka_ks.py \
	SEQS/Solanum_tuberosum.SolTub_3.0.pep.all.fa \
	SEQS/Solanum_tuberosum.SolTub_3.0.cds.all.fa  \
	RESULTS/SL_random_families.tabular
