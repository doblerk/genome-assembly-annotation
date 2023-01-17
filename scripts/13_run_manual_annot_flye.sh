#!/usr/bin/env bash

#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=16
#SBATCH --job-name=manual_annot
#SBATCH --output=/data/users/kdobler/assembly_annotation_course/outputs/output_manual_annot_%j.o
#SBATCH --error=/data/users/kdobler/assembly_annotation_course/outputs/error_manual_annot_%j.e
#SBATCH --partition=pcourseassembly

# Load required software
module load UHTS/Analysis/SeqKit/0.13.2;
module load Blast/ncbi-blast/2.9.0+;
module load SequenceAnalysis/MultipleSequenceAlignment/clustalw2/2.1;
module load Emboss/EMBOSS/6.6.0;

# Set the working directory
WORKDIR=/data/users/kdobler/assembly_annotation_course
OUTPUTDIR=${WORKDIR}/te_manual_annot/flye

cd ${OUTPUTDIR}

# 1.
seqkit fx2tab -l -n ${WORKDIR}/polished_assemblies/flye/pilon.fasta > ID.txt
sort -n -r -k 2 ID.txt > ID_sorted.txt
awk 'NR==1{print $1}' ID_sorted.txt > ID_max.txt

seqkit grep \
    -n -f ID_max.txt \
    ${WORKDIR}/polished_assemblies/flye/pilon.fasta \
    -o ${OUTPUTDIR}/contig.fasta

# 2.
seqkit sliding \
    -s 500 \
    -W 500 \
    -g contig.fasta \
    -o contig_windows.fasta

# 3.
mkdir blast_db

makeblastdb \
    -in contig_windows.fasta  \
    -dbtype nucl \
    -out blast_db/contig_windows

blastn \
    -query contig_windows.fasta \
    -db blast_db/contig_windows \
    -num_threads 10 \
    -outfmt 6 \
    -perc_identity 80 \
    -max_hsps 1 \
    > contig_windows.blastn

# 4.
awk '{print $1}' contig_windows.blastn > blast_ID.txt
sort blast_ID.txt > blast_ID_sorted.txt
uniq -c blast_ID_sorted.txt > blast_ID_cnt.txt
sort -n -r -k 1 blast_ID_cnt.txt > blast_ID_cnt_sorted.txt
awk '1<=NR && NR<=50{print $2}' blast_ID_cnt_sorted.txt > blast_50max.txt

seqkit \
    grep -n -f blast_50max.txt \
    contig_windows.fasta\
    -o contig_windows_TOP50.fa

# 5.
# do it locally with dotter

# 6.
makeblastdb \
    -in ${WORKDIR}/polished_assemblies/flye/pilon.fasta \
    -dbtype nucl \
    -out blast_db/genome

# Just chose the second contig of the file
awk -v seq='scaffold_453_pilon_sliding:2366501-2367000' \
 -v RS='>' '$1 == seq {print RS $0}' contig_windows_TOP50.fa > window.fa

blastn \
    -query window.fa \
    -db blast_db/genome \
    -outfmt 6 \
    -num_threads 10 \
    -perc_identity 80 \
    -qcov_hsp_perc 80 \
    -max_hsps 1 \
    > window.blastn

# 7.
awk '$10 < $9 {print($2"\t"$10-1"\t"$9)}' window.blastn > window.bed
awk '$10 > $9 {print($2"\t"$9-1"\t"$10)}' window.blastn >> window.bed

# 8.
seqkit subseq \
    --bed window.bed ${WORKDIR}/polished_assemblies/flye/pilon.fasta \
    -u 2000 \
    -d 2000 \
    > window.bed.fa

# 9.
clustalw2 -INFILE=window.bed.fa -OUTFILE=window.bed.aln

# 10.
# visualize it

# 11.
cons \
    -sequence window.bed.aln \
    -outseq TE1_family.cons \
    -name TE1_cons

# 12.
# blast it
makeblastdb \
    -in ../Brassicaceae_repbase_all_march2019.fasta \
    -dbtype nucl \
    -out blast_db/Brassicaceae

blastn \
    -query TE1_family.cons \
    -db blast_db/Brassicaceae \
    -outfmt 6 \
    -num_threads 10 \
    -max_hsps 1 \
    > cons_vs_genome.blastn
