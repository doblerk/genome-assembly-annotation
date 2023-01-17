#!/usr/bin/env bash

#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=16
#SBATCH --job-name=dating
#SBATCH --output=/data/users/kdobler/assembly_annotation_course/outputs/output_dating_canu_%j.o
#SBATCH --error=/data/users/kdobler/assembly_annotation_course/outputs/error_dating_canu_%j.e
#SBATCH --partition=pall

# Load required softwares
module load UHTS/Analysis/BEDTools/2.29.2;
module add Emboss/EMBOSS/6.6.0;

# Set the working directory
WORKDIR=/data/users/kdobler/assembly_annotation_course
INPUTDIR=${WORKDIR}/edta_outputs/canu
OUTPUTDIR=${WORKDIR}/te_dating/canu
PERLSCRIPTS=/data/courses/assembly-annotation-course/CDS_annotation/scripts

cd ${OUTPUTDIR}

# 1. Extract sequences of all intact LTR retrotransposons from the genome
# Create a GFF with only LTR retrotransposons
awk '$3~/retrotransposon/' ${INPUTDIR}/pilon.fasta.mod.EDTA.intact.gff3 > genome.mod.EDTA.intact.gff3_LTR

# Reformat gff
sed -i 's/ID.\+Name=//' genome.mod.EDTA.intact.gff3_LTR
sed -i 's/;.\+//' genome.mod.EDTA.intact.gff3_LTR
sed -i 's/_pi\s/_pilon\t/g' genome.mod.EDTA.intact.gff3_LTR
sed -i 's/_pil\s/_pilon\t/g' genome.mod.EDTA.intact.gff3_LTR
sed -i 's/_pilo\s/_pilon\t/g' genome.mod.EDTA.intact.gff3_LTR

awk '{print($1,$2,$9,$4,$5,$6,$7,$8,$3)}' genome.mod.EDTA.intact.gff3_LTR | sed 's/ /\t/g' > genome.mod.EDTA.intact.gff3_LTR_ref

# Extract sequences with bedtools getfasta
bedtools getfasta -fi ${WORKDIR}/polished_assemblies/canu/pilon.fasta -bed genome.mod.EDTA.intact.gff3_LTR_ref -name > genome.mod.EDTA.intact.LTR.fa

# 2. Rename LTR retrotransposons
# Replace : with _
sed -i 's/:/_/g' genome.mod.EDTA.intact.LTR.fa

# Put the abbreviation of the Arabidopsis accession in front of the TE names
sed -i 's/>/>Ler_/' genome.mod.EDTA.intact.LTR.fa

# 3. Copy fasta file containing the intact LTR retrotransposons into a new directory
mkdir intact_ltr
cp genome.mod.EDTA.intact.LTR.fa intact_ltr/
cd intact_ltr

# 4. Use split_flat to split the fasta file into multiple ones
${PERLSCRIPTS}/split_flat genome.mod.EDTA.intact.LTR.fa

# 5. Use the perl script LTR to align the two LTR sequences of each TE with Water
${PERLSCRIPTS}/LTR Ler_ N

# 6. use the perl script date_pair to estimate the divergence time between the LTR sequences of each TE
# Use the perl script date_pair to estimate the divergence time between the LTR sequences of each TE 
${PERLSCRIPTS}/date_pair > results.txt
