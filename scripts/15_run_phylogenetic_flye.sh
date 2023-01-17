#!/usr/bin/env bash

#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=16
#SBATCH --job-name=phylo
#SBATCH --output=/data/users/kdobler/assembly_annotation_course/outputs/output_phylo_flye_%j.o
#SBATCH --error=/data/users/kdobler/assembly_annotation_course/outputs/error_phylo_flye_%j.e
#SBATCH --partition=pall

# Load required software
module load UHTS/Analysis/SeqKit/0.13.2;
module load SequenceAnalysis/MultipleSequenceAlignment/clustal-omega/1.2.4;
module load Phylogeny/FastTree/2.1.10;

# Set the working directory
WORKDIR=/data/users/kdobler/assembly_annotation_course
INPUTDIR=${WORKDIR}/tesorter_outputs/flye
OUTPUTDIR=${WORKDIR}/phylogenetic_analysis/flye

cd ${OUTPUTDIR}

# Extract RT protein sequences ID for copia
more ${INPUTDIR}/Ler_flye_TEs.rexdb-plant.cls.pep | grep Ty1-RT | sed 's/>//'| sed 's/ .\+//' > copia/ID_ar.txt
more ${INPUTDIR}/brassicaceae.cls.pep | grep Ty1-RT | sed 's/>//' | sed 's/ .\+//' > copia/ID_brass.txt

# Extract RT protein sequences ID for gypsy
more ${INPUTDIR}/Ler_flye_TEs.rexdb-plant.cls.pep | grep Ty3-RT | sed 's/>//'| sed 's/ .\+//' > gypsy/ID_ar.txt
more ${INPUTDIR}/brassicaceae.cls.pep | grep Ty3-RT | sed 's/>//' | sed 's/ .\+//' > gypsy/ID_brass.txt

# Extract corresponding sequences from the pattern file for copia
seqkit grep \
    -f copia/ID_ar.txt \
    ${INPUTDIR}/Ler_flye_TEs.rexdb-plant.cls.pep \
    -o copia/ar.fa

seqkit grep \
    -f copia/ID_brass.txt \
    ${INPUTDIR}/brassicaceae.cls.pep \
    -o copia/brass.fa

# Extract corresponding sequences from the pattern file for gypsy
seqkit grep \
    -f gypsy/ID_ar.txt \
    ${INPUTDIR}/Ler_flye_TEs.rexdb-plant.cls.pep \
    -o gypsy/ar.fa

seqkit grep \
    -f gypsy/ID_brass.txt \
    ${INPUTDIR}/brassicaceae.cls.pep \
    -o gypsy/brass.fa

# Concatenate the files
cat copia/ar.fa copia/brass.fa > copia/copia.fa
cat gypsy/ar.fa gypsy/brass.fa > gypsy/gypsy.fa

# Shorten identifiers of RT sequences and replace : with _
sed -i 's/#.\+//' copia/copia.fa
sed -i 's/:/_/g' copia/copia.fa
sed -i 's/__/_/g' copia/copia.fa

sed -i 's/#.\+//' gypsy/gypsy.fa
sed -i 's/:/_/g' gypsy/gypsy.fa
sed -i 's/__/_/g' gypsy/gypsy.fa

# Align the sequences with clustal omega
clustalo -i copia/copia.fa -o copia/aligned_copia.fa
clustalo -i gypsy/gypsy.fa -o gypsy/aligned_gypsy.fa

# Infer approximately-maximum-likelihood phylogenetic tree
FastTree -out copia/copia.tre copia/aligned_copia.fa
FastTree -out gypsy/gypsy.tre gypsy/aligned_gypsy.fa

# Visualize and annotate the tree
# ...
