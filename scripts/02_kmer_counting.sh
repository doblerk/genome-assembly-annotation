#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=10G
#SBATCH --time=01:00:00
#SBATCH --job-name=kmer_counting
#SBATCH --output=/data/users/kdobler/assembly_annotation_course/outputs/output_kmer_%j.o
#SBATCH --error=/data/users/kdobler/assembly_annotation_course/outputs/error_kmer_%j.e
#SBATCH --partition=pall

# Load required packages
module load UHTS/Analysis/jellyfish/2.3.0;

# Define paths to directories
dir_output=/data/users/kdobler/assembly_annotation_course/ASSEMBLY/read_QC/kmer_counting
dir_input=/data/users/kdobler/assembly_annotation_course/participant_2

# Count kmers for Illumina
jellyfish count -C -m 21 -s 5G -t 4 <(zcat ${dir_input}/Illumina/*.fastq.gz) -o ${dir_output}/Illumina_reads.jf
# Export kmer count histogram for Illumina
jellyfish histo -t 4 ${dir_output}/Illumina_reads.jf > ${dir_output}/Illumina_reads.histo

# Count kmers for PacBio
jellyfish count -C -m 21 -s 5G -t 4 <(zcat ${dir_input}/pacbio/*.fastq.gz) -o ${dir_output}/PacBio_reads.jf
# Export kmer count histogram for PacBio
jellyfish histo -t 4 ${dir_output}/PacBio_reads.jf > ${dir_output}/PacBio_reads.histo

# Count kmers for RNAseq
jellyfish count -C -m 21 -s 5G -t 4 <(zcat ${dir_input}/RNAseq/*.fastq.gz) -o ${dir_output}/RNAseq_reads.jf
# Export kmer count histogram for RNAseq
jellyfish histo -t 4 ${dir_output}/RNAseq_reads.jf > ${dir_output}/RNAseq_reads.histo
