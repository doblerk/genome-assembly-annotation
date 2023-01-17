#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=12G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=flye_polishing
#SBATCH --output=/data/users/kdobler/assembly_annotation_course/outputs/output_flye_polishing_%j.o
#SBATCH --error=/data/users/kdobler/assembly_annotation_course/outputs/error_flye_polishing_%j.e
#SBATCH --partition=pcourseassembly

# Load the package
module add UHTS/Aligner/bowtie2/2.3.4.1
module load UHTS/Analysis/samtools/1.10
module load UHTS/Analysis/pilon/1.22

# Define paths to directories
dir_input_reads=/data/users/kdobler/assembly_annotation_course/participant_2/Illumina
dir_input_reference=/data/users/kdobler/assembly_annotation_course/assemblies/flye
dir_output=/data/users/kdobler/assembly_annotation_course/polished_assemblies/flye

# Go to output directory
cd ${dir_output}/

# Create index of the Illumina assemblies
bowtie2-build -f --threads 4 ${dir_input_reference}/assembly.fasta flye_index

# Run the alignment
bowtie2 -q --sensitive-local -p 4 -x flye_index \
 -1 ${dir_input_reads}/ERR3624574_1.fastq.gz -2 ${dir_input_reads}/ERR3624574_2.fastq.gz -S flye.sam

# Convert output SAM files to BAM
samtools sort -T $SCRATCH -@ $SLURM_CPUS_PER_TASK flye.sam -o flye_sorted.sam
samtools view -bS flye_sorted.sam > flye.bam
samtools index flye.bam

# Run pilon
java -Xmx45g -jar /mnt/software/UHTS/Analysis/pilon/1.22/bin/pilon-1.22.jar \
 --genome ${dir_input_reference}/assembly.fasta --bam flye.bam
