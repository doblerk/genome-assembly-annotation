#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=12G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=canu_polishing
#SBATCH --output=/data/users/kdobler/assembly_annotation_course/outputs/output_canu_polishing_%j.o
#SBATCH --error=/data/users/kdobler/assembly_annotation_course/outputs/error_canu_polishing_%j.e
#SBATCH --partition=pcourseassembly

# Load the package
module add UHTS/Aligner/bowtie2/2.3.4.1
module load UHTS/Analysis/samtools/1.10
module load UHTS/Analysis/pilon/1.22

# Define paths to directories
dir_input_reads=/data/users/kdobler/assembly_annotation_course/participant_2/Illumina
dir_input_reference=/data/users/kdobler/assembly_annotation_course/assemblies/canu
dir_output=/data/users/kdobler/assembly_annotation_course/polished_assemblies/canu

# Go to output directory
cd ${dir_output}/

# Create index of the Illumina assemblies
bowtie2-build -f --threads 4 ${dir_input_reference}/canu_pacbio.contigs.fasta canu_index

# Run the alignment
bowtie2 -q --sensitive-local -p 4 -x canu_index \
 -1 ${dir_input_reads}/ERR3624574_1.fastq.gz -2 ${dir_input_reads}/ERR3624574_2.fastq.gz -S canu.sam

# Convert output SAM files to BAM
samtools sort -T $SCRATCH -@ $SLURM_CPUS_PER_TASK canu.sam -o canu_sorted.sam
samtools view -bS canu_sorted.sam > canu.bam
samtools index canu.bam

# Run pilon
java -Xmx45g -jar /mnt/software/UHTS/Analysis/pilon/1.22/bin/pilon-1.22.jar \
 --genome ${dir_input_reference}/canu_pacbio.contigs.fasta --bam canu.bam
