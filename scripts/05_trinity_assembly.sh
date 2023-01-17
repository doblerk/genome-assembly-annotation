#!/usr/bin/env bash

#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=trinity_assembly
#SBATCH --output=/data/users/kdobler/assembly_annotation_course/outputs/output_trinity_%j.o
#SBATCH --error=/data/users/kdobler/assembly_annotation_course/outputs/error_trinity_%j.e
#SBATCH --partition=pcourseassembly

# Load the package
module load UHTS/Assembler/trinityrnaseq/2.5.1;

# Define paths to directories
dir_input=/data/users/kdobler/assembly_annotation_course/participant_2/RNAseq
dir_output=/data/users/kdobler/assembly_annotation_course/assemblies/trinity

# Run the trinity command
Trinity \
	--seqType fq \
	--left ${dir_input}/SRR1734309_1.fastq.gz \
	--right ${dir_input}/SRR1734309_2.fastq.gz \
	--CPU 12 \
	--max_memory 48G \
	--output ${dir_output}
