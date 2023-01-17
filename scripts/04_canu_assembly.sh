#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=01:00:00
#SBATCH --job-name=canu_assembly
#SBATCH --output=/data/users/kdobler/assembly_annotation_course/outputs/output_canu_%j.o
#SBATCH --error=/data/users/kdobler/assembly_annotation_course/outputs/error_canu_%j.e
#SBATCH --partition=pcourseassembly

# Load the package
module load UHTS/Assembler/canu/2.1.1;

# Define paths to directories 
dir_input=/data/users/kdobler/assembly_annotation_course/participant_2/pacbio/*
dir_output=/data/users/kdobler/assembly_annotation_course/assemblies/canu

# Run the canu command
canu \
    -p canu_pacbio \
    -d ${dir_output} \
    genomeSize=130m \
    -pacbio ${dir_input} \
    gridEngineResourceOption="--cpus-per-task=THREADS --mem-per-cpu=MEMORY" \
    gridOptions="--partition=pcourseassembly --mail-user=kalvin.dobler@students.unibe.ch"
