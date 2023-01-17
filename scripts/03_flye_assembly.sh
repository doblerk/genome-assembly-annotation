#!/usr/bin/env bash

#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --time=06:00:00
#SBATCH --job-name=flye_assembly
#SBATCH --output=/data/users/kdobler/assembly_annotation_course/outputs/output_flye_%j.o
#SBATCH --error=/data/users/kdobler/assembly_annotation_course/outputs/error_flye_%j.e
#SBATCH --partition=pall

# Load the package
module load UHTS/Assembler/flye/2.8.3;

# Define paths to directories
dir_input=/data/users/kdobler/assembly_annotation_course/participant_2/pacbio/
dir_output=/data/users/kdobler/assembly_annotation_course/assemblies/flye

# Run the assembly package for each file
flye \
    --pacbio-raw ${dir_input}/* \
    --out-dir ${dir_output} \
    --threads 16 \
    --genome-size 130M
