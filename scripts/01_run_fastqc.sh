#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000M
#SBATCH --time=01:00:00
#SBATCH --job-name=fastqc
#SBATCH --output=/data/users/kdobler/assembly_annotation_course/outputs/output_fastqc_%j.o
#SBATCH --error=/data/users/kdobler/assembly_annotation_course/outputs/error_fastqc_%j.e
#SBATCH --partition=pcourseassembly

module load UHTS/Quality_control/fastqc/0.11.9;

# choose the correct location of the output

fastqc -o ./read_QC/X --extract $1
