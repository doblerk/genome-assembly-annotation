#!/usr/bin/env bash

#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=16
#SBATCH --job-name=comparative
#SBATCH --output=/data/users/kdobler/assembly_annotation_course/outputs/output_flyeref_canu_comparing_%j.o
#SBATCH --error=/data/users/kdobler/assembly_annotation_course/outputs/error_flyeref_canu_comparing_%j.e
#SBATCH --partition=pcourseassembly

# Load the required package
module add UHTS/Analysis/mummer/4.0.0beta1
export PATH=/software/bin:$PATH

working_dir=/data/users/kdobler/assembly_annotation_course

cd ${working_dir}/compared_genomes/flyeref_vs_canu

nucmer \
 --prefix flyeref_vs_canu \
 --breaklen 1000 \
 --mincluster 1000 \
 ${working_dir}/polished_assemblies/flye/pilon.fasta \
 ${working_dir}/polished_assemblies/canu/pilon.fasta

mummerplot \
 -R ${working_dir}/polished_assemblies/flye/pilon.fasta \
 -Q ${working_dir}/polished_assemblies/canu/pilon.fasta \
 --filter \
 -t png \
 --large \
 --layout \
 --prefix flyeref_vs_canu \
 flyeref_vs_canu.delta
