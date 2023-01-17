#!/usr/bin/env bash

#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=16
#SBATCH --job-name=comparative
#SBATCH --output=/data/users/kdobler/assembly_annotation_course/outputs/output_canuref_flye_comparing_%j.o
#SBATCH --error=/data/users/kdobler/assembly_annotation_course/outputs/error_canuref_flye_comparing_%j.e
#SBATCH --partition=pcourseassembly

# Load the required package
module add UHTS/Analysis/mummer/4.0.0beta1
export PATH=/software/bin:$PATH

working_dir=/data/users/kdobler/assembly_annotation_course

cd ${working_dir}/compared_genomes/canuref_vs_flye

nucmer \
 --prefix canuref_vs_flye \
 --breaklen 1000 \
 --mincluster 1000 \
 ${working_dir}/polished_assemblies/canu/pilon.fasta \
 ${working_dir}/polished_assemblies/flye/pilon.fasta

mummerplot \
 -R ${working_dir}/polished_assemblies/canu/pilon.fasta \
 -Q ${working_dir}/polished_assemblies/flye/pilon.fasta \
 --filter \
 -t png \
 --large \
 --layout \
 --prefix canuref_vs_flye \
 canuref_vs_flye.delta
