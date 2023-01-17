#!/usr/bin/env bash

#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=16
#SBATCH --job-name=comparative
#SBATCH --output=/data/users/kdobler/assembly_annotation_course/outputs/output_canu_ref_comparing_%j.o
#SBATCH --error=/data/users/kdobler/assembly_annotation_course/outputs/error_canu_ref_comparing_%j.e
#SBATCH --partition=pcourseassembly

# Load the required package
module add UHTS/Analysis/mummer/4.0.0beta1
export PATH=/software/bin:$PATH

working_dir=/data/users/kdobler/assembly_annotation_course
dir_ref=/data/courses/assembly-annotation-course/references

cd ${working_dir}/compared_genomes/canu_vs_ref

nucmer \
 --prefix canu_vs_ref \
 --breaklen 1000 \
 --mincluster 1000 \
 ${dir_ref}/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa \
 ${working_dir}/polished_assemblies/canu/pilon.fasta

mummerplot \
 -R ${dir_ref}/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa \
 -Q ${working_dir}/polished_assemblies/canu/pilon.fasta \
 --filter \
 -t png \
 --large \
 --layout \
 --prefix canu_vs_ref \
 canu_vs_ref.delta
