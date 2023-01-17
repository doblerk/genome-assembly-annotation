#!/usr/bin/env bash

#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=16
#SBATCH --job-name=edta
#SBATCH --output=/data/users/kdobler/assembly_annotation_course/outputs/output_edta_%j.o
#SBATCH --error=/data/users/kdobler/assembly_annotation_course/outputs/error_edta_%j.e
#SBATCH --partition=pcourseassembly

WORKDIR=/data/users/kdobler/assembly_annotation_course
CONTAINER=/data/courses/assembly-annotation-course/containers2

cd ${WORKDIR}/edta_outputs/canu

singularity exec \
--bind ${CONTAINER} \
--bind ${WORKDIR} \
$CONTAINER/EDTA_v1.9.6.sif \
 EDTA.pl \
 --genome ${WORKDIR}/polished_assemblies/canu/pilon.fasta \
 --species others \
 --step all \
 --cds ${WORKDIR}/edta_outputs/TAIR10_cds_20110103_representative_gene_model \
 --anno 1 \
 --threads 16
