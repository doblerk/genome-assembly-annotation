#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=merqury
#SBATCH --output=/data/users/kdobler/assembly_annotation_course/outputs/output_merqury_polished_%j.o
#SBATCH --error=/data/users/kdobler/assembly_annotation_course/outputs/error_merqury_polished_%j.e
#SBATCH --partition=pcourseassembly

# Define the path variable
WORKDIR=/data/users/kdobler/assembly_annotation_course
dir_data=/data/courses/assembly-annotation-course
dir_input_reads=${WORKDIR}/participant_2/Illumina
dir_output_meryl=${WORKDIR}/evaluation_assemblies/merqury/meryl

cd ${dir_output_meryl}

# Step 1 - prepare mercury dbs
singularity exec \
--bind ${WORKDIR} \
--bind ${dir_data} \
/software/singularity/containers/Merqury-1.3-1.ubuntu20.sif \
meryl k=21 count output read1.meryl ${dir_input_reads}/ERR3624574_1.fastq.gz

singularity exec \
--bind ${WORKDIR} \
--bind ${dir_data} \
/software/singularity/containers/Merqury-1.3-1.ubuntu20.sif \
meryl k=21 count output read2.meryl ${dir_input_reads}/ERR3624574_2.fastq.gz

singularity exec \
--bind ${WORKDIR} \
--bind ${dir_data} \
/software/singularity/containers/Merqury-1.3-1.ubuntu20.sif \
meryl union-sum output Illumina.meryl read*.meryl

cd ${WORKDIR}/evaluation_assemblies/merqury/flye_polished

# Step 2 - merqury assembly evaluation for flye
# - k-mer counts of the read set
# - Assembly fasta file
# - Output prefix
singularity exec \
--bind ${WORKDIR} \
/software/singularity/containers/Merqury-1.3-1.ubuntu20.sif \
merqury.sh \
${dir_output_meryl}/Illumina.meryl \
${WORKDIR}/polished_assemblies/flye/pilon.fasta \
flye_polished

cd ${WORKDIR}/evaluation_assemblies/merqury/canu_polished

# Step 2 - merqury assembly evaluation for canu
singularity exec \
--bind ${WORKDIR} \
/software/singularity/containers/Merqury-1.3-1.ubuntu20.sif \
merqury.sh \
${dir_output_meryl}/Illumina.meryl \
${WORKDIR}/polished_assemblies/canu/pilon.fasta \
canu_polished
