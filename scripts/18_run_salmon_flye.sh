#!/usr/bin/env bash

#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=16
#SBATCH --job-name=salmon
#SBATCH --output=/data/users/kdobler/assembly_annotation_course/outputs/output_salmon_flye_%j.o
#SBATCH --error=/data/users/kdobler/assembly_annotation_course/outputs/error_salmon_flye_%j.e
#SBATCH --partition=pcourseassembly

# Load required softwares

# Set the working directory
WORKDIR=/data/users/kdobler/assembly_annotation_course
OUTPUTDIR=${WORKDIR}/salmon_outputs/flye
CONTAINERDIR=/software/singularity/containers
READSDIR=/data/courses/assembly-annotation-course/raw_data/Ler/participant_2/RNAseq

# 1. Quantify gene expression using Salmon
cd ${OUTPUTDIR}

# 2 Build index of your predicted genes
singularity exec \
--bind ${WORKDIR} \
--bind ${CONTAINERDIR} \
${CONTAINERDIR}/salmon-1.7.0-1.ubuntu18.sif \
salmon index \
 -t ${WORKDIR}/maker_outputs/flye/pilon.all.maker.transcripts.fasta \
 -i transcripts_index \
 -k 31

# 2.2 Quantify any set of reads directly against this index
singularity exec \
--bind ${WORKDIR} \
--bind ${CONTAINERDIR} \
--bind ${READSDIR} \
${CONTAINERDIR}/salmon-1.7.0-1.ubuntu18.sif \
salmon quant \
 -i transcripts_index \
 -l A \
 -1 <(zcat ${READSDIR}/SRR1734309_1.fastq.gz)\
 -2 <(zcat ${READSDIR}/SRR1734309_1.fastq.gz) \
 --validateMappings \
 -o transcripts_quant
