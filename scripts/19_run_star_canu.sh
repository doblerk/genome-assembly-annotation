#!/usr/bin/env bash

#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=16
#SBATCH --job-name=star
#SBATCH --output=/data/users/kdobler/assembly_annotation_course/outputs/output_star_canu_%j.o
#SBATCH --error=/data/users/kdobler/assembly_annotation_course/outputs/error_star_canu_%j.e
#SBATCH --partition=pcourseassembly

# Load required softwares
module load UHTS/Assembler/cufflinks/2.2.1;
module add UHTS/Aligner/STAR/2.7.9a;

# Set the working directory
WORKDIR=/data/users/kdobler/assembly_annotation_course
OUTPUTDIR=${WORKDIR}/star_outputs/canu
GFFDIR=${WORKDIR}/maker_outputs/canu/maker_functional.gff
DATADIR=${WORKDIR}/polished_assemblies/canu/pilon.fasta
READSDIR=${WORKDIR}/participant_2/RNAseq

# 0. Convert gff to gtf
gffread -E ${GFFDIR} -T -o ${OUTPUTDIR}/maker_functional.gtf

# 1. Generate index for STAR
STAR \
 --runThreadN 4 \
 --runMode genomeGenerate \
 --genomeDir ${OUTPUTDIR}/genomeDir \
 --genomeFastaFiles ${DATADIR} \
 --sjdbGTFfile ${OUTPUTDIR}/maker_functional.gtf

# 3. Run mapping of RNA-Seq reads
STAR \
 --runThreadN 4 \
 --genomeDir ${OUTPUTDIR}/genomeDir \
 --readFilesIn <(zcat ${READSDIR}/SRR1734309_1.fastq.gz) <(zcat ${READSDIR}/SRR1734309_2.fastq.gz) \
 --outFileNamePrefix ${OUTPUTDIR}/star \
 --outSAMtype BAM SortedByCoordinate \
 --outFilterMultimapNmax 10 \
 --outFilterMismatchNoverLmax 0.01 \
 --alignIntronMax 60000
