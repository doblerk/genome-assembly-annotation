#!/usr/bin/env bash

#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=quast
#SBATCH --output=/data/users/kdobler/assembly_annotation_course/outputs/output_quast_with_%j.o
#SBATCH --error=/data/users/kdobler/assembly_annotation_course/outputs/error_quast_with_%j.e
#SBATCH --partition=pcourseassembly

# Define paths to directories
PROJDIR=/data/users/kdobler/assembly_annotation_course
dir_ref=/data/courses/assembly-annotation-course

cd $PROJDIR/evaluation_assemblies/quast_with_ref_all

# Run QUAST without reference
singularity exec \
--bind $PROJDIR \
--bind $dir_ref \
/data/courses/assembly-annotation-course/containers/quast_5.1.0rc1.sif \
quast.py $PROJDIR/assemblies/flye/assembly.fasta \
 $PROJDIR/polished_assemblies/flye/pilon.fasta \
 $PROJDIR/assemblies/canu/canu_pacbio.contigs.fasta \
 $PROJDIR/polished_assemblies/canu/pilon.fasta \
 -r $dir_ref/references/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz \
 -g $dir_ref/references/TAIR10_GFF3_genes.gff \
 --large \
 --eukaryote \
 --est-ref-size 130000000 \
 --threads 4 \
 --labels "flye unpolished","flye polished","canu unpolished","canu polished"\
 -o $PROJDIR/evaluation_assemblies/quast_with_ref_all \
 --pacbio $PROJDIR/participant_2/pacbio/ERR3415825.fastq.gz \
 --pacbio $PROJDIR/participant_2/pacbio/ERR3415826.fastq.gz \
 --no-sv
