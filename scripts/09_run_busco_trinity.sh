#!/usr/bin/env bash

#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=busco
#SBATCH --output=/data/users/kdobler/assembly_annotation_course/outputs/output_busco_trinity_%j.o
#SBATCH --error=/data/users/kdobler/assembly_annotation_course/outputs/error_busco_trinity_%j.e
#SBATCH --partition=pall

# Define paths to directories
working_dir=/data/users/kdobler/assembly_annotation_course
dir_input=/data/users/kdobler/assembly_annotation_course/assemblies/trinity/Trinity.fasta
PROJDIR=/data/users/kdobler/assembly_annotation_course/evaluation_assemblies/busco_trinity

cd ${PROJDIR}

echo "busco \
    --in ${dir_input} \
    --mode transcriptome \
    --cpu 4 \
    --lineage brassicales_odb10 \
    --out trinity " > eval_trini
singularity exec \
 --bind ${working_dir} \
 /data/courses/assembly-annotation-course/containers2/busco_v5.1.2_cv1.sif \
 busco 

# Define paths to directories
working_dir=/data/users/kdobler/assembly_annotation_course
dir_input=/data/users/kdobler/assembly_annotation_course/assemblies/canu/
PROJDIR=/data/users/kdobler/assembly_annotation_course/evaluation_assemblies/busco_canu_unpolished

cd ${PROJDIR}

# UNPOLISHED
singularity exec \
--bind $working_dir \
/data/courses/assembly-annotation-course/containers2/busco_v5.1.2_cv1.sif \
busco -i $dir_input/canu_pacbio.contigs.fasta -l brassicales_odb10 -o unpolished_assembly -m genome --cpu 4 --out_path $PROJDIR

# POLISHED
singularity exec \
--bind $working_dir \
/data/courses/assembly-annotation-course/containers2/busco_v5.1.2_cv1.sif \
busco -i $working_dir/polished_assemblies/canu/pilon.fasta -l brassicales_odb10 -o polished_assembly -m genome --cpu 4 --out_path $PROJDIR
