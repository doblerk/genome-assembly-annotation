#!/usr/bin/env bash

#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=16
#SBATCH --job-name=te_sorter
#SBATCH --output=/data/users/kdobler/assembly_annotation_course/outputs/output_te_sorter_%j.o
#SBATCH --error=/data/users/kdobler/assembly_annotation_course/outputs/error_te_sorter_%j.e
#SBATCH --partition=pcourseassembly

# Load required software
module load UHTS/Analysis/BEDTools/2.29.2;

# Set the working directory
WORKDIR=/data/users/kdobler/assembly_annotation_course
INPUTDIR=${WORKDIR}/edta_outputs/canu
OUTPUTDIR=${WORKDIR}/tesorter_outputs/canu
CONTAINER=/data/courses/assembly-annotation-course/containers2


cd ${OUTPUTDIR}

# 1.
awk '$3~/retrotransposon/' ${INPUTDIR}/pilon.fasta.mod.EDTA.TEanno.gff3 > pilon.fasta.mod.EDTA.TEanno.gff3_edited

grep -v LTR ${INPUTDIR}/pilon.fasta.mod.EDTA.TEanno.gff3 >> pilon.fasta.mod.EDTA.TEanno.gff3_edited

sed 's/\(ID=\)\(\w*\)\(;.*\)/\2/g' pilon.fasta.mod.EDTA.TEanno.gff3_edited | \
 awk '{printf ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $9, $4, $5, $6, $7, $8, $3)}' \
 > parsed.tab

sed -i 's/_pi\s/_pilon\t/g' parsed.tab
sed -i 's/_pil\s/_pilon\t/g' parsed.tab
sed -i 's/_pilo\s/_pilon\t/g' parsed.tab

# Alternative
sed 's/\;Classification.*//' pilon.fasta.mod.EDTA.TEanno.gff3_edited > test.txt
sed 's/ID.*Name\=//' test.txt > test2.txt
awk -F'\t' -v OFS="\t" '{ print $1, $2, $9, $4, $5, $6, $7, $8, $3 }' test2.txt > pilon.fasta.mod.EDTA.TEanno.gff3_simplified
sed 's//contig/g' pilon.fasta.mod.EDTA.TEanno.gff3_simplified > correction
sed 's/\_pi.*\tEDTA/\_pilon\tEDTA/' pilon.fasta.mod.EDTA.TEanno.gff3_simplified > pilon.fasta.mod.EDTA.TEanno.gff3_final

bedtools getfasta \
 -s \
 -name \
 -fi ${WORKDIR}/polished_assemblies/canu/pilon.fasta \
 -bed pilon.fasta.mod.EDTA.TEanno.gff3_final \
 -fo Ler_canu_TEs

singularity exec \
 --bind $CONTAINER \
 --bind $WORKDIR \
${CONTAINER}/TEsorter_1.3.0.sif \
 TEsorter Ler_canu_TEs -db rexdb-plant

singularity exec \
 --bind $CONTAINER \
 --bind $WORKDIR \
 --bind /data/courses/assembly-annotation-course \
${CONTAINER}/TEsorter_1.3.0.sif \
 TEsorter /data/courses/assembly-annotation-course/Brassicaceae_repbase_all_march2019.fasta \
 -pre brassicaceae \
 -db rexdb-plant 
