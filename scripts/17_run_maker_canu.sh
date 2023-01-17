#!/usr/bin/env bash

#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=32
#SBATCH --job-name=maker
#SBATCH --output=/data/users/kdobler/assembly_annotation_course/outputs/output_maker_canu_%j.o
#SBATCH --error=/data/users/kdobler/assembly_annotation_course/outputs/error_maker_canu_%j.e
#SBATCH --partition=pcourseassembly

# Load required softwares
module load SequenceAnalysis/GenePrediction/maker/2.31.9;
module load Blast/ncbi-blast/2.9.0+;
module load UHTS/Analysis/busco/4.1.4;

# Set the working directory
WORKDIR=/data/users/kdobler/assembly_annotation_course
OUTPUTDIR=${WORKDIR}/maker_outputs/canu

cd ${OUTPUTDIR}

# 1. Create control files
maker -CTL

# 2. Edit maker_opts.ctl adapting following parameters
# -

# 3. Run MAKER
maker

# 4. Generate gff and fasta files
gff3_merge -d pilon.maker.output/pilon_master_datastore_index.log
fasta_merge -d pilon.maker.output/pilon_master_datastore_index.log

# 5. Replace gene IDs with shorter ones
# Build shorter IDs/Names for MAKER genes and transcripts
maker_map_ids --prefix LER_ pilon.all.gff > pilon.all.id.map

# Map short IDs/NAmes to MAKER fasta and gff files
map_fasta_ids pilon.all.id.map pilon.all.maker.proteins.fasta
map_fasta_ids pilon.all.id.map pilon.all.maker.transcripts.fasta
map_gff_ids pilon.all.id.map pilon.all.gff

# 7. Annotate proteins putative functions
# Blastp predicted proteins against the UniProt/SwissProt database
makeblastdb \
 -in /data/courses/assembly-annotation-course/CDS_annotation/uniprot-plant_reviewed.fasta \
 -dbtype prot \
 -out blast_db/uniprot_plant_reviewed

blastp \
 -query pilon.all.maker.proteins.fasta -db blast_db/uniprot_plant_reviewed \
 -num_threads 10 \
 -outfmt 6 \
 -evalue 1e-10 \
 > blastp_output.fa

# Map putative functions to the MAKER produced GFF3 and FASTA files
maker_functional_fasta \
 /data/courses/assembly-annotation-course/CDS_annotation/uniprot-plant_reviewed.fasta \
 blastp_output.fa \
 pilon.all.maker.proteins.fasta \
 > maker_functional.fa

maker_functional_gff \
 /data/courses/assembly-annotation-course/CDS_annotation/uniprot-plant_reviewed.fasta \
 blastp_output.fa \
 pilon.all.gff \
 > maker_functional.gff

# 8. Assess annotation quality with BUSCO
cd ${WORKDIR}/busco_anno_outputs/canu
busco \
 -i ${WORKDIR}/maker_outputs/canu/pilon.all.maker.proteins.fasta \
 -l brassicales_odb10 \
 -m proteins \
 -c 4 \
 --out busco_output
