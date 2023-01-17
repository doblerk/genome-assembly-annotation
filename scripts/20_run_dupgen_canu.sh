#!/usr/bin/env bash

#SBATCH --time=08:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=16
#SBATCH --job-name=dupgen
#SBATCH --output=/data/users/kdobler/assembly_annotation_course/outputs/output_dupgen_canu_%j.o
#SBATCH --error=/data/users/kdobler/assembly_annotation_course/outputs/error_dupgen_canu_%j.e
#SBATCH --partition=pcourseassembly

# Load required softwares
module load Blast/ncbi-blast/2.9.0+;
module load UHTS/Analysis/SeqKit/0.13.2;
module add Emboss/EMBOSS/6.6.0;
module load Phylogeny/paml/4.9j;

# Set the working directory
WORKDIR=/data/users/kdobler/assembly_annotation_course
OUTPUTDIR=${WORKDIR}/dupgen_outputs/canu
NNUDIR=/data/courses/assembly-annotation-course/CDS_annotation
PERLSCRIPTS=/data/courses/assembly-annotation-course/CDS_annotation/scripts

cd ${OUTPUTDIR}

# 1. Prepare input files
# Reformat some files
awk '{sub(/ID=/, "", $9); sub(/;.*/, "", $9); sub(/:.*/, "", $9); printf ("%s\t%s\t%s\t%s\n", $1, $9, $4, $5)}' ${WORKDIR}/maker_outputs/canu/pilon.all.gff | grep '\-RA' > Ath.gff
cat Ath.gff ${NNUDIR}/NNU_mRNA_single_model.gff > Ath_Nnu.gff 

# Blast fasta maker against itself, then protein against output
# Create a reference database
makeblastdb \
 -in ${WORKDIR}/maker_outputs/canu/pilon.all.maker.proteins.fasta \
 -dbtype prot \
 -out blast_db/Ath

# Align protein query sequences against the reference database
blastp \
 -query ${WORKDIR}/maker_outputs/canu/pilon.all.maker.proteins.fasta \
 -db blast_db/Ath \
 -evalue 1e-10 \
 -max_target_seqs 5 \
 -outfmt 6 \
 -num_threads 16 \
 -out Ath.blast

blastp \
 -query ${NNUDIR}/NNU.pep.fa.ref.single_model \
 -db blast_db/Ath \
 -evalue 1e-10 \
 -max_target_seqs 5 \
 -outfmt 6 \
 -num_threads 16 \
 -out Ath_Nnu.blast

# 2. Run DupGen_finder locally
DupGen_finder.pl -i ${OUTPUTDIR} -t Ath -c Nnu -o ${OUTPUTDIR}

# 3. Estimate genome duplication events

# Create two fasta files with CDS and protein sequences of WGD genes
cut -f 1 Ath.wgd.genes > ID.txt
seqkit grep -f ID.txt ${WORKDIR}/maker_outputs/canu/pilon.all.maker.transcripts.fasta -o Ler.wgd.genes.fa
seqkit translate Ler.wgd.genes.fa -o Ler.wgd.genes.fa.proteins
sed -i 's/_frame=1/_p/' Ler.wgd.genes.fa.proteins

# Create a list of WGD gene pairs
cut -f 1,3 Ath.wgd.pairs > wgd_pairs.csv
sed -i 's/RA/RA_p/g' wgd_pairs.csv

# Copy wgd_pairs.csv, Ler.wgd.genes.fa, and Ler.wgd.genes.fa.proteins to a new directory
mkdir LerDir
cd LerDir
cp ../wgd_pairs.csv .
cp ../Ler.wgd.genes.fa .
cp ../Ler.wgd.genes.fa.proteins .

# Split fasta files
${PERLSCRIPTS}/split_flat Ler.wgd.genes.fa
${PERLSCRIPTS}/split_flat Ler.wgd.genes.fa.proteins

# Align proteins with bestflash_from_list
${PERLSCRIPTS}/bestflash_from_list wgd_pairs.csv

# Generate codon-by-codon CDS alignments based on the protein alignments with pair_to_CDS_paml_pair
${PERLSCRIPTS}/pair_to_CDS_paml_pair LER pair

# Use yn00 from PAML (PAML_yn00_from_CDS_pair script) on all WGD pair alignments to calculate Ka/Ks (omega) values
${PERLSCRIPTS}/PAML_yn00_from_CDS_pair LER > PAML_yn00_results

# 4. Calculate average Ks value for each syntenic block

# Parse PAML_yn00_results
awk '{print($1,$1,$6,$7,$5)}' PAML_yn00_results | \
 sed 's/ /\t/g' | \
 sed 's/__x__/\t/' | \
 sed 's/_p//g' | \
 cut -f 1,2,4,5,6 | \
 sed 's/dN=//' | \
 sed 's/dS=//' | \
 sed 's/omega=//' | \
 awk '$4<5' > Ler.wgd.kaks

cd average_ks

# Combine the files using add_ka_ks_to_collinearity_file.pl
perl ${PERLSCRIPTS}/add_ka_ks_to_collinearity_file.pl Ler

# Calculate average Ks value for each syntenic block using compute_ks_for_synteny_blocks.pl
perl ${PERLSCRIPTS}/compute_ks_for_synteny_blocks.pl Ler.collinearity.kaks

# Estimate Ks peaks from Ks distribution using plot_syntenic_blocks_ks_distri.py
# do it locally python 2.7
python ${PERLSCRIPT}/plot_syntenic_blocks_ks_distri.py Ler.synteny.blocks.ks.info 2 Ler
