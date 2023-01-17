#!/usr/bin/env bash

#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=suppl
#SBATCH --output=/data/users/kdobler/assembly_annotation_course/outputs/output_suppl_%j.o
#SBATCH --error=/data/users/kdobler/assembly_annotation_course/outputs/error_suppl_%j.e
#SBATCH --partition=pall

# PART 1

WORKDIR=/data/users/kdobler/assembly_annotation_course/edta_outputs/canu
WORKDIR=/data/users/kdobler/assembly_annotation_course/edta_outputs/canu
DATA1=${WORKDIR}/pilon.fasta.mod.EDTA.TElib.fa
DATA2=${WORKDIR}/pilon.fasta.mod.EDTA.TEanno.gff3

cd ${WORKDIR}

# 1.
# Retrieve the TE family of each TE
awk '/TE_/' ${DATA1} | cut -d '#' -f 2 > TE_families.txt
awk '$9' ${DATA2} | cut -d ';' -f 3 | cut -d '=' -f 2 > TE_families2.txt

# Get the occurences of each family
cat TE_families2.txt | sort | uniq -c > TE_families_occurences2.txt

# 2.
# Compute the mean of active vs inactive TEs according to a chosen threshold.
awk '$9 ~ /Method=homology/' ${DATA2} | awk '$9' | cut -d ';' -f 5 | cut -d '=' -f 2 > Identity_perc.txt
awk 'BEGIN {active=0; inactive=0} {if ($1 >= 0.8) {++active} else {++inactive}} END {print "Active sites: ", active/NR, "Inactive sites: ", inactive/NR}' Identity_perc.txt > ActiveVsInactive.txt

# 3.
# Determine how many intact and fragmented copies are detected per TE family
# Get all fragmented TEs
awk '$9 ~ /Method=homology/' ${DATA2} | awk '$9' | cut -d ';' -f 3 | cut -d '=' -f 2 > fragmented_TEs.txt
cat fragmented_TEs.txt | sort | uniq -c > fragmented_TEs_occurences.txt

# Get all intact TEs
awk '$9 ~ /Method=structural/' ${DATA2} | awk '$9' | cut -d ';' -f 3 | cut -d '=' -f 2 > intact_TEs.txt
cat intact_TEs.txt | sort | uniq -c > intact_TEs_occurences.txt

# PART 2

# Define the directories of interest
WORKDIR=/data/users/kdobler/assembly_annotation_course/ANNOTATION
OUTPUTDIR=${WORKDIR}/phylogenetic_analysis/flye
CLADENAMES=${WORKDIR}/tesorter_outputs/flye/Ler_flye_TEs.rexdb-plant.cls.tsv
CLADENAMESBRASS=${WORKDIR}/tesorter_outputs/flye/brassicaceae.cls.tsv
COLORSETS=/data/courses/assembly-annotation-course/CDS_annotation/dataset_color_strip_template.txt

# 1.
# annotate main TE clades for copia
cd ${OUTPUTDIR}/gypsy

# Retrive all associated families and superclades
awk '$3 ~ /Gypsy/' ${CLADENAMES} | cut -f 4 | sort | uniq > gypsy_superclade.txt

# Annotate gypsy elements for the different superclades
for superclade in $(cat gypsy_superclade.txt)
do
 awk '$3 ~ /Gypsy/ && $4 ~ /"{$superclade}"/' ${CLADENAMES} | cut -f 1 | sed 's/:/_/g' | sed 's/$/ #9606/' > ${superclade}_ID.txt; 
done

# FOR GYPSY (CANU)
cd ${OUTPUTDIR}/gypsy

awk '$3 ~ /Gypsy/' ${CLADENAMES} | cut -f 4 | sort | uniq > gypsy_superclade.txt

grep -e "Retand" ${CLADENAMES} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #CD5C5C/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' > Retand_ID.txt
grep -e "Athila" ${CLADENAMES} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #FFB6C1/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' > Athila_ID.txt
grep -e "CRM" ${CLADENAMES} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #FFA07A/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' > CRM_ID.txt
grep -e "Chlamyvir" ${CLADENAMES} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #DDA0DD/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' > Chlamyvir_ID.txt
grep -e "Galadriel" ${CLADENAMES} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #9932CC/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' > Galadriel_ID.txt
grep -e "Reina" ${CLADENAMES} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #3CB371/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' > Reina_ID.txt
grep -e "TatII" ${CLADENAMES} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #4682B4/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' > TatII_ID.txt
grep -e "Tcn1" ${CLADENAMES} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #1E90FF/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' > Tcn1_ID.txt
grep -e "Tekay" ${CLADENAMES} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #87CEEB/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' > Tekay_ID.txt

grep -e "Retand" ${CLADENAMESBRASS} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #CD5C5C/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' >> Retand_ID.txt
grep -e "Athila" ${CLADENAMESBRASS} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #FFB6C1/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' >> Athila_ID.txt
grep -e "CRM" ${CLADENAMESBRASS} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #FFA07A/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' >> CRM_ID.txt
grep -e "Chlamyvir" ${CLADENAMESBRASS} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #DDA0DD/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' >> Chlamyvir_ID.txt
grep -e "Galadriel" ${CLADENAMESBRASS} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #9932CC/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' >> Galadriel_ID.txt
grep -e "Reina" ${CLADENAMESBRASS} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #3CB371/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' >> Reina_ID.txt
grep -e "TatII" ${CLADENAMESBRASS} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #4682B4/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' >> TatII_ID.txt
grep -e "Tcn1" ${CLADENAMESBRASS} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #1E90FF/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' >> Tcn1_ID.txt
grep -e "Tekay" ${CLADENAMESBRASS} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #87CEEB/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' >> Tekay_ID.txt

cat Retand_ID.txt Athila_ID.txt CRM_ID.txt Chlamyvir_ID.txt Galadriel_ID.txt Reina_ID.txt TatII_ID.txt Tcn1_ID.txt Tekay_ID.txt > all_ID.txt

# FOR COPIA (CANU)
cd ${OUTPUTDIR}/copia

awk '$3 ~ /Copia/' ${CLADENAMES} | cut -f 4 | sort | uniq > copia_superclade.txt

grep -e "Ale" ${CLADENAMES} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #CD5C5C/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' > Ale_ID.txt
grep -e "Tork" ${CLADENAMES} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #FFB6C1/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' > Tork_ID.txt
grep -e "Bianca" ${CLADENAMES} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #FFA07A/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' > Bianca_ID.txt
grep -e "Ikeros" ${CLADENAMES} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #DDA0DD/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' > Ikeros_ID.txt
grep -e "Ivana" ${CLADENAMES} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #9932CC/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' > Ivana_ID.txt
grep -e "SIRE" ${CLADENAMES} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #3CB371/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' > SIRE_ID.txt
grep -e "Angela" ${CLADENAMES} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #4682B4/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' > Angela_ID.txt
grep -e "Bryco" ${CLADENAMES} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #1E90FF/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' > Bryco_ID.txt
grep -e "Lyco" ${CLADENAMES} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #87CEEB/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' > Lyco_ID.txt
grep -e "TAR" ${CLADENAMES} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #FFA500/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' > TAR_ID.txt

grep -e "Ale" ${CLADENAMESBRASS} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #CD5C5C/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' >> Ale_ID.txt
grep -e "Tork" ${CLADENAMESBRASS} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #FFB6C1/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' >> Tork_ID.txt
grep -e "Bianca" ${CLADENAMESBRASS} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #FFA07A/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' >> Bianca_ID.txt
grep -e "Ikeros" ${CLADENAMESBRASS} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #DDA0DD/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' >> Ikeros_ID.txt
grep -e "Ivana" ${CLADENAMESBRASS} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #9932CC/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' >> Ivana_ID.txt
grep -e "SIRE" ${CLADENAMESBRASS} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #3CB371/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' >> SIRE_ID.txt
grep -e "Angela" ${CLADENAMESBRASS} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #4682B4/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' >> Angela_ID.txt
grep -e "Bryco" ${CLADENAMESBRASS} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #1E90FF/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' >> Bryco_ID.txt
grep -e "Lyco" ${CLADENAMESBRASS} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #87CEEB/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' >> Lyco_ID.txt
grep -e "TAR" ${CLADENAMESBRASS} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #FFA500/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' >> TAR_ID.txt

cat Ale_ID.txt Tork_ID.txt Bianca_ID.txt Ikeros_ID.txt Ivana_ID.txt SIRE_ID.txt Angela_ID.txt Bryco_ID.txt Lyco_ID.txt TAR_ID.txt > all_ID.txt

# FOR GYPSY (FLYE)
cd ${OUTPUTDIR}/gypsy

awk '$3 ~ /Gypsy/' ${CLADENAMES} | cut -f 4 | sort | uniq > gypsy_superclade.txt

grep -e "Retand" ${CLADENAMES} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #CD5C5C/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' > Retand_ID.txt
grep -e "Athila" ${CLADENAMES} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #FFB6C1/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' > Athila_ID.txt
grep -e "CRM" ${CLADENAMES} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #FFA07A/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' > CRM_ID.txt
grep -e "Chlamyvir" ${CLADENAMES} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #DDA0DD/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' > Chlamyvir_ID.txt
grep -e "Galadriel" ${CLADENAMES} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #9932CC/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' > Galadriel_ID.txt
grep -e "Reina" ${CLADENAMES} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #3CB371/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' > Reina_ID.txt
grep -e "TatII" ${CLADENAMES} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #4682B4/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' > TatII_ID.txt
grep -e "Tcn1" ${CLADENAMES} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #1E90FF/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' > Tcn1_ID.txt
grep -e "Tekay" ${CLADENAMES} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #87CEEB/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' > Tekay_ID.txt

grep -e "Retand" ${CLADENAMESBRASS} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #CD5C5C/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' >> Retand_ID.txt
grep -e "Athila" ${CLADENAMESBRASS} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #FFB6C1/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' >> Athila_ID.txt
grep -e "CRM" ${CLADENAMESBRASS} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #FFA07A/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' >> CRM_ID.txt
grep -e "Chlamyvir" ${CLADENAMESBRASS} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #DDA0DD/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' >> Chlamyvir_ID.txt
grep -e "Galadriel" ${CLADENAMESBRASS} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #9932CC/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' >> Galadriel_ID.txt
grep -e "Reina" ${CLADENAMESBRASS} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #3CB371/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' >> Reina_ID.txt
grep -e "TatII" ${CLADENAMESBRASS} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #4682B4/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' >> TatII_ID.txt
grep -e "Tcn1" ${CLADENAMESBRASS} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #1E90FF/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' >> Tcn1_ID.txt
grep -e "Tekay" ${CLADENAMESBRASS} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #87CEEB/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' >> Tekay_ID.txt

cat Retand_ID.txt Athila_ID.txt CRM_ID.txt Chlamyvir_ID.txt Galadriel_ID.txt Reina_ID.txt TatII_ID.txt Tcn1_ID.txt Tekay_ID.txt > all_ID.txt

cat ${COLORSETS} all_ID.txt > dataset_color_strip_template_gypsy.txt

# FOR COPIA (FLYE)
cd ${OUTPUTDIR}/copia

awk '$3 ~ /Copia/' ${CLADENAMES} | cut -f 4 | sort | uniq > copia_superclade.txt

grep -e "Ale" ${CLADENAMES} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #CD5C5C/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' > Ale_ID.txt
grep -e "Tork" ${CLADENAMES} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #FFB6C1/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' > Tork_ID.txt
grep -e "Bianca" ${CLADENAMES} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #FFA07A/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' > Bianca_ID.txt
grep -e "Ikeros" ${CLADENAMES} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #DDA0DD/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' > Ikeros_ID.txt
grep -e "Ivana" ${CLADENAMES} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #9932CC/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' > Ivana_ID.txt
grep -e "SIRE" ${CLADENAMES} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #3CB371/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' > SIRE_ID.txt
grep -e "Angela" ${CLADENAMES} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #4682B4/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' > Angela_ID.txt
grep -e "Bryco" ${CLADENAMES} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #1E90FF/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' > Bryco_ID.txt
grep -e "Lyco" ${CLADENAMES} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #87CEEB/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' > Lyco_ID.txt
grep -e "TAR" ${CLADENAMES} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #FFA500/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' > TAR_ID.txt

grep -e "Ale" ${CLADENAMESBRASS} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #CD5C5C/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' >> Ale_ID.txt
grep -e "Tork" ${CLADENAMESBRASS} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #FFB6C1/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' >> Tork_ID.txt
grep -e "Bianca" ${CLADENAMESBRASS} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #FFA07A/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' >> Bianca_ID.txt
grep -e "Ikeros" ${CLADENAMESBRASS} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #DDA0DD/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' >> Ikeros_ID.txt
grep -e "Ivana" ${CLADENAMESBRASS} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #9932CC/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' >> Ivana_ID.txt
grep -e "SIRE" ${CLADENAMESBRASS} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #3CB371/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' >> SIRE_ID.txt
grep -e "Angela" ${CLADENAMESBRASS} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #4682B4/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' >> Angela_ID.txt
grep -e "Bryco" ${CLADENAMESBRASS} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #1E90FF/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' >> Bryco_ID.txt
grep -e "Lyco" ${CLADENAMESBRASS} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #87CEEB/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' >> Lyco_ID.txt
grep -e "TAR" ${CLADENAMESBRASS} | cut -f 1 | sed 's/:/_/g' | sed 's/__/_/g' | sed 's/$/ #FFA500/' | sed 's/(+)//' | sed 's/(-)//' | sed 's/(?)//' >> TAR_ID.txt

cat Ale_ID.txt Tork_ID.txt Bianca_ID.txt Ikeros_ID.txt Ivana_ID.txt SIRE_ID.txt Angela_ID.txt Bryco_ID.txt Lyco_ID.txt TAR_ID.txt > all_ID.txt

cat ${COLORSETS} all_ID.txt > dataset_color_strip_template_copia.txt
