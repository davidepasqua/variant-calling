#!/bin/bash

# Set the base directory for the variant calling project
BASE_DIR="/Users/davide/Desktop/VARIANT_CALLING"


# run this into your terminal to find the path to the configuration file 
sudo find / -name "snpEff.config"

# Create a directory for the genome configuration files and copy the file 
mkdir -p config_file

cp -r /opt/homebrew/Cellar/snpeff/4.3t/share/snpeff/ snpEff.config $BASE_DIR"/config_file"

# SnpEff allows downloading many pre-set genomes in the "snpEff.config" configuration file
# Example for downloading the human genome
#  snpEff download -c config_file/snpEff.config -v hg38
#  snpEff download -c config_file/snpEff.config -v GRCh37.75

# If the genome of interest is not present, as in our case, the configuration file must be modified
# and parameters must be inserted to "build" snpEff files for variant annotation

# We will need the names of our "chromosomes" present in the genome:
#   NC_003197.2
#   NC_003277.2


# Copy the gbf file to a dedicated folder in the configuration file
mkdir -p $BASE_DIR"/config_file/data/ASM694v2"
cp -r $BASE_DIR/Genome/GCF_000006945.2_ASM694v2_genomic.gbff $BASE_DIR/config_file/data/ASM694v2/genes.gbk

cd $BASE_DIR
# Modify the configuration file as follows
# More information on modifying the configuration file can be found here: https://pcingola.github.io/SnpEff/se_buildingdb/
# Paragraph "Step 2 - building a database from GBK files"
echo "" >> config_file/snpEff.config
echo "# S.enterica thypinurium" >> config_file/snpEff.config
echo "ASM694v2.genome: Salmonella enterica thypinurium" >> config_file/snpEff.config
echo "   ASM694v2.chromosomes: NC_003197.2, NC_003277.2" >> config_file/snpEff.config
echo "   ASM694v2.NC_003197.2.codonTable: Bacterial_and_Plant_Plastid" >> config_file/snpEff.config
echo "   ASM694v2.NC_003277.2.codonTable: Bacterial_and_Plant_Plastid" >> config_file/snpEff.config
echo "" >> config_file/snpEff.config

# Create the file for annotation
snpEff build -genbank -c config_file/snpEff.config -v ASM694v2

# perform annotation
snpEff -c config_file/snpEff.config ASM694v2 VCF/ISZ-2_S5.VARIANTS.vcf > VCF/ISZ-2_S5.ANN.vcf
