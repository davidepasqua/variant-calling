#!/bin/bash

# Set the base directory for the variant calling project
BASE_DIR="/Users/davide/Desktop/VARIANT_CALLING"

# Create the Genome directory inside the base directory
mkdir -p $BASE_DIR Genome

# Download the genome files from NCBI
wget -P $BASE_DIR/Genome/ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.fna.gz
wget -P $BASE_DIR/Genome/ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.gff.gz
wget -P $BASE_DIR/Genome/ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.gbff.gz

# Unzip the downloaded genome files
gunzip $BASE_DIR/Genome/*

# Set the input directory containing the input files
INP_DIR="$BASE_DIR/InputData/Trim_out_dir"

# Set the output directory for the alignments
OUT_DIR="$BASE_DIR/Alignments"

# Create the output directory
mkdir -p "$OUT_DIR"

# Set the reference genome file name
GENOME_FILE="$BASE_DIR/Genome/GCF_000006945.2_ASM694v2_genomic.fna"

# Index the reference genome
bwa index "$GENOME_FILE"

# Process all fastq files in the input directory
for FILE in "$INP_DIR"/*.fastq
do
  # Check if the file is a fastq file
  if [[ $FILE == *.fastq ]]; then
    # Get the file name without extension
    FILENAME=$(basename -- "$FILE")
    FILENAME="${FILENAME%.*}"

    # Align paired-end reads
    if [[ "$FILENAME" == *_R1_pai ]]; then
      # Get the R1 paired-end file name without extension
      FILENAME_R1="${FILENAME%_R1_pai}"
      
      # Perform the alignment
      bwa mem -t 4 "$GENOME_FILE" \
        "$INP_DIR/$FILENAME_R1""_R1_pai.fastq" "$INP_DIR/$FILENAME_R1""_R2_pai.fastq" > "$OUT_DIR/$FILENAME_R1"".paired.sam"
    fi

    # Align unpaired reads
    if [[ "$FILENAME" == *_R1_unp ]]; then
      # Get the R1 unpaired file name without extension
      FILENAME_R1="${FILENAME%_R1_unp}"
      
      # Perform the alignment
      bwa mem -t 4 "$GENOME_FILE" \
        "$INP_DIR/$FILENAME_R1""_R1_unp.fastq" > "$OUT_DIR/$FILENAME_R1"".unpaired_1.sam"
      bwa mem -t 4 "$GENOME_FILE" \
        "$INP_DIR/$FILENAME_R1""_R2_unp.fastq" > "$OUT_DIR/$FILENAME_R1"".unpaired_2.sam"
    fi
  fi
done