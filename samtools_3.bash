#!/bin/bash

BASE_DIR="/Users/davide/Desktop/VARIANT_CALLING"
ALIGNMENTS_DIR="${BASE_DIR}/Alignments"

FILES=(
    "ISZ-20_S50.paired"
    "ISZ-20_S50.unpaired_2"
    "ISZ-2_S5.unpaired_1"
    "ISZ-7_S10.paired"
    "ISZ-7_S10.unpaired_2"
    "ISZ-20_S50.unpaired_1"
    "ISZ-2_S5.paired"
    "ISZ-2_S5.unpaired_2"
    "ISZ-7_S10.unpaired_1"
)

for FILE in "${FILES[@]}"; do
    SAM_FILE="${ALIGNMENTS_DIR}/${FILE}.sam"
    BAM_FILE="${ALIGNMENTS_DIR}/${FILE}.bam"
    SORTED_BAM_FILE="${ALIGNMENTS_DIR}/${FILE}.SORTED.bam"

    # Convert SAM to BAM
    samtools view -S -b "${SAM_FILE}" > "${BAM_FILE}"

    # Sort BAM files
    samtools sort "${BAM_FILE}" > "${SORTED_BAM_FILE}"

    # Merge sorted BAM files
    samtools merge "${ALIGNMENTS_DIR}/MERGED.bam" "${ALIGNMENTS_DIR}"/*.SORTED.bam

    # Sort merged BAM file
    samtools sort "${ALIGNMENTS_DIR}/MERGED.bam" > "${ALIGNMENTS_DIR}/MERGED.SORTED.bam"

    # Index merged and sorted BAM file
    samtools index "${ALIGNMENTS_DIR}/MERGED.SORTED.bam" > "${ALIGNMENTS_DIR}/MERGED.SORTED.bam.bai"

    # Calculate coverage position by position (excluding positions where no reads map)
    samtools depth "${ALIGNMENTS_DIR}/MERGED.SORTED.bam" > "${ALIGNMENTS_DIR}/MERGED.DEPTH"

    # Calculate coverage position by position (including positions where no reads map)
    samtools depth -a "${ALIGNMENTS_DIR}/MERGED.SORTED.bam" > "${ALIGNMENTS_DIR}/MERGED.DEPTH.ALL"

    # Subset BAM file with only one alignment in specific positions of a specific chromosome
    samtools view -h -b "${ALIGNMENTS_DIR}/MERGED.SORTED.bam" "NC_003197.2:10000-12000" > "${ALIGNMENTS_DIR}/subset.bam"

    # Get detailed information on all alignments using flags
    samtools flagstat "${ALIGNMENTS_DIR}/MERGED.SORTED.bam"

    # Calculate phred score of nucleotides of reads covering each position
    samtools mpileup "${ALIGNMENTS_DIR}/MERGED.SORTED.bam" > "${ALIGNMENTS_DIR}/MERGED.SORTED.MPILEUP"

    # Obtain consensus sequence
    samtools consensus "${ALIGNMENTS_DIR}/MERGED.SORTED.bam" > "${ALIGNMENTS_DIR}/consensus.fasta"

    
done

# Visualize with tview
    samtools tview "${ALIGNMENTS_DIR}/MERGED.SORTED.bam" "${BASE_DIR}/Genome/GCF_000006945.2_ASM694v2_genomic.fna"
    samtools tview "${ALIGNMENTS_DIR}/MERGED.SORTED.bam" "${BASE_DIR}/Genome/GCF_000006945.2_ASM694v2_genomic.fna" -p "NC_003197.2:831"