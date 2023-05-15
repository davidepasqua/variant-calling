# Set the base directory for the variant calling project
BASE_DIR="/Users/davide/Desktop/VARIANT_CALLING"

# Create the VCF directory inside the base directory
mkdir -p $BASE_DIR/VCF

# Perform mpileup using bcftools
# -o: output file
# -f: reference genome file
# --annotate: specify the annotations to include in the output

bcftools mpileup -o VCF/ISZ-2_S5.vcf \
	-f $BASE_DIR/Genome/GCF_000006945.2_ASM694v2_genomic.fna \
	--annotate FORMAT/DP,FORMAT/AD \
	$BASE_DIR/Alignments/MERGED.SORTED.bam 

# Call the variants using bcftools
# -o: output file
# --ploidy: set ploidy to 1 since it's a bacterium
# -m: call rare variants
# -v: report only positions with at least one difference from the reference
# In the output file, the most probable genotype will be called based on the PL column

bcftools call --ploidy 1 -m -v -o VCF/ISZ-2_S5.VARIANTS.vcf VCF/ISZ-2_S5.vcf