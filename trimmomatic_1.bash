source activate trimmo_env

INP_DIR="/Users/davide/Desktop/VARIANT_CALLING/InputData"
########### TRIMMOMATIC ###########

# Create output directory for Trimmomatic output
Trim_out_dir="$INP_DIR/Trimmomatic_output"
mkdir -p $Trim_out_dir

# Loop through each forward read file in the directory
for forward in $INP_DIR/*R1*.fastq; do
	
	printf "\n#############################################\n"
	
	reverse=${forward/R1/R2}
	
	printf "analisi dei campioni: %s\t%s\n" $forward $reverse
	ls $forward
	ls $reverse
	
	basename_file=$(basename ${forward})
	
	output_trimlog=${Trim_out_dir}/${basename_file/_L001_R1_001.fastq/_trimlog.txt}
	output_for_pai=${Trim_out_dir}/${basename_file/_L001_R1_001.fastq/_R1_pai.fastq}
	output_rev_pai=${Trim_out_dir}/${basename_file/_L001_R1_001.fastq/_R2_pai.fastq}
	output_for_unp=${Trim_out_dir}/${basename_file/_L001_R1_001.fastq/_R1_unp.fastq}
	output_rev_unp=${Trim_out_dir}/${basename_file/_L001_R1_001.fastq/_R2_unp.fastq}
	
	trimmomatic PE -phred33 -threads 10 -trimlog $output_trimlog $forward $reverse $output_for_pai $output_for_unp $output_rev_pai $output_rev_unp HEADCROP:15 TRAILING:27 LEADING:27 MINLEN:250 AVGQUAL:27
	
done
