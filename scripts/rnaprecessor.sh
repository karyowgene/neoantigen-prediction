#!/bin/bash

set -e  

#Source directory where your files are located
source_dir="/user/HDD/rna_seq"

#Loop through files in the source directory
for file in "$source_dir"/*_R1_001.fastq.gz; do
    if [ -f "$file" ]; then
        #extract the sample ID from the filename (e.g., "3_Ambs_203_S65_R1_001.fastq.gz" -> "203")
        sample_id=$(basename "$file" | cut -d '_' -f 3)
        
        #Rename R1 file
        mv "$file" "$source_dir/${sample_id}_R1.fastq.gz"
        
        #Rename R2 file
        mv "$source_dir/$(basename "$file" | sed 's/R1/R2/')" "$source_dir/${sample_id}_R2.fastq.gz"
        
        echo "Renamed files for sample $sample_id."
    fi
done

echo "File renaming completed."

echo "checking reads quality"

fastqc /user/HDD/rna_seq/*fastq.gz -t 4

echo "Done fastqc!"

echo "performing multiqc...."

multiqc /user/HDD/rna_seq -o /user/HDD/rna_seq

echo "Done multiqc!"

echo "Trimming the reads"

# Define the input and output directories
input_dir="/user/HDD/rna_seq"
output_dir="/user/HDD2/rna_seq/trimmed"

#Create output directory if it doesn't exist
mkdir -p "$output_dir"

#Iterate over your sample pairs
for R1_input in "${input_dir}"/*_R1.fastq.gz; do
    #Check if there are matching R2 reads
    R2_input="${R1_input/_R1.fastq.gz/_R2.fastq.gz}"
    if [ -e "$R2_input" ]; then
        #Define sample name
        sample_name=$(basename "$R1_input" _R1.fastq.gz)

        #Define output file paths
        R1_output="${output_dir}/${sample_name}_R1_paired.fastq.gz"
        R1_unpaired_output="${output_dir}/${sample_name}_R1_unpaired.fastq.gz"
        R2_output="${output_dir}/${sample_name}_R2_paired.fastq.gz"
        R2_unpaired_output="${output_dir}/${sample_name}_R2_unpaired.fastq.gz"

        #Check if the output files already exist (indicating that the sample has been trimmed)
        if [ -e "$R1_output" ] && [ -e "$R1_unpaired_output" ] && [ -e "$R2_output" ] && [ -e "$R2_unpaired_output" ]; then
            echo "Sample $sample_name has already been trimmed. Skipping..."
        else
            #Run Trimmomatic
            trimmomatic PE -threads 4 "$R1_input" "$R2_input" \
                "$R1_output" "$R1_unpaired_output" "$R2_output" "$R2_unpaired_output" \
                ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:151
        fi
    else
        echo "Warning: No matching R2 read found for $R1_input"
    fi
done

echo "trimming completed successfully"

echo "fastqc for trimmed reads...."

fastqc /user/HDD2/rna_seq/trimmed/*_paired.fastq.gz -t 4

echo "fastqc trimmed reads Done!"

echo "multiqc trimmed reads...."

multiqc /user/HDD2/rna_seq/trimmed -o /user/HDD2/rna_seq/trimmed

echo "multiqc trimmed reads Done!"

echo "kallisto alignment starts.."

cd /user/HDD2/RNA/rnaref

curl -O ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

kallisto index -i Homo_sapiens.GRCh38.cdna.all.release-110.idx -t 4 Homo_sapiens.GRCh38.cdna.all.fa.gz

#reference transcriptome index
index="/user/HDD2/RNA/rnaref/Homo_sapiens.GRCh38.cdna.all.release-110.idx"

#FASTQ files
input_dir="/user/HDD2/rna_seq/trimmed"

#output directory
output_dir="/user/HDD2/rna_seq/kallisto_quants"

#if out doesn't exist
mkdir -p "$output_dir"

#Loop through FASTQs in the input_dir
for R1_file in "$input_dir"/*_R1_paired.fastq.gz; do
    if [ -f "$R1_file" ]; then
        #Extract the sample name
        sample_name=$(basename "$R1_file" | sed 's/_R1_paired.fastq.gz//')
        
        #Create output subdirectory for current sample
        sample_output_dir="$output_dir/$sample_name"
        mkdir -p "$sample_output_dir"
        
        #Define paired R2 file
        R2_file="$input_dir/${sample_name}_R2_paired.fastq.gz"
        
        #Run Kallisto for the current sample
        kallisto quant -i "$index" -o "$sample_output_dir" "$R1_file" "$R2_file"
        
        echo "Kallisto pseudo-alignment and abundances quantification completed for $sample_name."
    fi
done

echo "Kallisto analysis for all samples completed."