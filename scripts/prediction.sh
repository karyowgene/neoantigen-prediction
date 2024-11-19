#!/bin/bash

set -e
set -u

echo "predicting neoantigens with pvacseq"

#STEP 2 predict neoantigens
# Define other parameters
#OTHER_PARAMS="-e1 7,8,9,10,11 --iedb-install-directory /opt/iedb -t 4"
OTHER_PARAMS="-e1 9,10,11 --iedb-install-directory /opt/iedb -t 1"
MORE_PARAMS="--normal-cov 0 --tdna-cov 0 --tdna-vaf 0 --normal-vaf 0.99 --trna-cov 0 --trna-vaf 0"

# Iterate over VCF files in the specified directory
ANNO_DIR="/media/user/Expansion/exome_data/annotated"
OUTPUT_BASE_DIR="/home/user/neoant3"
VCF_DIR="/pvacseq_mydataall_data"

#echo "VCF files in $VCF_DIR:"
#ls "$VCF_DIR"

for vcf_file in "$ANNO_DIR"/*.gx.vcf; do
    if [ -f "$vcf_file" ]; then
        # Extract sample name from the VCF file name (e.g., "203_annotated.vcf" -> "203")
        sample_name=$(basename "$vcf_file" | sed 's/_annotated.gx.vcf//')
        
        # Extract HLA alleles from hla-hd results
        normal_sample_name=$((sample_name + 1))
        
        # Define the input file path
        input_file="$VCF_DIR/$(basename "$vcf_file")"
        
        # Define the output directory for the current sample
        OUTPUT_DIR="$OUTPUT_BASE_DIR/$sample_name"
        
        # create if it doest exist
        mkdir -p "$OUTPUT_DIR"

        # Use the correct syntax for running the awk command and storing the result in HLA_TYPES
        HLA_TYPES="HLA-A*02:01,HLA-A*30:01,HLA-B*58:01,HLA-C*06:01"
        #HLA_TYPES=$(awk -F'\t' '{for (i = 2; i <= NF; i++) if ($i != "-" && $i ~ /HLA-[A-C]/) { gsub(/:[0-9][0-9]:/, ":", $i); printf "%s%s", sep, $i; sep = ","; } } END { print ""; }' /media/user/Expansion/exome_data/estimation/"$normal_sample_name"/result/"$normal_sample_name"_final.result.txt | sed 's/,$//')

        echo "executing pvacseq in docker container"
        # Execute the pvacseq run command for the current sample
        #docker run -v /media/user/Expansion/exome_data/estimation/"$normal_sample_name"/result:/result -v /media/user/Expansion/exome_data/annotated:/pvacseq_mydataall_data -v /home/user/neoant3:/pvacseq_outputall_mydata 981a3b8b2bd8 pvacseq run "$input_file" "$sample_name" "$HLA_TYPES" MHCflurry MHCnuggetsI SMM SMMPMBEC /pvacseq_outputall_mydata $OTHER_PARAMS -s 500 -d 500 --normal-sample-name "$normal_sample_name" $MORE_PARAMS
        
        # Execute the pvacseq run command for the current sample
        docker run -v /media/user/Expansion/exome_data/annotated:/pvacseq_mydataall_data -v "$OUTPUT_DIR":/pvacseq_outputall_mydata 981a3b8b2bd8 pvacseq run "$input_file" "$sample_name" "$HLA_TYPES" MHCflurry MHCnuggetsI SMM SMMPMBEC /pvacseq_outputall_mydata $OTHER_PARAMS -s 500 -d 500 --normal-sample-name "$normal_sample_name" $MORE_PARAMS

        echo "Processed sample: $sample_name"
    fi
done