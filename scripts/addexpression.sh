#!/bin/bash

set -e

echo "preparing expression files"
#Define paths
output_dir="/user/rna_seq"
reference_fa="/user/home/freedom/snakemake/rnaref/Homo_sapiens.GRCh38.cdna.all.fa.gz"
tx2gene_csv="$output_dir/tx2gene.Homo_sapiens.GRCh38.cdna.csv"

#Create tx2gene file
zcat "$reference_fa" | grep '>' | awk '{FS= " "}BEGIN{ print "TXNAME,GENEID"};{print substr($1,2) "," substr($4,6)};' > "$tx2gene_csv"

# Run Rscript for tximport
Rscript - <<RSCRIPT
library(tximport)
library(rhdf5)

#Read the annotation file
tx2gene <- read.csv("$tx2gene_csv")

#Path to abundance files
setwd("/user/rna_seq")
files <- file.path("kallisto_quants", list.files("kallisto_quants"), "abundance.tsv")
names(files) <- list.files("kallisto_quants")

#Gene level counts/abundance
txi.kallisto.g <- tximport(files, type = "kallisto", tx2gene = tx2gene)

#Extract
abundance_values <- (txi.kallisto.g$abundance)

#Add a header for the first column
#abundance_values <- cbind("gene_name" = rownames(abundance_values), abundance_values)

#Write only the abundance values to the output file
write.table(abundance_values, file = "$output_dir/gene_expression.tsv", quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)
RSCRIPT

awk 'BEGIN {OFS="\t"; print "gene_name\t201\t202\t203\t204\t205\t206\t207\t208\t209\t210\t211\t212\t213\t214\t215\t216\t217\t218\t219\t220\t221\t222\t223\t224\t225\t226\t227\t228\t229\t230\t231\t232\t233\t234\t235\t236\t237\t238\t239\t240\t241\t242\t243\t244\t245\t246" } NR>1 {gsub("abundance.", "", $1); print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13 "\t" $14 "\t" $15 "\t" $16 "\t" $17 "\t" $18 "\t" $19 "\t" $20 "\t" $21 "\t" $22 "\t" $23 "\t" $24 "\t" $25 "\t" $26 "\t" $27 "\t" $28 "\t" $29 "\t" $30 "\t" $31 "\t" $32 "\t" $33 "\t" $34 "\t" $35 "\t" $36 "\t" $37 "\t" $38 "\t" $39 "\t" $40 "\t" $41 "\t" $42 "\t" $43 "\t" $44 "\t" $45 "\t" $46 "\t" $47}' $output_dir/gene_expression.tsv > $output_dir/modified_gene_expression.tsv

echo "Done preparing expression file" 

echo "adding expression to vcf...."

#Define the folder containing your *.vcf.gz files
vcf_folder="/user/exome_data/annotated"
#Define the paths to your VCF files folder and expression data folder
expression_file="/user/rna_seq/modified_gene_expression.tsv"
#Define the output folder
output_folder="/user/exome_data/annotated"

#Iterate over all *.vcf.gz files in the folder
for gz_file in "$vcf_folder"/*.vcf.gz; do
    if [ -f "$gz_file" ]; then
        #Extract the base name without the extension
        base_name=$(basename "$gz_file" .vcf.gz)

        #Use cat to decompress and concatenate to a new .vcf file
        cat "$gz_file" > "$vcf_folder/$base_name.vcf"

        echo "Processed: $gz_file"
    fi
done

#Iterate over VCF files in the specified folder
for vcf_file in "$vcf_folder"/*.vcf; do
    if [ -f "$vcf_file" ]; then
        #Extract sample name from VCF file name
        sample_name=$(basename "$vcf_file" | sed 's/_annotated.vcf//')
        
        #Path to the output annotated VCF file
        output_vcf="$output_folder/${sample_name}_annotated.gx.vcf"

        #Run vcf-expression-annotator
        vcf-expression-annotator "$vcf_file" "$expression_file" -s "$sample_name" custom gene --id-column gene_name --expression-column "$sample_name" -o "$output_vcf" --ignore-ensembl-id-version

        echo "Expression added to the VCF for $sample_name"
    fi
done