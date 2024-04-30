#!/bin/bash

set -e

#STEP 1 CHECK QUALITY OF RAW READS

echo "checking for raw reads quality"

fastqc /user/exome_data/all/*fastq.gz -t 4

echo "Done raw reads fastqc!"

echo "performing multiqc on raw reads ...."

multiqc /user/exome_data/all -o /user/exome_data/all

echo "Done raw reads multiqc!"


#STEP 2 TRIMMING

# Define the input and output directories
input_dir="/user/exome_data/all"
output_dir="/user/exome_data/trimmed"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Iterate over your sample pairs
for R1_input in "${input_dir}"/*_R1.fastq.gz; do
    # Check if there are matching R2 reads
    R2_input="${R1_input/_R1.fastq.gz/_R2.fastq.gz}"
    if [ -e "$R2_input" ]; then
        # Define sample name
        sample_name=$(basename "$R1_input" _R1.fastq.gz)

        # Define output file paths
        R1_output="${output_dir}/${sample_name}_R1_paired.fastq.gz"
        R1_unpaired_output="${output_dir}/${sample_name}_R1_unpaired.fastq.gz"
        R2_output="${output_dir}/${sample_name}_R2_paired.fastq.gz"
        R2_unpaired_output="${output_dir}/${sample_name}_R2_unpaired.fastq.gz"

        # Run Trimmomatic
        trimmomatic PE -threads 4 "$R1_input" "$R2_input" \
            "$R1_output" "$R1_unpaired_output" "$R2_output" "$R2_unpaired_output" \
            ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:151
    else
        echo "Warning: No matching R2 read found for $R1_input"
    fi
    echo "Trimming is done"
done


# STEP 3 CHECK QUALITY OF TRIMMED READS

echo "fastqc trimmed reads...."

fastqc /user/exome_data/trimmed/*_paired.fastq.gz -t 4

echo "fastqc trimmed Done!"

echo "multiqc trimmed reads...."

multiqc /user/exome_data/trimmed -o /user/exome_data/trimmed

echo "multiqc for timmed reads Done!"


# STEP 4 MAP READS TO REF
echo "bwa mem index..."
bwa index /user/home/bcancer/ref/hg38.fasta

echo "bwa mem reads mapping..."
find /user/exome_data/trimmed/ -name "*_paired.fastq.gz" | grep -v _R1_paired.fastq.gz | sed 's/_R2_paired.fastq.gz//' | parallel bwa mem -t 4 /user/home/bcancer/ref/hg38.fasta {}_R1_paired.fastq.gz {}_R2_paired.fastq.gz '>' '{}'.sam

echo "reads mapping done!"

echo "going inside trimmed folder..."
cd /user/exome_data/trimmed

echo "samtools convert sam to bam...."

# Iterate over SAM files and convert to BAM in parallel
ls *.sam | parallel '
    # Extract the base name without the file extension
    base_name={.}

    # Check if the corresponding BAM file exists
    if [ ! -f "${base_name}.bam" ]; then
        # Convert SAM to BAM using samtools
        samtools view -@ 4 -b {} > "${base_name}.bam"
    else
        echo "BAM file already exists for ${base_name}. Skipping conversion."
    fi
'

echo "cleaning up the sam files"
rm *.sam

echo "make dir mapped"
mkdir -p /user/exome_data/mapped

echo "move bams to mapped"
mv *.bam /user/exome_data/mapped

echo "get in mapped "
cd /user/exome_data/mapped

echo "samtools sort...."
ls *.bam | parallel 'samtools sort -o {.}_sorted.bam {} -@ 4'

mkdir -p /user/exome_data/sortedbams

mv *_sorted.bam /user/exome_data/sortedbams

echo "samtools index...."
ls *_sorted.bam | parallel 'samtools index {} -@ 2'

echo "we are done...."


#STEP 5 MARK DUPLICATE 
echo "Entering the sortedbams folder"
cd /user/exome_data/sortedbams

echo "Running Picard MarkDuplicates"

# run Picard MarkDuplicates for each BAM file
ls *_sorted.bam | parallel 'picard MarkDuplicates -I {} -O {.}_dupmarked.bam -M {.}_metrics.txt --VALIDATION_STRINGENCY LENIENT 2>&1 | tee {.}_markduplicates.log' 

echo "Duplicate marking completed..."

rm *_sorted.bam

echo "Done and cleaned!"


#STEP 6 ADD READ GROUP
# Define the read group information
RGID="HV3HWDSXY.4"
RGLB="library1"
RGPL="illumina"
RGPU="unit1"

# Specify the input and output directory
input_dir="/user/exome_data/sortedbams"
output_dir="/user/exome_data/rgbams"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Iterate over BAM files in the input directory
for bam_file in "$input_dir"/*_sorted_dupmarked.bam; do
    # Check if the file is a BAM file
    if [ -f "$bam_file" ]; then
        # Extract the sample name from the BAM file name
        sample=$(basename "$bam_file" | cut -d'_' -f1)
        
        # Add read groups using Picard
        picard AddOrReplaceReadGroups \
            I="$bam_file" \
            O="$output_dir/$sample"_sorted_dupmarked_rg.bam \
            RGID="$RGID" \
            RGLB="$RGLB" \
            RGPL="$RGPL" \
            RGPU="$RGPU" \
            RGSM="$sample"
    fi
done

cd ..
rm -rf /user/exome_data/sortedbams

# STEP 7 BASE QUALITY RECALIBRATION
echo "Entering the rgbams folder"
cd /user/exome_data/rgbams

# Define reference genome and known sites for BQSR
#reference=/user/home/bcancer/ref/hg38.fasta
known_sites=/user/home/bcancer/hg38.dbsnp138.vcf

# index the reference using samtools
echo "indexing the reference with samtools"
samtools faidx /user/home/bcancer/ref/hg38.fasta 

# create ref dic using picard
echo "creating reference dic using picard"
picard CreateSequenceDictionary -R /user/home/bcancer/ref/GCF_000001405.40_GRCh38.p14_genomic.fna -O /user/home/bcancer/ref/GCF_000001405.40_GRCh38.p14_genomic.dict 

# Define GATK java options
java_options="-DGATK_STACKTRACE_ON_USER_EXCEPTION=true"

# Step 1: Base Quality Score Recalibration (BQSR)
echo "Running GATK BaseRecalibrator for BQSR"
ls *.bam | parallel "gatk --java-options '$java_options' BaseRecalibrator -R $reference -I {} --known-sites $known_sites -O {.}_recal_data.table" 

# Step 2: Apply BQSR
echo "Applying BQSR with GATK ApplyBQSR"
ls *_sorted_dupmarked_rg.bam | parallel "gatk --java-options '$java_options' ApplyBQSR -R $reference -I {} --bqsr-recal-file {.}_recal_data.table -O {.}_recal.bam" 


echo "Base Quality Score Recalibration completed..."
echo "Done!"

rm *_rg.bam

# STEP 8 METRICS COLLECTION FROM BAM FILES
echo "Entering the rgbams folder"
cd /user/exome_data/rgbams

echo "Running CollectAlignmentSummaryMetrics"

# Start a background process to run Picard MarkDuplicates for each BAM file
ls *_sorted_dupmarked_rg_recal.bam | parallel 'gatk CollectAlignmentSummaryMetrics -I {} -R /user/home/bcancer/ref/hg38.fasta -O {.}_alignmetrics.txt' 

echo "collecting alignmentmetrics completed..."

echo "Running CollectInsertSizeMetrics"

# Start a background process to run Picard MarkDuplicates for each BAM file
ls *_sorted_dupmarked_rg_recal.bam | parallel 'gatk CollectInsertSizeMetrics -I {} -O {.}_insertmetrics.txt -H /user/exome_data/rgbams/insert_size_histogram.pdf' 

echo "collecting insert metrics completed..."

echo "Performing Multi quality check"

# Create the output directory if it doesn't exist
mkdir -p /user/exome_data/all/trimmed/mapped/rgbams/multiqc_report/

# Run MultiQC
multiqc /user/exome_data/rgbams/ -o /user/exome_data/rgbams

echo "Multiqc successfully completed"

echo "Done!"


# STEP 9 VARIANT CALLING
#STEP 9-1 CALL NORMAL VARIANTS
cd /user/exome_data/rgbams

#reference_genome="/user/home/bcancer/ref/hg38.fasta"

# Output directory for VCF files
output_dir="/user/exome_data/nomvariants"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Process all BAM files in the current directory
for normal_bam in *_sorted_dupmarked_rg_recal.bam; do
    # Extract the sample name from the BAM file (e.g., "204_sorted_dupmarked_rg_recal.bam" -> "204")
    sample_name="${normal_bam%%_*}"

    # Check if the sample name is an even number composed of three digits
    if [[ $sample_name =~ ^[0-9]{3}$ && $((sample_name % 2)) -eq 0 ]]; then
        # Normal sample found
        normal_sample="$normal_bam"
        
        # Run Mutect2 for the normal sample
        output_vcf="${output_dir}/${sample_name}_normal.vcf.gz"
        gatk --java-options "-Xms16G -Xmx20G" Mutect2 \
             -R "$reference" \
             -I "$normal_sample" \
             -max-mnp-distance 0 \
             -O "$output_vcf"
    fi
done

echo "Somatic variant calling for normal samples completed."


#STEP 9-2 CREATE GENOMICDB FOR NORMAL BAM

#prep exome interval file
curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz" |\
    gunzip -c | cut -f 3,5,6 | sort -t $'\t' -k1,1 -k2,2n | bedtools merge -i - > /user/home/bcancer/exome.bed 
grep '>' hg38.fasta | awk '{print $1}' > chroms
awk 'NR==FNR{a[$1]; next} $1 in a' chroms exome.bed > exome_new.bed

# Define the paths and filenames
#reference="/user/home/bcancer/ref/hg38.fasta"
intervals_file="/user/home/bcancer/exome_new.bed"
workspace_path="/user/exome_data/pon_db"

    
# Path to the directory containing normal VCF files
normal_vcf_dir="/user/exome_data/nomvariants"

echo "Creating GenomicsDB.."

# Create the GenomicsDB workspace
gatk --java-options "-Xms10G -Xmx15G" GenomicsDBImport \
    -R "$reference" \
    -L "$intervals_file" \
    --merge-input-intervals true \
    --genomicsdb-workspace-path "${workspace_path}" \
    -V "${normal_vcf_dir}/202_normal.vcf.gz" \
    -V "${normal_vcf_dir}/204_normal.vcf.gz" \
    -V "${normal_vcf_dir}/206_normal.vcf.gz" \
    -V "${normal_vcf_dir}/208_normal.vcf.gz" \
    -V "${normal_vcf_dir}/210_normal.vcf.gz" \
    -V "${normal_vcf_dir}/212_normal.vcf.gz" \
    -V "${normal_vcf_dir}/214_normal.vcf.gz" \
    -V "${normal_vcf_dir}/216_normal.vcf.gz" \
    -V "${normal_vcf_dir}/218_normal.vcf.gz" \
    -V "${normal_vcf_dir}/220_normal.vcf.gz" \
    -V "${normal_vcf_dir}/222_normal.vcf.gz" \
    -V "${normal_vcf_dir}/224_normal.vcf.gz" \
    -V "${normal_vcf_dir}/226_normal.vcf.gz" \
    -V "${normal_vcf_dir}/228_normal.vcf.gz" \
    -V "${normal_vcf_dir}/230_normal.vcf.gz" \
    -V "${normal_vcf_dir}/232_normal.vcf.gz" \
    -V "${normal_vcf_dir}/234_normal.vcf.gz" \
    -V "${normal_vcf_dir}/236_normal.vcf.gz" \
    -V "${normal_vcf_dir}/238_normal.vcf.gz" \
    -V "${normal_vcf_dir}/240_normal.vcf.gz" \
    -V "${normal_vcf_dir}/242_normal.vcf.gz" \
    -V "${normal_vcf_dir}/244_normal.vcf.gz" \
    -V "${normal_vcf_dir}/246_normal.vcf.gz" 

echo "GenomicsDBImport completed."


#STEP 9-3 CREATE PANEL OF NORMALS
# Reference genome file
#reference_genome="/user/home/bcancer/ref/hg38.fasta"

# GenomicsDB workspace path (path to your GenomicsDB workspace for normal samples)
genomicsdb_workspace_path="/user/exome_data/all/trimmed/mapped/rgbams/variants/pon_db"

# Output VCF file for the Panel of Normals (PoN)

output_pon_vcf="/user/exome_data/nomvariants/pon.vcf.gz"

echo "Creating somatic panel of normals"

# Run GATK CreateSomaticPanelOfNormals
gatk --java-options "-Xms16G -Xmx20G" CreateSomaticPanelOfNormals \
     -R "$reference" \
     -V gendb://"${workspace_path}" \
     -O "$output_pon_vcf"

echo "Panel of Normals (PoN) creation completed."


#STEP 9-4 VARIANT CALLING
echo "Calling variants from tumor-normal samples using mutect2.."
cd /user/exome_data/rgbams

# Reference genome file (replace with your reference genome)
#reference_genome="/user/home/bcancer/ref/hg38.fasta"

# Output directory for VCF files
output_dir="/user/exome_data/variants/"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

#germline_resource="/user/home/bcancer/af-only-gnomad.vcf.gz"

panel_normal="/user/exome_data/nomvariants/pon.vcf.gz"

# Process all BAM files in the current directory
for tumor_bam in *_sorted_dupmarked_rg_recal.bam; do
    # Extract the sample name from the BAM file (e.g., "203_sorted_dupmarked_rg_recal.bam" -> "203")
    sample_name="${tumor_bam%%_*}"

    # Check if the sample name is an odd number composed of three digits
    if [[ $sample_name =~ ^[0-9]{3}$ && $((sample_name % 2)) -ne 0 ]]; then
        # Tumor sample found
        tumor_sample="$tumor_bam"
        
        # Deduce the corresponding normal sample name (e.g., "203" -> "204")
        normal_sample="$((sample_name + 1))_sorted_dupmarked_rg_recal.bam"
        
        # Check if the deduced normal sample file exists
        if [ -e "$normal_sample" ]; then
            # Run Mutect2 for the tumor-normal pair
            output_vcf="${output_dir}/${sample_name}_somatic.vcf.gz"
            gatk --java-options "-Xms16G -Xmx20G" Mutect2 \
                 -R "$reference" \
                 -I "$tumor_sample" \
                 -I "$normal_sample" \
                 --panel-of-normals "$panel_normal" \
                 -O "$output_vcf"
        else
            echo "Corresponding normal sample not found for tumor sample $sample_name."
        fi
    fi
done

echo "Somatic variant calling completed."


# STEP 10 VARIANT NORMALIZATION

echo "normalizing variants.."
# Specify the directory containing the input VCF files
input_dir="/user/exome_data/variants/"

# Specify the reference genome in FASTA format
#reference_genome="/user/home/bcancer/ref/hg38.fasta"

# Create an output directory for normalized VCF files
output_dir="/user/exome_data/normalized/"
mkdir -p "$output_dir"

# Iterate over all VCF files in the input directory
for input_vcf in "$input_dir"/*somatic.vcf.gz; do
    if [ -f "$input_vcf" ]; then
        # Extract the sample name from the VCF file (e.g., "sample_somatic.vcf.gz" -> "sample_somatic")
        sample_name=$(basename "$input_vcf" | sed 's/\.vcf\.gz//')

        # Specify the output VCF file name
        output_vcf="$output_dir/${sample_name}_normalized.vcf.gz"

        # Normalize the VCF using vt
        vt normalize "$input_vcf" -r "$reference" -n -o "$output_vcf"

        # Check if the normalization process was successful
        if [ $? -eq 0 ]; then
            echo "Normalization completed for $input_vcf. Output saved as $output_vcf"
        else
            echo "Normalization failed for $input_vcf."
        fi
    fi
done

#echo "Normalization of VCF files in $input_dir completed."


#STEP 11 FILTER SNPs

echo "filtering varants for SNP and INDELs and for depth"
# Define paths and filenames
#reference_genome="/user/home/bcancer/ref/hg38.fasta"
input_dir="/user/exome_data/normalized"  # Directory containing your VCF files
output_dir="/user/exome_data/filteredvars"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through all VCF files in the input directory
for input_vcf in "$input_dir"/*_somatic_normalized.vcf.gz; do
    if [ -f "$input_vcf" ]; then
        # Extract the sample name from the VCF file (e.g., "203_somatic_normalized.vcf.gz" -> "203")
        sample_name=$(basename "$input_vcf" | sed 's/_somatic_normalized.vcf.gz//')

        # Index the original normalized VCF file
        gatk IndexFeatureFile -I "$input_vcf"

        # Create a filtered VCF for SNPs
        gatk SelectVariants \
            -R "$reference" \
            -V "$input_vcf" \
            -select-type SNP \
            -select-type INDEL \
            -O "$output_dir/${sample_name}_variants.vcf.gz"

        # Index the SNP VCF
        gatk IndexFeatureFile -I "$output_dir/${sample_name}_variants.vcf.gz"

        # Perform VariantFiltration on the SNP VCF
        gatk VariantFiltration \
            -R "$reference" \
            -V "$output_dir/${sample_name}_variants.vcf.gz" \
            --filter-expression "DP < 10" \
            --filter-name "LowQD" \
            -O "$output_dir/${sample_name}_filtered_variants.vcf.gz"
    fi
done

echo "SNP filtration completed."


#STEP 12 VCF ANNOTATION

echo "variant annotation using vep.."

#ANNOATE THE VCFs
# Define paths and parameters
gff_file="/user/home/bcancer/ref/Homo_sapiens.GRCh38.110.gff3.gz"
fasta_file="/user/home/bcancer/ref/hg38.fasta.gz"
input_dir="/user/exome_data/filteredvars"
output_dir="/user/exome_data/annotated"
plugin_dir="/user/home/bcancer/software/VEP_plugins"

# Ensure the output directory exists
mkdir -p "$output_dir"

# Loop through VCF files in the input directory
for input_vcf in "$input_dir"/*_filtered_variants.vcf.gz; do
    # Extract the sample name from the VCF file (e.g., "203_filtered_snps.vcf.gz" -> "203")
    sample_name=$(basename "$input_vcf" | sed 's/_filtered_variants.vcf.gz//')

    # Define the output VCF path for this sample
    output_vcf="$output_dir/${sample_name}_annotated.vcf.gz"

    # Run VEP with the provided options
    vep --hgvs --fasta "$fasta_file" \
    	--force_overwrite \
    	--gff "$gff_file" \
    	-i "$input_vcf" \
    	--format vcf \
        -o "$output_vcf" \
        --vcf \
        --symbol --terms SO --tsl --biotype \
        --plugin Downstream --plugin Frameshift --plugin Wildtype \
        --coding_only \
        --no_intergenic \
        --pick --pick_allele \
        --dir_plugins "$plugin_dir"
        
        
    echo "VEP annotation completed for $sample_name. Annotated VCF is saved as $output_vcf."
done


#STEP 13 ID HLA ALLELES FROM NORMAL SAMPLES

echo "preparing files for hla typing using hla-hd.."
#Decompress fastq.gz
# Source directory where your files are located
source_dir="/user/exome_data/trimmed"

# Target directory where you want to copy and decompress the files
target_dir="/user/home/bcancer/software/hlahd.1.7.0/data"
# Loop through files in the source directory
for r1_file in "$source_dir"/*_R1_paired.fastq.gz; do
    if [ -f "$r1_file" ]; then
        # Extract the sample ID from the filename (e.g., "203_R1.fastq.zip" -> "203")
        sample_id=$(basename "$r1_file" | cut -d '_' -f 1)
        
        # Check if the sample ID starts with an even number (first three digits)
        if [[ "$sample_id" =~ ^[0-9]{3}$ && $((sample_id % 2)) -eq 0 ]]; then
            # Find the corresponding R2 file
            r2_file="$source_dir/${sample_id}_R2_paired.fastq.gz"

            # Check if the corresponding R2 file exists
            if [ -f "$r2_file" ]; then
                # Copy both R1 and R2 files to the target directory
                cp "$r1_file" "$target_dir"
                cp "$r2_file" "$target_dir"
            fi
        fi
    fi
done

# Decompress all files in the target directory
for compressed_file in "$target_dir"/*.gz; do
    if [ -f "$compressed_file" ]; then
        zcat "$compressed_file" > "${compressed_file%.gz}"  # Decompress and remove .gz extension
    fi
done

echo "File copying, decompression completed."


echo "hlahd starting.."
# hlahd just because I had already saved the fastq files in its data folder after it not working
#cd /user/home/bcancer/software/hlahd.1.7.0

# Define paths and parameters
#hlahd_script="./bin/hlahd.sh"
freq_data="/user/home/bcancer/software/hlahd.1.7.0/freq_data"
data_dir="/user/home/bcancer/software/hlahd.1.7.0/data"
gene_split_file="/user/home/bcancer/software/hlahd.1.7.0/HLA_gene.split.txt"
dictionary_dir="/user/home/bcancer/software/hlahd.1.7.0/dictionary"
estimation_dir="/user/exome_data/estimation"

#create extimation directory if it doesnt exist

mkdir -p "$estimation_dir"

# Number of threads to use
threads=4

# Loop through the sample data directory
for sample_R1 in "$data_dir"/*_R1_paired.fastq; do
    if [ -f "$sample_R1" ]; then
        # Extract the sample name from the R1 file name
        sample_name=$(basename "$sample_R1" | sed 's/_R1_paired.fastq//')
        
        # Define input R1 and R2 file paths
        input_R1="$sample_R1"
        input_R2="$data_dir/${sample_name}_R1_paired.fastq"

        # Run hlahd.sh v1.2.1 for the current sample
        hlahd.sh -t "$threads" -m 100 -f "$freq_data" "$input_R1" "$input_R2" \
            "$gene_split_file" "$dictionary_dir" "$sample_name" "$estimation_dir"

        echo "HLA-HD processing completed for $sample_name."
    fi
done

echo "hla_hd completed."

echo "wesprocessor completed successfully"