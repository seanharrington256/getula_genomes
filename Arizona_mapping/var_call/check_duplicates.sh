#!/bin/bash

# Check if a VCF file was provided as an argument
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <vcf_file>"
    exit 1
fi

# Define the input VCF file from the command line argument
VCF_FILE_INDELS=$1

VCF_FILE="NoIndel$VCF_FILE_INDELS"

# Create a file with no indels
vcftools --gzvcf $VCF_FILE_INDELS --remove-indels --recode --stdout | bgzip -c > $VCF_FILE


# index the file:
bcftools index $VCF_FILE

# Extract the base name of the VCF file (without extension)
PREFIX=$(basename $VCF_FILE .vcf.gz)

# # Define output file prefixes
# PREFIX="DUPCHK_${BASE_NAME}"

# Step 1: Extract snps and invariant sites, but not indels
# bcftools query -f '%CHROM\t%POS\n' $VCF_FILE > ${PREFIX}_all_positions.txt

bcftools query -f '%CHROM\t%POS\n' -e 'TYPE="indel"' $VCF_FILE > ${PREFIX}_all_positions.txt



# Step 2: Find duplicates
sort ${PREFIX}_all_positions.txt | uniq -d > ${PREFIX}_duplicates_positions.txt

# # Step 3: Create a list of duplicate positions
# awk '{print $1" "$2}' ${PREFIX}_duplicates_positions.txt > ${PREFIX}_duplicates_regions.txt

# Step 4: Extract and print duplicate positions
bcftools view -R ${PREFIX}_duplicates_positions.txt $VCF_FILE -Oz -o ${PREFIX}_extracted_duplicates.vcf.gz
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' ${PREFIX}_extracted_duplicates.vcf.gz > ${PREFIX}_extracted_positions.txt
