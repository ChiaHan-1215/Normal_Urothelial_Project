
# =============================================================================
# FGFR3 VARIANT ANNOTATION PIPELINE (hg38)
# ==============================================================================

# --- PHASE 1: REFERENCE & DBSNP PREPARATION ---

# 1. Download hg38 Reference Genome
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

# 2. Download dbSNP Raw Data and Index
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz.tbi

# 3. Download the Assembly Report
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_report.txt

# 4. Generate the Chromosome Mapping file (local_rename.txt)
# This maps NCBI RefSeq IDs (Column 7) to UCSC Names (Column 10)
grep -v "^#" GCF_000001405.40_GRCh38.p14_assembly_report.txt | \
awk -F'\t' '$2 == "assembled-molecule" {print $7, $10}' > local_rename.txt

# 5. Standardize dbSNP Chromosome Names
bcftools annotate --rename-chrs local_rename.txt \
    GCF_000001405.40.gz \
    -Oz -o dbsnp_hg38_chr.vcf.gz

bcftools index -t dbsnp_hg38_chr.vcf.gz


# --- PHASE 2: SAMPLE VCF CONVERSION & FORMATTING ---

# 6. Convert PLINK binary files to VCF
plink2 --bfile ../samples \
    --fa hg38.fa.gz \
    --chr 1-22, X, Y \
    --set-all-var-ids '@:#' \
    --new-id-max-allele-len 50 \
    --export vcf bgz \
    --out step1_raw

# 7. Rename Sample Chromosomes to UCSC style (chr1, chr2...)
bcftools annotate --rename-chrs local_rename.txt \
    -Oz -o step2_chr.vcf.gz \
    step1_raw.vcf.gz

# 8. Sort and Index the Sample VCF
bcftools sort step2_chr.vcf.gz -Oz -o step3_sorted.vcf.gz
bcftools index -t step3_sorted.vcf.gz


# --- PHASE 3: ANNOTATION ---

# 9. Annotate IDs from dbSNP
bcftools annotate -a dbsnp_hg38_chr.vcf.gz \
    -c ID \
    -Oz -o step4_annotated.vcf.gz \
    step3_sorted.vcf.gz

# 10. Set custom ID format (rsID:REF:ALT)
bcftools annotate --set-id '%ID:%REF:%ALT' \
    -Oz -o final_confluence_hg38.vcf.gz \
    step4_annotated.vcf.gz

bcftools index -t final_confluence_hg38.vcf.gz


# --- PHASE 4: FGFR3 GENE EXTRACTION ---

# 11. Extract FGFR3 Coordinates
bcftools view -r chr4:1707662-1816636 \
    -Oz -o FGFR3.sub.vcf.gz \
    final_confluence_hg38.vcf.gz

# 12. Normalize and De-duplicate Variants
bcftools norm -d all \
    -Oz -o FGFR3.unique.vcf.gz \
    FGFR3.sub.vcf.gz

bcftools index -t FGFR3.unique.vcf.gz

echo "Pipeline Complete. Final file: FGFR3.unique.vcf.gz"
