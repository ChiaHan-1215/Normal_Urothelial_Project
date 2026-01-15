# FGFR3 Variant Annotation Pipeline (hg38)

This pipeline handles the downloading of references, standardization of chromosome naming, annotation of sample VCFs using dbSNP, and the extraction of the FGFR3 gene region.

---

### Full Pipeline Code

```bash
# ==========================================
# 1. REFERENCE & DBSNP PREPARATION
# ==========================================

# Download hg38 reference genome
wget [https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz)

# Download dbSNP VCF and index from NCBI
wget [https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz](https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz)
wget [https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz.tbi](https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz.tbi)

# Download the assembly report to create a chromosome mapping file
wget [https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_report.txt](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_report.txt)

# Create the mapping file (local_rename.txt) to convert NCBI names to UCSC 'chr' names
grep -v "^#" GCF_000001405.40_GRCh38.p14_assembly_report.txt | \
awk -F'\t' '$2 == "assembled-molecule" {print $7, $10}' > local_rename.txt

# Rename chromosomes in the dbSNP VCF and index it
bcftools annotate --rename-chrs local_rename.txt GCF_000001405.40.gz -Oz -o dbsnp_hg38_chr.vcf.gz
bcftools index -t dbsnp_hg38_chr.vcf.gz

# ==========================================
# 2. SAMPLE VCF PROCESSING
# ==========================================

# Export PLINK binary files to VCF
plink2 --bfile ../samples \
  --fa hg38.fa.gz \
  --chr 1-22, X, Y \
  --set-all-var-ids '@:#' \
  --new-id-max-allele-len 50 \
  --export vcf bgz \
  --out step1_raw

# Rename sample chromosomes to match UCSC 'chr' format
bcftools annotate --rename-chrs local_rename.txt -Oz -o step2_chr.vcf.gz step1_raw.vcf.gz
bcftools index -t step2_chr.vcf.gz

# Coordinate-sort the sample VCF
bcftools sort step2_chr.vcf.gz -Oz -o step3_sorted.vcf.gz
bcftools index -t step3_sorted.vcf.gz

# ==========================================
# 3. ANNOTATION & REFINEMENT
# ==========================================

# Annotate with rsIDs from the standardized dbSNP file
bcftools annotate -a dbsnp_hg38_chr.vcf.gz -c ID -Oz -o step4_annotated.vcf.gz step3_sorted.vcf.gz

# Set custom IDs in the format 'rsID:REF:ALT'
bcftools annotate --set-id '%ID:%REF:%ALT' -Oz -o final_confluence_hg38.vcf.gz step4_annotated.vcf.gz
bcftools index -t final_confluence_hg38.vcf.gz

# ==========================================
# 4. GENE EXTRACTION (FGFR3)
# ==========================================

# Extract the FGFR3 region based on coordinates
bcftools view -r chr4:1707662-1816636 -Oz -o FGFR3.sub.vcf.gz final_confluence_hg38.vcf.gz

# Normalize and remove duplicate variants
bcftools norm -d all -Oz -o FGFR3.unique.vcf.gz FGFR3.sub.vcf.gz
bcftools index -t FGFR3.unique.vcf.gz

echo "Process Complete. Final output: FGFR3.unique.vcf.gz"
