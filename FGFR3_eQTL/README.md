### Date: 01142025

- For FGFR3, mark wuch isofom is main one and which is skipping exons

- GTEx:

  - Run GWAS marker for FGFR3 first

  - To explore FGFR3 bulk epression across tissue, plot is below:  https://gtexportal.org/home/gene/ZNF320/geneExpressionTab

  - To explre FGFR3 iso-expression for 77 blader vs smoth muscle, make a plot, extract iso TPM for these group, make plot

  - Finalized table of GTEx 77 bladder tissue isoform normalized value with sex age race with GT FGFR3-tacc3 region



- Parse Normal Urth:
Run GWAS marker for FGFR3 first

		FGFR3:
		- GTEx compare FGFR3 isoform expression plot to see if similiar GTEx 
		- Finalized table of bladder tissue isoform normalized value with sex race with GT FGFR3-tacc3 region


		Alos for GWAS wide gene, consult Oscar to select TAD region and GWAS








### Date: 12232025

### Goal: The FGFR3 expression with SNPs eQTL

### Data and steps: 
Use the sample RNA-seq TPM with Genotype vcf file

1. generate TPM files from RNA-seq using salmon/rsem etc

2. Do quantile normalization of TPM file

```R

# using Deseq2 to output normlaized count
library(DESeq2)
cts <- as.matrix(raw_genecount)
# Optional: ensure integers
storage.mode(cts) <- "integer"
coldata <- DataFrame(row.names = colnames(cts))  # zero-column DataFrame with sample IDs

colData <- data.frame(
  sample_id = colnames(cts),
  row.names = colnames(cts)
)

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design = ~ 1) # ~ 1 means no variables in the design


keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

dds <- estimateSizeFactors(dds)                 # or DESeq(dds) if you'll do DE later
norm_counts <- counts(dds, normalized = TRUE) 
norm_counts
rownames(norm_counts) <- geneNames[rownames(norm_counts)]
norm_counts <- as.data.frame(norm_counts)

norm_counts$gene_symbol <- rownames(norm_counts)
norm_counts <- norm_counts[,c(349,1:348)]


#####################Do z score on Oscar's code ##############################

norm_counts <- norm_counts[,-1]
str(norm_counts)

# The file format should look like this:

#TCGA.2F.A9KO TCGA.2F.A9KP TCGA.2F.A9KQ TCGA.2F.A9KR TCGA.2F.A9KT
# TACC3         4.9714       5.7981       4.9421       7.0476       6.9486
# FGFR3         6.0588       7.7760       7.4456       8.8113       8.9363
# SLBP          5.3292       6.1125       5.6256       5.8066       6.8386
# TMEM129       3.7971       5.4871       5.4430       6.0605       5.9451
# FAM53A        0.2522       1.7702       0.2400       0.8246       1.8036

# Need to matched with data format 

head(norm_counts[1:5,1:5])

qnorm.DeseqNC.set <- t(apply(norm_counts, 1, rank, ties.method = "average"))
qnorm.DeseqNC.set <- qnorm(qnorm.DeseqNC.set / (ncol(qnorm.DeseqNC.set) + 1))
qnorm.DeseqNC.set[1:5,1:5]
qnorm.DeseqNC.set <- as.data.frame(qnorm.DeseqNC.set)
qnorm.DeseqNC.set$Gene_symbol <- rownames(qnorm.DeseqNC.set)
qnorm.DeseqNC.set <- qnorm.DeseqNC.set[,c(349,1:348)]
#write.table(qnorm.DeseqNC.set,'IMvigor_gene_zscore_from_Deseq2_normalizedcount_usingOscarMethod.csv',col.names = T,row.names = F,quote = F,sep = ',')


```

3. vcf file from Genotype folder using plink

4. extract FGFR3-TACC3 region GT of interest.

```
# In CCAD2

ml slurm/
ml bcftools/
bcftools view -r chr4:1718666-1810880 NU_sample.vcf.gz -Oz -o FGFR3_sub_NU_sample.vcf.gz
```
   
6. use lm model in R to test the unadjust (SNP~TPM) and adjust (SNP~TPM + RACE + SEX)


### The eQTL analysis can be found also in Note "model analysis code and writing up for Github readme"

   
