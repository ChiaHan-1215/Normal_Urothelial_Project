# Goal: Counting Splicing Junction Reads for FGFR3 Exon Skipping and Normalization Methods

## Description:

### Reference Method from CHRNA5 Smoking Paper:

The quantification of splice junctions for all CHRNA5 exon 5 isoforms was performed using FeatureCounts (v4.5.0, in R) and a custom .bed file, which provided coordinates for each isoform-specific junction. Only samples with ≥30 total exon 5 splice junction reads were retained for analysis (A549: n=3, UMUC3: n=2; 67–211 junctions per sample); this excluded one UMUC3 replicate with insufficient junction read counts. Isoform ratios were calculated as the fraction of specific junctions over total junction reads. For each cell line, differences in isoform ratios between LT-CSCs and LT-DMSO were compared by unpaired two-sided t tests; p < 0.05 was considered significant.

**FeatureCounts command used by Max:**
```bash
##### The pET01_SJ.saf format: ##########
# Tab-delimited file:
# GeneID	Chr	Start	End	Strand
pET5p	pET01_CHRNA5_Exon5	1290	3593	+
##############################################
# featureCounts v2.0.6 
featureCounts -a pET01_SJ.saf -F SAF -o output input.bam --splitOnly -J
# The count output will have multiple columns that need to be annotated based on exon location
```

**Reference data:**
- Excel file location: `/LTG/Prokunina Lab/Paper Files/2025 CHRNA5 splicing/Submission 111025_Human Genomics/Hogshead Supplementary Tables_HG.xlsx`
- See **Table S4** for rationale and example workflow

### Additional Method: SH-SY5Y Junction Counting in CHRNA5 Region

```bash
# Using FeatureCounts to count junction reads (similar to Sashimi plot in IGV)
# System: ccad2
# Input: SH-SY5Y RNA-seq BAMs
# SAF file format (tab-separated):
# CHRNA5_iso	chr15	78588357	78593275	+

ml slurm/
ml subread/
for i in *.bam
do 
  featureCounts -t 12 -a Featurecount_CHRNA5_isoform/C5_iso.saf \
    -o Featurecount_CHRNA5_isoform/${i%.hg38.md.bam}_fct ${i} \
    --splitOnly -J -F SAF -p --countReadPairs
done
```

**Analysis code:**  
GitHub repository: [SHSY5Y_featurecount_CHRNA5.R](https://github.com/ChiaHan-1215/SH-SY5Y_cell_line_RNA_seq_with_long_read/blob/main/SHSY5Y_featurecount_CHRNA5.R)



---------
---------


## Date: 12/16/2025
## Test RT4 BAM file

```
# region of interest,hg38, FGFR3 exon skipping region
chr4:1801409-1804912

# The saf file:
# FGFR3_region.saf

FGFR3_skip  chr4  1801409  1804912  +

```
The featurecount command: 

```
ml slurm/
ml subread/

 featureCounts -t 12 -a FGFR3_region.saf  -F SAF -o output input.bam --splitOnly -J 

```

The output named `xxx.jcounts` is our input for further analysis in Rstudio

```R
library(dplyr)
library(tidyr)
library(car)


setwd('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/Parse_scRNA-seq/RT4_2D3D_featurecunt/')

df <- read.delim('output.jcounts',header = T)
head(df)

# chr4 and FGFR3 
# filter chr4:1801409-1804912

df.sub <- df %>%
  filter(
    (Site1_chr == 'chr4' & Site1_location >= 1801409 & Site1_location <= 1804912) |
      (Site2_chr == 'chr4' & Site2_location >= 1801409 & Site2_location <= 1804912)
  )


df.sub <- df.sub[,c(1,3,4,7,9:ncol(df.sub))]
#names(df.sub)[5] <- "RT4_2D"
# filter > 1 read 
#df.sub <- df.sub %>% filter(RT4_2D > 1)
# Now have to determine which region is skipling 

# make tag
df.sub$jct_region <- paste0(df.sub$Site1_chr,":",df.sub$Site1_location,"-",df.sub$Site2_location) 

# select to keep
# chr4:1728531-1803335 <- FGFR3_TACC3 fusion?
# chr4:1801536-1801620 <- ex5-ex6

df.sub$jct_tag <- car::recode("")


### Get the exon corrdine ######### 
# Should follow this https://github.com/ChiaHan-1215/Project_mapBC/blob/main/porject_TCGA_BLCA_rs2236786_%20TMEM129/Extract_junction_count_from_Recount3_with_exon_region_Gene_forward_strand.R

library(rtracklayer)
library(dplyr)

gtf <- import('gencode.v39.annotation.gtf')     

# c("FGFR3")

# 5 to 3 ->
genes_of_interest <-c("FGFR3")
# 3 to 5 <-
# genes_of_interest <-c("FAM53A","SLBP","TMEM129")


# Since there's difference of 5' to 3' gene and 3' to 5' gene 
# so need to seperated by gene direct, the +1 -1 base for this code is for FGFR3,TACC3

# Subset exons for your genes
exons <- subset(gtf, type == "exon" & gene_name %in% genes_of_interest)
#exon_df <- as.data.frame(exons)[, c("seqnames", "start", "end", "gene_name", "exon_id", "exon_number")]
exon_df <- as.data.frame(exons)
# just get protein coding exon only
#exon_df <- exon_df[!is.na(exon_df$ccdsid),]
exon_df <- exon_df[, c("seqnames", "start", "end", "gene_name","transcript_name", "exon_id", "exon_number")]
exon_df$location <- paste0(exon_df$seqnames,':',exon_df$start,"-",exon_df$end)

# exon_df.u <- exon_df[!duplicated(exon_df[, c("exon_id", "location")]), ]

iso_list <- exon_df$transcript_name %>% unique() %>% sort()

jct_per_iso <- data.frame()

for (k in iso_list){
  
  # k <- iso_list[2]
  exon_df.u <- exon_df %>% filter(transcript_name == k)
  
  jct_df <- exon_df.u %>%
    arrange(transcript_name, as.numeric(exon_number)) %>%
    group_by(gene_name, transcript_name, seqnames) %>%
    mutate(
      next_start = lead(start),
      next_exon  = lead(exon_number)
    ) %>%
    filter(!is.na(next_start), as.numeric(next_exon) == as.numeric(exon_number) + 1) %>%
    ungroup() %>%    
    transmute(
      gene_name,
      transcript_name,
      jct_label = paste0("jct_ex", exon_number, "_ex", next_exon),
      location = paste0(seqnames, ":", end ,"-", next_start )
    ) 
  jct_per_iso <- rbind(jct_per_iso,jct_df)
  
  
  
  
}


jct_per_iso <- jct_per_iso[!duplicated(jct_per_iso$location),]

# Since the featurecount output kind of -1 with all start and end. need to md based on it 



```

