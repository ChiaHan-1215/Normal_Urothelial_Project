### Copy from CGR ChernobylThyroidCancer-RNAseq snakemake file 

```
rule star_index:
    input:
          "/fdb/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa",
          "/fdb/igenomes/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf"
    output:
          "star_index/complete.txt"
    threads: 24
    shell:
          """
          STAR --runThreadN 24 --runMode genomeGenerate --genomeDir star_index --sjdbGTFfile {input[1]} --sjdbOverhang 149 --genomeFastaFiles {input[0]} 2>log/star_index.err
          touch {output} 
          """

rule star_align:
    input:
          "/fdb/igenomes/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf",
          "trimmed/{sample}_filtered_1P.fq.gz",
          "trimmed/{sample}_filtered_2P.fq.gz",
          "star_index/complete.txt"
    output:
          "star_align/{sample}/{sample}Aligned.sortedByCoord.out.bam",
          "star_align/{sample}/{sample}ReadsPerGene.out.tab",
          "star_align/{sample}/{sample}Log.final.out"
    threads: 24
    params:
          index="star_index"
    shell:
          """
          STAR --runThreadN 24 --genomeDir {params.index} --readFilesIn {input[1]} {input[2]} --outFileNamePrefix star_align/{wildcards.sample}/{wildcards.sample} --readFilesCommand zcat --sjdbGTFfile {input[0]} --sjdbOverhang 149 --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate --twopassMode Basic 2>log/{wildcards.sample}_star_align.err
          """

```

#### ChiaHan's script for doing indexing and aligning in CCAD2

```
#!/bin/bash

ml star/2.7.10b 


ANNO=gencode.v39.annotation.gtf
REF=references_Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta 

#echo start STAR
#
# STAR --runThreadN 24 --runMode genomeGenerate\
# --genomeDir star_index --sjdbGTFfile ${ANNO} --sjdbOverhang 149 --genomeFastaFiles ${REF}
#
#echo fin generateding index



# L1 to L8, get the fastq files and save as txt of list of fastq files

for i in {1..8}
do
ls /DCEG/Projects/DataDelivery/Prokunina/TP0325-RS7-Urothelial-Samples-RNA-seq/250410_LH00324_0032_A22NLVTLT3/L${i}/Project_TP0325-RS7/Sample*/*R1*_001.fastq.gz >> fastqR1.txt
done


# Use the txt above to get the file of _R2 and do star alignment

cat fastqR1.txt | while read i
do
sample_name=$(basename "$(dirname "$i")" | cut -d- -f1)
echo "start $sample_name"

fastqR1=`echo $i`
fastqR2=`echo $i | sed 's/_R1_/_R2_/g'`

STAR --runThreadN 24 --genomeDir star_index_hg38_150bp --readFilesIn ${fastqR1} ${fastqR2} --outFileNamePrefix Align_results_hg38/${sample_name}\
 --readFilesCommand zcat --sjdbGTFfile ${ANNO} --sjdbOverhang 149\
 --quantMode TranscriptomeSAM GeneCounts\ # add TranscriptomeSAM to get the trasncriptome bam file for RSEM
 --outSAMtype BAM SortedByCoordinate --twopassMode Basic &&  echo finfin ${sample_name}

done


```
Once the alignment is complete, do the remove duplication and RSEM to calcuate TPM

```
# remove duplication

MarkDuplicates \
    -I /data/star_out/${sample_id}.Aligned.sortedByCoord.out.patched.bam \
    -O Aligned.sortedByCoord.out.patched.md.md.bam \
    -PROGRAM_RECORD_ID null \
    -MAX_RECORDS_IN_RAM 500000 \
    -SORTING_COLLECTION_SIZE_RATIO 0.25 \
    -M ${sample_id}.Aligned.sortedByCoord.out.patched.md.marked_dup_metrics.txt \
    -ASSUME_SORT_ORDER coordinate \
    -TAGGING_POLICY DontTag \
    -OPTICAL_DUPLICATE_PIXEL_DISTANCE 100


# making RSEM index, default is only bowtie2
# can add --star to alos include star index 

rsem-prepare-reference \
  --gtf  /data/gencode.v39.GRCh38.annotation.gtf \
  --star \                       # build STAR index as well
  --num-threads 4 \
  /data/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta \
  /data/rsem_reference



# generate TPM
rsem-calculate-expression \
  --num-threads 4 \
  --fragment-length-max 1000 \
  --no-bam-output \
  --paired-end \
  --estimate-rspd \
  --bam /data/star_out/${sample_id}.Aligned.toTranscriptome.out.bam \
  /data/rsem_reference/rsem_reference \
  /data/${sample_id}.rsem


### This is the example of code:
# NOTE: the reference folder direction need to indicated to ref/reference_name

rsem-calculate-expression \
  --num-threads 4   --fragment-length-max 1000 \
  --no-bam-output   --estimate-rspd \
  --bam SHSY5Y_S03_U3Aligned.toTranscriptome.out.bam \
--paired-end  ../RSEM_Ref/rsem_reference \
 Output.rsem

### This is the example of code START from fastq files:
rsem-calculate-expression --star --paired-end --strandedness reverse --no-bam-output -p 6\
 --star-gzipped-read-file \
 fastqs/SHSY5Y_S01_U1_R1.fastq.gz fastqs/SHSY5Y_S01_U1_R2.fastq.gz \
 ./RSEM_REF/rsem_ref_star ./SY5Y_rsem/test


```
Since using for loop takes longer time to complete whole sample, can use `swarm`. below is the way to create the swarm file based on list of bam file

```
for i in *.bam; do OUTPUT=$(echo $i | sed 's/Aligned.sortedByCoord.out.bam//g'); echo "ml picard/2.26.11 ; picard MarkDuplicates -I ${i} -O ${OUTPUT}.md.bam -PROGRAM_RECORD_ID null -MAX_RECORDS_IN_RAM 500000 -SORTING_COLLECTION_SIZE_RATIO 0.25 -M ${OUTPUT}_metrics.txt -ASSUME_SORT_ORDER coordinate -TAGGING_POLICY DontTag -OPTICAL_DUPLICATE_PIXEL_DISTANCE 100"; done >> mkdup.swarm

# once the mkdup.swarm is created, use swarm command to run it

swarm -t [number cpus] --time [4:00:00] -g [gb for memory] mkdup.swarm

```

### RSEM can start from fastqs or BAM files
### We can also use salmon tool to generate TPM

https://salmon.readthedocs.io/en/latest/salmon.html

#### Salmon instruction
ref:
https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/

- download the fastq file of transcripts and genome sequence, I used the GENCODE V39 HG38

```
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.transcripts.fa.gz

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/GRCh38.primary_assembly.genome.fa.gz

```

- get the chromsome id and remove `>`

```

grep "^>" <(gunzip -c GRCh38.primary_assembly.genome.fa.gz) | cut -d " " -f 1 > decoys.txt

# .bak is generated the backup file for decoy.txt, but it jsut backup  
sed -i.bak -e 's/>//g' decoys.txt

```

- concat the transcripts and genome

```
cat gencode.v39.transcripts.fa.gz GRCh38.primary_assembly.genome.fa.gz > gentrome.fa.gz

```

- now first do `salmon index` and then `salmon quant`

```
salmon index -t gentrome.fa.gz -d decoys.txt -p 12 -i salmon_index --gencode

# quant
salmon quant -i transcripts_index -l salmon_index -1 reads1.fq -2 reads2.fq --validateMappings -o transcripts_quant

```



Once we get the quant.sf file, to make gene TPM we can use `tximport` 


https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html


******

### Below is the Rscipt for handleing mergeing TPMs in isoform and gene level 
### For salmon and RSEM output


```R
library(dplyr)
library(tximport)
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86)

setwd('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/CCLE and other RNA-seq Bam files/SH-SY5Y_short_read_RNA_seq_bams_hg38/salmon_quant_result_for_TPM/')

edb <- EnsDb.Hsapiens.v86

tx2gene <- AnnotationDbi::select(
  edb,
  keys(edb, keytype="TXNAME"),
  columns = c("TXNAME", "GENEID"),
  keytype = "TXNAME"
)  

names(tx2gene) <- c('transcript_id','gene_id','TXID')

lff <- list.files('.',pattern = ".sf")
# like SetName()
names(lff) <- gsub("_salmon.sf$", "", basename(lff))
# check lff to see name match with quant.sf or not

# Just output isoform TPM 
txi_iso <- tximport(
  lff,
  type    = "salmon",
  txIn    = TRUE,      # you’re supplying Salmon outputs
  txOut   = TRUE       # **do not** summarize—keep transcript estimates
)

isoform_tpm <- txi_iso$abundance

# output summary of isofrom TPM into gene TPM

txi <- tximport(
  lff,
  type    = "salmon",
  tx2gene = tx2gene,
  ignoreTxVersion = TRUE
)

gene_tpm <- txi$abundance

# 1. Extract the unique Ensembl gene IDs you have
gene_ids <- rownames(gene_tpm)

gene_map <- ensembldb::select(
  edb,
  keys    = gene_ids,
  keytype = "GENEID",
  columns = c("GENEID","GENENAME")
)

# 3. Turn that into a named vector for lookup
symb <- setNames(gene_map$GENENAME, gene_map$GENEID)

df <- read.delim('SC917163_salmon_quant.sf')

# set name 
rownames(gene_tpm) <- ifelse(is.na(symb[gene_ids]),
                             gene_ids,symb[gene_ids])

gene_tpm <- gene_tpm %>% mutate(GENE=rownames(gene_tpm), .before = SHSY5Y_S01_U1 )
 
```


******


```R


library(dplyr)
library(tximport)
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86)
library(xlsx)


setwd('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/RNA-seq/RNAseq_Normal_Urothelial_05152025/Noraml_rsem/')

edb <- EnsDb.Hsapiens.v86

tx2gene <- AnnotationDbi::select(
  edb,
  keys(edb, keytype="TXNAME"),
  columns = c("TXNAME", "GENEID"),
  keytype = "TXNAME"
)  

names(tx2gene) <- c('transcript_id','gene_id','TXID')

lff <- list.files('.',pattern = ".genes.results")
# like SetName()
names(lff) <- gsub(".genes.results", "", basename(lff))
# check lff to see name match with quant.sf or not

# output summary of isofrom TPM into gene TPM

txi <- tximport(
  lff,
  type    = "rsem",
  txIn = F,
  txOut = F)

gene_tpm <- txi$abundance

# 1. Extract the unique Ensembl gene IDs you have
gene_ids <- rownames(gene_tpm)
gene_ids <- gsub('\\.[0-9]+','',gene_ids)


gene_map <- ensembldb::select(
  edb,
  keys    = gene_ids,
  keytype = "GENEID",
  columns = c("GENEID","GENENAME")
)


# 3. Turn that into a named vector for lookup
symb <- setNames(gene_map$GENENAME, gene_map$GENEID)

gene_tpm <- as.data.frame(gene_tpm)
rownames(gene_tpm) <- gene_ids 

gene_tpm <- gene_tpm %>% mutate(GENEID=rownames(gene_tpm), .before = names(gene_tpm)[1] )

# Keep original Ensembl ID row names
# Add gene symbols as a separate column
gene_tpm <- gene_tpm %>% mutate(GENESYMBOL=ifelse(is.na(symb[gene_ids]), gene_ids, symb[gene_ids]), .before = names(gene_tpm)[1] )


# match names with sampleID

id <- readxl::read_xls('/Volumes/ifs/DCEG/Projects/DataDelivery/Prokunina/TP0325-RS7-Urothelial-Samples-RNA-seq/TP0325-RS7_QC-SUMMARY.xls')
id <- data.frame(id)
id$Vial.Label <- gsub(' ','_',id$Vial.Label)
id$Vial.Label <- gsub("_RNA$","",id$Vial.Label)

# get the name match 

names(gene_tpm) <- gsub('Sample_','',names(gene_tpm))

# Make a named vector for mapping
id_map <- setNames(id$Vial.Label, id$CGR.Sample.ID)

# Replace using mapping
colnames(gene_tpm) <- ifelse(
  colnames(gene_tpm) %in% names(id_map),
  id_map[colnames(gene_tpm)],
  colnames(gene_tpm)   # keep unchanged if not in mapping
)


# Now we save table
write.csv(gene_tpm,'../Sample_rsem_TPM.csv',row.names = F,quote = F)

############################################################
############################################################
############################################################

# For salmon result 
setwd('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/RNA-seq/RNAseq_Normal_Urothelial_05152025/SALMON_Quant_TPM/')

lff <- list.files('.',pattern = '.sf')

names(lff) <- gsub("_salmon_quant.sf", "", lff)


# output summary of isofrom TPM into gene TPM

txi <- tximport(
  lff,
  type    = "salmon",
  tx2gene = tx2gene,
  ignoreTxVersion = TRUE
)

gene_tpm_salmon <- txi$abundance %>% as.data.frame()

# 1. Extract the unique Ensembl gene IDs you have
gene_ids <- rownames(gene_tpm_salmon)
gene_ids <- gsub('\\.[0-9]+','',gene_ids)


gene_map <- ensembldb::select(
  edb,
  keys    = gene_ids,
  keytype = "GENEID",
  columns = c("GENEID","GENENAME")
)

# 3. Turn that into a named vector for lookup
symb <- setNames(gene_map$GENENAME, gene_map$GENEID)

gene_tpm_salmon <- as.data.frame(gene_tpm_salmon)
gene_tpm_salmon <- gene_tpm_salmon %>% mutate(GENEID=rownames(gene_tpm_salmon), .before = names(gene_tpm_salmon)[1] )

# Keep original Ensembl ID row names
# Add gene symbols as a separate column
gene_tpm_salmon <- gene_tpm_salmon %>% mutate(GENESYMBOL=ifelse(is.na(symb[gene_ids]), gene_ids, symb[gene_ids]), .before = names(gene_tpm_salmon)[1] )


# match names with sampleID

id <- readxl::read_xls('/Volumes/ifs/DCEG/Projects/DataDelivery/Prokunina/TP0325-RS7-Urothelial-Samples-RNA-seq/TP0325-RS7_QC-SUMMARY.xls')
id <- data.frame(id)
id$Vial.Label <- gsub(' ','_',id$Vial.Label)
id$Vial.Label <- gsub("_RNA$","",id$Vial.Label)

# Make a named vector for mapping
id_map <- setNames(id$Vial.Label, id$CGR.Sample.ID)

# Replace using mapping
colnames(gene_tpm_salmon) <- ifelse(
  colnames(gene_tpm_salmon) %in% names(id_map),
  id_map[colnames(gene_tpm_salmon)],
  colnames(gene_tpm_salmon)   # keep unchanged if not in mapping
)

# save 
write.csv(gene_tpm_salmon,'../Sample_salmon_TPM.csv',row.names = F,quote = F)




```
