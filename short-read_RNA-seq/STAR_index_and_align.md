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


# making RSEM index

rsem-prepare-reference \
        /data/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta \
        /data/rsem_reference/rsem_reference \
        --gtf /data/gencode.v39.GRCh38.annotation.gtf \
        --num-threads 4"


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
rsem-calculate-expression \
  --num-threads 4   --fragment-length-max 1000 \
  --no-bam-output   --estimate-rspd \
  --bam SHSY5Y_S03_U3Aligned.toTranscriptome.out.bam \
--paired-end  ../RSEM_Ref/rsem_reference Output.rsem

```
Since using for loop takes longer time to complete whole sample, can use `swarm`. below is the way to create the swarm file based on list of bam file

```
for i in *.bam; do OUTPUT=$(echo $i | sed 's/Aligned.sortedByCoord.out.bam//g'); echo "ml picard/2.26.11 ; picard MarkDuplicates -I ${i} -O ${OUTPUT}.md.bam -PROGRAM_RECORD_ID null -MAX_RECORDS_IN_RAM 500000 -SORTING_COLLECTION_SIZE_RATIO 0.25 -M ${OUTPUT}_metrics.txt -ASSUME_SORT_ORDER coordinate -TAGGING_POLICY DontTag -OPTICAL_DUPLICATE_PIXEL_DISTANCE 100"; done >> mkdup.swarm

# once the mkdup.swarm is created, use swarm command to run it

swarm -t [number cpus] --time [4:00:00] -g [gb for memory] mkdup.swarm

```


### Since the RSEM reuqired aligned to transcriptome sequence, can use salmon tool to generate TPM

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

- now first do `salmon index`

```
salmon index -t gentrome.fa.gz -d decoys.txt -p 12 -i salmon_index --gencode

```

Once we get the quant.sf file, to make gene TPM we can use `tximport` 


https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html


