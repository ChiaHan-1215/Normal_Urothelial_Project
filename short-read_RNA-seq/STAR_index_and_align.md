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
 --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate --twopassMode Basic &&  echo finfin ${sample_name}

done


```

