### Copy from CGR RANseq snakemake

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
