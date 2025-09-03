### Goal: 

Using sc long-read data to detect potential isofrom and intergrated with short-read sc-RNA-seq data



### Data: 

Fastq files: 

    - since fastq files is too large for split-pipe, we use `seqkit` to split fastqs
    - script folder: `/Volumes/data/parse_single_cell/Long_read_pacbio_with_parse_test/ALL_32SMRT_cell_pacbio_files/Fastqs`

    
------

#### Approachs:

- Using Pasre own script to clean, align fastqs

- the source code example is in the folder `LR_generate_pairs.py`

```
# running example

cd B02.part_110 ; source myconda ; conda activate parse_v130 ; python /data/leec20/parse_single_cell/Long_read_pacbio_with_parse_test/Parse_Pacbio_PIPELINE/LR_generate_pairs.py --out_dir ./B02.part_110_OUT/ --chemistry v2 --fastq ./B02.part_110.fastq.gz --multiple_fq --new_fname P02split --l1dist 3 --l2dist 2 && pigz -p 12 ./B02.part_110_OUT/*.fastq && echo finned

```

- files, as we have B0{1..4} sub-library fastqs:
      - The cleaned fastqs of B01,B02, `/Volumes/data/parse_single_cell/Long_read_pacbio_with_parse_test/ALL_32SMRT_cell_pacbio_files/Original_PacBio_BAM/For_alignment_ParsePac/B0{1..2}-out`

      - The cleaned fastqs of B03,B04, NOT YET START Align, the folder `/Volumes/data/parse_single_cell/Long_read_pacbio_with_parse_test/ALL_32SMRT_cell_pacbio_files/Original_PacBio_BAM/Lib2_fastqs_and_aligned_BAMs/B0{3..4}_barcode_head.fastq.gz`


------

- using TALON to create and annotated chr15 ONLY, can be done only 5 sample per run somehow...
    
  - the result and script running folder: `/Volumes/data/parse_single_cell/Long_read_pacbio_with_parse_test/ALL_32SMRT_cell_pacbio_files/aligned_BAM/TALON_chr15_database_output`
    
  - chr15 BAM slice files: `/data/leec20/parse_single_cell/Long_read_pacbio_with_parse_test/ALL_32SMRT_cell_pacbio_files/aligned_BAM/TALON_chr15_database_output/Parse_pacbio_Merged_chr15_slice/`

------

- Using TALON to cerate and annotated CHRNA5 only:

- the result folder: `/Volumes/data/parse_single_cell/Long_read_pacbio_with_parse_test/ALL_32SMRT_cell_pacbio_files/Original_PacBio_BAM/TALON_Anno_place/Merged_Parse_PacBIo_chr15_slice`

------

- Using PacBio long-read tools, pbmm2 -> iso-seq -> pigen or SQANTI3 
- Files:
      - the folder: `/Volumes/data/parse_single_cell/Long_read_pacbio_with_parse_test/ALL_32SMRT_cell_pacbio_files/ISO_seq_collpase/Pigeon_Pac`
   
      - fastqs files is after parse script preocessed above [files, as we have B0{1..4} sub-library fastqs]
  
      - BAM files are aligned using `pbmm2`, scripts in `Lib1_ParseWay_processing` and `Lib2_ParseWay_processing`; the BAMs in ``

```
#!/bin/bash
export TMPDIR=/lscratch/$SLURM_JOB_ID
module load samtools
source myconda
conda activate PacBio

# command line
hg38=/fdb/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa
pbmm2 align $hg38 myfiles.fofn Merged_lib1_Parseway.PacBio.sorted.bam -j 24 -m 10G --sort --preset ISOSEQ && echo finned

```


  - Next step:

  - merge lib1 and lib2 as single BAM
  - 
```
#!/bin/bash
# Date: 07312025
# loc: Biowulf

ml samtools

samtools merge -@ 12 Lib12_Parseway.PacBio.sorted.bam Merged_lib1_Parseway.PacBio.sorted.bam Merged_lib2_Parseway.PacBio.sorted.bam && echo fined


```
   
  - Do `collapse`
      
      - `isoseq collapse -j 36 Parse_Pac_Merge_isoseq_pbmm2_hg38_sorted.bam Parse_Merged_collapse_isoseq_pbmm2_collapse.gff && echo finndededed`
   
        

#### Peigne workthrough

https://isoseq.how/classification/workflow.html

- the CAGE file and intropolis can be download form:

https://github.com/ConesaLab/SQANTI3/wiki/Running-SQANTI3-Quality-Control#additional-inputs-optional

or just download from Iso-seq website, which is faster:

https://downloads.pacbcloud.com/public/dataset/MAS-Seq/REF-pigeon_ref_sets/Human_hg38_Gencode_v39/

The file are for CAGE: `refTSS_v3.3_human_coordinate.hg38.sorted.bed`,  intropolis: `intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified2.sorted.tsv`, polyA: `polyA.list.txt`


- If you want to generate a filtered GFF, you need to also provide the GFF that was used as input to pigeon classify

- Current all required ref are in:

`/data/leec20/parse_single_cell/Long_read_pacbio_with_parse_test/ALL_32SMRT_cell_pacbio_files/ISO_seq_collpase/Pigeon_Pac/Pigeon_files`

- gtf file was from GENECODE V47 [Comprehensive gene annotation]

https://www.gencodegenes.org/human/release_47.html

- Download all the files need describe above, first generated sorted gtf files

`pigeon prepare gencode.v47.basic.annotation.gtf /data/leec20/hg38.fa refTSS_v3.3_human_coordinate.hg38.bed intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified2.tsv --log-level INFO`

- Once generated sorted files, can remove old unsorted file

- then do pigen classify

`pigeon classify  ../aligned_bams/TEST_TT/ff.gff gencode.v47.basic.annotation.sorted.gtf /data/leec20/hg38.fa --cage-peak refTSS_v3.3_human_coordinate.hg38.sorted.bed --poly-a ployAlist.txt -o SET1 --fl ../aligned_bams/TEST_TT/ff.flnc_count.txt  --log-version INFO`

- Note: the `--fl` can be set for single cell or bulk analysis if info provided

- After build up file, do filter

`pigeon filter SET1_classification.txt -i ../aligned_bams/TEST_TT/ff.gff --mono-exon --log-level  INFO`

- Note: `--mono-exon` remove only single exon transcipts.
  
- The reuslt will generated `xxx.filtered.gff`, can be view on IGV

  


#### SQANTI3 is also a good choice, Note: Piegen is based on this tool

https://github.com/ConesaLab/SQANTI3/wiki/Introduction-to-SQANTI3#workflow




----
----
# CUT OFF 

### Code:


### Result:

This is the script focus on Parse PacBio Long read and long read
application tools

Date: 02062025

### FLAIR tools for analyzing Nanopore reads on CHRNA5 isoform changes in smoking (CSC) vs DMSO in cells

-   Package ref: FLAIR
    <https://flair.readthedocs.io/en/latest/index.html>

-   **Installation**

FLAIR can be conda installed

for using conda on Biowulf, follow the link below
<https://hpc.nih.gov/docs/diy_installation/conda.html>

```         
# create env name called flair and install required packages
conda create -n flair -c conda-forge -c bioconda flair

# activate the env
conda activate flair
```

-   **FLAIR modules**

# Making markdown
- Currently at the step of making plots

## processing pipeline
  - LR-split pipe

## Annotation tools
  - FLAIR: https://flair.readthedocs.io/en/latest/index.html
  - TALON: https://github.com/mortazavilab/TALON/tree/master/example
  - SWAN: https://freese.gitbook.io/swan/tutorials/data_processing
  - ISO-seq with pigen piepline: https://isoseq.how/classification/pigeon.html

## Good pipline refernece for guidence 
- ENCODE Long Read RNA-Seq Analysis Protocol for Human v3.0 (from Minimap2 with TALON)

https://www.encodeproject.org/documents/81af563b-5134-4f78-9bc4-41cb42cc6a48/@@download/attachment/ENCODE%20Long%20Read%20RNA-Seq%20Analysis%20Pipeline%20v3%20%28Human%29.pdf

### FLAIR tools for analyzing Nanopore reads on CHRNA5 isoform changes in smoking (CSC) vs DMSO in cells

-   Package ref: FLAIR
    <https://flair.readthedocs.io/en/latest/index.html>

-   **Installation**

FLAIR can be conda installed

for using conda on Biowulf, follow the link below
<https://hpc.nih.gov/docs/diy_installation/conda.html>

```         
# create env name called flair and install required packages
conda create -n flair -c conda-forge -c bioconda flair

# activate the env
conda activate flair
```

-   **FLAIR modules**
