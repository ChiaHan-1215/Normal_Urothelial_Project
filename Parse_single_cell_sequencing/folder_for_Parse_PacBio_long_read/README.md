### Goal: 

Using sc long-read data to detect potential isofrom and intergrated with short-read sc-RNA-seq data



### Data: 

Fastq files: 

    - since fastq files is too large for split-pipe, we use `seqkit` to split fastqs
    - script folder: `/Volumes/data/parse_single_cell/Long_read_pacbio_with_parse_test/ALL_32SMRT_cell_pacbio_files/Fastqs`

    
------

Approachs:

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

- Using PacBio long-read tools, pbmm2 -> iso-seq -> pigen 





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
