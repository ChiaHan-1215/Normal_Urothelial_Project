### Goal: Detecting FGFR3 expression in Single cell RAN-seq of Normal urothelial tissue sample and sep by sex

### Data: 
- The Whole 75 sample of scRNA-seq RDS file are loacted in my Biowulf account: `/data/leec20/parse_single_cell/2024_output/MEGA_cb_all_75_bladder_sample/all-sample/Per_sample_filtering_and_cellAnnted/BP_FINAL_MERGEed_1M_parse_withcellanno_Intered.rds`

- The original sample sex age info: `/LTG/Prokunina Lab/All about samples/Urothelial samples Victor Romanov`
  - Note: the sample id "sample_xx" is the same in this sheet sample number "xx" as the tool add "sample_" before the number "xx"
 
 
- The project folder in Biowulf: `/data/leec20/project_FGFR3_GADD4_parse_single_cell`, In the folder: 
  -   `VR_sample_info.csv` : the original sample info saved as csv file
  -   `project_scRNA_Urth_tissue_FGFR3.R`: The working progress script 

### Code: 
- Scripts are run in Biowulf
