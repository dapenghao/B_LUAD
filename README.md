# B_LUAD

This repository contains matadata and codes necessary for analysis of the scRNA and scBCR data of human LUAD presented in Hao et al. Cancer Discovery, in revision (2022). 

Cells with low complexity libraries or likely cell debris have been removed. Batch effect was corrected by Harmony and verified by k-BET. The Change-O repertoire clonal assignment toolkit was used to define B cell clones.

Please contact LWang22@mdanderson.org and hkadara@mdanderson.org with any questions.



## Downloading the data

All the required input data have been organized into R objects that can be downloaded from /Input_data/.

Cell level metadata is available in the provided /Input_data/Cell metaData.rda, which contains QC statistics of cells, clinical information of associated samples and cell types. 

Cell level BCR clones and SHM are available in the provided /Input_data/10Xclone_mutation.rda, which contains all the required BCR data of each cell, including VDJ calls, sequence and alignments, clones, SHM and associated information of corresponding samples.

Cell_ID to cell type association for all the TME cells are available in the provided /input_data/cellType_All.rda.

The Seurat object of all the B lineage cells is available using the following dropbox link:

https://www.dropbox.com/s/duliy7xawf0e351/Bcell.rds?dl=0

The Seurat object contains the count matrices, cell metadata and UMAP clustering results. Normalized data matrix is removed for the sake of file size. Please run the seurate function NormalizeData() before using it. 

## Data visualization

### Requirements

Tested on macOS Big Sur

1. R version: 4.1.2
2. R packages
   - ggplot2
   - data.table
   - Seurat
   - dplyr
   - tidyr
   - ggpubr
   - RColorBrewer
   - monocle
   - smoother
   - pheatmap
   - Hmisc
   - monocle3
   - ggrepel
   - CytoTRACE
   - shazam

3. igblast_1.17.1
4. Change-O toolkit
