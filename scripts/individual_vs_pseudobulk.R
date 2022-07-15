################################################################################
# individual_vs_pseudobulk.R
#
# Compare the results of stellarscope pseudobulk and stellarscope individual 
#
# (1) compare detected and undetected TEs
# (2) compare read count values to determine differences
# (3) compare different features cell by cell 
################################################################################
# load libraries
#library(Matrix)
#library(magrittr)
#library(UpSetR)
#library(ggplot2)
#library(ggVennDiagram)
#library(grid)
#library(ggrepel)
#library(gplots)
#library(Seurat)

######################### declare functions
# function to 


######################### read in data
#sample.name = "10k_PBMC_3p_nextgem_Chromium_X"
#sample.name = "20k_PBMC_3p_HT_nextgem_Chromium_X"

# these matrices have the full features dataset
file.path("data", "seurat_raw", sample.name, paste0(sample.name, "_pseudobulk_seurat_qc_raw_exclude.Rds")) %>%
  readRDS() -> pseudobulk

file.path("data", "seurat_raw", sample.name, paste0(sample.name, "_pseudobulk_seurat_qc_raw_exclude.Rds")) %>%
  readRDS() -> individual

######################### keep only TE features