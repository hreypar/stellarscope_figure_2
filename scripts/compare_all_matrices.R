################################################################################
# compare_all_matrices.R
#
# Compare the results of stellarscope pseudobulk, individual and 10x clusters matrices 
#
# (1) compare detected and undetected TEs
# (2) compare filtered cells

# set theory and lists of lists are my jam here
################################################################################
# load libraries
#library(Matrix)
library(magrittr)
#library(UpSetR)
library(ggplot2)
library(ggVennDiagram)
#library(grid)
#library(ggrepel)
#library(gplots)
library(Seurat)

######################### declare functions
# function to subset TEs from a list of feature names names
te_subset <- function(fnames) {
  r <- fnames[!grepl("^ENSG", fnames)]
  return(r)
}

# function to plot a venn diagram of the TE features and cells in 
# (1) counts matrix, (2) raw qc seurat, (3) normalized seurat
compare_matrices_fun <- function(sample.name, stellarscope.mode, reassign.method = "exclude") {
  # read in files
  file.path("data", "counts_matrix_R",
            sample.name,
            reassign.method,
            paste0(sample.name, "_", stellarscope.mode, "_matrix_counts_", reassign.method, ".Rds")) %>%
    readRDS() -> counts.mat
  
  file.path("data", "seurat_raw",
            sample.name,
            reassign.method,
            paste0(sample.name, "_", stellarscope.mode, "_seurat_qc_raw_", reassign.method, ".Rds")) %>%
    readRDS() -> raw.seurat
  
  file.path("data", "seurat_analysis",
            sample.name,
            reassign.method,
            paste0(sample.name, "_", stellarscope.mode, "_seurat_analysis_", reassign.method, ".Rds")) %>%
    readRDS() -> an.seurat
  
  # seurat changes the feature names gggg
  # substitutes underscore for dash
  # at least there aren't TE names with dash so it's safe to substitute back
  table(grepl("-", te_subset(rownames(counts.mat))))
  
  ##### Venn diagram of TE features 
  png(file.path("results", "comparison_features_cells", "features",
                paste0(sample.name, "_", stellarscope.mode, "_", reassign.method, "_venn_features.png")),
      width = 16, height = 10, units = "in", res=300)
  
  print(ggVennDiagram(list(MatCounts=te_subset(rownames(counts.mat)), 
                     SeuratRawQC=gsub("-", "_", te_subset(rownames(raw.seurat))), 
                     SeuratNorm=gsub("-", "_", te_subset(rownames(an.seurat)))
  )
  ) + 
    labs(title=paste0(sample.name, " (", stellarscope.mode, "; ", reassign.method, ")"),
         subtitle = "Features (only TE transcripts) in my objects", caption = Sys.Date()) +
    scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")
  )
  
  dev.off()
  
  ##### Venn diagram of Cells
  png(file.path("results", "comparison_features_cells", "cells",
                paste0(sample.name, "_", stellarscope.mode, "_", reassign.method, "_venn_cells.png")),
      width = 16, height = 10, units = "in", res=300)
  
  print(ggVennDiagram(list(MatCounts=colnames(counts.mat),
                     SeuratRaw=colnames(raw.seurat),
                     SeuratNorm=colnames(an.seurat))
  ) +
    labs(title=paste0(sample.name, " (", stellarscope.mode, "; ", reassign.method, ")"),
         subtitle = "Cells in my objects", caption = Sys.Date()) +
    scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")
  )
  
  dev.off()
  
  ##### save an object with the sets 
  TE_features=te_subset(rownames(counts.mat))

  rawSeurat_features = gsub("-", "_", te_subset(rownames(raw.seurat)))
  normSeurat_features = gsub("-", "_", te_subset(rownames(an.seurat)))
  
  if(length(setdiff(rawSeurat_features[rawSeurat_features %in% TE_features],
          normSeurat_features[normSeurat_features %in% TE_features])) > 0) {
    stop("rawSeurat and normSeurat don't have the same features\n")
  }
  
  # since we checked and the counts matrices all have the full TE annotation set
  # and both raw seurat and analysed seurat have the same TE sets
  detection <- data.frame(TE_features, TE_features %in% rawSeurat_features)
  colnames(detection) <- c("TE_transcripts", paste(sample.name, stellarscope.mode, sep="_"))
  
  return(detection)
}

################################################################################
######################### compare all matrices by dataset
#sample.name = "500_PBMC_3p_LT_Chromium_X"
#reassign.method = "exclude"
#stellarscope.mode = "pseudobulk"
pbmc.datasets <- c("500_PBMC_3p_LT_Chromium_X", "10k_PBMC_3p_nextgem_Chromium_X", "20k_PBMC_3p_HT_nextgem_Chromium_X")

detected_pseudobulk <- lapply(pbmc.datasets, compare_matrices_fun, stellarscope.mode = "pseudobulk")
detected_individual <- lapply(pbmc.datasets, compare_matrices_fun, stellarscope.mode = "individual")
detected_clusters <- lapply(pbmc.datasets, compare_matrices_fun, stellarscope.mode = "clusters")


################################################################################
######################### compare all matrices













