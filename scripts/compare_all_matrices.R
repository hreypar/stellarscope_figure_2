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
library(UpSetR)
library(ggplot2)
library(ggVennDiagram)
library(grid)
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
compare_matrices_fun <- function(sample.name, stellarscope.mode, reassign.method = "exclude", clusters.source="") {
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
  ## features
  TE_features=te_subset(rownames(counts.mat))

  rawSeurat_features = gsub("-", "_", te_subset(rownames(raw.seurat)))
  normSeurat_features = gsub("-", "_", te_subset(rownames(an.seurat)))
  
  if(length(setdiff(rawSeurat_features[rawSeurat_features %in% TE_features],
          normSeurat_features[normSeurat_features %in% TE_features])) > 0) {
    stop("rawSeurat and normSeurat don't have the same features!\n")
  }
  
  # since we checked and the counts matrices all have the full TE annotation set
  # and both raw seurat and analysed seurat have the same TE sets
  detection <- data.frame(TE_features, TE_features %in% normSeurat_features)
  if(stellarscope.mode == "clusters") {
    colnames(detection) <- c("TE_transcripts", paste(sample.name, "celltypes", clusters.source, sep="_"))
  } else {
    colnames(detection) <- c("TE_transcripts", paste(sample.name, stellarscope.mode, sep="_"))
  }
  
  ## cells
  initial_cells=colnames(counts.mat)
  
  if(length(setdiff(colnames(raw.seurat), colnames(an.seurat))) > 0) {
    stop("rawSeurat and normSeurat don't have the same cells!\n")
  }
  
  true_cells = data.frame(initial_cells, initial_cells %in% colnames(an.seurat))
  
  if(stellarscope.mode == "clusters") {
    colnames(true_cells) <- c(paste0("initial_cells_", sample.name),
                              paste0("filtered_cells_", stellarscope.mode, "_", clusters.source))
  } else {
    colnames(true_cells) <- c(paste0("initial_cells_", sample.name),
                              paste0("filtered_cells_", stellarscope.mode))
  }
  
  out <- list(detection, true_cells)
  names(out) <- c("te", "cells") # DO I HAVE TO NAME THESE SPECIFICALLY?
  
  return(out)
}

################################################################################
######################### compare all matrices by dataset
#sample.name = "500_PBMC_3p_LT_Chromium_X"
#reassign.method = "exclude"
#stellarscope.mode = "pseudobulk"
#clusters.source = "10X"

pbmc.datasets <- c("500_PBMC_3p_LT_Chromium_X", "10k_PBMC_3p_nextgem_Chromium_X", "20k_PBMC_3p_HT_nextgem_Chromium_X")

detected_pseudobulk <- lapply(pbmc.datasets, compare_matrices_fun, stellarscope.mode = "pseudobulk")
#names(detected_pseudobulk) <- pbmc.datasets
detected_individual <- lapply(pbmc.datasets, compare_matrices_fun, stellarscope.mode = "individual")
detected_clusters <- lapply(pbmc.datasets, compare_matrices_fun, stellarscope.mode = "clusters", clusters.source="10X")

##### check that the features are in the same order
pse <- do.call(cbind, lapply(detected_pseudobulk, function(x) {
  x[["te"]]$TE_transcripts
}))

ind <- do.call(cbind, lapply(detected_individual, function(x) {
  x[["te"]]$TE_transcripts
}))

clu <- do.call(cbind, lapply(detected_clusters, function(x) {
  x[["te"]]$TE_transcripts
}))

all <- cbind(pse, ind, clu)
all <- apply(all, 1, unique)

if(length(all) != 27812) {
  stop("The order of features doesn't match!\n")
} else {
  rm(pse, ind, clu, all)
}

# create features and cells tables and save
pbmc_features_table <- cbind(detected_pseudobulk[[1]]$te,
                             detected_individual[[1]]$te,
                             detected_clusters[[1]]$te,
                             detected_pseudobulk[[2]]$te,
                             detected_individual[[2]]$te,
                             detected_clusters[[2]]$te,
                             detected_pseudobulk[[3]]$te,
                             detected_individual[[3]]$te,
                             detected_clusters[[3]]$te)

pbmc_features_table <- pbmc_features_table[, !duplicated(colnames(pbmc_features_table))]
write.csv(x = pbmc_features_table, file = "results/comparison_features_cells/pbmc_datasets_detected_features.csv", quote = F, row.names = F)


################################################################################
######################### compare all matrices
##### plot upsetR

pbmc_upset <- lapply(colnames(pbmc_features_table), function(x) {
  if(x == "TE_transcripts") {
    return(pbmc_features_table[,`x`])
  } else{
    ind <- pbmc_features_table[,`x`]
    return(pbmc_features_table[ind, "TE_transcripts"])
  }
})

names(pbmc_upset) <- colnames(pbmc_features_table)

png("results/comparison_features_cells/pbmc_datasets_detected_features.png",
    width = 26, height = 12, units = "in", res=300)

upset(fromList(pbmc_upset),
      order.by = c("freq"), point.size = 3.5, nsets = 20, #mb.ratio = c(0.45, 0.55),
      mainbar.y.label = "Detected Transposable Elements",  
      text.scale = c(2, 1.35, 1.35, 1.25, 1.65, 1.95)
) 
grid.text("Comparison of TE Features", x = 0.65, y=0.95, gp=gpar(fontsize=20))

dev.off()













