################################################################################
# individual_vs_pseudobulk.R
#
# Compare the results of stellarscope pseudobulk and stellarscope individual 
#
# (1) compare detected and undetected TEs
# (2) compare read count diffs to determine differences
# (3) compare different features cell by cell 
################################################################################
#################### load libraries
library(Matrix)
library(scopetools)
library(Seurat)
library(patchwork)
library(magrittr)
# library(UpSetR)
library(ggplot2)
# library(ggVennDiagram)
# library(VennDiagram)
#
################################################################################
#################### declare functions
harmonize_matrix <- function(matrix_list) {
  featnames <- Reduce(intersect, lapply(matrix_list, rownames))
  cellnames <- Reduce(intersect, lapply(matrix_list, colnames))
  
  lapply(matrix_list, function(mat) mat[featnames, cellnames])
}
################################################################################
#################### declare variables
sample.name = "20k_PBMC_3p_HT_nextgem_Chromium_X"
reassign.method = "exclude"
#stellarscope.mode = "pseudobulk"
# clusters.source = "10X"
#
#################### read in stellarscope results and apply QC to them
pseudobulk_20k <- scopetools::load_stellarscope_seurat(stellarscope_dir = "data/telescope_pseudobulk/20k_PBMC_3p_HT_nextgem_Chromium_X/",
                                      TE_count_file = "data/telescope_pseudobulk/20k_PBMC_3p_HT_nextgem_Chromium_X/20k_PBMC_3p_HT_nextgem_Chromium_X_pseudobulk-TE_counts_exclude.mtx", 
                                     starsolo_dir = "data/starsolo_alignment/20k_PBMC_3p_HT_nextgem_Chromium_X/20k_PBMC_3p_HT_nextgem_Chromium_X.Solo.out/Gene/filtered/")

pseudobulk_20k.qc <- scopetools::stellarscope_cell_qc(pseudobulk_20k)

individual_20k <- scopetools::load_stellarscope_seurat(stellarscope_dir = "data/telescope_individual/20k_PBMC_3p_HT_nextgem_Chromium_X/", 
                                                       TE_count_file = "data/telescope_individual/20k_PBMC_3p_HT_nextgem_Chromium_X/20k_PBMC_3p_HT_nextgem_Chromium_X_individual-TE_counts_exclude.mtx",
                                                       starsolo_dir = "data/starsolo_alignment/20k_PBMC_3p_HT_nextgem_Chromium_X/20k_PBMC_3p_HT_nextgem_Chromium_X.Solo.out/Gene/filtered/")

individual_20k.qc <- scopetools::stellarscope_cell_qc(individual_20k)

rm(pseudobulk_20k, individual_20k)
################################################################################
################ Find the differences between the two matrices #################
######## extract counts matrices
counts.list <- list(individual_20k.qc@assays$RNA@counts, pseudobulk_20k.qc@assays$RNA@counts)
names(counts.list) <- c("ind.counts", "psb.counts") 

counts.list.harmonized <- harmonize_matrix(counts.list)
rm(counts.list)

# remove protein coding transcripts
# psb.counts <- psb.counts[!grepl("^ENSG", rownames(psb.counts)), ]
# ind.counts <- ind.counts[!grepl("^ENSG", rownames(ind.counts)), ]
# # the barcodes and elements match
# setdiff(colnames(psb.counts), colnames(ind.counts))
# setdiff(rownames(psb.counts), rownames(ind.counts))
# # the elements are in the same order but the barcodes are not
# table(rownames(ind.counts) == rownames(psb.counts))
# table(colnames(ind.counts) == colnames(psb.counts))
# # now they are
# psb.counts <- psb.counts[,colnames(ind.counts)]
# table(colnames(ind.counts) == colnames(psb.counts))

######## subtract pseudobulk from individual
# zero are the same counts in both
# positive numbers are higher counts in individual
# negative numbers are higher counts in pseudobulk
ind_minus_psb <- counts.list.harmonized$ind.counts - counts.list.harmonized$psb.counts

########  denominator for barplot
# total number of elements in subtracted matrix
nelems <- Reduce('*', dim(ind_minus_psb))
#https://blog.zhaw.ch/datascience/r-reduce-applys-lesser-known-brother/

# number of elements that are nonzero in pseudobulk OR individual, before drop0
nelems.nz <- length(ind_minus_psb@x)

# logical matrices that indicate nonzero elements
psb.counts.nz.logical <- counts.list.harmonized$psb.counts != 0
ind.counts.nz.logical <- counts.list.harmonized$ind.counts != 0

x <- mat.pseudo.logical | mat.indiv.logical
sum(x)
length(x@x)







# sparse matrices are the best
diffs <- table(ind_minus_psb@x) 

#o = c("-1", "0", "1")
o = ""

png(file.path("results", "barplot_individual_pseudobulk_difference.png"),
    width = 15, height = 11, units = "in", res=300)

par(oma=c(3,6,3,3))
barplot(diffs[!names(diffs) %in% o], horiz = TRUE, las=1, border="grey", col = "cornflowerblue", 
        xlim = c(0,250), space = 1, xaxt="n", 
        main = paste("Individual - Pseudobulk\n", sample.name),
        names.arg = paste0(names(diffs[!names(diffs) %in% o]),
                           "  (n=", prettyNum(diffs[!names(diffs) %in% o], big.mark=","), ")")
)
axis(side = 1, at = seq(0,250,25))
lapply(seq(0,250,10), function(y) {
  abline(v = y, lty=2, lwd=0.5, col="grey")
})

dev.off()


# which are they?

ind_minus_psb %<>% drop0()

#which(ind_minus_psb < -5, arr.ind = TRUE) %>% rownames() %>% unique() -> my.super.dif.elements

which(ind_minus_psb < 0, arr.ind = TRUE) %>% rownames() %>% table() %>% sort(decreasing = T) -> minus
which(ind_minus_psb > 0, arr.ind = TRUE) %>% rownames() %>% table() %>% sort(decreasing = T) -> plus

my.super.dif.elements <- c(names(minus), names(plus)[1:6]) # arbitrary selection of six

elements.counts.df <- do.call(rbind, lapply(my.super.dif.elements, function(e) {
  
  data.frame(Individual=ind.counts[e, ], Pseudobulk=psb.counts[e,], Element=e)
}))

elements.counts.df$Element <- factor(elements.counts.df$Element, levels = my.super.dif.elements)

png(file.path("results", "scatterplots_individual_pseudobulk_difference.png"),
    width = 18, height = 16, units = "in", res=300)

ggplot(elements.counts.df, aes(Pseudobulk, Individual)) + geom_point() + theme_linedraw() + 
  geom_abline(slope = 1) + ggtitle(sample.name) + facet_wrap(~Element, ncol = 4) +
  theme(strip.text.x = element_text(size = 12, colour = "white", face="bold"),
        axis.title = element_text(size=13)) +
  ylab("Individual counts per cell") + xlab("Pseudobulk counts per cell") + 
  xlim(0,20) + ylim(0,20)

dev.off()

# bubble plot  at the family level (plotting all the differently counted loci from each family)

different.loci.pse.ind <- sort(unique(c(names(minus), names(plus)))) 




different.loci.pse.ind <- do.call(rbind, lapply(different.loci.pse.ind, function(e) {
  separated <- strsplit(e, "-")
  data.frame(Family=separated[1], Element=e,
             Individual=ind.counts[e, ], Pseudobulk=psb.counts[e,])
}))

# 
different.loci.pse.ind <- do.call(rbind, lapply(strsplit(different.loci.pse.ind, "-"), function(x) {
  data.frame(family=x[1], element=paste(x[1], x[2], sep = "-"))
}))

sort(table(different.loci.pse.ind$family), decreasing = TRUE)




#table(pseudobulk_20k.qc@assays$RNA@meta.features$te_family)


# https://r-graph-gallery.com/320-the-basis-of-bubble-plot.html

#https://www.rdocumentation.org/packages/Matrix/versions/1.4-1/topics/dgCMatrix-class

#grid.newpage()
#grid::grid.draw(VennDiagram::venn.diagram(lapply(te.total.counts, names), NULL))

# # check that features match
# table(rownames(pseudobulk_20k.qc) == rownames(individual_20k.qc))
# 
# # check the cells
# setdiff(colnames(pseudobulk_20k.qc), colnames(individual_20k.qc))
# table(colnames(pseudobulk_20k.qc) == colnames(individual_20k.qc))

# harmonize seurats
# pbmc_20k.qc.list <- list(pseudobulk_20k.qc, individual_20k.qc)
# rm(pseudobulk_20k.qc, individual_20k.qc)
# 
# pbmc_20k.qc.cells <- Reduce(intersect, lapply(pbmc_20k.qc.list, colnames))
# 
# new <- lapply(pbmc_20k.qc.list, function(s) {
#   subset(s, cells= pbmc_20k.qc.cells)
# })


################################################################################
# check distributions 
#
# v <- Seurat::VlnPlot(pseudobulk_20k.qc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# t <- Seurat::VlnPlot(pseudobulk_20k.qc, features = c("percent.TE"))
# # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
# plot1 <- FeatureScatter(pseudobulk_20k.qc, feature1 = "nFeature_RNA", feature2 = "percent.TE")
# plot2 <- FeatureScatter(pseudobulk_20k.qc, feature1 = "nCount_RNA", feature2 = "percent.TE")
# 
# pdf(filename=file.path("individual_vs_pseudobulk", paste0(sample.name, "_qc.pdf")),
#     width=11, height=8.5, units="in")
# print(v /  (plot1 + plot2) + 
#         plot_layout( guides = "collect") + 
#         plot_annotation(title = paste0("\n", sample.name," (Pseudobulk)\n"), 
#                         theme = theme(plot.title = element_text(face = "bold", size = 14, hjust=0.5))) 
# )
# dev.off()
#
#rm(plot1, plot2, v, t)
################################################################################
########################### Detected features ##################################
#####  function to get the sum of counts per row (feature) from a seurat object
# only delivers TE features counts
# get_feature_counts_sum <- function(seurat.object) {
#   # get matrix
#   seurat.object@assays$RNA@counts %>% Matrix::rowSums() -> features.total.counts
#   
#   # keep only TE features
#   #table(seurat.object@assays$RNA@meta.features$feattype)
#   #grepl("^ENSG", row.names(seurat.object.counts)) %>% table()
#   features.total.counts <- features.total.counts[!grepl("^ENSG", names(features.total.counts))]
#   
#   # only return the ones that have counts
#   features.total.counts <- features.total.counts[features.total.counts > 0]
#   
#   return(features.total.counts)
# }
##########################
# we need the clusters matrix for the venn diagram
# how did Matthew make the Seurat objects keep all the features? min cells and features equal to zero
# clusters_20k.qc <- readRDS("data/seurat_raw/20k_PBMC_3p_HT_nextgem_Chromium_X/exclude/20k_PBMC_3p_HT_nextgem_Chromium_X_clusters_seurat_qc_raw_exclude.Rds")
# 
# # venn diagram of detected features 
# te.total.counts <- lapply(list(pseudobulk_20k.qc, individual_20k.qc, clusters_20k.qc),
#                           get_feature_counts_sum)
# 
# names(te.total.counts) <- c("pseudobulk", "individual", paste0("clusters_", clusters.source))
# 
# png(file.path("results", paste0(sample.name, "_", reassign.method, "_venn_features.png")),
#     width = 16, height = 10, units = "in", res=300)
# 
# ggVennDiagram(lapply(te.total.counts, names)) + 
#   labs(title=paste0(sample.name, " (", reassign.method, ")"),
#        subtitle = "Detected TE transcripts", caption = Sys.Date()) +
#   scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")
# 
# dev.off()
# 
# # upsetR diagram
# upset(fromList(lapply(te.total.counts, names)),
#       order.by = c("freq"), point.size = 3.5, 
#       matrix.color = "black",
#       mainbar.y.label = "Detected Transposable Elements",  
#       sets.x.label = "Total TEs (set size)",
#       text.scale = c(2, 1.35, 1.35, 1.25, 1.65, 1.95)
# ) 
# grid.text(sample.name,x = 0.65, y=0.95, gp=gpar(fontsize=20))
# 
# 
# # get the sets of transcripts
# TE.sets <- VennDiagram::get.venn.partitions(lapply(te.total.counts, names))
# 
# # with the different azimuth clusters
################################################################################