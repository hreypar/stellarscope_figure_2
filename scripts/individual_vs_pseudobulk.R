################################################################################
# individual_vs_pseudobulk.R
#
# Compare the results of stellarscope pseudobulk and stellarscope individual 
#
# (1) compare detected and undetected TEs
# (2) compare read count diffs to determine differences
# (3) compare different features cell by cell 
################################################################################
# Detected TEs (venn or upset diagram 3 stellarscope methods)

# Detected TEs (output lists of elements e.g. as a table)

# Read counts differences (barplot of table of subtracted matrix) individual vs pseudobulk
### TE that was counted most differently across cells (table of "is it equal to zero?", or not equal to zero)
### For a given TE, number of reads that are different (abs of rowsum)

# Read counts differences (output which are the most differently counted elements)
# Read counts diffetences (scatterplot of a given feature)

######
# Use read in objects function
# Apply QC function
# Harmonize function when it's neccesary 


# load libraries
library(Matrix)
library(scopetools)
library(magrittr)
library(UpSetR)
library(ggplot2)
library(ggVennDiagram)
library(VennDiagram)
#library(grid)
#library(ggrepel)
#library(gplots)
#library(Seurat)
library(scater)

######################### declare functions
##### QC seurat function
# Scater library dependency
# maybe add argument to keep/discard non-qc-passing-cells 
stellarscope_cell_qc <- function(the.seurat) {
  fmeta <- the.seurat[['RNA']]@meta.features
  mt_feats <- grepl('^MT-', fmeta$symbol)
  the.seurat[['percent.mt']] <- Seurat::PercentageFeatureSet(the.seurat, features=fmeta[mt_feats, 'id'])
  herv_feats <- !is.na(fmeta$te_class) & fmeta$te_class == 'LTR'
  the.seurat[['percent.HERV']] <- Seurat::PercentageFeatureSet(the.seurat, features=fmeta[herv_feats, 'id'])
  l1_feats <- !is.na(fmeta$te_class) & fmeta$te_class == 'LINE'
  the.seurat[['percent.L1']] <- Seurat::PercentageFeatureSet(the.seurat, features=fmeta[l1_feats, 'id'])
  te_feats <- fmeta$feattype == 'TE'
  the.seurat[['percent.TE']] <- Seurat::PercentageFeatureSet(the.seurat, features=fmeta[te_feats, 'id'])
  
  qc.ncount_rna <- scater::isOutlier(the.seurat$nCount_RNA, log = TRUE, type = "both")
  qc.nfeature_rna <- scater::isOutlier(the.seurat$nFeature_RNA, log = TRUE, type = "both")
  qc.percent_mt <- scater::isOutlier(the.seurat$percent.mt,  type="higher")
  
  thresh <- data.frame(ncount = attr(qc.ncount_rna, "thresholds"),
                       nfeature = attr(qc.nfeature_rna, "thresholds"),
                       mt = attr(qc.percent_mt, "thresholds"))
  
  the.seurat <- subset(the.seurat, subset = nCount_RNA > thresh["lower", "ncount"] &
                         nCount_RNA < thresh["higher", "ncount"] &
                         nFeature_RNA >  thresh["lower", "nfeature"] & 
                         nFeature_RNA < thresh["higher", "nfeature"] & 
                         percent.mt < thresh["higher", "mt"])
  the.seurat
}

#####  function to get the sum of counts per row (feature) from a seurat object
# only delivers TE features counts
get_feature_counts_sum <- function(seurat.object) {
  # get matrix
  seurat.object@assays$RNA@counts %>% Matrix::rowSums() -> features.total.counts
  
  # keep only TE features
  #table(seurat.object@assays$RNA@meta.features$feattype)
  #grepl("^ENSG", row.names(seurat.object.counts)) %>% table()
  features.total.counts <- features.total.counts[!grepl("^ENSG", names(features.total.counts))]
  
  # only return the ones that have counts
  features.total.counts <- features.total.counts[features.total.counts > 0]
  
  return(features.total.counts)
}

######################### read in data
sample.name = "20k_PBMC_3p_HT_nextgem_Chromium_X"
reassign.method = "exclude"
#stellarscope.mode = "pseudobulk"
clusters.source = "10X"

######################### read in stellarscope results
pseudobulk_20k <- scopetools::load_stellarscope_seurat(stellarscope_dir = "data/telescope_pseudobulk/20k_PBMC_3p_HT_nextgem_Chromium_X/",
                                      TE_count_file = "data/telescope_pseudobulk/20k_PBMC_3p_HT_nextgem_Chromium_X/20k_PBMC_3p_HT_nextgem_Chromium_X_pseudobulk-TE_counts_exclude.mtx", 
                                     starsolo_dir = "data/starsolo_alignment/20k_PBMC_3p_HT_nextgem_Chromium_X/20k_PBMC_3p_HT_nextgem_Chromium_X.Solo.out/Gene/filtered/")

individual_20k <- scopetools::load_stellarscope_seurat(stellarscope_dir = "data/telescope_individual/20k_PBMC_3p_HT_nextgem_Chromium_X/", 
                                                       TE_count_file = "data/telescope_individual/20k_PBMC_3p_HT_nextgem_Chromium_X/20k_PBMC_3p_HT_nextgem_Chromium_X_individual-TE_counts_exclude.mtx",
                                                       starsolo_dir = "data/starsolo_alignment/20k_PBMC_3p_HT_nextgem_Chromium_X/20k_PBMC_3p_HT_nextgem_Chromium_X.Solo.out/Gene/filtered/")

# QC seurat object
pseudobulk_20k.qc <- stellarscope_cell_qc(pseudobulk_20k)

individual_20k.qc <- stellarscope_cell_qc(individual_20k)

################################################################################
################### Detected features
# we need the clusters matrix for the venn diagram
# how did Matthew make the Seurat objects keep all the features? min cells and features equal to zero
clusters_20k.qc <- readRDS("data/seurat_raw/20k_PBMC_3p_HT_nextgem_Chromium_X/exclude/20k_PBMC_3p_HT_nextgem_Chromium_X_clusters_seurat_qc_raw_exclude.Rds")

# venn diagram of detected features 
te.total.counts <- lapply(list(pseudobulk_20k.qc, individual_20k.qc, clusters_20k.qc),
                          get_feature_counts_sum)

names(te.total.counts) <- c("pseudobulk", "individual", paste0("clusters_", clusters.source))

png(file.path("results", paste0(sample.name, "_", reassign.method, "_venn_features.png")),
    width = 16, height = 10, units = "in", res=300)

ggVennDiagram(lapply(te.total.counts, names)) + 
  labs(title=paste0(sample.name, " (", reassign.method, ")"),
       subtitle = "Detected TE transcripts", caption = Sys.Date()) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")

dev.off()

# upsetR diagram
upset(fromList(lapply(te.total.counts, names)),
      order.by = c("freq"), point.size = 3.5, 
      matrix.color = "black",
      mainbar.y.label = "Detected Transposable Elements",  
      sets.x.label = "Total TEs (set size)",
      text.scale = c(2, 1.35, 1.35, 1.25, 1.65, 1.95)
) 
grid.text(sample.name,x = 0.65, y=0.95, gp=gpar(fontsize=20))



# get the sets of transcripts
#Reduce(intersect, lapply(te.total.counts, names))

TE.sets <- VennDiagram::get.venn.partitions(lapply(te.total.counts, names))

# with the different azimuth clusters


# you really have to harmonize to get scatterplots
# t <- unlist(TE.sets[2,5])
# 
# plot(individual_20k.qc[which(rownames(individual_20k.qc) == t)]@assays$RNA@counts,
#      clusters_20k.qc[which(rownames(individual_20k.qc) == t)]@assays$RNA@counts)


# subtract matrices
individual.counts <- individual_20k.qc@assays$RNA@counts
pseudobulk.counts <- pseudobulk_20k.qc@assays$RNA@counts

table(rownames(individual.counts) == rownames(pseudobulk.counts))
table(colnames(individual.counts) == colnames(pseudobulk.counts))

pseudobulk.counts <- pseudobulk.counts[,colnames(individual.counts)]
table(colnames(individual.counts) == colnames(pseudobulk.counts))


ind_minus_pse <- individual.counts - pseudobulk.counts 
ind_minus_pse <- ind_minus_pse[!grepl("^ENSG", rownames(ind_minus_pse)), ]
#https://stat.ethz.ch/R-manual/R-devel/library/Matrix/html/drop0.html


# sparse matrices are the best
diffs <- table(ind_minus_pse@x) 

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



#which(ind_minus_pse < -5, arr.ind = TRUE) %>% rownames() %>% unique() -> my.super.dif.elements

which(ind_minus_pse < -5, arr.ind = TRUE) %>% rownames() %>% table() %>% sort(decreasing = T) -> minus
which(ind_minus_pse > 2, arr.ind = TRUE) %>% rownames() %>% table() %>% sort(decreasing = T) -> plus

my.super.dif.elements <- c(names(minus), names(plus)[1:6]) # arbitrary selection of six

test <- do.call(rbind, lapply(my.super.dif.elements, function(e) {
  
  data.frame(Individual=individual.counts[e, ], Pseudobulk=pseudobulk.counts[e,], Element=e)
}))

test$Element <- factor(test$Element, levels = my.super.dif.elements)

png(file.path("results", "scatterplots_individual_pseudobulk_difference.png"),
    width = 18, height = 16, units = "in", res=300)

ggplot(test, aes(Pseudobulk, Individual)) + geom_point() + theme_linedraw() + 
  geom_abline(slope = 1) + ggtitle(sample.name) + facet_wrap(~Element, ncol = 4) +
  theme(strip.text.x = element_text(size = 12, colour = "white", face="bold"),
        axis.title = element_text(size=13)) +
  ylab("Individual counts per cell") + xlab("Pseudobulk counts per cell") + 
  xlim(0,20) + ylim(0,20)

dev.off()

# 

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





