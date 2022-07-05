################################################################################
# plot_upset.R
#
# plot an upsetR diagram of the TE features recovered by
# each stellarscope method
################################################################################
# load libraries
library(Matrix)
library(magrittr)
library(UpSetR)
library(ggplot2)
library(ggVennDiagram)
library(grid)
library(ggrepel)
library(gplots)

######################### declare functions
keepTEs_function <- function(fullmat) {
  fullmat %>% rownames() %>%
    grepl("^ENSG", .) %>%
    `!` %>%
    fullmat[., ] -> fullmat
  return(fullmat)
}

######################### read in data
#sample.name = "10k_PBMC_3p_nextgem_Chromium_X"
sample.name = "20k_PBMC_3p_HT_nextgem_Chromium_X"

# these matrices have the full features dataset
file.path("data", "counts_matrix_R", sample.name, paste0(sample.name, "_pseudobulk_matrix_counts_exclude.Rds")) %>%
  readRDS() -> pseudobulk

file.path("data", "counts_matrix_R", sample.name, paste0(sample.name, "_individual_matrix_counts_exclude.Rds")) %>%
  readRDS() -> individual

file.path("data", "counts_matrix_R", sample.name, paste0(sample.name, "_clusters_matrix_counts_exclude.Rds")) %>%
  readRDS() -> clusters

######################### keep only TE features
pseudobulk %<>% keepTEs_function()
individual %<>% keepTEs_function()
clusters %<>% keepTEs_function()

# up to this point the matrices have all the cells and all the TEs
if(sum(dim(pseudobulk) == dim(individual) & dim(pseudobulk) == dim(clusters)) != 2) {
  stop("The matrices don't have the same elements\n")
}

# they also must have the same TEs and in the same order
if(sum(rownames(pseudobulk) != rownames(individual) | rownames(pseudobulk) != rownames(clusters)) != 0) {
  stop("The matrices don't have the same TEs")
}

#c(rownames(pseudobulk), rownames(individual), rownames(clusters)) %>%
#  unique() %>% length()

TEs.universe <- rownames(pseudobulk)
################################################################################
##### plot counts distribution 

te.counts.sum <- rbind(data.frame(stellarscope.method="Pseudobulk", te.counts.sum=rowSums(pseudobulk), element.name=names(rowSums(pseudobulk))),
                       data.frame(stellarscope.method="Cell-by-Cell", te.counts.sum=rowSums(individual), element.name=names(rowSums(individual))),
                       data.frame(stellarscope.method="By-Cluster", te.counts.sum=rowSums(clusters), element.name=names(rowSums(clusters)))
                       )
rownames(te.counts.sum) <- NULL

png(file.path("results", sample.name, paste0(sample.name, "_jitter_rowsums.png")),
   width = 17, height = 10, units = "in", res=300)
# png(file.path("results", sample.name, paste0(sample.name, "_jitter_rowsums_subset.png")),
#     width = 17, height = 10, units = "in", res=300)

ggplot(te.counts.sum, aes(x=stellarscope.method, y=te.counts.sum)) + 
  geom_violin() + geom_jitter(position = position_jitter(seed = 1)) + theme_minimal() +
  geom_text_repel(aes(label=ifelse(te.counts.sum>11000, as.character(element.name),''), fontface="bold"),
           position = position_jitter(seed = 1), color="blue") +
  # geom_text_repel(aes(label=ifelse(te.counts.sum>4700, as.character(element.name),''), fontface="bold"), 
  #                 position = position_jitter(seed = 1), color="blue") +
    ggtitle(sample.name) + xlab("\nStellarscope Method") + ylab("Sum of TE counts (across all cells)\n") +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14))
  #+ ylim(0,5000)
  
dev.off()

rm(te.counts.sum)
# rowsums_all <- list(rowSums(pseudobulk),
#                    rowSums(individual),
#                    rowSums(clusters))
# names(rowsums_all) <- c("Pseudobulk", "Cell-by-Cell", "By-Cluster")

# png(file.path("results", sample.name, paste0(sample.name, "_readcounts_sum_boxplot.png")), 
#     width = 15, height = 9, units = "in", res = 300)
# par(mar = c(4, 8, 4, 2))
# boxplot(rowsums_all, outline=FALSE, horizontal = TRUE, notch = TRUE, las = 1,
#         xlab = "Sum of read counts by feature",
#         main = paste0("Counts by feature distribution ", sample.name))
# dev.off()
# 
# rm(rowsums_all)
################################################################################
###### get TEs in the different matrices
which(rowSums(pseudobulk)>0) %>% names() -> nonzero.TEs.pseudobulk
which(rowSums(individual)>0) %>% names() -> nonzero.TEs.individual
which(rowSums(clusters)>0) %>% names() -> nonzero.TEs.clusters

TE.sets <- list(TEs.universe, nonzero.TEs.individual, nonzero.TEs.pseudobulk, nonzero.TEs.clusters)
names(TE.sets) <- c("annotation",  "Cell-by-Cell", "Pseudobulk", "By-Cluster")

rm(TEs.universe, nonzero.TEs.pseudobulk, nonzero.TEs.individual, nonzero.TEs.clusters)

##### plot venn diagram
# png(file.path("results", sample.name, paste0(sample.name, "_venn_annotation.png")),
#     width = 16, height = 10, units = "in", res=300)
# ggVennDiagram(TE.sets) + labs(title=sample.name, subtitle = "Number of detected TE transcripts",
#                               caption = Sys.Date())
# dev.off()

##### plot upset
png(file.path("results", sample.name, paste0(sample.name, "_upsetr_annotation.png")),
    width = 16, height = 10, units = "in", res=300)

upset(fromList(TE.sets),
      order.by = c("degree"), point.size = 3.5, 
      matrix.color = "black", main.bar.color = c(rep("grey25",7), "steelblue"),
      mainbar.y.label = "Detected Transposable Elements",  
      sets.x.label = "Total TEs (set size)", sets.bar.color = c("steelblue","maroon","#93af8b","orange"),
      text.scale = c(2, 1.35, 1.35, 1.25, 1.65, 1.95)
      ) 
grid.text(sample.name,x = 0.65, y=0.95, gp=gpar(fontsize=20))

dev.off()

##### plot venn diagram
png(file.path("results", sample.name, paste0(sample.name, "_venn_methods.png")),
    width = 16, height = 10, units = "in", res=300)

ggVennDiagram(TE.sets[c("Pseudobulk", "Cell-by-Cell", "By-Cluster")]) + 
  labs(title=sample.name, subtitle = "Number of detected TE transcripts", caption = Sys.Date()) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")

dev.off()

##### which are the intersections (save the table)
vdata <- venn(TE.sets[c("Pseudobulk", "Cell-by-Cell", "By-Cluster")], show.plot = FALSE)
te.intersections <- attr(x = vdata, "intersections")

te.intersections <- do.call(rbind, lapply(names(te.intersections), function(x){
  data.frame(x, te.intersections[[x]])
  }))
colnames(te.intersections) <- c("set", "TE")
te.intersections$set <- gsub(pattern = "-", replacement = "", te.intersections$set)

te.intersections <- data.frame(te.intersections,
                               rowsums.pseudobulk=rowSums(pseudobulk)[te.intersections$TE],
                               rowsums.cellbycell=rowSums(individual)[te.intersections$TE],
                               rowsums.bycluster=rowSums(clusters)[te.intersections$TE]
)

cells.pseudobulk = unlist(lapply(te.intersections$TE, function(x) { sum(pseudobulk[x, ] != 0) }))
cells.individual = unlist(lapply(te.intersections$TE, function(x) { sum(individual[x, ] != 0) }))
cells.clusters = unlist(lapply(te.intersections$TE, function(x) { sum(clusters[x, ] != 0) }))

te.intersections <- data.frame(te.intersections, cells.pseudobulk, cells.individual, cells.clusters)

te.intersections <- data.frame(te.intersections, sample.name)

write.csv(x = te.intersections,
          file = file.path("results", sample.name, paste0(sample.name, "_te_intersections_sets.csv")),
          quote = FALSE, row.names = FALSE)

rm(TE.sets, vdata, cells.pseudobulk, cells.individual, cells.clusters)#, te.intersections)
###############################################################################
# differentiate cell barcodes from each matrix because we'll merge them
colnames(pseudobulk) %<>% paste0("_pse")
colnames(individual) %<>% paste0("_ind") 
colnames(clusters) %<>% paste0("_clu") 

# bind them
all <- cbind(pseudobulk, individual, clusters)
# remove rows (TEs) that have zero counts in all the cells in all the matrices
all <- all[rowSums(all) > 0, ]


##### plot differences

pse.ind.diff <- all[, grepl("_pse$", colnames(all)) ] - all[, grepl("_ind$", colnames(all)) ]

pse.clu.diff <- all[, grepl("_pse$", colnames(all)) ] - all[, grepl("_clu$", colnames(all)) ]

ind.clu.diff <- all[, grepl("_ind$", colnames(all)) ] - all[, grepl("_clu$", colnames(all)) ]


as.matrix(ind.clu.diff) %>% table() -> ind.clu.table


######## this needs to be a function
png(file.path("results", sample.name, paste0(sample.name, "_barplot_IND_CLU.png")),
    width = 30, height = 10, units = "in", res=300)

barplot(ind.clu.table[names(ind.clu.table) != "0"], las=1,
        main = paste0("Difference in Read Count between Cell-by-Cell and By-Cluster (",sample.name, ")"),
        xlab = "subtraction value (number of elements)", 
        ylab = "",
        names.arg = paste0(names(ind.clu.table[names(ind.clu.table) != "0"]),
                           "\n(", prettyNum(ind.clu.table[names(ind.clu.table) != "0"], big.mark=","), ")")
        )

dev.off()



