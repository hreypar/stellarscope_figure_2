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

######################### declare functions
keepTEs_function <- function(fullmat) {
  fullmat %>% rownames() %>%
    grepl("^ENSG", .) %>%
    `!` %>%
    fullmat[., ] -> fullmat
  return(fullmat)
}

######################### read in data
sample.name = "10k_PBMC_3p_nextgem_Chromium_X"

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

c(rownames(pseudobulk), rownames(individual), rownames(clusters)) %>%
  unique() %>% length()

TEs.universe <- rownames(pseudobulk)

##### plot counts distribution 
rowsums_all <- list(rowSums(pseudobulk),
                   rowSums(individual),
                   rowSums(clusters))
names(rowsums_all) <- c("Pseudobulk", "Cell-by-Cell", "By-Cluster")

png(file.path("results", sample.name, paste0(sample.name, "_readcounts_sum_boxplot.png")), 
    width = 15, height = 9, units = "in", res = 300)
par(mar = c(4, 8, 4, 2))
boxplot(rowsums_all, outline=FALSE, horizontal = TRUE, notch = TRUE, las = 1,
        xlab = "Sum of read counts by feature",
        main = paste0("Counts by feature distribution ", sample.name))
dev.off()

rm(rowsums_all)
######################### get TEs in the different matrices
which(rowSums(pseudobulk)>0) %>% names() -> nonzero.TEs.pseudobulk
which(rowSums(individual)>0) %>% names() -> nonzero.TEs.individual
which(rowSums(clusters)>0) %>% names() -> nonzero.TEs.clusters

TEs.sets <- list(TEs.universe, nonzero.TEs.pseudobulk, nonzero.TEs.individual, nonzero.TEs.clusters)
names(TEs.sets) <- c("TE annotation", "Pseudobulk", "Cell-by-Cell", "By-Cluster")

rm(TEs.universe, nonzero.TEs.pseudobulk, nonzero.TEs.individual, nonzero.TEs.clusters)

##### plot upset

png(file.path("results", sample.name, paste0(sample.name, "_upsetr_all.png")),
    width = 16, height = 10, units = "in", res=300)

upset(fromList(TEs.sets[2:4]), order.by = "freq", point.size = 3.5,
      mainbar.y.label = "Detected Transposable Elements", main.bar.color = "cornflowerblue",
      sets.x.label = "Total TEs (set size)", sets.bar.color = "cornflowerblue",
      text.scale = c(2, 1.35, 1.35, 1.25, 1.65, 1.95)
      )

dev.off()


# which are the intersections (save the table)

rm(TEs.sets)
###############################################################################

# differentiate cell barcodes from each matrix because we'll merge them
colnames(pseudobulk) %<>% paste0("_pse")
colnames(individual) %<>% paste0("_ind") 
colnames(clusters) %<>% paste0("_clu") 

# bind them
all <- cbind(pseudobulk, individual, clusters)

all <- all[rowSums(all) > 0, ]

# plot differences

ind.pse.diff <- all[, grepl("_ind$", colnames(all)) ] - all[, grepl("_pse$", colnames(all)) ]

pse.clu.diff <- all[, grepl("_pse$", colnames(all)) ] - all[, grepl("_clu$", colnames(all)) ]

ind.clu.diff <- all[, grepl("_ind$", colnames(all)) ] - all[, grepl("_clu$", colnames(all)) ]



#table(as.matrix(ind.clu.diff))

table(ind.clu.diff)

######## this needs to be a function
png(paste0("barplot_subtraction_", s, ".png"), width = 28, height = 11, units = "in", res=300)

barplot(subtraction, las=1,
        main = paste("Difference in Read Count between Cell-by-Cell and Pseudobulk\n",sample.name, "matrices", size),
        xlab = "subtraction value (number of elements)", 
        ylab = "",
        names.arg = paste0(names(te.diff), "\n(", prettyNum(te.diff, big.mark=","), ")"))
#ylab = paste0("Transposable Elements  (total=", length(te.diff),")"))

# PBMC_500
#abline(h = seq(500,6500, 500), lty=2, col="grey79")

# PBMC_10k and PBMC_20k
#abline(h = seq(500,10000, 500), lty=2, col="grey79")
#axis(side = 2, at = seq(0,10000,1000), las=1)

dev.off()







