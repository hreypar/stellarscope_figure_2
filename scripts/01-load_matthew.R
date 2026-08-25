library(scopetools)
library(scater)
library(Matrix)

source('scripts/functions.R')

pbmc20K.pseudo.orig <- load_stellarscope_seurat(
  'data/telescope_pseudobulk/20k_PBMC_3p_HT_nextgem_Chromium_X',
  'data/starsolo_alignment/20k_PBMC_3p_HT_nextgem_Chromium_X/20k_PBMC_3p_HT_nextgem_Chromium_X.Solo.out/Gene/filtered',
  TE_count_file = 'data/telescope_pseudobulk/20k_PBMC_3p_HT_nextgem_Chromium_X/20k_PBMC_3p_HT_nextgem_Chromium_X_pseudobulk-TE_counts_exclude.mtx',
  project='pbmc20K.pseudo'
)
pbmc20K.pseudo <- stellarscope_cell_qc(pbmc20K.pseudo.orig)



pbmc20K.indiv.orig <- load_stellarscope_seurat(
  'data/telescope_individual/20k_PBMC_3p_HT_nextgem_Chromium_X',
  'data/starsolo_alignment/20k_PBMC_3p_HT_nextgem_Chromium_X/20k_PBMC_3p_HT_nextgem_Chromium_X.Solo.out/Gene/filtered',
  TE_count_file = 'data/telescope_individual/20k_PBMC_3p_HT_nextgem_Chromium_X/20k_PBMC_3p_HT_nextgem_Chromium_X_individual-TE_counts_exclude.mtx',
  project='pbmc20K.indiv'
)
pbmc20K.indiv <- stellarscope_cell_qc(pbmc20K.indiv.orig)


# clusters

# seurat_list <- list(pbmc20K.pseudo, pbmc20K.indiv)
# harmonize_seurat <- function(seurat_list) {
#   features <- Reduce(intersect, lapply(seurat_list, rownames))
#   cellnames <- Reduce(intersect, lapply(seurat_list, colnames))
#
#   new_seurat_list <- lapply(seurat_list, function(sobj) {
#     tmp <- sobj[features, ]
#
#     tmp <- tmp[ , cellnames[1:(length(cellnames)-1)]]
#     tmp
#   })
# }


mat.pseudo <- pbmc20K.pseudo[['RNA']]@counts
mat.indiv <- pbmc20K.indiv[['RNA']]@counts
# mat.clust <- pbmc20K.clust[['RNA']]@counts

harmonized <- harmonize_matrix(list(mat.pseudo, mat.indiv))
# harmonized <- harmonize_matrix(list(mat.pseudo, mat.indiv, mat.clust))

mat.pseudo <- harmonized[[1]]
mat.indiv <- harmonized[[2]]
# mat.clust <- harmonized[[3]]

stopifnot(all(colnames(mat.pseudo) == colnames(mat.indiv)))
stopifnot(all(rownames(mat.pseudo) == rownames(mat.indiv)))
# stopifnot(all(colnames(mat.pseudo) == colnames(mat.clust)))
# stopifnot(all(rownames(mat.pseudo) == rownames(mat.clust)))


# All the matrices are harmonized
# Remove genes

diff.pseudo.indiv <- mat.pseudo - mat.indiv

# denominator for barplot
# total number of elements in matrix
nelems <- Reduce('*', dim(diff.pseudo.indiv))
# number of elem that are nonzero in pseudo OR indiv, before drop0
nelems.nz <- length(diff.pseudo.indiv@x)

mat.pseudo.logical <- mat.pseudo != 0
mat.indiv.logical <- mat.indiv != 0
x <- mat.pseudo.logical | mat.indiv.logical
sum(x)
length(x@x)




locus <- 'MER101-16p12.2a'

fdf <- pbmc20K.indiv[['RNA']]@meta.features
fdf <- fdf[complete.cases(fdf),]
alllocs <- fdf[fdf$te_family == 'HERVH', 'id']


by_locus <- lapply(alllocs, function(locus) {
  tmp <-
    data.frame(locus = locus,
               pseudo = mat.pseudo[locus,],
               indiv = mat.indiv[locus,]) %>%
    tibble::rownames_to_column('cell_bc') %>%
    dplyr::group_by(pseudo, indiv) %>%
    dplyr::summarise(count = n(), .groups = 'keep') %>%
    data.frame
  tmp$locus <- locus
  tmp
}) %>% bind_rows

# plot only one locus
by_locus %>%
  dplyr::filter(locus=='HML2-1q22') %>%
  dplyr::filter(pseudo != 0 | indiv != 0) %>%
  ggplot(aes(x=pseudo, y=indiv)) +
  geom_point(aes(size=count)) +
  xlim(0,11) + ylim(0,11) +
  geom_abline(slope=1, color='#99000033')

# plot all loci
by_locus %>%
  dplyr::filter(pseudo != 0 | indiv != 0) %>%
  dplyr::group_by(pseudo, indiv) %>%
  dplyr::summarise(count=sum(count)) %>%
  ggplot(aes(x=pseudo, y=indiv)) +
      geom_point(aes(size=count)) +
      xlim(0,11) + ylim(0,11) +
      geom_abline(slope=1, color='#99000033')



# png(file.path("results", "barplot_individual_pseudobulk_difference.png"),
#     width = 15, height = 11, units = "in", res=300)


o <- ''
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


library(Matrix)
# this rowSums gives length of 27,812, but should be 88,461
# z <- rowSums(diff.pseudo.indiv)
z <- Matrix::rowSums(diff.pseudo.indiv)

z.TE <- z[pbmc20K.indiv[['RNA']]@meta.features$feattype == 'TE']
head(z.TE)

most_diff_loc <- names(z[which(z == max(z))])




# scale_x_discrete(breaks=seq(0,10,1)) + scale_y_discrete(breaks=seq(0,10,1)) +



# CTTCTCTCAAGAATAC





