library(scopetools)
library(scater)

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


harmonize_matrix <- function(matrix_list) {
  featnames <- Reduce(intersect, lapply(matrix_list, rownames))
  cellnames <- Reduce(intersect, lapply(matrix_list, colnames))

  lapply(matrix_list, function(mat) mat[featnames, cellnames])
}

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
library(Matrix)
# this rowSums gives length of 27,812, but should be 88,461
# z <- rowSums(diff.pseudo.indiv)
z <- Matrix::rowSums(diff.pseudo.indiv)

z.TE <- z[pbmc20K.indiv[['RNA']]@meta.features$feattype == 'TE']
head(z.TE)

most_diff_loc <- names(z[which(z == max(z))])


locus <- 'MER101-16p12.2a'
data.frame(pseudo=mat.pseudo[locus, ], indiv=mat.indiv[locus, ]) %>%
  ggplot(aes(x=pseudo, y=indiv)) +
      geom_jitter(size=1, shape=12, height=0.01, width=0.1) +
      xlim(0,11) + ylim(0,11) +
      geom_abline(slope=1, color='#99000033')

      # scale_x_discrete(breaks=seq(0,10,1)) + scale_y_discrete(breaks=seq(0,10,1)) +
  


# CTTCTCTCAAGAATAC





