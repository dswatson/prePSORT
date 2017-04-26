### HEATMAPS ###

# Load libraries
library(limma)
library(RColorBrewer)
library(ggsci)
library(NMF)
library(dplyr)

# Import data
mat <- readRDS('./Data/mat_mRNA.rds')
fit <- readRDS('./Data/fit_mRNA.rds')
clin <- read.csv('./Data/Clinical.csv')
sup <- read.csv('./Results/Clusters/Supervised/PAM.csv') 
sup <- sup[, grep('mRNA', colnames(sup))]
unsup <- read.csv('./Results/Clusters/Unsupervised/PAM.csv')
unsup <- unsup[, grep('mRNA', colnames(unsup))]

# Prep data, colors
coefs <- colnames(fit)[grep('wk00.Response', colnames(fit))]
n_genes <- round(0.01 * nrow(mat))
rb <- colorRampPalette(rev(brewer.pal(10, 'RdBu')))(n = 256)
greys <- colorRampPalette(brewer.pal(9, 'Greys'))(n = 256)
d3 <- pal_d3()(6)
cols <- list(greys, d3[1:2], d3[3:4], d3[5:6])

# Create heatmaps
for (coef in coefs) {
  top <- topTable(fit, coef = coef, number = n_genes, sort.by = 'p')
  genes <- rownames(top)
  tissue <- gsub('.wk00.Response', '', coef)
  pheno <- clin %>% filter(Tissue == tissue, Time == 'wk00')
  top <- mat[genes, match(pheno$Sample, colnames(mat))]
  colnames(top) <- gsub('\\..*', '', colnames(top))
  pdf(paste0('./Results/Heatmaps/', tissue, '_wk00.pdf'), height = 8, width = 5)
  aheatmap(top, distfun = 'pearson', scale = 'row', col = rb, hclustfun = 'average',
           main = paste0('Top 1% of Genes by Response:\n', tissue, ', wk00'),
           annCol = list('Delta PASI' = pheno$DeltaPASI,
                            'PASI 75' = pheno$PASI_75,
                         'Supervised' = sup[, grep(tissue, colnames(sup))],
                       'Unsupervised' = unsup[, grep(tissue, colnames(unsup))]),
           annColors = cols, border_color = 'black')
  dev.off()
}



