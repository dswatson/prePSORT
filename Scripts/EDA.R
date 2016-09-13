# Load libraries
library(data.table)
library(tximport)
library(DESeq2)
library(ggplot2)
library(wordspace)
library(RColorBrewer)
library(NMF)
library(tidyr)
library(dplyr)

# Import data
pheno <- read.csv('Adalimumab_LS_clin.csv', stringsAsFactors=FALSE)
t2g <- fread('Ensembl.Hs79.Tx.csv')
files <- file.path(getwd(), pheno$Count_file, 'MB.abundance.tsv')
txi <- tximport(files, type='kallisto', tx2gene=t2g, reader=fread)

# Remove underexpressed genes
mat <- txi$counts
mat <- mat[rowMeans(mat) > 1, ]

# Density plot
lmat <- log2(mat + 0.05)
X <- Y <- matrix(0, 512, ncol(lmat))
  for (j in 1:ncol(lmat)) {
    d <- density(lmat[, j], na.rm=TRUE)
    X[, j] <- d$x
    Y[, j] <- d$y
  }
molten <- data.frame(Batch = rep(pheno$batch, each=512), gather(data.frame(X), Sample, lCounts),
                     Density = gather(data.frame(Y), Sample, Density) %>% select(Density))
ggplot(molten, aes(lCounts, Density, colour=Batch)) + geom_line(size=0.5) + 
  theme_bw() + theme(legend.justification=c(1, 1), legend.position=c(1, 1)) + 
  labs(title='Density Plot', x=expression('log'[2]*'(Counts + 0.05)'))

# PCA
mat <- rlog(round(mat))
pca <- prcomp(t(mat), center=TRUE, scale.=TRUE)
dat <- data.frame(pca$x[, 1:2], Batch=pheno$batch)
imp1 <- round(summary(pca)$importance['Proportion of Variance', 1] * 100, 2)
imp2 <- round(summary(pca)$importance['Proportion of Variance', 2] * 100, 2)
ggplot(dat, aes(PC1, PC2, colour=Batch, shape=Batch)) + 
  geom_point() + geom_hline(yintercept=0, size=.2) + geom_vline(xintercept=0, size=.2) + 
  labs(title='PCA', x=paste0('PC1 (', imp1, '%)'), y=paste0('PC2 (', imp2, '%)')) + theme_bw()

# Subject similarity
dm <- dist.matrix(scale(t(mat)), method='euclidean')
rb <- colorRampPalette(brewer.pal(10, 'RdBu'))(n=256)
aheatmap(dm, col=rb, Rowv=FALSE, annCol=list(Batch = pheno$batch), 
  distfun=function(x) as.dist(x), hclustfun='average', main='Subject Similarity Matrix')


