### LESIONAL VS. NONLESIONAL SKIN ###

# Load libraries
library(readr)
library(tximport)
library(sva)
library(DESeq2)
library(dplyr)

# Import data
setwd('/Users/David/Documents/QMUL/PSORT')
pheno <- read.csv('Clinical_factors_Pre_PSORT_Amy_v7_15_8_16.csv')
t2g <- read_tsv('Transcripts.tsv')
dir_L <- paste(getwd(), 'LesionalSkin', sep='/')
files_L <- file.path(dir_L, pheno$sample, 'abundance.tsv')
txi_L <- tximport(files_L, type='kallisto', tx2gene=t2g, reader=read_csv)
dir_N <- paste(getwd(), 'NonlesionalSkin', sep='/')
files_N <- file.path(dir_N, pheno$sample, 'abundance.tsv')
txi_N <- tximport(files_N, type='kallisto', tx2gene=t2g, reader=read_csv)
a <- cbind(txi_L$abundance, txi_N$abundance)
c <- cbind(txi_L$counts, txi_N$counts)
l <- cbind(txi_L$length, txi_N$length)
txi <- list(abundance=a, counts=c, length=l, countsFromAbundance='no')
pheno <- rbind(pheno, pheno)
pheno$site <- c(rep('L', 30), rep('N', 30))
pheno$sample <- paste(pheno$sample, pheno$site, sep='.')
pheno$resp <- paste(pheno$time, pheno$site, sep='.')

# Normalise and filter data before running SVA
dds <- DESeqDataSetFromTximport(txi, colData=pheno, design= ~ subject + resp)
dds <- estimateSizeFactors(dds)
dat <- counts(dds, normalized=TRUE)
dat <- dat[rowVars(dat) > 1, ]

# Run SVA
mod <- model.matrix(~ subject + resp, data=pheno)
mod0 <- model.matrix(~ subject, data=pheno)
svobj <- svaseq(dat, mod, mod0)

# Add surrogate variables to clinical factors
pheno <- cbind(pheno, svobj$sv)
colnames(pheno)[49:61] <- paste0('SV', 1:13)

# Differential expression
dds <- DESeqDataSetFromTximport(txi, colData=pheno, design= ~ SV1 + SV2 + SV3 + 
                                SV4 + SV5 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + 
                                SV11 + SV12 + SV13 + subject + site)
dds <- DESeq(dds)
res <- na.omit(data.frame(results(dds, filterfun=ihw, 
               contrast=c('resp', 'wk0.L', 'wk0.N'))))


