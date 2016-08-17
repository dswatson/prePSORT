### LESIONAL VS. NONLESIONAL SKIN ###

# Load libraries
library(readr)
library(tximport)
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

# Differential expression
dds <- DESeqDataSetFromTximport(txi, colData=pheno, design= ~ subject + resp)
dds <- DESeq(dds)
res <- na.omit(data.frame(results(dds, filterfun=ihw, 
               contrast=c('resp', 'wk0.L', 'wk0.N'))))


