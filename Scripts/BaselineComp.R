# Load libraries
library(readr)
library(tximport)
library(sva)
library(DESeq2)
library(dplyr)

# Phenotypic data
pheno <- read_csv('Clinical_factors_Pre_PSORT_Amy_v7_15_8_16.csv')
pheno <- subset(pheno, time != 'wk12')

# Genotypic data
t2g <- read_tsv('Transcripts.tsv')
dir <- paste(getwd(), 'LesionalSkin', sep='/')  # Or whatever
files <- file.path(dir, pheno$sample, 'abundance.tsv')
txi <- tximport(files, type='kallisto', tx2gene=t2g, reader=read_csv)

# Normalise and filter data before running SVA
dds <- DESeqDataSetFromTximport(txi, colData=pheno, design= ~ 1)
dds <- estimateSizeFactors(dds)
counts <- counts(dds, normalized=TRUE)
counts <- counts[rowVars(counts) > 1, ]

# Run SVA
mod <- model.matrix(~ 0 + subject + time + time:Delta_PASI, data=pheno)
mod <- mod[, colnames(mod) != 'timewk0:Delta_PASI']
mod0 <- mod[, colnames(mod) != 'timewk1:Delta_PASI']
svobj <- svaseq(counts, mod, mod0)
mod <- cbind(mod, svobj$sv)

# Differential expression
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds, modelMatrix=mod, betaPrior=FALSE)
res <- na.omit(data.frame(results(dds, filterfun=ihw, name='timewk1.Delta_PASI')))
sum(res$padj < 0.05)