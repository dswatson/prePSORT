# Load libraries
library(readr)
library(tximport)
library(sva)
library(DESeq2)
library(dplyr)

# Phenotypic data
pheno <- read_csv('Clinical_factors_Pre_PSORT_Amy_v7_15_8_16.csv')
pheno <- subset(pheno, time != 'wk12')
pheno <- pheno %>%
  mutate(wk0.Delta_PASI = ifelse(time == 'wk0', Delta_PASI, 0),
         wk1.Delta_PASI = ifelse(time == 'wk1', Delta_PASI, 0))

# Genotypic data
t2g <- read_tsv('Transcripts.tsv')
dir <- paste(getwd(), 'LesionalSkin', sep='/')
files <- file.path(dir, pheno$sample, 'abundance.tsv')
txi <- tximport(files, type='kallisto', tx2gene=t2g, reader=read_csv)

# Normalise and filter data before running SVA
dat <- DESeqDataSetFromTximport(txi, colData=pheno, design= 
                                ~ 0 + subject + time + wk1.Delta_PASI)
dat <- estimateSizeFactors(dat)
dat <- counts(dat, normalized=TRUE)
dat <- dat[rowVars(dat) > 1, ]

# Run SVA
mod <- model.matrix(~ 0 + subject + time + wk1.Delta_PASI, data=pheno)
mod0 <- model.matrix(~ 0 + subject + time, data=pheno)
svobj <- svaseq(dat, mod, mod0)
pheno <- cbind(pheno, svobj$sv)
colnames(pheno)[51:53] <- paste0('SV', 1:3)

# Differential expression
dds <- DESeqDataSetFromTximport(txi, colData=pheno, design= ~ 0 + subject + time + 
                                SV1 + SV2 + SV3 + wk1.Delta_PASI)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds, betaPrior=FALSE, maxit=150)
res <- na.omit(data.frame(results(dds, filterfun=ihw)))


