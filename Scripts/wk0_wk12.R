# Load libraries
library(data.table)
library(tximport)
library(sva)
library(DESeq2)
library(dplyr)

# Prep data
pheno <- read.csv('Clinical.csv', stringsAsFactors=FALSE)
pheno <- filter(pheno, time != 'wk1')
pheno <- pheno %>%
  mutate(wk0.Delta_PASI = ifelse(time == 'wk0', Delta_PASI, 0),
         wk12.Delta_PASI = ifelse(time == 'wk12', Delta_PASI, 0))
t2g <- fread('Ensembl.Hs79.Tx.csv')
e2g <- fread('Ensembl.Hs79.GeneSymbols.csv')

# Pipeline
for (tissue in c('Blood', 'LesionalSkin', 'NonlesionalSkin')) {

  # TxImport
  dir <- paste(getwd(), tissue, sep='/')
  files <- file.path(dir, pheno$sample, 'MB.abundance.tsv')
  txi <- tximport(files, type='kallisto', tx2gene=t2g, reader=fread)

  # Normalise and filter data for SVA
  dds <- DESeqDataSetFromTximport(txi, colData=pheno, design= ~ 1)
  dds <- estimateSizeFactors(dds)
  mat <- counts(dds, normalized=TRUE)
  mat <- mat[rowMeans(mat) > 1, ] 

  # Run SVA
  mod <- model.matrix(~ 0 + subject + time + wk12.Delta_PASI, data=pheno)
  mod0 <- model.matrix(~ 0 + subject + time, data=pheno)
  svobj <- svaseq(mat, mod, mod0)
  des <- cbind(mod, svobj$sv)

  # Differential expression
  dds <- estimateDispersions(dds, modelMatrix=des, maxit=1000)
  dds <- nbinomWaldTest(dds, betaPrior=FALSE, modelMatrix=des, maxit=1000)
  res <- data.frame(results(dds, name='wk12.Delta_PASI', filterfun=ihw))
  eid <- rownames(res)
  res <- res %>%
    mutate(gene_id = eid) %>%
    inner_join(e2g, by='gene_id') %>%
    arrange(pvalue) %>%
    select(gene_name, baseMean:padj)
  write.csv(res, paste0('wk0_wk12,', tissue, '.csv'), row.names=FALSE)

}


