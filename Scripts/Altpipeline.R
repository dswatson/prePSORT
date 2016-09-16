# Load libraries
library(data.table)
library(tximport)
library(sva)
library(DESeq2)
library(dplyr)

# Prep data
pheno <- read.csv(paste0(getwd(), '/Data/Clinical.csv'), stringsAsFactors=FALSE)
t2g <- fread(paste0(getwd(), '/Data/Ensembl.Hs79.Tx.csv'))
e2g <- fread(paste0(getwd(), '/Data/Ensembl.Hs79.GeneSymbols.csv'))

# Loop
for (tissue in c('Blood', 'LesionalSkin', 'NonlesionalSkin')) {

  # TxImport
  dir <- paste(getwd(), 'Data', tissue, sep='/')
  files <- file.path(dir, pheno$sample, 'MB.abundance.tsv')
  txi <- tximport(files, type='kallisto', tx2gene=t2g, reader=fread)

  # Normalise and filter matrix for SVA
  dds <- DESeqDataSetFromTximport(txi, colData=pheno, design= ~ 1)
  dds <- estimateSizeFactors(dds)
  mat <- counts(dds, normalized=TRUE)
  mat <- mat[rowMeans(mat) > 1, ] 

  # Run SVA
  mod <- model.matrix(~ sex + age + bmi + time:Delta_PASI, data=pheno)
  mod0 <- model.matrix(~ sex + age + bmi, data=pheno)
  svobj <- svaseq(mat, mod, mod0)
  des <- cbind(mod, svobj$sv)
  colnames(des)[5:(7 + svobj$n.sv)] <- c('wk0.DeltaPASI', 'wk1.DeltaPASI', 'wk12.DeltaPASI', 
                                         paste0('SV', 1:svobj$n.sv))

  # DESeq
  dds <- estimateDispersions(dds, modelMatrix=des, maxit=1000)
  dds <- nbinomWaldTest(dds, betaPrior=FALSE, modelMatrix=des, maxit=1000)

  # Baseline
  res <- data.frame(results(dds, name='wk0.DeltaPASI', filterfun=ihw))
  id <- rownames(res)
  res <- res %>%
    mutate(gene_id = id) %>%
    inner_join(e2g, by='gene_id') %>%
    arrange(pvalue) %>%
    select(gene_name, baseMean:padj)
  fwrite(res, paste0('Baseline,', tissue, '.csv'))

  # One week change
  res <- data.frame(results(dds, contrast=list('wk0.DeltaPASI', 'wk1.DeltaPASI'), filterfun=ihw))
  id <- rownames(res)
  res <- res %>%
    mutate(gene_id = id) %>%
    inner_join(e2g, by='gene_id') %>%
    arrange(pvalue) %>%
    select(gene_name, baseMean:padj)
  fwrite(res, paste0('wk0_wk1,', tissue, '.csv'))

  # Twelve week change
  res <- data.frame(results(dds, contrast=list('wk0.DeltaPASI', 'wk12.DeltaPASI'), filterfun=ihw))
  id <- rownames(res)
  res <- res %>%
    mutate(gene_id = id) %>%
    inner_join(e2g, by='gene_id') %>%
    arrange(pvalue) %>%
    select(gene_name, baseMean:padj)
  fwrite(res, paste0('wk0_wk12,', tissue, '.csv'))

}


