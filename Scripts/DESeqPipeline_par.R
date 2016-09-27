# Load libraries
library(data.table)
library(tximport)
library(sva)
library(DESeq2)
library(dplyr)

# Prep data
pheno <- fread(paste0(getwd(), '/Data/PrePSORT_Clinical.csv')) %>%
  mutate(wk00 = ifelse(Time == 'wk00', 1, 0),
         wk01 = ifelse(Time == 'wk01', 1, 0),
         wk12 = ifelse(Time == 'wk12', 1, 0),
         wk00.Response = ifelse(Time == 'wk00', DeltaPASI, 0),
         wk01.Response = ifelse(Time == 'wk01', DeltaPASI, 0),
         wk12.Response = ifelse(Time == 'wk12', DeltaPASI, 0))
t2g <- fread(paste0(getwd(), '/Data/Ensembl.Hs79.Tx.csv'))
e2g <- fread(paste0(getwd(), '/Data/Ensembl.Hs79.GeneSymbols.csv'))

# Define loop
loop <- function(tissue) {

  # Subparallelize
  library(BiocParallel)
  register(MulticoreParam(4))

  ### Week Zero ###

  # TxImport
  dir <- paste0(getwd(), '/Data/', tissue)
  files <- file.path(dir, pheno$Sample, 'abundance.tsv')
  txi <- tximport(files, type='kallisto', tx2gene=t2g, reader=fread)
  dds <- DESeqDataSetFromTximport(txi, colData=pheno, design= ~ 1)
  dds <- estimateSizeFactors(dds)
  mat <- counts(dds, normalized=TRUE)
  keep <- rowMeans(mat) > 1
  mat <- mat[keep, ] 
  dds <- dds[keep, ]

  # SVA
  mod <- model.matrix(~ 0 + Time + Sex + Age + BMI + PASI_A + 
                      wk00.Response + wk01.Response + wk12.Response, data=pheno)
  mod0 <- model.matrix(~ 0 + Time + Sex + Age + BMI + PASI_A, data=pheno)
  svobj <- svaseq(mat, mod, mod0)
  des <- cbind(mod, svobj$sv)

  # DEseq
  dds1 <- estimateDispersions(dds, modelMatrix=des, maxit=10000)
  dds1 <- nbinomWaldTest(dds1, betaPrior=FALSE, modelMatrix=des, maxit=10000)

  # Extract
  res <- data.frame(results(dds1, name='wk00.Response', filterfun=ihw))
  id <- rownames(res)
  res <- res %>%
    mutate(gene_id = id) %>%
    inner_join(e2g, by='gene_id') %>%
    rename(EnsemblID  = gene_id,
           GeneSymbol = gene_name, 
           AvgExpr    = baseMean,
           logFC      = log2FoldChange, 
           p.value    = pvalue,
           q.value    = padj) %>%
    arrange(p.value) %>%
    select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
  fwrite(res, paste0(getwd(), '/Results/DESeq/', tissue, ',wk00.txt'), sep='\t')

  ### Week One ###

  # Extract
  res <- data.frame(results(dds1, name='wk01.Response', filterfun=ihw))
  id <- rownames(res)
  res <- res %>%
    mutate(gene_id = id) %>%
    inner_join(e2g, by='gene_id') %>%
    rename(EnsemblID  = gene_id,
           GeneSymbol = gene_name, 
           AvgExpr    = baseMean,
           logFC      = log2FoldChange, 
           p.value    = pvalue,
           q.value    = padj) %>%
    arrange(p.value) %>%
    select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
  fwrite(res, paste0(getwd(), '/Results/DESeq/', tissue, ',wk01.txt'), sep='\t')

  ### Week Twelve ###

  # Extract
  res <- data.frame(results(dds1, name='wk12.Response', filterfun=ihw))
  id <- rownames(res)
  res <- res %>%
    mutate(gene_id = id) %>%
    inner_join(e2g, by='gene_id') %>%
    rename(EnsemblID  = gene_id,
           GeneSymbol = gene_name, 
           AvgExpr    = baseMean,
           logFC      = log2FoldChange, 
           p.value    = pvalue,
           q.value    = padj) %>%
    arrange(p.value) %>%
    select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
  fwrite(res, paste0(getwd(), '/Results/DESeq/', tissue, ',wk12.txt'), sep='\t')

  ### One week change ###

  # SVA
  mod <- model.matrix(~ 0 + Subject + wk01 + wk12 + wk01.Response + wk12.Response, data=pheno)
  mod0 <- model.matrix(~ 0 + Subject + wk01 + wk12, data=pheno)
  svobj <- svaseq(mat, mod, mod0)
  des <- cbind(mod, svobj$sv)

  # DESeq
  dds2 <- estimateDispersions(dds, modelMatrix=des, maxit=10000)
  dds2 <- nbinomWaldTest(dds2, betaPrior=FALSE, modelMatrix=des, maxit=10000)

  # Extract
  res <- data.frame(results(dds2, name='wk01.Response', filterfun=ihw))
  id <- rownames(res)
  res <- res %>%
    mutate(gene_id = id) %>%
    inner_join(e2g, by='gene_id') %>%
    rename(EnsemblID  = gene_id,
           GeneSymbol = gene_name, 
           AvgExpr    = baseMean,
           logFC      = log2FoldChange, 
           p.value    = pvalue,
           q.value    = padj) %>%
    arrange(p.value) %>%
    select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
  fwrite(res, paste0(getwd(), '/Results/DESeq/', tissue, ',wk00-wk01.txt'), sep='\t')

  ### Twelve Week Change ###

  # Extract
  res <- data.frame(results(dds2, name='wk12.Response', filterfun=ihw))
  id <- rownames(res)
  res <- res %>%
    mutate(gene_id = id) %>%
    inner_join(e2g, by='gene_id') %>%
    rename(EnsemblID  = gene_id,
           GeneSymbol = gene_name, 
           AvgExpr    = baseMean,
           logFC      = log2FoldChange, 
           p.value    = pvalue,
           q.value    = padj) %>%
    arrange(p.value) %>%
    select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
  fwrite(res, paste0(getwd(), '/Results/DESeq/', tissue, ',wk00-wk12.txt'), sep='\t')

  ### Eleven Week Change ###

  # SVA
  mod <- model.matrix(~ 0 + Subject + wk00 + wk12 + wk00.Response + wk12.Response, data=pheno)
  mod0 <- model.matrix(~ 0 + Subject + wk00 + wk12, data=pheno)
  svobj <- svaseq(mat, mod, mod0)
  des <- cbind(mod, svobj$sv)

  # DESeq
  dds3 <- estimateDispersions(dds, modelMatrix=des, maxit=10000)
  dds3 <- nbinomWaldTest(dds3, betaPrior=FALSE, modelMatrix=des, maxit=10000)

  # Extract
  res <- data.frame(results(dds3, name='wk12.Response', filterfun=ihw))
  id <- rownames(res)
  res <- res %>%
    mutate(gene_id = id) %>%
    inner_join(e2g, by='gene_id') %>%
    rename(EnsemblID  = gene_id,
           GeneSymbol = gene_name, 
           AvgExpr    = baseMean,
           logFC      = log2FoldChange, 
           p.value    = pvalue,
           q.value    = padj) %>%
    arrange(p.value) %>%
    select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
  fwrite(res, paste0(getwd(), '/Results/DESeq/', tissue, ',wk01-wk12.txt'), sep='\t')

}

# Compute in parallel
library(doParallel)
registerDoParallel(cores=3)
foreach(i=c('Blood', 'LesionalSkin', 'NonlesionalSkin')) %dopar% loop(i)


