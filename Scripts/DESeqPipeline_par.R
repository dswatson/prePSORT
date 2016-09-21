# Load libraries
library(data.table)
library(tximport)
library(sva)
library(DESeq2)
library(dplyr)

# Prep data
pheno <- fread(paste0(getwd(), '/Data/PrePSORT_Clinical.csv')) %>%
  mutate(wk0.DeltaPASI = ifelse(time == 'wk0', Delta_PASI, 0),
         wk1.DeltaPASI = ifelse(time == 'wk1', Delta_PASI, 0),
         wk12.DeltaPASI = ifelse(time == 'wk12', Delta_PASI, 0))
t2g <- fread(paste0(getwd(), '/Data/Ensembl.Hs79.Tx.csv'))
e2g <- fread(paste0(getwd(), '/Data/Ensembl.Hs79.GeneSymbols.csv'))

# Define loop
loop <- function(tissue) {

  # Subparallelize
  library(BiocParallel)
  register(MulticoreParam(4))

  ### Week Zero ###

  # TxImport
  dir <- paste(getwd(), 'Data', tissue, sep='/')
  files <- file.path(dir, pheno$sample, 'MB.abundance.tsv')
  txi <- tximport(files, type='kallisto', tx2gene=t2g, reader=fread)

  # Normalise, filter
  dds <- DESeqDataSetFromTximport(txi, colData=pheno, design= ~ 1)
  dds <- estimateSizeFactors(dds)
  mat <- counts(dds, normalized=TRUE)
  mat <- mat[rowMeans(mat) > 1, ] 

  # SVA
  mod <- model.matrix(~ 0 + time + sex + age + bmi +
                      wk0.DeltaPASI + wk1.DeltaPASI + wk12.DeltaPASI, data=pheno)
  mod0 <- model.matrix(~ 0 + time + sex + age + bmi, data=pheno)
  svobj <- svaseq(mat, mod, mod0)
  des <- cbind(mod, svobj$sv)

  # DEseq
  dds <- estimateDispersions(dds, modelMatrix=des, maxit=1000)
  dds <- nbinomWaldTest(dds, betaPrior=FALSE, modelMatrix=des, maxit=1000)
  res <- data.frame(results(dds, name='wk0.DeltaPASI', filterfun=ihw))
  id <- rownames(res)
  res <- res %>%
    mutate(gene_id = id) %>%
    inner_join(e2g, by='gene_id') %>%
    rename(EnsemblID = gene_id,
           GeneSymbol = gene_name, 
           AvgExpr = baseMean,
           logFC = log2FoldChange, 
           p.value = pvalue,
           q.value = padj) %>%
    arrange(p.value) %>%
    select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
  fwrite(res, paste0(getwd(), '/Results/DESeq/Response/', tissue, ',wk0.txt'), sep='\t')

  ### Week One ###

  # DESeq
  res <- data.frame(results(dds, name='wk1.DeltaPASI', filterfun=ihw))
  id <- rownames(res)
  res <- res %>%
    mutate(gene_id = id) %>%
    inner_join(e2g, by='gene_id') %>%
    rename(EnsemblID = gene_id,
           GeneSymbol = gene_name, 
           AvgExpr = baseMean,
           logFC = log2FoldChange, 
           p.value = pvalue,
           q.value = padj) %>%
    arrange(p.value) %>%
    select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
  fwrite(res, paste0(getwd(), '/Results/DESeq/Response/', tissue, ',wk1.txt'), sep='\t')

  ### Week Twelve ###

  # DESeq
  res <- data.frame(results(dds, name='wk12.DeltaPASI', filterfun=ihw))
  id <- rownames(res)
  res <- res %>%
    mutate(gene_id = id) %>%
    inner_join(e2g, by='gene_id') %>%
    rename(EnsemblID = gene_id,
           GeneSymbol = gene_name, 
           AvgExpr = baseMean,
           logFC = log2FoldChange, 
           p.value = pvalue,
           q.value = padj) %>%
    arrange(p.value) %>%
    select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
  fwrite(res, paste0(getwd(), '/Results/DESeq/Response/', tissue, ',wk12.txt'), sep='\t')

  ### One week change ###

  # TxImport
  pheno2 <- filter(pheno, time != 'wk12')
  files <- file.path(dir, pheno2$sample, 'MB.abundance.tsv')
  txi <- tximport(files, type='kallisto', tx2gene=t2g, reader=fread) 

  # Normalise, filter
  dds <- DESeqDataSetFromTximport(txi, colData=pheno2, design= ~ 1)
  dds <- estimateSizeFactors(dds)
  mat <- counts(dds, normalized=TRUE)
  mat <- mat[rowMeans(mat) > 1, ] 

  # SVA
  mod <- model.matrix(~ 0 + time + subject + wk1.DeltaPASI, data=pheno2)
  mod0 <- model.matrix(~ 0 + time + subject, data=pheno2)
  svobj <- svaseq(mat, mod, mod0)
  des <- cbind(mod, svobj$sv)

  # DESeq
  dds <- estimateDispersions(dds, modelMatrix=des, maxit=1000)
  dds <- nbinomWaldTest(dds, betaPrior=FALSE, modelMatrix=des, maxit=1000)
  res <- data.frame(results(dds, name='wk1.DeltaPASI', filterfun=ihw))
  id <- rownames(res)
  res <- res %>%
    mutate(gene_id = id) %>%
    inner_join(e2g, by='gene_id') %>%
    rename(EnsemblID = gene_id,
           GeneSymbol = gene_name, 
           AvgExpr = baseMean,
           logFC = log2FoldChange, 
           p.value = pvalue,
           q.value = padj) %>%
    arrange(p.value) %>%
    select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
  fwrite(res, paste0(getwd(), '/Results/DESeq/Response/', tissue, ',wk0-wk1.txt'), sep='\t')

  ### Eleven Week Change ###

  # TxImport
  pheno2 <- filter(pheno, time != 'wk0')
  files <- file.path(dir, pheno2$sample, 'MB.abundance.tsv')
  txi <- tximport(files, type='kallisto', tx2gene=t2g, reader=fread)

  # Normalise, filter
  dds <- DESeqDataSetFromTximport(txi, colData=pheno2, design= ~ 1)
  dds <- estimateSizeFactors(dds)
  mat <- counts(dds, normalized=TRUE)
  mat <- mat[rowMeans(mat) > 1, ] 

  # SVA 
  mod <- model.matrix(~ 0 + time + subject + wk12.DeltaPASI, data=pheno2)
  mod0 <- model.matrix(~ 0 + time + subject, data=pheno2)
  svobj <- svaseq(mat, mod, mod0)
  des <- cbind(mod, svobj$sv)

  # DESeq
  dds <- estimateDispersions(dds, modelMatrix=des, maxit=1000)
  dds <- nbinomWaldTest(dds, betaPrior=FALSE, modelMatrix=des, maxit=1000)
  res <- data.frame(results(dds, name='wk12.DeltaPASI', filterfun=ihw))
  id <- rownames(res)
  res <- res %>%
    mutate(gene_id = id) %>%
    inner_join(e2g, by='gene_id') %>%
    rename(EnsemblID = gene_id,
           GeneSymbol = gene_name, 
           AvgExpr = baseMean,
           logFC = log2FoldChange, 
           p.value = pvalue,
           q.value = padj) %>%
    arrange(p.value) %>%
    select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
  fwrite(res, paste0(getwd(), '/Results/DESeq/Response/', tissue, ',wk1-wk12.txt'), sep='\t')

  ### Twelve Week Change ###

  # TxImport
  pheno2 <- filter(pheno, time != 'wk1')
  files <- file.path(dir, pheno2$sample, 'MB.abundance.tsv')
  txi <- tximport(files, type='kallisto', tx2gene=t2g, reader=fread)

  # Normalise, filter
  dds <- DESeqDataSetFromTximport(txi, colData=pheno2, design= ~ 1)
  dds <- estimateSizeFactors(dds)
  mat <- counts(dds, normalized=TRUE)
  mat <- mat[rowMeans(mat) > 1, ] 

  # SVA
  mod <- model.matrix(~ 0 + time + subject + wk12.DeltaPASI, data=pheno2)
  mod0 <- model.matrix(~ 0 + time + subject, data=pheno2)
  svobj <- svaseq(mat, mod, mod0)
  des <- cbind(mod, svobj$sv)

  # DESeq
  dds <- estimateDispersions(dds, modelMatrix=des, maxit=1000)
  dds <- nbinomWaldTest(dds, betaPrior=FALSE, modelMatrix=des, maxit=1000)
  res <- data.frame(results(dds, name='wk12.DeltaPASI', filterfun=ihw))
  id <- rownames(res)
  res <- res %>%
    mutate(gene_id = id) %>%
    inner_join(e2g, by='gene_id') %>%
    rename(EnsemblID = gene_id,
           GeneSymbol = gene_name, 
           AvgExpr = baseMean,
           logFC = log2FoldChange, 
           p.value = pvalue,
           q.value = padj) %>%
    arrange(p.value) %>%
    select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
  fwrite(res, paste0(getwd(), '/Results/DESeq/Response/', tissue, ',wk0-wk12.txt'), sep='\t')

}

# Compute in parallel
library(doParallel)
registerDoParallel(cores=3)
foreach(i=c('Blood', 'LesionalSkin', 'NonlesionalSkin')) %dopar% loop(i)


