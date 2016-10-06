# Load libraries
library(data.table)
library(tximport)
library(sva)
library(DESeq2)
library(dplyr)

# Prep data
pheno <- fread(paste0(getwd(), '/Data/PrePSORT_Clinical.csv'))
t2g <- fread(paste0(getwd(), '/Data/Ensembl.Hs79.Tx.csv'))
e2g <- fread(paste0(getwd(), '/Data/Ensembl.Hs79.GeneSymbols.csv'))

# Define loop
loop <- function(tissue, resp, cov) {

  # TxImport
  dir <- paste0(getwd(), '/Data/', tissue)
  files <- file.path(dir, pheno$Sample, 'abundance.tsv')
  txi <- tximport(files, type = 'kallisto', tx2gene = t2g, reader = fread)
  dds <- DESeqDataSetFromTximport(txi, colData = pheno, design = ~ 1)
  dds <- estimateSizeFactors(dds)
  mat <- counts(dds, normalized = TRUE)
  keep <- rowMeans(mat) > 1
  mat <- mat[keep, ] 
  dds <- dds[keep, ]

  # SVA
  if (resp == 'PASI_75') {
    if (cov == 'noCov') {
      mod <- model.matrix(~ 0 + Time + Time:PASI_75)
    }
    else {
      mod <- model.matrix(~ 0 + Time + Sex + Age + BMI + HLACW6 + 
                          PASI_wk00 + Time:PASI_75, data = pheno)
    }
  } else if (resp == 'DeltaPASI') {
    if (cov == 'noCov') {
      mod <- model.matrix(~ 0 + Time + Time:DeltaPASI, data = pheno)
      } else {
        mod <- model.matrix(~ 0 + Time + Sex + Age + BMI + 
                            PASI_wk00 + Time:DeltaPASI, data = pheno)
      }
  }
  if (cov == 'clinAndSV') {
    mod0 <- model.matrix(~ 0 + Time + Sex + Age + BMI + PASI_wk00, 
                         data = pheno)
    svobj <- svaseq(mat, mod, mod0)
    des <- cbind(mod, svobj$sv)
    colnames(des)[8:ncol(des)] <- c('wk00.Response', 'wk01.Response', 'wk12.Response',
                                    paste0('SV', 1:svobj$n.sv))
  } else {
    des <- mod
    colnames(des)[8:ncol(des)] <- c('wk00.Response', 'wk01.Response', 'wk12.Response')
  }

  # DEseq
  dds <- estimateDispersions(dds, modelMatrix = des, maxit = 10000)
  dds <- nbinomWaldTest(dds, modelMatrix = des, maxit = 10000, betaPrior = FALSE)

  ### Week Zero ###

  # Extract
  res <- data.frame(results(dds, 
                            name = 'wk00.Response', 
                            filterfun = ihw))
  id <- rownames(res)
  res <- res %>%
    mutate(gene_id = id) %>%
    inner_join(e2g, by = 'gene_id') %>%
    rename(EnsemblID  = gene_id,
           GeneSymbol = gene_name, 
           AvgExpr    = baseMean,
           logFC      = log2FoldChange, 
           p.value    = pvalue,
           q.value    = padj) %>%
    arrange(p.value) %>%
    select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
  fwrite(res, paste(getwd(), 'Results/WithinTissue/DESeq2', resp, tissue, cov,
         paste('DESeq2', resp, tissue, cov, 'wk00.txt', sep = '.'), sep = '/'),
         sep = '\t')

  ### Week One ###

  # Extract
  res <- data.frame(results(dds, 
                            name = 'wk01.Response', 
                            filterfun = ihw))
  id <- rownames(res)
  res <- res %>%
    mutate(gene_id = id) %>%
    inner_join(e2g, by = 'gene_id') %>%
    rename(EnsemblID  = gene_id,
           GeneSymbol = gene_name, 
           AvgExpr    = baseMean,
           logFC      = log2FoldChange, 
           p.value    = pvalue,
           q.value    = padj) %>%
    arrange(p.value) %>%
    select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
  fwrite(res, paste(getwd(), 'Results/WithinTissue/DESeq2', resp, tissue, cov,
         paste('DESeq2', resp, tissue, cov, 'wk01.txt', sep = '.'), sep = '/'),
         sep = '\t')

  ### Week Twelve ###

  # Extract
  res <- data.frame(results(dds, 
                            name = 'wk12.Response', 
                            filterfun = ihw))
  id <- rownames(res)
  res <- res %>%
    mutate(gene_id = id) %>%
    inner_join(e2g, by = 'gene_id') %>%
    rename(EnsemblID  = gene_id,
           GeneSymbol = gene_name, 
           AvgExpr    = baseMean,
           logFC      = log2FoldChange, 
           p.value    = pvalue,
           q.value    = padj) %>%
    arrange(p.value) %>%
    select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
  fwrite(res, paste(getwd(), 'Results/WithinTissue/DESeq2', resp, tissue, cov,
         paste('DESeq2', resp, tissue, cov, 'wk12.txt', sep = '.'), sep = '/'),
         sep = '\t')

  ### One week change ###

  # Extract
  res <- data.frame(results(dds, 
                            contrast  = list('wk01.Response', 'wk00.Response'), 
                            filterfun = ihw))
  id <- rownames(res)
  res <- res %>%
    mutate(gene_id = id) %>%
    inner_join(e2g, by = 'gene_id') %>%
    rename(EnsemblID  = gene_id,
           GeneSymbol = gene_name, 
           AvgExpr    = baseMean,
           logFC      = log2FoldChange, 
           p.value    = pvalue,
           q.value    = padj) %>%
    arrange(p.value) %>%
    select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
  fwrite(res, paste(getwd(), 'Results/WithinTissue/DESeq2', resp, tissue, cov,
         paste('DESeq2', resp, tissue, cov, 'wk00-wk01.txt', sep = '.'), sep = '/'),
         sep = '\t')

  ### Eleven Week Change ###

  # Extract
  res <- data.frame(results(dds, 
                            contrast = list('wk12.Response', 'wk01.Response'), 
                            filterfun = ihw))
  id <- rownames(res)
  res <- res %>%
    mutate(gene_id = id) %>%
    inner_join(e2g, by = 'gene_id') %>%
    rename(EnsemblID  = gene_id,
           GeneSymbol = gene_name, 
           AvgExpr    = baseMean,
           logFC      = log2FoldChange, 
           p.value    = pvalue,
           q.value    = padj) %>%
    arrange(p.value) %>%
    select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
  fwrite(res, paste(getwd(), 'Results/WithinTissue/DESeq2', resp, tissue, cov,
         paste('DESeq2', resp, tissue, cov, 'wk01-wk12.txt', sep = '.'), sep = '/'),
         sep = '\t')

  ### Twelve Week Change ###

  # Extract
  res <- data.frame(results(dds, 
                            contrast = list('wk12.Response', 'wk00.Response'), 
                            filterfun = ihw))
  id <- rownames(res)
  res <- res %>%
    mutate(gene_id = id) %>%
    inner_join(e2g, by = 'gene_id') %>%
    rename(EnsemblID  = gene_id,
           GeneSymbol = gene_name, 
           AvgExpr    = baseMean,
           logFC      = log2FoldChange, 
           p.value    = pvalue,
           q.value    = padj) %>%
    arrange(p.value) %>%
    select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
  fwrite(res, paste(getwd(), 'Results/WithinTissue/DESeq2', resp, tissue, cov,
         paste('DESeq2', resp, tissue, cov, 'wk00-wk12.txt', sep = '.'), sep = '/'),
         sep = '\t')

}

# Compute in parallel
library(doParallel)
registerDoParallel(cores = 12)
foreach(t = c('Blood', 'LesionalSkin', 'NonlesionalSkin')) %:%
  foreach(r = c('PASI_75', 'DeltaPASI')) %:%
    foreach(c = c('noCov', 'clinCov', 'clinAndSV')) %dopar% 
      loop(tissue = t, resp = r, cov = c)


