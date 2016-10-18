# Load libraries
library(data.table)
library(tximport)
library(sva)
library(DESeq2)
library(dplyr)

# Prep data
pheno <- fread('./Data/Clinical.csv') %>%
  mutate(Time.Tissue = paste(Time, Tissue, sep = '.'))
t2g <- fread('./Data/Ensembl.Hs79.Tx.csv')
e2g <- fread('./Data/Ensembl.Hs79.GeneSymbols.csv')

# TxImport
files <- file.path('./Data/RawCounts', pheno$Sample, 'abundance.tsv')
txi <- tximport(files, type = 'kallisto', tx2gene = t2g, reader = fread)
dds <- DESeqDataSetFromTximport(txi, colData = pheno, design = ~ 1)
dds <- estimateSizeFactors(dds)
mat <- counts(dds, normalized = TRUE)
keep <- rowMeans(mat) > 1
mat <- mat[keep, ] 
dds <- dds[keep, ]
rld <- assay(rlog(dds))
id <- rownames(rld)
means <- rowMeans(rld)

# Define loop
loop <- function(resp, cov) {
  
  # Parallelize
  library(BiocParallel)
  register(MulticoreParam(2))

  # Design
  if (resp == 'Dichotomous') {
    if (cov == 'None') {
      mod <- model.matrix(~ 0 + Time.Tissue + Time.Tissue:PASI_75, data = pheno)
    } else {
      mod <- model.matrix(~ 0 + Time.Tissue + Sex + Age + BMI + HLACW6 + PASI_wk00 +  
                          Time.Tissue:PASI_75, data = pheno)
    }
  } else if (resp == 'Continuous') {
    if (cov == 'None') {
      mod <- model.matrix(~ 0 + Time.Tissue + Time.Tissue:DeltaPASI, data = pheno)
    } else {
      mod <- model.matrix(~ 0 + Time.Tissue + Sex + Age + BMI + HLACW6 + PASI_wk00 +  
                          Time.Tissue:DeltaPASI, data = pheno)
    }
  }
  if (cov == 'All') {
    mod0 <- model.matrix(~ 0 + Time.Tissue + Sex + Age + BMI + HLACW6 + PASI_wk00, 
                         data = pheno)
    svobj <- svaseq(mat, mod, mod0)
    des <- cbind(mod, svobj$sv)
    colnames(des)[15:ncol(des)] <- c(paste(rep(unique(pheno$Time), each = 3), 
                                           unique(pheno$Tissue), 'Response', sep = '.'),
                                     paste0('SV', 1:svobj$n.sv))
  } else {
    des <- mod
    colnames(des)[(ncol(des) - 8):ncol(des)] <- 
      c(paste(rep(unique(pheno$Time), each = 3), 
              unique(pheno$Tissue), 'Response', sep = '.'))
  }

  # DEseq
  dds <- estimateDispersions(dds, modelMatrix = des, maxit = 10000)
  dds <- nbinomWaldTest(dds, modelMatrix = des, maxit = 10000, betaPrior = FALSE)
  
  for (tissue in unique(pheno$Tissue)) {
    
    ### AT TIME ###

    # wk0
    res <- data.frame(results(dds, filterfun = ihw,
                      name = paste0('wk00.', tissue, '.Response'))) %>%
      mutate(gene_id = id,
             AvgExpr = means) %>%
      inner_join(e2g, by = 'gene_id') %>%
      rename(EnsemblID  = gene_id,
             GeneSymbol = gene_name, 
             logFC      = log2FoldChange, 
             p.value    = pvalue,
             q.value    = padj) %>%
      arrange(p.value) %>%
      select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
    fwrite(res, paste0('./Results/Response/DESeq2/',
           paste('DESeq2', tissue, resp, cov, 'wk00.txt', 
           sep = '.')), sep = '\t')
  
    # wk01
    res <- data.frame(results(dds, filterfun = ihw,
                              name = paste0('wk01.', tissue, '.Response'))) %>%
      mutate(gene_id = id,
             AvgExpr = means) %>%
      inner_join(e2g, by = 'gene_id') %>%
      rename(EnsemblID  = gene_id,
             GeneSymbol = gene_name, 
             logFC      = log2FoldChange, 
             p.value    = pvalue,
             q.value    = padj) %>%
      arrange(p.value) %>%
      select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
    fwrite(res, paste0('./Results/Response/DESeq2/',
           paste('DESeq2', tissue, resp, cov, 'wk01.txt', 
           sep = '.')), sep = '\t')
  
    # wk12
    res <- data.frame(results(dds, filterfun = ihw,
                              name = paste0('wk12.', tissue, '.Response'))) %>%
      mutate(gene_id = id,
             AvgExpr = means) %>%
      inner_join(e2g, by = 'gene_id') %>%
      rename(EnsemblID  = gene_id,
             GeneSymbol = gene_name, 
             logFC      = log2FoldChange, 
             p.value    = pvalue,
             q.value    = padj) %>%
      arrange(p.value) %>%
      select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
    fwrite(res, paste0('./Results/Response/DESeq2/',
           paste('DESeq2', tissue, resp, cov, 'wk12.txt', 
           sep = '.')), sep = '\t')
    
    ### OVER TIME ###
    
    # Delta01
    res <- data.frame(results(dds, filterfun = ihw, contrast = 
                              list(paste0('wk01.', tissue, '.Response'), 
                                   paste0('wk00.', tissue, '.Response')))) %>%
      mutate(gene_id = id,
             AvgExpr = means) %>%
      inner_join(e2g, by = 'gene_id') %>%
      rename(EnsemblID  = gene_id,
             GeneSymbol = gene_name, 
             logFC      = log2FoldChange, 
             p.value    = pvalue,
             q.value    = padj) %>%
      arrange(p.value) %>%
      select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
    fwrite(res, paste0('./Results/Response/DESeq2/',
                       paste('DESeq2', tissue, resp, cov, 'wk00-wk01.txt', 
                             sep = '.')), sep = '\t')
    
    # Delta11
    res <- data.frame(results(dds, filterfun = ihw, contrast = 
                              list(paste0('wk12.', tissue, '.Response'), 
                                   paste0('wk01.', tissue, '.Response')))) %>%
      mutate(gene_id = id,
             AvgExpr = means) %>%
      inner_join(e2g, by = 'gene_id') %>%
      rename(EnsemblID  = gene_id,
             GeneSymbol = gene_name, 
             logFC      = log2FoldChange, 
             p.value    = pvalue,
             q.value    = padj) %>%
      arrange(p.value) %>%
      select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
    fwrite(res, paste0('./Results/Response/DESeq2/',
                       paste('DESeq2.Blood', resp, cov, 'wk01-wk12.txt', 
                             sep = '.')), sep = '\t')
    
    # Delta12
    res <- data.frame(results(dds, filterfun = ihw, contrast = 
                              list(paste0('wk12.', tissue, '.Response'), 
                                   paste0('wk00.', tissue, '.Response')))) %>%
      mutate(gene_id = id,
             AvgExpr = means) %>%
      inner_join(e2g, by = 'gene_id') %>%
      rename(EnsemblID  = gene_id,
             GeneSymbol = gene_name, 
             logFC      = log2FoldChange, 
             p.value    = pvalue,
             q.value    = padj) %>%
      arrange(p.value) %>%
      select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
    fwrite(res, paste0('./Results/Response/DESeq2/',
                       paste('DESeq2', tissue, resp, cov, 'wk00-wk12.txt', 
                             sep = '.')), sep = '\t')
  
  }

}

# Compute in parallel
library(doParallel)
registerDoParallel(6)
foreach(r = c('Continuous', 'Dichotomous')) %:%
  foreach(c = c('None', 'Some', 'All')) %dopar% 
    loop(resp = r, cov = c)


