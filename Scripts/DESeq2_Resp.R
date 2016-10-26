# Load libraries
library(data.table)
library(tximport)
library(DESeq2)
library(sva)
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
dds <- dds[rowSums(mat > 1) >= 9, ]
mat <- mat[rowSums(mat > 5) >= 9, ]
idx <- rownames(dds)

# Define loop
loop <- function(resp, cov) {
  
  # Parallelize
  library(BiocParallel)
  register(MulticoreParam(3))

  # Design
  if (resp == 'Dichotomous') {
    mod <- model.matrix(~ 0 + Time.Tissue + Sex + Age + BMI + HLACW6 + PASI_wk00 +  
                        Time.Tissue:PASI_75, data = pheno)
  } else {
    mod <- model.matrix(~ 0 + Time.Tissue + Sex + Age + BMI + HLACW6 + PASI_wk00 +  
                        Time.Tissue:DeltaPASI, data = pheno)
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
    colnames(des)[15:ncol(des)] <- c(paste(rep(unique(pheno$Time), each = 3), 
                                           unique(pheno$Tissue), 'Response', sep = '.'))
  }

  # DEseq
  dds <- estimateDispersions(dds, modelMatrix = des, maxit = 10000)
  dds <- nbinomWaldTest(dds, modelMatrix = des, maxit = 10000, betaPrior = FALSE)
  
  # Extract
  for (tissue in unique(pheno$Tissue)) {
    
    ### AT TIME ###
    
    for (time in unique(pheno$Time)) {
      
      as.data.frame(results(dds, filterfun = ihw, alpha = 0.05, name = 
                            paste(time, tissue, 'Response', sep = '.'))) %>%
        mutate(gene_id = idx,
               AvgExpr = log2(baseMean)) %>%
        inner_join(e2g, by = 'gene_id') %>%
        rename(EnsemblID  = gene_id,
               GeneSymbol = gene_name, 
               logFC      = log2FoldChange, 
               p.value    = pvalue,
               q.value    = padj) %>%
        arrange(p.value) %>%
        select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value) %>%
        fwrite(paste0('./Results/Response/DESeq2/', 
                      paste('DESeq2', tissue, resp, cov, time, 'txt', 
                            sep = '.')), sep = '\t')
      
    }
    
    ### OVER TIME ###
    
    deltas <- data.frame(Delta01 = c('wk01', 'wk00'), 
                         Delta11 = c('wk12', 'wk01'), 
                         Delta12 = c('wk12', 'wk00'))
    
    for (j in 1:ncol(deltas)) {
      
      as.data.frame(results(dds, filterfun = ihw, alpha = 0.05, contrast = 
                            list(paste(deltas[1, j], tissue, 'Response', sep = '.'),
                                 paste(deltas[2, j], tissue, 'Response', sep = '.')))) %>%
        mutate(gene_id = idx,
               AvgExpr = log2(baseMean)) %>%
        inner_join(e2g, by = 'gene_id') %>%
        rename(EnsemblID  = gene_id,
               GeneSymbol = gene_name, 
               logFC      = log2FoldChange, 
               p.value    = pvalue,
               q.value    = padj) %>%
        arrange(p.value) %>%
        select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value) %>%
        fwrite(paste0('./Results/Response/DESeq2/',
                      paste('DESeq2', tissue, resp, cov, 
                            paste(delta[2], delta[1], sep = '-'), 
                            'txt', sep = '.')), sep = '\t')
      
    }
  
  }

}

# Compute in parallel
library(doParallel)
registerDoParallel(4)
foreach(r = c('Continuous', 'Dichotomous')) %:%
  foreach(c = c('Some', 'All')) %dopar% 
    loop(resp = r, cov = c)



