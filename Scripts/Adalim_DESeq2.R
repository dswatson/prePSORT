# Load libraries
library(data.table)
library(tximport)
library(DESeq2)
library(sva)
library(dplyr)

# Prep data
pheno <- fread('./Data/Adalimumab/Clinical.csv') %>%
  mutate(Time.Tissue = paste(Time, Tissue, sep = '.'))
t2g <- fread('./Data/Ensembl.Hs79.Tx.csv')
e2g <- fread('./Data/Ensembl.Hs79.GeneSymbols.csv')

# TxImport
files <- file.path('./Data/Adalimumab/Counts', pheno$Sample, 'MB.abundance.tsv')
txi <- tximport(files, type = 'kallisto', tx2gene = t2g, reader = fread)
dds <- DESeqDataSetFromTximport(txi, colData = pheno, design = ~ 1)
dds <- estimateSizeFactors(dds)
mat <- counts(dds, normalized = TRUE)
dds <- dds[rowSums(mat > 1) >= 2, ]
mat <- mat[rowSums(mat > 5) >= 2, ]
idx <- rownames(dds)

# SVA
mod <- model.matrix(~ 0 + Time.Tissue + PASI_wk00 + Batch + 
                      Time.Tissue:DeltaPASI, data = pheno)
mod <- mod[, -c(13:15)]  # Remove response for wks 1, 4, 7
mod0 <- mod[, !grepl('DeltaPASI', colnames(mod))]
svobj <- svaseq(mat, mod, mod0)
des <- cbind(mod, svobj$sv)
colnames(des)[11:ncol(des)] <- c(paste(rep(c('wk00', 'wk12'), each = 2),
                                       rep(c('Lesional', 'Nonlesional'), times = 2), 
                                       'Response', sep = '.'),
                                 paste0('SV', 1:svobj$n.sv))

# DEseq
dds <- estimateDispersions(dds, modelMatrix = des, maxit = 10000)
dds <- nbinomWaldTest(dds, modelMatrix = des, maxit = 10000, betaPrior = FALSE)

# Define loop
loop <- function(tissue, time) {
  
  ### AT TIME ###
    
  data.frame(results(dds, filterfun = ihw, name = 
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
               fwrite(paste0('./Results/Adalimumab/', 
                             paste('DESeq2', tissue, time, 'txt', 
                                   sep = '.')), sep = '\t')
  
  ### OVER TIME ###
    
  data.frame(results(dds, filterfun = ihw, contrast = 
                       list(paste('wk12', tissue, 'Response', sep = '.'),
                            paste('wk00', tissue, 'Response', sep = '.')))) %>%
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
    fwrite(paste0('./Results/Adalimumab/',
                  paste0('DESeq2.', tissue, 'Delta12.txt')), sep = '\t')
  
}

# Compute in parallel
library(doParallel)
registerDoParallel(4)
foreach(i = c('Lesional', 'Nonlesional')) %:%
  foreach(j = c('wk00', 'wk12')) %dopar% loop(i, j)


