# Load libraries
library(data.table)
library(tximport)
library(edgeR)
library(limma)
library(sva)
library(qvalue)
library(dplyr)

# Prep data
pheno <- fread('./Data/Adalimumab/Clinical.csv') %>%
  mutate(Time.Tissue = paste(Time, Tissue, sep = '.'))
t2g <- fread('./Data/Ensembl.Hs79.Tx.csv')
e2g <- fread('./Data/Ensembl.Hs79.GeneSymbols.csv')

# TxImport
files <- file.path('./Data/Adalimumab/Counts', pheno$Sample, 'MB.abundance.tsv')
txi <- tximport(files, type = 'kallisto', tx2gene = t2g, reader = fread, 
                countsFromAbundance = 'lengthScaledTPM')
keep <- rowSums(cpm(txi$counts) > 1) >= 2
y <- DGEList(txi$counts[keep, ])
y <- calcNormFactors(y)

# SVA
mod <- model.matrix(~ 0 + Time.Tissue + PASI_wk00 + Batch + 
                      Time.Tissue:DeltaPASI, data = pheno)
mod <- mod[, -c(13:15)]  # Remove response for wks 1, 4, 7
mod0 <- mod[, !grepl('DeltaPASI', colnames(mod))]
svobj <- svaseq(cpm(y), mod, mod0)
des <- cbind(mod, svobj$sv)
colnames(des)[11:ncol(des)] <- c(paste(rep(c('wk00', 'wk12'), each = 2),
                                       rep(c('Lesional', 'Nonlesional'), times = 2), 
                                       'Response', sep = '.'),
                                 paste0('SV', 1:svobj$n.sv))

# Voom
v <- voom(y, des)
idx <- rownames(v)
urFit <- lmFit(v, des)

# Define loop
loop <- function(tissue, time)  {
  
  ### AT TIME ###
  
  fit <- eBayes(urFit)
  topTable(fit, number = Inf, sort.by = 'none',
           coef = paste(time, tissue, 'Response', sep = '.')) %>%
    mutate(q.value = qvalue(P.Value)$qvalues, 
           gene_id = idx) %>%
    inner_join(e2g, by = 'gene_id') %>%
    rename(EnsemblID  = gene_id,
           GeneSymbol = gene_name,
           p.value    = P.Value, 
           AvgExpr    = AveExpr) %>%
    arrange(p.value) %>%
    select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value) %>%
    fwrite(paste0('./Results/Adalimumab/', 
                  paste('voom', tissue, time, 'txt', sep = '.')), sep = '\t')
  
  ### OVER TIME ###
  
  cm <- makeContrasts('Lesional.Delta12' = wk12.Lesional.Response - wk00.Lesional.Response,
                      'Nonlesional.Delta12' = wk12.Nonlesional.Response - wk00.Nonlesional.Response, 
                      levels = des)
  fit <- eBayes(contrasts.fit(urFit, cm))
  topTable(fit, number = Inf, sort.by = 'none',
           coef = paste0(tissue, '.Delta12')) %>%
    mutate(q.value = qvalue(P.Value)$qvalues, 
           gene_id = idx) %>%
    inner_join(e2g, by = 'gene_id') %>%
    rename(EnsemblID  = gene_id,
           GeneSymbol = gene_name,
           p.value    = P.Value, 
           AvgExpr    = AveExpr) %>%
    arrange(p.value) %>%
    select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value) %>%
    fwrite(paste0('./Results/Adalimumab/', 
                  paste0('voom.', tissue, '.Delta12.txt')), sep = '\t')
  
}

# Compute in parallel
library(doParallel)
registerDoParallel(4)
foreach(i = c('Lesional', 'Nonlesional')) %:%
  foreach(j = c('wk00', 'wk12')) %dopar% loop(i, j)


