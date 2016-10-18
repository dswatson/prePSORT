# Load libraries
library(data.table)
library(tximport)
library(edgeR)
library(sva)
library(qvalue)
library(dplyr)

# Prep data
pheno <- fread('./Data/Clinical.csv') %>%
  mutate(Time.Tissue = paste(Time, Tissue, sep = '.'))
t2g <- fread('./Data/Ensembl.Hs79.Tx.csv')
e2g <- fread('./Data/Ensembl.Hs79.GeneSymbols.csv')

# TxImport
files <- file.path('./Data/RawCounts', pheno$Sample, 'abundance.tsv')
txi <- tximport(files, type = 'kallisto', tx2gene = t2g, reader = fread, 
                countsFromAbundance = 'lengthScaledTPM')
keep <- rowSums(cpm(txi$counts) > 1) >= 3
y <- DGEList(txi$counts[keep, ])
y <- calcNormFactors(y)

# Define loop
loop <- function(resp, cov) {

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
    svobj <- svaseq(cpm(y), mod, mod0)
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

  # Voom
  v <- voom(y, des)
  fit <- eBayes(lmFit(v, des))

  ### AT TIME ###

  for (tissue in unique(pheno$Tissue)) {

    # wk0
    top <- topTable(fit, number = Inf, sort.by = 'p',
                    coef = paste0('wk00.', tissue, '.Response'))
    id <- rownames(top)
    top <- top %>%
      mutate(q.value = qvalue(P.Value)$qvalues, gene_id = id) %>%
      inner_join(e2g, by = 'gene_id') %>%
      rename(EnsemblID  = gene_id,
             GeneSymbol = gene_name,
             p.value    = P.Value, 
             AvgExpr    = AveExpr) %>%
      select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
    fwrite(top, paste0('./Results/FullRun/voom/', 
           paste('voom', tissue, resp, cov, 'wk00.txt', 
           sep = '.')), sep = '\t')

    # wk1
    top <- topTable(fit, number = Inf, sort.by = 'p',
                    coef = paste0('wk01.', tissue, '.Response'))
    id <- rownames(top)
    top <- top %>%
      mutate(q.value = qvalue(P.Value)$qvalues, gene_id = id) %>%
      inner_join(e2g, by = 'gene_id') %>%
      rename(EnsemblID  = gene_id,
             GeneSymbol = gene_name,
             p.value    = P.Value, 
             AvgExpr    = AveExpr) %>%
      select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
    fwrite(top, paste0('./Results/FullRun/voom/', 
           paste('voom', tissue, resp, cov, 'wk01.txt', 
           sep = '.')), sep = '\t')

    # wk12
    top <- topTable(fit, number = Inf, sort.by = 'p',
                    coef = paste0('wk12.', tissue, '.Response'))
    id <- rownames(top)
    top <- top %>%
      mutate(q.value = qvalue(P.Value)$qvalues, gene_id = id) %>%
      inner_join(e2g, by = 'gene_id') %>%
      rename(EnsemblID  = gene_id,
             GeneSymbol = gene_name,
             p.value    = P.Value, 
             AvgExpr    = AveExpr) %>%
      select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
    fwrite(top, paste0('./Results/FullRun/voom/', 
           paste('voom', tissue, resp, cov, 'wk12.txt', 
           sep = '.')), sep = '\t')

  }

  ### OVER TIME ###

  fit <- lmFit(v, des)
  cm <- makeContrasts('Blood.Delta01' = wk01.Blood.Response - wk00.Blood.Response,
                      'Blood.Delta11' = wk12.Blood.Response - wk01.Blood.Response,
                      'Blood.Delta12' = wk12.Blood.Response - wk00.Blood.Response, 
                      'Lesional.Delta01' = wk01.Lesional.Response - wk00.Lesional.Response,
                      'Lesional.Delta11' = wk12.Lesional.Response - wk01.Lesional.Response,
                      'Lesional.Delta12' = wk12.Lesional.Response - wk00.Lesional.Response,
                      'Nonlesional.Delta01' = wk01.Nonlesional.Response - wk00.Nonlesional.Response,
                      'Nonlesional.Delta11' = wk12.Nonlesional.Response - wk01.Nonlesional.Response,
                      'Nonlesional.Delta12' = wk12.Nonlesional.Response - wk00.Nonlesional.Response, 
                      levels = des)
  fit <- eBayes(contrasts.fit(fit, cm))

  for (tissue in unique(pheno$Tissue)) {

    # Delta01
    top <- topTable(fit, number = Inf, sort.by = 'p',
                    coef = paste0(tissue, '.Delta01'))
    id <- rownames(top)
    top <- top %>%
      mutate(q.value = qvalue(P.Value)$qvalues, gene_id = id) %>%
      inner_join(e2g, by = 'gene_id') %>%
      rename(EnsemblID  = gene_id,
             GeneSymbol = gene_name,
             p.value    = P.Value, 
             AvgExpr    = AveExpr) %>%
      select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
    fwrite(top, paste0('./Results/WithinTissue/voom/', 
           paste('voom', tissue, resp, cov, 'wk00-wk01.txt', 
           sep = '.')), sep = '\t')

    # Delta11
    top <- topTable(fit, number = Inf, sort.by = 'p',
                    coef = paste0(tissue, '.Delta11'))
    id <- rownames(top)
    top <- top %>%
      mutate(q.value = qvalue(P.Value)$qvalues, gene_id = id) %>%
      inner_join(e2g, by = 'gene_id') %>%
      rename(EnsemblID  = gene_id,
             GeneSymbol = gene_name,
             p.value    = P.Value, 
             AvgExpr    = AveExpr) %>%
      select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
    fwrite(top, paste0('./Results/WithinTissue/voom/', 
           paste('voom', tissue, resp, cov, 'wk01-wk12.txt', 
           sep = '.')), sep = '\t')

    # Delta12
    top <- topTable(fit, number = Inf, sort.by = 'p',
                    coef = paste0(tissue, '.Delta12'))
    id <- rownames(top)
    top <- top %>%
      mutate(q.value = qvalue(P.Value)$qvalues, gene_id = id) %>%
      inner_join(e2g, by = 'gene_id') %>%
      rename(EnsemblID  = gene_id,
             GeneSymbol = gene_name,
             p.value    = P.Value, 
             AvgExpr    = AveExpr) %>%
      select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
    fwrite(top, paste0('./Results/WithinTissue/voom/', 
           paste('voom', tissue, resp, cov, 'wk00-wk12.txt', 
           sep = '.')), sep = '\t')

  }
  
}

# Compute in parallel
library(doParallel)
registerDoParallel(12)
foreach(r = c('Continuous', 'Dichotomous')) %:%
  foreach(c = c('None', 'Some', 'All')) %dopar% 
    loop(resp = r, cov = c)


