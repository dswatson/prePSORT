# Load libraries
library(data.table)
library(tximport)
library(edgeR)
library(sva)
library(limma)
library(qvalue)
library(dplyr)

# Prep data
pheno <- fread(paste0(getwd(), '/Data/PrePSORT_Clinical.csv'))
t2g <- fread(paste0(getwd(), '/Data/Ensembl.Hs79.Tx.csv'))
e2g <- fread(paste0(getwd(), '/Data/Ensembl.Hs79.GeneSymbols.csv'))

# Define loop
loop <- function(resp) {

  # TxImport
  dir <- paste0(getwd(), '/Data/Blood')
  files <- file.path(dir, pheno$Sample, 'abundance.tsv')
  txi <- tximport(files, type = 'kallisto', tx2gene = t2g, reader = fread, 
                  countsFromAbundance = 'lengthScaledTPM')
  keep <- rowSums(cpm(txi$counts) > 1) >= 3
  y <- DGEList(txi$counts[keep, ])
  y <- calcNormFactors(y)

  # SVA
  if (resp == 'PASI_50') { 
    mod <- model.matrix(~ 0 + Time + Sex + Age + BMI + PASI_wk00 + Time:PASI_50,
                        data = pheno)
  } else if (resp == 'PASI_75') {
    mod <- model.matrix(~ 0 + Time + Sex + Age + BMI + PASI_wk00 + Time:PASI_75,
                        data = pheno)
  } else if (resp == 'DeltaPASI') {
    mod <- model.matrix(~ 0 + Time + Sex + Age + BMI + PASI_wk00 + Time:DeltaPASI,
                        data = pheno)
  }
  mod0 <- model.matrix(~ 0 + Time + Sex + Age + BMI + PASI_wk00, 
                       data = pheno)
  svobj <- svaseq(cpm(y), mod, mod0)
  des <- cbind(mod, svobj$sv)
  colnames(des)[8:ncol(des)] <- c('wk00.Response', 'wk01.Response', 'wk12.Response',
                                  paste0('SV', 1:svobj$n.sv))

  ### Week Zero ###

  # Voom
  v <- voom(y, des)
  fit <- lmFit(v, des)
  fit <- eBayes(fit, robust = TRUE)

  # Extract
  top <- topTable(fit, coef = 'wk00.Response', number = Inf, sort.by = 'p')
  id <- rownames(top)
  top <- top %>%
    mutate(q.value = qvalue(P.Value)$qvalues, gene_id = id) %>%
    inner_join(e2g, by = 'gene_id') %>%
    rename(EnsemblID  = gene_id,
           GeneSymbol = gene_name,
           p.value    = P.Value, 
           AvgExpr    = AveExpr) %>%
    select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
  fwrite(top, paste0(getwd(), '/Results/ResponseComparisons/voom/Blood/', 
    'voom.Blood.', resp, ',wk00.txt'), sep = '\t')

  ### Week One ###

  # Extract
  top <- topTable(fit, coef = 'wk01.Response', number = Inf, sort.by = 'p')
  id <- rownames(top)
  top <- top %>%
    mutate(q.value = qvalue(P.Value)$qvalues, gene_id = id) %>%
    inner_join(e2g, by = 'gene_id') %>%
    rename(EnsemblID  = gene_id,
           GeneSymbol = gene_name,
           p.value    = P.Value, 
           AvgExpr    = AveExpr) %>%
    select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
  fwrite(top, paste0(getwd(), '/Results/ResponseComparisons/voom/Blood/', 
    'voom.Blood.', resp, ',wk01.txt'), sep = '\t')

  ### Week Twelve ###

  # Extract
  top <- topTable(fit, coef = 'wk12.Response', number = Inf, sort.by = 'p')
  id <- rownames(top)
  top <- top %>%
    mutate(q.value = qvalue(P.Value)$qvalues, gene_id = id) %>%
    inner_join(e2g, by = 'gene_id') %>%
    rename(EnsemblID  = gene_id,
           GeneSymbol = gene_name,
           p.value    = P.Value, 
           AvgExpr    = AveExpr) %>%
    select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
  fwrite(top, paste0(getwd(), '/Results/ResponseComparisons/voom/Blood/', 
    'voom.Blood.', resp, ',wk12.txt'), sep = '\t')

  ### One Week Change ###

  # Build contrasts
  cm <- makeContrasts('Delta01' = wk01.Response - wk00.Response,
                      'Delta11' = wk12.Response - wk01.Response,
                      'Delta12' = wk12.Response - wk00.Response, levels=des)
  fit <- lmFit(v, des)
  fit <- contrasts.fit(fit, cm)
  fit <- eBayes(fit, robust = TRUE)

  # Extract
  top <- topTable(fit, coef = 'Delta01', number = Inf, sort.by = 'p')
  id <- rownames(top)
  top <- top %>%
    mutate(q.value = qvalue(P.Value)$qvalues, gene_id = id) %>%
    inner_join(e2g, by = 'gene_id') %>%
    rename(EnsemblID  = gene_id,
           GeneSymbol = gene_name,
           p.value    = P.Value, 
           AvgExpr    = AveExpr) %>%
    select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
  fwrite(top, paste0(getwd(), '/Results/ResponseComparisons/voom/Blood/', 
    'voom.Blood.', resp, ',wk00-wk01.txt'), sep = '\t')

  ### Eleven Week Change ###

  # Extract
  top <- topTable(fit, coef = 'Delta11', number = Inf, sort.by = 'p')
  id <- rownames(top)
  top <- top %>%
    mutate(q.value = qvalue(P.Value)$qvalues, gene_id = id) %>%
    inner_join(e2g, by = 'gene_id') %>%
    rename(EnsemblID  = gene_id,
           GeneSymbol = gene_name,
           p.value    = P.Value, 
           AvgExpr    = AveExpr) %>%
    select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
  fwrite(top, paste0(getwd(), '/Results/ResponseComparisons/voom/Blood/', 
    'voom.Blood.', resp, ',wk01-wk12.txt'), sep = '\t')

  ### Twelve Week Change ###

  # Extract
  top <- topTable(fit, coef = 'Delta12', number = Inf, sort.by = 'p')
  id <- rownames(top)
  top <- top %>%
    mutate(q.value = qvalue(P.Value)$qvalues, gene_id = id) %>%
    inner_join(e2g, by = 'gene_id') %>%
    rename(EnsemblID  = gene_id,
           GeneSymbol = gene_name,
           p.value    = P.Value, 
           AvgExpr    = AveExpr) %>%
    select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
  fwrite(top, paste0(getwd(), '/Results/ResponseComparisons/voom/Blood/', 
    'voom.Blood.', resp, '.wk00-wk12.txt'), sep = '\t')

}

# Compute in parallel
library(doParallel)
registerDoParallel(cores = 3)
foreach(i = c('PASI_50', 'PASI_75', 'DeltaPASI')) %dopar% loop(i)


