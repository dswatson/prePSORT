# Load libraries
library(data.table)
library(tximport)
library(edgeR)
library(sva)
library(limma)
library(qvalue)
library(dplyr)

# Prep data
pheno <- fread(paste0(getwd(), '/Data/PrePSORT_Clinical,time.csv')) %>%
  mutate(wk00 = ifelse(Time == 'wk00', 1, 0),
         wk01 = ifelse(Time == 'wk01', 1, 0),
         wk12 = ifelse(Time == 'wk12', 1, 0),
         wk00.Response = ifelse(Time == 'wk00', Delta_PASI, 0),
         wk01.Response = ifelse(Time == 'wk01', Delta_PASI, 0),
         wk12.Response = ifelse(Time == 'wk12', Delta_PASI, 0))
t2g <- fread(paste0(getwd(), '/Data/Ensembl.Hs79.Tx.csv'))
e2g <- fread(paste0(getwd(), '/Data/Ensembl.Hs79.GeneSymbols.csv'))

# Loop
for (tissue in c('Blood', 'LesionalSkin', 'NonlesionalSkin')) {

  # TxImport
  dir <- paste0(getwd(), '/Data/', tissue)
  files <- file.path(dir, pheno$Sample, 'abundance.tsv')
  txi <- tximport(files, type='kallisto', tx2gene=t2g, reader=fread, 
                  countsFromAbundance='lengthScaledTPM')
  keep <- rowSums(cpm(txi$counts) > 1) >= 3
  y <- DGEList(txi$counts[keep, ])
  y <- calcNormFactors(y)

  ### Week Zero ###

  # SVA
  mod <- model.matrix(~ 0 + Time + Sex + Age + BMI + PASI_A + 
                      wk00.Response + wk01.Response + wk12.Response, data=pheno)
  mod0 <- model.matrix(~ Sex + Age + BMI + PASI_A, data=pheno)
  svobj <- svaseq(cpm(y), mod, mod0)
  des <- cbind(mod, svobj$sv)
  colnames(des)[11:ncol(des)] <- paste0('SV', 1:svobj$n.sv)

  # Voom
  v <- voom(y, des)
  fit <- lmFit(v, des)
  fit <- eBayes(fit, robust=TRUE)

  # Extract
  top <- topTable(fit, coef='wk00.Response', number=Inf, sort.by='p')
  id <- rownames(top)
  top <- top %>%
    mutate(q.value = qvalue(P.Value)$qvalues, gene_id = id) %>%
    inner_join(e2g, by='gene_id') %>%
    rename(EnsemblID  = gene_id,
           GeneSymbol = gene_name,
           p.value    = P.Value, 
           AvgExpr    = AveExpr) %>%
    select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
  fwrite(top, paste0(getwd(), '/Results/voom/', tissue, ',wk00.txt'), sep='\t')

  ### Week One ###

  # Extract
  top <- topTable(fit, coef='wk01.Response', number=Inf, sort.by='p')
  id <- rownames(top)
  top <- top %>%
    mutate(q.value = qvalue(P.Value)$qvalues, gene_id = id) %>%
    inner_join(e2g, by='gene_id') %>%
    rename(EnsemblID  = gene_id,
           GeneSymbol = gene_name,
           p.value    = P.Value, 
           AvgExpr    = AveExpr) %>%
    select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
  fwrite(top, paste0(getwd(), '/Results/voom/', tissue, ',wk01.txt'), sep='\t')
  
  ### Week Twelve ###

  # Extract
  top <- topTable(fit, coef='wk12.Response', number=Inf, sort.by='p')
  id <- rownames(top)
  top <- top %>%
    mutate(q.value = qvalue(P.Value)$qvalues, gene_id = id) %>%
    inner_join(e2g, by='gene_id') %>%
    rename(EnsemblID  = gene_id,
           GeneSymbol = gene_name,
           p.value    = P.Value, 
           AvgExpr    = AveExpr) %>%
    select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
  fwrite(top, paste0(getwd(), '/Results/voom/', tissue, ',wk12.txt'), sep='\t')

  ### One Week Change ###

  # SVA
  mod <- model.matrix(~ 0 + Subject + wk01 + wk12 + wk01.Response + wk12.Response, data=pheno)
  mod0 <- model.matrix(~ Subject, data=pheno)
  svobj <- svaseq(cpm(y), mod, mod0)
  des <- cbind(mod, svobj$sv)                                                           
  colnames(des)[15:ncol(des)] <- paste0('SV', 1:svobj$n.sv) 

  # Voom
  v <- voom(y, des)
  fit <- lmFit(v, des)
  fit <- eBayes(fit, robust=TRUE)

  # Extract
  top <- topTable(fit, coef='wk01.Response', number=Inf, sort.by='p')
  id <- rownames(top)
  top <- top %>%
    mutate(q.value = qvalue(P.Value)$qvalues, gene_id = id) %>%
    inner_join(e2g, by='gene_id') %>%
    rename(EnsemblID  = gene_id,
           GeneSymbol = gene_name,
           p.value    = P.Value, 
           AvgExpr    = AveExpr) %>%
    select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
  fwrite(top, paste0(getwd(), '/Results/voom/', tissue, ',wk00-wk01.txt'), sep='\t')

  ### Twelve Week Change ###

  # Extract
  top <- topTable(fit, coef='wk12.Response', number=Inf, sort.by='p')
  id <- rownames(top)
  top <- top %>%
    mutate(q.value = qvalue(P.Value)$qvalues, gene_id = id) %>%
    inner_join(e2g, by='gene_id') %>%
    rename(EnsemblID  = gene_id,
           GeneSymbol = gene_name,
           p.value    = P.Value, 
           AvgExpr    = AveExpr) %>%
    select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
  fwrite(top, paste0(getwd(), '/Results/voom/', tissue, ',wk00-wk12.txt'), sep='\t')

  ### Eleven Week Change ###

  # SVA
  mod <- model.matrix(~ 0 + Subject + wk00 + wk12 + wk00.Response + wk12.Response, data=pheno)
  svobj <- svaseq(cpm(y), mod, mod0)
  des <- cbind(mod, svobj$sv)                                                           
  colnames(des)[15:ncol(des)] <- paste0('SV', 1:svobj$n.sv) 

  # Voom
  v <- voom(y, des)
  fit <- lmFit(v, des)
  fit <- eBayes(fit, robust=TRUE)

  # Extract
  top <- topTable(fit, coef='wk12.Response', number=Inf, sort.by='p')
  id <- rownames(top)
  top <- top %>%
    mutate(q.value = qvalue(P.Value)$qvalues, gene_id = id) %>%
    inner_join(e2g, by='gene_id') %>%
    rename(EnsemblID  = gene_id,
           GeneSymbol = gene_name,
           p.value    = P.Value, 
           AvgExpr    = AveExpr) %>%
    select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
  fwrite(top, paste0(getwd(), '/Results/voom/', tissue, ',wk01-wk12.txt'), sep='\t')

}


