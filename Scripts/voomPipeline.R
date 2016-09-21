# Load libraries
library(data.table)
library(tximport)
library(edgeR)
library(sva)
library(limma)
library(qvalue)
library(dplyr)

# Prep data
pheno <- fread(paste0(getwd(), '/Data/PrePSORT_Clinical.csv')) %>%
  mutate(wk0.DeltaPASI = ifelse(time == 'wk0', Delta_PASI, 0),
         wk1.DeltaPASI = ifelse(time == 'wk1', Delta_PASI, 0),
         wk12.DeltaPASI = ifelse(time == 'wk12', Delta_PASI, 0))
t2g <- fread(paste0(getwd(), '/Data/Ensembl.Hs79.Tx.csv'))
e2g <- fread(paste0(getwd(), '/Data/Ensembl.Hs79.GeneSymbols.csv'))

# Loop
for (tissue in c('Blood', 'LesionalSkin', 'NonlesionalSkin')) {

  # TxImport
  dir <- paste(getwd(), 'Data', tissue, sep='/')
  files <- file.path(dir, pheno$sample, 'MB.abundance.tsv')
  txi <- tximport(files, type='kallisto', tx2gene=t2g, reader=fread, 
                  countsFromAbundance='lengthScaledTPM')
  keep <- rowSums(cpm(txi$counts) > 1) >= 3
  y <- DGEList(txi$counts[keep, ])
  y <- calcNormFactors(y)

  ### Baseline ###

  # SVA
  mod <- model.matrix(~ 0 + time + sex + age + bmi +
                      wk0.DeltaPASI + wk1.DeltaPASI + wk12.DeltaPASI, data=pheno)
  mod0 <- model.matrix(~ 0 + time + sex + age + bmi, data=pheno)
  svobj <- svaseq(cpm(y), mod, mod0)
  des <- cbind(mod, svobj$sv)
  colnames(des)[10:ncol(des)] <- paste0('SV', 1:svobj$n.sv)

  # Voom
  v <- voom(y, des)
  fit <- lmFit(v, des)
  fit <- eBayes(fit, robust=TRUE)

  # Extract
  top <- topTable(fit, coef='wk0.DeltaPASI', number=Inf, sort.by='p')
  id <- rownames(top)
  top <- top %>%
    mutate(q.value = qvalue(P.Value)$qvalues, gene_id = id) %>%
    inner_join(e2g, by='gene_id') %>%
    rename(p.value = P.Value, Gene = gene_name) %>%
    select(Gene, logFC:p.value, q.value, B)
  fwrite(top, paste0(getwd(), '/Results/voom_', tissue, ',wk0.csv'))

  ### One week change ###

  # SVA
  mod <- model.matrix(~ 0 + time + subject + wk0.DeltaPASI + wk1.DeltaPASI, data=pheno)
  mod0 <- model.matrix(~ 0 + time + subject, data=pheno)
  svobj <- svaseq(cpm(y), mod, mod0)
  mod2 <- model.matrix(~ 0 + time + sex + age + bmi + 
                      wk0.DeltaPASI + wk1.DeltaPASI + wk12.DeltaPASI, data=pheno) 
  des <- cbind(mod2, svobj$sv)                                                           
  colnames(des)[10:ncol(des)] <- paste0('SV', 1:svobj$n.sv) 

  # Voom
  v <- voom(y, des)
  corfit <- duplicateCorrelation(v, des, block=pheno$subject)
  v <- voom(y, des, correlation=corfit$consensus, block=pheno$subject)
  corfit <- duplicateCorrelation(v, des, block=pheno$subject)
  fit <- lmFit(v, des, correlation=corfit$consensus, block=pheno$subject)
  cm <- makeContrasts('wk0-wk1' = wk1.DeltaPASI - wk0.DeltaPASI, levels=des)
  fit <- contrasts.fit(fit, cm)
  fit <- eBayes(fit, robust=TRUE)

  # Extract
  top <- topTable(fit, coef='wk0-wk1', number=Inf, sort.by='p')
  id <- rownames(top)
  top <- top %>%
    mutate(q.value = qvalue(P.Value)$qvalues, gene_id = id) %>%
    inner_join(e2g, by='gene_id') %>%
    rename(p.value = P.Value, Gene = gene_name) %>%
    select(Gene, logFC:p.value, q.value, B)
  fwrite(top, paste0(getwd(), '/Results/voom_', tissue, ',wk0-wk1.csv'))

  ### Twelve week change ###

  # SVA
  mod <- model.matrix(~ 0 + time + subject + wk0.DeltaPASI + wk12.DeltaPASI, data=pheno)
  svobj <- svaseq(cpm(y), mod, mod0)
  des <- cbind(mod2, svobj$sv)
  colnames(des)[10:ncol(des)] <- paste0('SV', 1:svobj$n.sv) 

  # Voom
  v <- voom(y, des)
  corfit <- duplicateCorrelation(v, des, block=pheno$subject)
  v <- voom(y, des, correlation=corfit$consensus, block=pheno$subject)
  corfit <- duplicateCorrelation(v, des, block=pheno$subject)
  fit <- lmFit(v, des, correlation=corfit$consensus, block=pheno$subject)
  cm <- makeContrasts('wk0-wk12' = wk12.DeltaPASI - wk0.DeltaPASI, levels=des)
  fit <- contrasts.fit(fit, cm)
  fit <- eBayes(fit, robust=TRUE)

  # Extract
  top <- topTable(fit, coef='wk0-wk12', number=Inf, sort.by='p')
  id <- rownames(top)
  top <- top %>%
    mutate(q.value = qvalue(P.Value)$qvalues, gene_id = id) %>%
    inner_join(e2g, by='gene_id') %>%
    rename(p.value = P.Value, Gene = gene_name) %>%
    select(Gene, logFC:p.value, q.value, B)
  fwrite(top, paste0(getwd(), '/Results/voom_', tissue, ',wk0-wk12.csv'))

}


