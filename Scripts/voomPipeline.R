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

  # Run SVA
  mod <- model.matrix(~ 0 + sex + age + bmi + time:Delta_PASI, data=pheno)
  mod0 <- model.matrix(~ 0 + sex + age + bmi, data=pheno)
  svobj <- svaseq(cpm(y), mod, mod0)
  des <- cbind(mod, svobj$sv)
  colnames(des)[5:(7 + svobj$n.sv)] <- c('wk0.DeltaPASI', 'wk1.DeltaPASI', 'wk12.DeltaPASI', 
                                         paste0('SV', 1:svobj$n.sv))

  # Build linear voom model
  v <- voom(y, des)
  corfit <- duplicateCorrelation(v, des, block=pheno$subject)
  v <- voom(y, des, correlation=corfit$consensus, block=pheno$subject)
  corfit <- duplicateCorrelation(v, des, block=pheno$subject)
  fit <- lmFit(v, des, correlation=corfit$consensus, block=pheno$subject)
  fit2 <- eBayes(fit, robust=TRUE)

  # Baseline
  top <- topTable(fit2, coef='wk0.DeltaPASI', number=Inf, sort.by='p')
  id <- rownames(top)
  top <- top %>%
    mutate(q.value = qvalue(P.Value)$qvalues, gene_id = id) %>%
    inner_join(e2g, by='gene_id') %>%
    rename(p.value = P.Value, Gene = gene_name) %>%
    select(Gene, logFC:p.value, q.value, B)
  fwrite(top, paste0(getwd(), '/Results/Baseline,', tissue, '_voom.csv'))

  # One week change
  cm <- makeContrasts('wk0_wk1' = wk1.DeltaPASI - wk0.DeltaPASI, 
                      'wk0_wk12' = wk12.DeltaPASI - wk0.DeltaPASI, levels=des)
  fit2 <- contrasts.fit(fit, cm)
  fit3 <- eBayes(fit2, robust=TRUE)
  top <- topTable(fit3, coef='wk0_wk1', number=Inf, sort.by='p')
  id <- rownames(top)
  top <- top %>%
    mutate(q.value = qvalue(P.Value)$qvalues, gene_id = id) %>%
    inner_join(e2g, by='gene_id') %>%
    rename(p.value = P.Value, Gene = gene_name) %>%
    select(Gene, logFC:p.value, q.value, B)
  fwrite(top, paste0(getwd(), '/Results/wk0_wk1,', tissue, '_voom.csv'))

  # Twelve week change  
  top <- topTable(fit3, coef='wk0_wk12', number=Inf, sort.by='p')
  id <- rownames(top)
  top <- top %>%
    mutate(q.value = qvalue(P.Value)$qvalues, gene_id = id) %>%
    inner_join(e2g, by='gene_id') %>%
    rename(p.value = P.Value, Gene = gene_name) %>%
    select(Gene, logFC:p.value, q.value, B)
  fwrite(top, paste0(getwd(), '/Results/wk0_wk12,', tissue, '_voom.csv'))

}


