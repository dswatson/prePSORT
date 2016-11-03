# Load libraries
library(data.table)
library(tximport)
library(edgeR)
library(limma)
library(qvalue)
library(dplyr)

# Prep data
pheno <- fread('./Data/Clinical.csv')
t2g <- fread('./Data/Ensembl.Hs79.Tx.csv')
e2g <- fread('./Data/Ensembl.Hs79.GeneSymbols.csv')

# TxImport
files <- file.path('./Data/RawCounts', pheno$Sample, 'abundance.tsv')
txi <- tximport(files, type = 'kallisto', tx2gene = t2g, reader = fread, 
                countsFromAbundance = 'lengthScaledTPM')
keep <- rowSums(cpm(txi$counts) > 1) >= 9
y <- DGEList(txi$counts[keep, ])
y <- calcNormFactors(y)

# Fit model
des <- model.matrix(~ 0 + Time:Tissue + Time:Tissue:DeltaPASI, data = pheno)
colnames(des)[10:18] <- c(paste(rep(unique(pheno$Time), each = 3), 
                                unique(pheno$Tissue), 
                                'Response', sep = '.'))
v <- voomWithQualityWeights(y, des, method = 'reml', maxiter = 1000)
corfit <- duplicateCorrelation(v, des, block = pheno$Subject)
v <- voomWithQualityWeights(y, des, method = 'reml', maxiter = 1000,
                            correlation = corfit$consensus, block = pheno$Subject)
corfit <- duplicateCorrelation(v, des, block = pheno$Subject)
idx <- rownames(v)
urFit <- lmFit(v, des, correlation = corfit$consensus, block = pheno$Subject)

# At time
fit <- eBayes(urFit)
for (tissue in unique(pheno$Tissue))  {
  for (time in unique(pheno$Time)) {
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
      fwrite(paste0('./Results/Response/MixedModel/', 
                    paste(tissue, time, 'txt', sep = '.')), sep = '\t')
  }
}
  
# Over time
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
fit <- eBayes(contrasts.fit(urFit, cm))
for (tissue in unique(pheno$Tissue))  {
  for (delta in c('Delta01', 'Delta11', 'Delta12')) {
    topTable(fit, number = Inf, sort.by = 'none',
             coef = paste(tissue, delta, sep = '.')) %>%
      mutate(q.value = qvalue(P.Value)$qvalues, 
             gene_id = idx) %>%
      inner_join(e2g, by = 'gene_id') %>%
      rename(EnsemblID  = gene_id,
             GeneSymbol = gene_name,
             p.value    = P.Value, 
             AvgExpr    = AveExpr) %>%
      arrange(p.value) %>%
      select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value) %>%
      fwrite(paste0('./Results/Response/MixedModel/', 
                    paste(tissue, delta, 'txt', sep = '.')), sep = '\t')
  }
}


