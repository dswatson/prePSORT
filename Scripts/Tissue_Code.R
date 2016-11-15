# Load libraries
library(data.table)
library(tximport)
library(edgeR)
library(limma)
library(qvalue)
library(dplyr)

# Prep data
pheno <- fread('./Data/Clinical.csv') %>%
  mutate(Tissue = relevel(as.factor(Tissue), ref = 'Nonlesional'))
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
des <- model.matrix(~ 0 + Subject:Time + Tissue:Time, data = pheno)
v <- voomWithQualityWeights(y, des)
idx <- rownames(v)
fit <- eBayes(lmFit(v, des), robust = TRUE)

# Export results
topTable(fit, coef = 34, number = Inf, sort.by = 'none') %>%
  mutate(gene_id = idx) %>%
  inner_join(e2g, by = 'gene_id') %>%
  rename(EnsemblID  = gene_id,
         GeneSymbol = gene_name,
         p.value    = P.Value,
         q.value    = adj.P.Value,
         AvgExpr    = AveExpr) %>%
  arrange(p.value) %>%
  select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value) %>%
  fwrite('./Results/Tissue/LvsN_wk00.txt', sep = '\t')


