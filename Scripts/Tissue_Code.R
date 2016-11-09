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
des <- model.matrix(~ Subject + Tissue:Time, data = pheno)
des <- des[, c(1:10, 12:13)]
colnames(des)[11:12] <- c('wk00.Lesional', 'wk00.Nonlesional')
v <- voomWithQualityWeights(y, des)
idx <- rownames(v)
fit <- lmFit(v, des)
cm <- makeContrasts('Tissue' = wk00.Lesional - wk00.Nonlesional, levels = des)
fit <- eBayes(contrasts.fit(fit, cm))

# Export results
topTable(fit, coef = 'Tissue', number = Inf, sort.by = 'none') %>%
  mutate(q.value = qvalue(P.Value)$qvalues, 
         gene_id = idx) %>%
  inner_join(e2g, by = 'gene_id') %>%
  rename(EnsemblID  = gene_id,
         GeneSymbol = gene_name,
         p.value    = P.Value, 
         AvgExpr    = AveExpr) %>%
  arrange(p.value) %>%
  select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value) %>%
  fwrite('./Results/Tissue/LvsN_wk00.txt', sep = '\t')


