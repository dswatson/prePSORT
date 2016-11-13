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

# Define results function
res <- function(contrast) {
  topTable(fit, number = Inf, sort.by = 'none',
           coef = contrast) %>%
    mutate(q.value = qvalue(P.Value)$qvalues, 
           gene_id = idx) %>%
    inner_join(e2g, by = 'gene_id') %>%
    rename(EnsemblID  = gene_id,
           GeneSymbol = gene_name,
           p.value    = P.Value, 
           AvgExpr    = AveExpr) %>%
    arrange(p.value) %>%
    select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value) %>%
    fwrite(paste0('./Results/Response/', 
                  paste0(contrast, '.txt')), sep = '\t')
}

# Voom
des <- model.matrix(~ 0 + Subject:Tissue + Tissue:Time, data = pheno)
colnames(des) <- c(paste(paste0('S', 1:10), 
                         rep(unique(pheno$Tissue), each = 10), sep = '.'),
                   paste(rep(unique(pheno$Tissue), times = 2), 
                         rep(c('Delta01', 'Delta12'), each = 3), sep = '.'))
v <- voomWithQualityWeights(y, des)
urFit <- lmFit(v, des)
fit <- eBayes(urFit, robust = TRUE)
for (i in colnames(des)[31:36]) res(i)
cm <- makeContrasts('Blood.Delta11' = Blood.Delta12 - Blood.Delta01,
                    'Lesional.Delta11' = Lesional.Delta12 - Lesional.Delta01,
                    'Nonlesional.Delta11' = Nonlesional.Delta12 - Nonlesional.Delta01,
                    levels = des)
fit <- eBayes(contrasts.fit(urFit, cm), robust = TRUE)
for (i in colnames(cm)) res(i)

