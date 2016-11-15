Differential Expression Over Time
================

-   [Load Data](#load-data)
-   [Fit model](#fit-model)

All analysis was conducted in R version 3.3.2 using the following script. Computations were performed on a MacBook Pro with 16GB of RAM.

``` r
# Load libraries
library(data.table)
library(tximport)
library(edgeR)
library(limma)
library(qvalue)
library(dplyr)
```

Load Data
=========

We begin by importing clinical and gene expression data. Genes are filtered and normalised using the guidelines laid out in Response\_Script.md.

``` r
# Load data
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
```

We create a custom function to export results.

``` r
# Define results function
res <- function(contrast) {
  topTable(fit, number = Inf, sort.by = 'none',
           coef = contrast) %>%
    mutate(gene_id = idx) %>%
    inner_join(e2g, by = 'gene_id') %>%
    rename(EnsemblID  = gene_id,
           GeneSymbol = gene_name,
           p.value    = P.Value,
           q.value    = adj.P.Val,
           AvgExpr    = AveExpr) %>%
    arrange(p.value) %>%
    select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value) %>%
    fwrite(paste0('./Results/Response/', 
                  paste0(contrast, '.txt')), sep = '\t')
}
```

Fit model
=========

We fit a model using the following formula.

``` r
# Voom
des <- model.matrix(~ 0 + Subject:Tissue + Tissue:Time, data = pheno)
colnames(des) <- c(paste(paste0('S', 1:10), 
                         rep(unique(pheno$Tissue), each = 10), sep = '.'),
                   paste(rep(unique(pheno$Tissue), times = 2), 
                         rep(c('Delta01', 'Delta12'), each = 3), sep = '.'))
```

The first 30 coefficients pick up baseline expression for each subject-tissue interaction. The final 6 coefficients represent the changes from baseline to week 1 and 12 in each tissue type. We use sample weights to account for variable library quality in our model.

``` r
v <- voomWithQualityWeights(y, des)
urFit <- lmFit(v, des)
fit <- eBayes(lmFit(v, des), robust = TRUE)
for (i in colnames(des)[31:36]) res(i)
```

We now create a contrast matrix to extract 11 week changes by contrasting 12 and 1 week changes.

``` r
cm <- makeContrasts('Blood.Delta11' = 
                      Blood.Delta12 - Blood.Delta01,
                    'Lesional.Delta11' = 
                      Lesional.Delta12 - Lesional.Delta01,
                    'Nonlesional.Delta11' = 
                      Nonlesional.Delta12 - Nonlesional.Delta01,
                    levels = des)
fit <- eBayes(contrasts.fit(urFit, cm), robust = TRUE)
for (i in colnames(cm)) res(i)
```
