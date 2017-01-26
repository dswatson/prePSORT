# Load libraries
library(data.table)
library(tximport)
library(edgeR)
library(limma)
library(qvalue)
library(dplyr)

# Prep data
pheno <- fread('./Data/Clinical.csv') %>%
  mutate(Subject.Tissue = paste(Subject, Tissue, sep = '.'))
t2g <- fread('./Data/Ensembl.Hs79.Tx.csv')
e2g <- fread('./Data/Ensembl.Hs79.GeneSymbols.csv')

# TxImport
files <- file.path('./Data/RawCounts', pheno$Sample, 'abundance.tsv')
txi <- tximport(files, type = 'kallisto', tx2gene = t2g, reader = fread, 
                countsFromAbundance = 'lengthScaledTPM')
keep <- rowSums(cpm(txi$counts) > 1) >= 9
y <- DGEList(txi$counts[keep, ])
y <- calcNormFactors(y)
e2g <- e2g %>% filter(gene_id %in% rownames(y)) 

# Winsorise delta PASI distribution
winsorise <- function(x, multiple = 2) {
  y <- x - median(x)
  lim <- mad(y, center = 0) * multiple
  y[y > lim] <- lim
  y[y < -lim] <- -lim
  y <- y + median(x)
  return(y)
}
pheno <- pheno %>% mutate(DeltaPASI = winsorise(DeltaPASI))

# Define results function
res <- function(coef) {
  topTable(fit, coef = coef, number = Inf, sort.by = 'none') %>%
    rename(AvgExpr = AveExpr,
           p.value = P.Value) %>%
    mutate(q.value = qvalue(p.value)$qvalues,
              Gene = rownames(v)) %>%
    arrange(p.value) %>%
    mutate(Idx = row_number()) %>%
    select(Idx, Gene, AvgExpr, logFC, p.value, q.value) %>%
    fwrite(paste0('./Results/Response/RNAseq/', 
                  paste0(coef, '.txt')), sep = '\t')
}

### AT TIME ###
des <- model.matrix(~ 0 + Tissue:Time + Tissue:Time:DeltaPASI, data = pheno)
colnames(des)[10:18] <- paste(unique(pheno$Tissue), rep(unique(pheno$Time), each = 3),
                              'Response', sep = '.')
v <- voomWithQualityWeights(y, des)
icc <- duplicateCorrelation(v, des, block = pheno$Subject.Tissue)
v <- voomWithQualityWeights(y, des, correlation = icc$cor, 
                            block = pheno$Subject.Tissue)
icc <- duplicateCorrelation(v, des, block = pheno$Subject.Tissue)  
v <- avereps(v, e2g$gene_name)
fit <- lmFit(v, des, correlation = icc$cor, block = pheno$Subject.Tissue)
fit <- eBayes(fit, robust = TRUE)
for (i in colnames(des)[10:18]) res(i)
  
### OVER TIME ###
des <- model.matrix(~ 0 + Subject:Tissue + Tissue:Time + Tissue:Time:DeltaPASI, 
                    data = pheno)
des <- des[, !grepl('wk00', colnames(des))]
colnames(des) <- c(paste(paste0('S', 1:10), 
                         rep(unique(pheno$Tissue), each = 10), sep = '.'),
                   paste(unique(pheno$Tissue),
                         rep(c('wk01', 'wk12'), each = 3), sep = '.'),
                   paste(unique(pheno$Tissue), 
                         rep(c('Delta01', 'Delta12'), each = 3), 
                         'Response', sep = '.'))
v <- voomWithQualityWeights(y, des)
v <- avereps(v, e2g$gene_name)
urFit <- lmFit(v, des)
fit <- eBayes(urFit, robust = TRUE)
for (i in colnames(des)[37:42]) res(i)
cm <- makeContrasts('Blood.Delta11.Response' = 
                      Blood.Delta12.Response - Blood.Delta01.Response,
                    'Lesional.Delta11.Response' = 
                      Lesional.Delta12.Response - Lesional.Delta01.Response,
                    'Nonlesional.Delta11.Response' = 
                      Nonlesional.Delta12.Response - Nonlesional.Delta01.Response,
                    levels = des)
fit <- eBayes(contrasts.fit(urFit, cm), robust = TRUE)
for (i in colnames(cm)) res(i)


