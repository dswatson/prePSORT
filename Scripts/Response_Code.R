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

# Winsorise delta PASI distribution
winsorise <- function(x, multiple = 2) {
  y <- x - median(x)
  lim <- mad(y, center = 0) * multiple
  y[y > lim] <- lim
  y[y < -lim] <- -lim
  y <- y + median(x)
  return(y)
}
dt <- pheno %>%
  distinct(Subject, DeltaPASI) %>%
  mutate(Winsorised = winsorise(DeltaPASI)) %>%
  data.table()
pheno[Subject == 'S09', DeltaPASI := dt[Subject == 'S09', Winsorised]]

# Fit model
des <- model.matrix(~ 0 + Time:Tissue + Time:Tissue:DeltaPASI, data = pheno)
colnames(des) <- c(paste(rep(unique(pheno$Time), times = 3), 
                         rep(unique(pheno$Tissue), each = 3), sep = '.'),
                   paste(rep(unique(pheno$Time), times = 3),
                         rep(unique(pheno$Tissue), each = 3),
                         'Response', sep = '.'))
v <- voomWithQualityWeights(y, des)
corfit <- duplicateCorrelation(v, des, block = pheno$Subject)
v <- voomWithQualityWeights(y, des, 
                            correlation = corfit$consensus, block = pheno$Subject)
corfit <- duplicateCorrelation(v, des, block = pheno$Subject)
idx <- rownames(v)
fit <- lmFit(v, des, correlation = corfit$consensus, block = pheno$Subject)
fit <- eBayes(fit, robust = TRUE)

### AT TIME ###
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
      fwrite(paste0('./Results/Response/', 
                    paste(tissue, time, 'txt', sep = '.')), sep = '\t')
  }
}
  
### OVER TIME ###

# One week change
y1 <- y[, !grepl('wk12', pheno$Sample)]
int <- model.matrix(~ 0 + Subject:Tissue:Time, data = pheno)
int1 <- int[, grepl('wk00', colnames(int))]
slope <- model.matrix(~ 0 + Tissue:Time:DeltaPASI, data = pheno)
slope1 <- slope[, grepl('wk01', colnames(slope))]
des <- cbind(int1, slope1)
colnames(des)[28:30] <- paste0(unique(pheno$Tissue), '.Delta01.Response')
v <- voomWithQualityWeights(y1, des)
fit <- eBayes(lmFit(v, des), robust = TRUE)
for (tissue in unique(pheno$Tissue)) {
  topTable(fit, number = Inf, sort.by = 'none',
           coef = paste0(tissue, '.Delta01.Response')) %>%
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
                  paste0(tissue, '.Delta01.txt')), sep = '\t')
}   

# Eleven week change
y11 <- y[, !grepl('wk00', pheno$Sample)]
int11 <- int[, grepl('wk01', colnames(int))]
slope11 <- slope[, grepl('wk12', colnames(slope))]
des <- cbind(int11, slope11)
colnames(des)[28:30] <- paste0(unique(pheno$Tissue), '.Delta11.Response')
v <- voomWithQualityWeights(y11, des)
fit <- eBayes(lmFit(v, des), robust = TRUE)
for (tissue in unique(pheno$Tissue)) {
  topTable(fit, number = Inf, sort.by = 'none',
           coef = paste0(tissue, '.Delta11.Response')) %>%
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
                  paste0(tissue, '.Delta11.txt')), sep = '\t')
} 

# Twelve week change
y12 <- y[, !grepl('wk01', pheno$Sample)]
int12 <- int[, grepl('wk00', colnames(int))]
slope12 <- slope[, grepl('wk12', colnames(slope))]
des <- cbind(int12, slope12)
colnames(des)[28:30] <- paste0(unique(pheno$Tissue), '.Delta12.Response')
v <- voomWithQualityWeights(y12, des)
fit <- eBayes(lmFit(v, des), robust = TRUE)
for (tissue in unique(pheno$Tissue)) {
  topTable(fit, number = Inf, sort.by = 'none',
           coef = paste0(tissue, '.Delta12.Response')) %>%
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
                  paste0(tissue, '.Delta12.txt')), sep = '\t')
} 

