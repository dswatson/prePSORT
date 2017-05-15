### mRNA ###

# Load libraries
library(data.table)
library(tximport)
library(edgeR)
library(limma)
library(qvalue)
library(dplyr)

# Import data
clin <- fread('./Data/Clinical.csv') %>%
  mutate(Subject.Tissue = paste(Subject, Tissue, sep = '.'))
t2g <- fread('./Data/Hs79.t2g.csv')
files <- file.path('./Data/RawCounts', clin$Sample, 'abundance.tsv')
txi <- tximport(files, type = 'kallisto', tx2gene = t2g, importer = fread, 
                countsFromAbundance = 'lengthScaledTPM')

# Collapse, filter counts
keep <- rowSums(cpm(txi$counts) > 1) >= 9
y <- DGEList(txi$counts[keep, ])
y <- calcNormFactors(y)
colnames(y) <- clin$Sample

# Winsorise delta PASI distribution
winsorise <- function(x, multiple = 2) {
  y <- x - median(x)
  lim <- mad(y, center = 0) * multiple
  y[y > lim] <- lim
  y[y < -lim] <- -lim
  y <- y + median(x)
  return(y)
}
clin <- clin %>%
  group_by(Tissue, Time) %>%
  mutate(DeltaPASI = winsorise(DeltaPASI)) %>%
  ungroup()

# Results function
res <- function(coef) {
  
  # Genes
  topTable(fit, coef = coef, number = Inf, sort.by = 'none') %>%
    rename(AvgExpr = AveExpr,
           p.value = P.Value) %>%
    mutate(q.value = qvalue(p.value)$qvalues,
           Gene = rownames(v)) %>%
    arrange(p.value) %>%
    select(Gene, AvgExpr, logFC, p.value, q.value) %>%
    fwrite(paste0('./Results/Response/mRNA/', coef, '.txt'), sep = '\t')
  
}

# Fit model
des <- model.matrix(~ 0 + Tissue:Time + Tissue:Time:DeltaPASI, data = clin)
colnames(des)[10:18] <- paste(unique(clin$Tissue), 
                              rep(unique(clin$Time), each = 3), sep = '_')
v <- voomWithQualityWeights(y, des)
icc <- duplicateCorrelation(v, des, block = clin$Subject.Tissue)
v <- voomWithQualityWeights(y, des, correlation = icc$cor, 
                            block = clin$Subject.Tissue)
icc <- duplicateCorrelation(v, des, block = clin$Subject.Tissue)  
fit <- lmFit(v, des, correlation = icc$cor, block = clin$Subject.Tissue)
fit <- eBayes(fit, robust = TRUE)

# Export
for (j in colnames(des)[10:18]) res(j)
saveRDS(v$E, './Data/mat_mRNA.rds')
saveRDS(fit, './Data/fit_mRNA.rds')


