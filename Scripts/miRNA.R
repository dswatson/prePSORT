# Load libraries
library(data.table)
library(edgeR)
library(limma)
library(qvalue)
library(dplyr)

# Prep data
clin <- fread('./Data/Clinical.csv') %>%
  filter(Tissue == 'Blood') %>%
  mutate(Sample = gsub('.Blood.', '_', Sample))
mat <- read.csv('./Data/miRNA.csv', row.names = 1) %>% 
  as.matrix()

# Collapse, filter counts
keep <- rowSums(cpm(mat) > 1) >= 3
y <- DGEList(mat[keep, ])
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
clin <- clin %>%
  group_by(Tissue, Time) %>%
  mutate(DeltaPASI = winsorise(DeltaPASI))

# Define results function
res <- function(contrast) {
  topTable(fit, number = Inf, sort.by = 'none', coef = contrast) %>%
    mutate(Gene = rownames(y),
        q.value = qvalue(P.Value)$qvalues) %>%
    rename(p.value = P.Value,
           AvgExpr = AveExpr) %>%
    arrange(p.value) %>%
    select(Gene, AvgExpr, logFC, p.value, q.value) %>%
    fwrite(paste0('./Results/Response/miRNA/', 
                  paste0(contrast, '.txt')), sep = '\t')
}

# Fit model
des <- model.matrix(~ 0 + Time + Time:DeltaPASI, data = clin)
colnames(des)[4:6] <- c(paste(unique(clin$Time), 'Response', sep = '.'))
v <- voomWithQualityWeights(y, des)
icc <- duplicateCorrelation(v, des, block = clin$Subject)
v <- voomWithQualityWeights(y, des, correlation = icc$cor,
                            block = clin$Subject)
icc <- duplicateCorrelation(v, des, block = clin$Subject)
fit <- lmFit(v, des, correlation = icc$cor, block = clin$Subject)
fit <- eBayes(fit)
for (j in colnames(des)[4:6]) res(j)
saveRDS(v$E, './Data/mat_miRNA.rds')
saveRDS(fit, './Data/fit_miRNA.rds')
