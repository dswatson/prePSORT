# Load libraries
library(data.table)
library(limma)
library(qvalue)
library(dplyr)

# Prep data
clin <- fread('./Data/Clinical.csv') %>%
  filter(Tissue == 'Blood') %>%
  mutate(Sample = gsub('.Blood.', '_', Sample))
mat <- read.csv('./Data/SomaScan.csv', row.names = 1) %>% 
  as.matrix()
mat <- log2(mat)
mat <- mat[, match(clin$Sample, colnames(mat))]

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
    mutate(Analyte = rownames(mat),
           q.value = qvalue(P.Value)$qvalues) %>%
    rename(p.value = P.Value,
           AvgExpr = AveExpr) %>%
    arrange(p.value) %>%
    select(Analyte, AvgExpr, logFC, p.value, q.value) %>%
    fwrite(paste0('./Results/Response/Proteomics/', 
                  paste0(contrast, '.txt')), sep = '\t')
}

# Fit model
des <- model.matrix(~ 0 + Time + Time:DeltaPASI, data = clin)
colnames(des)[4:6] <- c(paste(unique(clin$Time), 'Response', sep = '.'))
icc <- duplicateCorrelation(mat, des, block = clin$Subject)
fit <- lmFit(mat, des, correlation = icc$cor, block = clin$Subject)
fit <- eBayes(fit)
for (j in colnames(des)[4:6]) res(j)
saveRDS(mat, './Data/mat_prot.rds')
saveRDS(fit, './Data/fit_prot.rds')

