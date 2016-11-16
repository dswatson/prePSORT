# Load libraries
library(data.table)
library(limma)
library(tidyverse)

# Prep data
pheno <- fread('./Data/Clinical.csv') %>%
  filter(Tissue == 'Blood') %>%
  mutate(Sample = paste0('S', 1:30))
mat <- read.csv('./Data/SomaScan.csv', row.names = 1) %>% 
  as.matrix()
mat <- log2(mat)
idx <- rownames(mat)

# Winsorise delta PASI distribution
winsorise <- function(x, multiple = 2) {
  y <- x - median(x)
  lim <- mad(y, center = 0) * multiple
  y[y > lim] <- lim
  y[y < -lim] <- -lim
  y <- y + median(x)
  return(y)
}
df <- pheno %>%
  distinct(Subject, DeltaPASI) %>%
  mutate(Winsorised = winsorise(DeltaPASI)) 
pheno$DeltaPASI[pheno$Subject == 'S09'] <- df %>%
  filter(Subject == 'S09') %>%
  select(Winsorised) %>%
  as.numeric()

# Define results function
res <- function(contrast) {
  topTable(fit, number = Inf, sort.by = 'none',
           coef = contrast) %>%
    mutate(Analyte = idx) %>%
    rename(p.value    = P.Value,
           q.value    = adj.P.Val,
           AvgExpr    = AveExpr) %>%
    arrange(p.value) %>%
    select(Analyte, AvgExpr, logFC, p.value, q.value) %>%
    fwrite(paste0('./Results/Response/Proteomics/', 
                  paste0(contrast, '.txt')), sep = '\t')
}

### AT TIME ###
des <- model.matrix(~ 0 + Time + Time:DeltaPASI, data = pheno)
colnames(des)[4:6] <- c(paste(unique(pheno$Time), 'Response', sep = '.'))
corfit <- duplicateCorrelation(mat, des, block = pheno$Subject)
fit <- lmFit(mat, des, correlation = corfit$consensus)
fit <- eBayes(fit, robust = TRUE)
for (i in colnames(des)[4:6]) res(i)

### OVER TIME ###
des <- model.matrix(~ 0 + Subject + Time + Time:DeltaPASI, data = pheno)
des <- des[, !grepl('wk00', colnames(des))]
colnames(des)[13:14] <- paste(c('Delta01', 'Delta12'), 'Response', sep = '.')
urFit <- lmFit(mat, des)
fit <- eBayes(urFit, robust = TRUE)
for (i in colnames(des)[13:14]) res(i)
cm <- makeContrasts('Delta11.Response' = Delta12.Response - Delta01.Response,
                    levels = des)
fit <- eBayes(contrasts.fit(urFit, cm), robust = TRUE)
res('Delta11.Response')

