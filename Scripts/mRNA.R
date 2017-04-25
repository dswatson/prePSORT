# Load libraries, register cores, set seed
library(data.table)
library(tximport)
library(edgeR)
library(limma)
library(qvalue)
library(qusage)
library(dplyr)
library(doMC)
registerDoMC(10)
set.seed(123)

# Import data
clin <- fread('./Data/Clinical.csv') %>%
  mutate(Subject.Tissue = paste(Subject, Tissue, sep = '.'))
t2g <- fread('./Data/Hs79.t2g.csv')
mods <- readRDS('./Data/LiModules.rds')
files <- file.path('./Data/RawCounts', clin$Sample, 'abundance.tsv')
txi <- tximport(files, type = 'kallisto', tx2gene = t2g, reader = fread, 
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
clin <- clin %>% mutate(DeltaPASI = winsorise(DeltaPASI))

# Results function
res <- function(coef) {
  
  # Genes
  topTable(fit, coef = coef, number = Inf, sort.by = 'none') %>%
    rename(AvgExpr = AveExpr,
           p.value = P.Value) %>%
    mutate(q.value = qvalue(p.value)$qvalues,
           Gene = rownames(v)) %>%
    arrange(p.value) %>%
    mutate(Idx = row_number()) %>%
    select(Idx, Gene, AvgExpr, logFC, p.value, q.value) %>%
    fwrite(paste0('./Results/Response/RNAseq/', 
                  paste0(coef, '.Genes.txt')), sep = '\t')
  
  # Modules
  se <- sqrt(fit$s2.post) * fit$stdev.unscaled[, coef]
  sd.a <- se / (fit$sigma * fit$stdev.unscaled[, coef])
  sd.a[is.infinite(sd.a)] <- 1
  resid_mat <- residuals(fit, v)
  overlap <- sapply(mods$m2g, function(g) sum(g %in% rownames(v)))
  mods$m2g <- mods$m2g[overlap > 1L]
  res <- newQSarray(mean = fit$coefficients[, coef], SD = se,
                    sd.alpha = sd.a, dof = fit$df.total, 
                    labels = rep('resid', ncol(v)))        # Create QSarray object
  res <- aggregateGeneSet(res, mods$m2g, n.points = 2^14)  # PDF per gene set
  res <- calcVIF(resid_mat, res, useCAMERA = FALSE)        # VIF on resid_mat
  qsTable(res, number = Inf, sort.by = 'p') %>%
    rename(p.value = p.Value,
            Module = pathway.name,
             logFC = log.fold.change) %>%
    inner_join(mods$anno, by = 'Module') %>%
    mutate(q.value = qvalue(p.value)$qvalues,
               Idx = row_number()) %>%
    select(Idx, Module, Title, Category, logFC, p.value, q.value) %>%
    fwrite(paste0('./Results/Response/RNAseq/',
                  paste0(coef, '.Modules.txt')), sep = '\t')
  
}

# Fit model
des <- model.matrix(~ 0 + Tissue:Time + Tissue:Time:DeltaPASI, data = clin)
colnames(des)[10:18] <- paste(unique(clin$Tissue), rep(unique(clin$Time), each = 3),
                              'Response', sep = '_')
v <- voomWithQualityWeights(y, des)
icc <- duplicateCorrelation(v, des, block = clin$Subject.Tissue)
v <- voomWithQualityWeights(y, des, correlation = icc$cor, 
                            block = clin$Subject.Tissue)
icc <- duplicateCorrelation(v, des, block = clin$Subject.Tissue)  
fit <- lmFit(v, des, correlation = icc$cor, block = clin$Subject.Tissue)
fit <- eBayes(fit)
foreach(j = colnames(des)[10:18]) %dopar% res(j)
saveRDS(v$E, './Data/mat_mRNA.rds')
saveRDS(fit, './Data/fit_mRNA.rds')

