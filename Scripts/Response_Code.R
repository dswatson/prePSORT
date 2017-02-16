# Load libraries, register cores, set seed
library(data.table)
library(tximport)
library(edgeR)
library(limma)
library(qvalue)
library(qusage)
library(dplyr)
library(doParallel)
registerDoParallel(20)
set.seed(123)

# Import data
clin <- fread('./Data/Clinical.csv') %>%
  mutate(Subject.Tissue = paste(Subject, Tissue, sep = '.'))
t2g <- fread('./Data/Ensembl.Hs79.Tx.csv')
e2g <- fread('./Data/Ensembl.Hs79.GeneSymbols.csv')
mods <- fread('./Data/ChaussabelModules.csv')
mod_list <- lapply(unique(mods$Module), function(m) mods[Module == m, gene_name])
names(mod_list) <- unique(mods$Module)
files <- file.path('./Data/RawCounts', clin$Sample, 'abundance.tsv')
txi <- tximport(files, type = 'kallisto', tx2gene = t2g, reader = fread, 
                countsFromAbundance = 'lengthScaledTPM')

# Collapse, filter counts
y <- avereps(txi$counts, ID = e2g[gene_id %in% rownames(txi$counts), gene_name])
keep <- rowSums(cpm(y) > 1) >= 9
y <- DGEList(y[keep, ])
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
  overlap <- sapply(mod_list, function(p) sum(p %in% rownames(v)))
  mod_list <- mod_list[overlap > 0]
  res <- newQSarray(mean = fit$coefficients[, coef],       # Create QSarray obj
                      SD = se,
                sd.alpha = sd.a,
                     dof = fit$df.total,
                  labels = rep('resid', ncol(v)))
  res <- aggregateGeneSet(res, mod_list, n.points = 2^14)  # PDF per gene set
  res <- calcVIF(resid_mat, res, useCAMERA = FALSE)        # VIF on resid_mat
  qsTable(res, number = Inf, sort.by = 'p') %>%
    rename(p.value = p.Value,
            Module = pathway.name,
             logFC = log.fold.change) %>%
    mutate(q.value = qvalue(p.value)$qvalues,
               Idx = row_number()) %>%
    select(Idx, Module:p.value, q.value) %>%
    fwrite(paste0('./Results/Response/RNAseq/',
                  paste0(coef, 'Modules.txt')), sep = '\t')
  
}

### AT TIME ###
des <- model.matrix(~ 0 + Tissue:Time + Tissue:Time:DeltaPASI, data = clin)
colnames(des)[10:18] <- paste(unique(clin$Tissue), rep(unique(clin$Time), each = 3),
                              'Response', sep = '.')
v <- voomWithQualityWeights(y, des)
icc <- duplicateCorrelation(v, des, block = clin$Subject.Tissue)
v <- voomWithQualityWeights(y, des, correlation = icc$cor, 
                            block = clin$Subject.Tissue)
icc <- duplicateCorrelation(v, des, block = clin$Subject.Tissue)  
fit <- lmFit(v, des, correlation = icc$cor, block = clin$Subject.Tissue)
fit <- eBayes(fit)
foreach(j = colnames(des)[10:18]) %dopar% res(j)
  
### OVER TIME ###
des <- model.matrix(~ 0 + Subject:Tissue + Tissue:Time + Tissue:Time:DeltaPASI, 
                    data = clin)
des <- des[, !grepl('wk00', colnames(des))]
colnames(des) <- c(paste(paste0('S', 1:10), 
                         rep(unique(clin$Tissue), each = 10), sep = '.'),
                   paste(unique(clin$Tissue),
                         rep(c('wk01', 'wk12'), each = 3), sep = '.'),
                   paste(unique(clin$Tissue), 
                         rep(c('Delta01', 'Delta12'), each = 3), 
                         'Response', sep = '.')) 
v <- voomWithQualityWeights(y, des)
urFit <- lmFit(v, des)
fit <- eBayes(urFit)
foreach(j = colnames(des)[37:42]) %dopar% res(j)
cm <- makeContrasts('Blood.Delta11.Response' = 
                      Blood.Delta12.Response - Blood.Delta01.Response,
                    'Lesional.Delta11.Response' = 
                      Lesional.Delta12.Response - Lesional.Delta01.Response,
                    'Nonlesional.Delta11.Response' = 
                      Nonlesional.Delta12.Response - Nonlesional.Delta01.Response,
                    levels = des)
fit <- eBayes(contrasts.fit(urFit, cm))
foreach(j = colnames(cm)) %dopar% res(j)


