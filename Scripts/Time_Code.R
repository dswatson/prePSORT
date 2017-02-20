# Load libraries, register cores, set seed
library(data.table)
library(tximport)
library(edgeR)
library(limma)
library(qvalue)
library(qusage)
library(dplyr)
library(doParallel)
registerDoParallel(6)
set.seed(123)

# Import data
clin <- fread('./Data/Clinical.csv') 
t2g <- fread('./Data/Ensembl.Hs79.Tx.csv')
e2g <- fread('./Data/Ensembl.Hs79.GeneSymbols.csv')
mods <- readRDS('./Data/LiModules.rds')
files <- file.path('./Data/RawCounts', clin$Sample, 'abundance.tsv')
txi <- tximport(files, type = 'kallisto', tx2gene = t2g, reader = fread, 
                countsFromAbundance = 'lengthScaledTPM')

# Collapse, filter counts
y <- avereps(txi$counts, ID = e2g[gene_id %in% rownames(txi$counts), gene_name])
keep <- rowSums(cpm(y) > 1) >= 9
y <- DGEList(y[keep, ])
y <- calcNormFactors(y)

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
    fwrite(paste0('./Results/Time/', paste0(coef, '.Genes.txt')), sep = '\t')
  
  # Modules
  se <- sqrt(fit$s2.post) * fit$stdev.unscaled[, coef]
  sd.a <- se / (fit$sigma * fit$stdev.unscaled[, coef])
  sd.a[is.infinite(sd.a)] <- 1
  resid_mat <- residuals(fit, v)
  overlap <- sapply(mods$m2g, function(p) sum(p %in% rownames(v)))
  mods$m2g <- mods$m2g[overlap > 0]
  res <- newQSarray(mean = fit$coefficients[, coef],       # Create QSarray obj
                    SD = se,
                    sd.alpha = sd.a,
                    dof = fit$df.total,
                    labels = rep('resid', ncol(v)))
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
    fwrite(paste0('./Results/Time/', paste0(coef, '.Modules.txt')), sep = '\t')
  
}

# Fit model
des <- model.matrix(~ 0 + Subject:Tissue + Tissue:Time, data = clin)
colnames(des)[31:36] <- paste0(rep(c('wk00_vs_wk01_', 'wk00_vs_wk12_'), each = 3),
                               unique(clin$Tissue))
v <- voomWithQualityWeights(y, des)
fit <- eBayes(lmFit(v, des)) 

# Execute in parallel
foreach(j = colnames(des)[31:36]) %dopar% res(j)


