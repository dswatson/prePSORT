# Load libraries, register cores, set seed
library(data.table)
library(tximport)
library(edgeR)
library(limma)
library(qvalue)
library(qusage)
library(dplyr)
library(doParallel)
registerDoParallel(10)
set.seed(123)

# Import data
clin <- fread('./Data/Clinical.csv') 
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
    fwrite(paste0('./Results/Tissue', paste0(coef, '.Genes.txt')), sep = '\t')
  
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
    fwrite(paste0('./Results/Tissue/', paste0(coef, '.Modules.txt')), sep = '/t')
  
}

# Blood vs. skin
des <- model.matrix(~ 0 + Subject:Time + Tissue:Time, data = clin)
colnames(des)[31:36] <- c(paste0('Blood_vs_Lesional_', unique(clin$Time)), 
                          paste0('Blood_vs_Nonlesional_', unique(clin$Time)))
v <- voomWithQualityWeights(y, des)
fit <- eBayes(lmFit(v, des))
foreach(j = colnames(des)[31:36]) %dopar% res(j)

# Nonlesional vs. lesional
clin <- clin %>%
  mutate(Tissue = relevel(as.factor(Tissue), ref = 'Nonlesional'))
des <- model.matrix(~ 0 + Subject:Time + Tissue:Time, data = clin)
colnames(des)[34:36] <- paste0('Nonlesional_vs_Lesional_', unique(clin$Time))
v <- voomWithQualityWeights(y, des)
fit <- eBayes(lmFit(v, des))
foreach(j = colnames(des)[34:36]) %dopar% res(j)


