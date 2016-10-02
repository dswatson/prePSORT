# Load libraries
library(data.table)
library(tximport)
library(sva)
library(DESeq2)
library(dplyr)

# If parallel
library(BiocParallel)
register(MulticoreParam(12))

# Prep data
pheno <- fread(paste0(getwd(), '/Data/PrePSORT_Clinical.csv'))
t2g <- fread(paste0(getwd(), '/Data/Ensembl.Hs79.Tx.csv'))
e2g <- fread(paste0(getwd(), '/Data/Ensembl.Hs79.GeneSymbols.csv'))

# TxImport
dir.B <- paste0(getwd(), '/Data/Blood')
files.B <- file.path(dir.B, pheno$Sample, 'abundance.tsv')
txi.B <- tximport(files.B, type = 'kallisto', tx2gene = t2g, reader = fread)
dir.L <- paste0(getwd(), '/Data/LesionalSkin')
files.L <- file.path(dir.L, pheno$Sample, 'abundance.tsv')
txi.L <- tximport(files.L, type = 'kallisto', tx2gene = t2g, reader = fread)
dir.N <- paste0(getwd(), '/Data/NonlesionalSkin')
files.N <- file.path(dir.N, pheno$Sample, 'abundance.tsv')
txi.N <- tximport(files.N, type = 'kallisto', tx2gene = t2g, reader = fread)
txi <- list(abundance = cbind(txi.B$abundance, txi.L$abundance, txi.N$abundance),
            counts = cbind(txi.B$counts, txi.L$counts, txi.N$counts),
	          length = cbind(txi.B$length, txi.L$length, txi.N$length),
	          countsFromAbundance = 'no')
rm(dir.B, files.B, txi.B, dir.L, files.L, txi.L, dir.N, files.N, txi.N)

# Big pheno
pheno <- rbind(pheno, pheno, pheno) %>%
  mutate(Site = c(rep('B', 30), rep('L', 30), rep('N', 30)),
         Sample = paste(Sample, Site, sep = '.'),
         Time.Site = paste(Time, Site, sep = '.'))

# Normalise, filter
dds <- DESeqDataSetFromTximport(txi, colData=pheno, design = ~ 1)
dds <- estimateSizeFactors(dds)
mat <- counts(dds, normalized = TRUE)
keep <- rowMeans(mat) > 1
mat <- mat[keep, ]
dds <- dds[keep, ]

# SVA
mod <- model.matrix(~ 0 + Site + Subject + Time, data = pheno)
mod0 <- model.matrix(~ Subject + Time, data = pheno)
svobj <- svaseq(mat, mod, mod0)
des <- cbind(mod, svobj$sv)

# DEseq
dds <- estimateDispersions(dds, modelMatrix = des, maxit = 10000)
dds <- nbinomWaldTest(dds, modelMatrix = des, maxit = 10000, betaPrior = FALSE)

# Blood vs. LesionalSkin
res <- data.frame(results(dds, 
                          contrast  = list('SiteL', 'SiteB'), 
                          filterfun = ihw))
id <- rownames(res)
res <- res %>%
  mutate(gene_id = id) %>%
  inner_join(e2g, by = 'gene_id') %>%
  rename(EnsemblID  = gene_id,
         GeneSymbol = gene_name, 
         AvgExpr    = baseMean,
         logFC      = log2FoldChange, 
         p.value    = pvalue,
         q.value    = padj) %>%
  arrange(p.value) %>%
  select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
fwrite(res, paste0(getwd(), '/Results/BvsL.txt'), sep = '\t')

# Blood vs. NonlesionalSkin
res <- data.frame(results(dds, 
                          contrast  = list('SiteN', 'SiteB'), 
                          filterfun = ihw))
id <- rownames(res)
res <- res %>%
  mutate(gene_id = id) %>%
  inner_join(e2g, by = 'gene_id') %>%
  rename(EnsemblID  = gene_id,
         GeneSymbol = gene_name, 
         AvgExpr    = baseMean,
         logFC      = log2FoldChange, 
         p.value    = pvalue,
         q.value    = padj) %>%
  arrange(p.value) %>%
  select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
fwrite(res, paste0(getwd(), '/Results/BvsN.txt'), sep = '\t')

# LesionalSkin vs NonlesionalSkin
res <- data.frame(results(dds, 
                          contrast  = list('SiteL', 'SiteN'), 
                          filterfun = ihw))
id <- rownames(res)
res <- res %>%
  mutate(gene_id = id) %>%
  inner_join(e2g, by = 'gene_id') %>%
  rename(EnsemblID  = gene_id,
         GeneSymbol = gene_name, 
         AvgExpr    = baseMean,
         logFC      = log2FoldChange, 
         p.value    = pvalue,
         q.value    = padj) %>%
  arrange(p.value) %>%
  select(EnsemblID, GeneSymbol, AvgExpr, logFC, p.value, q.value)
fwrite(res, paste0(getwd(), '/Results/LvsN.txt'), sep = '\t')



