# Load libraries, register cores, set seed
library(data.table)
library(tximport)
library(edgeR)
library(limma)
library(fastcluster)
library(ConsensusClusterPlus)
library(ggplot2)
library(ggsci)
library(dplyr)
library(doMC)
registerDoMC(10)
set.seed(123)

# Import data
clin <- fread('./Data/Clinical.csv') %>%
  mutate(Subject.Tissue = paste(Subject, Tissue, sep = '.'))
t2g <- fread('./Data/Hs79.t2g.csv')
files <- file.path('./Data/RawCounts', clin$Sample, 'abundance.tsv')
txi <- tximport(files, type = 'kallisto', tx2gene = t2g, reader = fread, 
                countsFromAbundance = 'lengthScaledTPM')

# Collapse, filter counts
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
clin <- clin %>% mutate(DeltaPASI = winsorise(DeltaPASI))

# Fit model
des <- model.matrix(~ 0 + Tissue:Time + Tissue:Time:DeltaPASI, data = clin)
coefs <- paste(unique(clin$Tissue), rep(unique(clin$Time), each = 3),
               'Response', sep = '.')
colnames(des)[10:18] <- coefs
v <- voomWithQualityWeights(y, des)
icc <- duplicateCorrelation(v, des, block = clin$Subject.Tissue)
v <- voomWithQualityWeights(y, des, correlation = icc$cor, 
                            block = clin$Subject.Tissue)
icc <- duplicateCorrelation(v, des, block = clin$Subject.Tissue)  
fit <- lmFit(v, des, correlation = icc$cor, block = clin$Subject.Tissue)
fit <- eBayes(fit)

# Consensus Clusters
clusters <- clin %>% 
  select(Subject) %>%
  unique()
consensus <- function(coef) {
  
  # Consensus cluster on top 1% of genes
  top <- topTable(fit, coef = coef, number = Inf, sort.by = 'p') 
  n_genes <- round(0.01 * nrow(top))
  genes <- rownames(top)[seq_len(n_genes)]
  tissue.time <- gsub('.Response', '', coef)
  samples <- grep(tissue.time, clin$Sample)
  mat <- v$E[genes, samples]
  cc <- ConsensusClusterPlus(mat, maxK = 5, reps = 1e4, pItem = 0.9, pFeature = 0.9,
                             plot = 'png')
  
  # Determine optimal k from PAC index
  suppressWarnings(
    PAC <- expand.grid(k = 2:5, Idx = c(0.1, 0.9)) %>%
      rowwise() %>%
      mutate(CDF = ecdf(unlist(cc[[k]][['consensusMatrix']]))(Idx),
               k = as.factor(k)) %>%
      group_by(k) %>%
      mutate(PAC = diff(CDF)) %>%
      select(k, PAC) %>%
      unique() %>%
      as.data.frame()
  )
  k <- PAC %>%
    mutate(k = as.numeric(k) + 1) %>%
    filter(PAC == min(PAC)) %>%
    filter(k == min(k)) %>%
    select(k) %>%
    as.numeric() 
  
  # Fill clusters data frame
  dm <- as.dist(1 - cor(mat))
  hc <- hclust(dm, method = 'average')
  assignment <- cutree(hc, k = k)
  if (tissue.time == 'Lesional.wk12') {
    assignment <- c(assignment, NA)
  }
  clusters$new <- assignment
  colnames(clusters)[colnames(clusters) == 'new'] <- coef
  return(clusters)
  
}

# Execute in parallel
combo <- function(x, y) inner_join(x, y, by = 'Subject')
clusters <- foreach(coef = coefs, .combine = combo) %dopar% consensus(coef)
fwrite(clusters, './Results/ConsensusClusters/clusters.csv')

