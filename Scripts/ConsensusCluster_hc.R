# Load libraries, register cores, set seed
library(data.table)
library(limma)
library(fastcluster)
library(ConsensusClusterPlus)
library(dplyr)
library(doMC)
registerDoMC(5)
set.seed(123)

# Create clusters data frame
clusters <- data.frame(Subject = c(paste0('S0', seq_len(9)),
                                   paste0('S', c(10, 11)))) %>%
  filter(Subject != 'S06')

# Define consensus function
consensus <- function(plat) {
  
  # Load data
  if (grepl('mRNA', plat)) {
    mat <- readRDS('./Data/mat_mRNA.rds')
    fit <- readRDS('./Data/fit_mRNA.rds')
    clin <- fread('./Data/Clinical.csv')
    if (plat == 'mRNA_Blood') {
      coef <- colnames(fit)[grepl('Blood.wk00', colnames(fit))]
      samples <- grep('Blood.wk00', colnames(mat))
    } else if (plat == 'mRNA_Lesional') {
      coef <- colnames(fit)[grepl('Lesional.wk00', colnames(fit))]
      samples <- grep('Lesional.wk00', colnames(mat))
    } else if (plat == 'mRNA_Nonlesional') {
      coef <- colnames(fit)[grepl('Nonlesional.wk00', colnames(fit))]
      samples <- grep('Nonlesional.wk00', colnames(mat))
    }
  } else if (plat == 'miRNA_Blood') {
    mat <- readRDS('./Data/mat_miRNA.rds')
    fit <- readRDS('./Data/fit_miRNA.rds')
    coef <- colnames(fit)[grepl('wk00.Response', colnames(fit))]
    samples <- grep('wk00', colnames(mat))
  } else if (plat == 'Prot_Blood') {
    mat <- readRDS('./Data/mat_prot.rds')
    fit <- readRDS('./Data/fit_prot.rds')
    coef <- colnames(fit)[grepl('wk00.Response', colnames(fit))]
    samples <- grep('wk00', colnames(mat))
  }
  
  # Consensus cluster on top 1% of probes
  n_probes <- round(0.01 * nrow(mat))
  top <- topTable(fit, coef = coef, number = n_probes, sort.by = 'p')
  probes <- rownames(top)
  mat <- mat[probes, samples]
  mat <- sweep(mat, 1, apply(mat, 1, median))
  cc <- ConsensusClusterPlus(mat, maxK = 4, reps = 1e4, pItem = 0.9, pFeature = 0.9,
                             plot = 'png')
  
  # Determine optimal k from PAC index
  suppressWarnings(
    PAC <- expand.grid(k = 2:4, Idx = c(0.1, 0.9)) %>%
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
  clusters$new <- cutree(hc, k = k)
  colnames(clusters)[colnames(clusters) == 'new'] <- plat
  return(clusters)
  
}

# Execute in parallel
combo <- function(x, y) inner_join(x, y, by = 'Subject')
plats <- c('mRNA_Blood', 'mRNA_Lesional', 'mRNA_Nonlesional',
           'miRNA_Blood', 'Prot_Blood')
clusters <- foreach(plat = plats, .combine = combo) %dopar% consensus(plat)
fwrite(clusters, './Results/ConsensusClusters/clusters.csv')


