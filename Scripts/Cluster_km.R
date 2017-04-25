# Load libraries, register cores, set seed
library(data.table)
library(limma)
library(Rtsne)
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
    n_probes <- round(0.1 * nrow(mat))
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
    n_probes <- round(0.1 * nrow(mat))
  } else if (plat == 'Prot_Blood') {
    mat <- readRDS('./Data/mat_prot.rds')
    fit <- readRDS('./Data/fit_prot.rds')
    coef <- colnames(fit)[grepl('wk00.Response', colnames(fit))]
    samples <- grep('wk00', colnames(mat))
    n_probes <- nrow(mat)
  }
  
  # Consensus cluster on t-SNE projection of top 10% of probes
  top <- topTable(fit, coef = coef, number = n_probes, sort.by = 'p')
  probes <- rownames(top)
  mat <- mat[probes, samples]
  dm <- dist(t(mat))
  tsne <- Rtsne(dm, perplexity = 2, theta = 0.1, check_duplicates = FALSE, 
                is_distance = TRUE)
  
  # Fill clusters data frame
  clusters$new <- kmeans(tsne$Y, 2)$cluster
  colnames(clusters)[colnames(clusters) == 'new'] <- plat
  return(clusters)
  
}

# Execute in parallel
combo <- function(x, y) inner_join(x, y, by = 'Subject')
plats <- c('mRNA_Blood', 'mRNA_Lesional', 'mRNA_Nonlesional',
           'miRNA_Blood', 'Prot_Blood')
clusters <- foreach(plat = plats, .combine = combo) %dopar% consensus(plat)
fwrite(clusters, './Results/ConsensusClusters/clusters_km.csv')


