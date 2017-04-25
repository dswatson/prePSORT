# Load libraries, register cores, set seed
library(data.table)
library(limma)
library(wordspace)
library(Rtsne)
library(cluster)
library(dplyr)
library(doMC)
registerDoMC(5)
set.seed(123)

# Create summary data frame
out <- fread('./Data/Clinical.csv') %>%
  distinct(Subject) 

# Define clusters function
clusters <- function(plat, supervised = TRUE) {
  
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
  
  # Median center data, filter 
  mat <- mat[, samples]
  mat <- sweep(mat, 1, apply(mat, 1, median))
  top <- round(0.5 * nrow(mat))
  if (supervised) {
    hits <- topTable(fit, coef = coef, number = top, sort.by = 'p')
    genes <- rownames(hits)
    dm <- dist.matrix(t(mat[genes, ]), method = 'euclidean')
  } else {
    dm <- matrix(nrow = ncol(mat), ncol = ncol(mat))
    top_idx <- nrow(mat) - top + 1
    for (i in 2:ncol(mat)) {
      for (j in 1:(i - 1)) {
        dm[i, j] <- sqrt(sum(sort.int((mat[, i] - mat[, j])^2,
                                      partial = top_idx)[top_idx:nrow(mat)]))
      }
    }
  }
  
  # Embed samples
  tsne <- Rtsne(as.dist(dm), perplexity = 2, theta = 0.1, 
                check_duplicates = FALSE, is_distance = TRUE)
  
  # Run PAM
  out$new <- pam(tsne$Y, k = 2, cluster.only = TRUE)
  out$new <- paste0('C', out$new)
  colnames(out)[colnames(out) == 'new'] <- plat
  return(out)
  
}

# Execute in parallel
combo <- function(x, y) inner_join(x, y, by = 'Subject')
plats <- c('mRNA_Blood', 'mRNA_Lesional', 'mRNA_Nonlesional',
           'miRNA_Blood', 'Prot_Blood')
foreach(plat = plats, .combine = combo) %dopar% 
  clusters(plat, supervised = TRUE) %>%
  fwrite('./Results/Clusters/Supervised/PAM.csv')
foreach(plat = plats, .combine = combo) %dopar% 
  clusters(plat, supervised = FALSE) %>%
  fwrite('./Results/Clusters/Unsupervised/PAM.csv')


