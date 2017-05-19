### CLUSTERS ###

# Load libraries, register cores, set seed
library(data.table)
library(limma)
library(edgeR)
library(wordspace)
library(Rtsne)
library(cluster)
library(dplyr)
library(doMC)
registerDoMC(4)
set.seed(123)

# Define continuous response
clin <- fread('./Data/Clin_Baseline.csv')
winsorise <- function(x, multiple = 2) {
  y <- x - median(x)
  lim <- mad(y, center = 0) * multiple
  y[y > lim] <- lim
  y[y < -lim] <- -lim
  y <- y + median(x)
  return(y)
}
clin <- clin %>% mutate(DeltaPASI = winsorise(DeltaPASI))
y <- clin$DeltaPASI

# Create output data frame
out <- clin %>% select(Subject)

# Import data
mRNA <- readRDS('./Data/mRNA_RawCounts.rds')
mRNA <- mRNA[, grepl('wk00', colnames(mRNA))]
miRNA <- readRDS('./Data/miRNA_RawCounts.rds')
miRNA <- miRNA[, grepl('wk00', colnames(miRNA))]
prot <- readRDS('./Data/Prot_RawDat.rds')
prot <- prot[, grepl('wk00', colnames(prot))]

# Define clusters function
clusters <- function(data_type, supervised = TRUE) {
  
  # Normalize, transform data
  if (data_type == 'Blood_mRNA') {
    mat <- mRNA[, grepl('Blood', colnames(mRNA))]
    keep <- rowSums(cpm(mat) > 1) >= 3
    mat <- DGEList(mat[keep, ])
    mat <- calcNormFactors(mat)
  } else if (data_type == 'Lesional_mRNA') {
    mat <- mRNA[, grepl('Lesional', colnames(mRNA))]
    keep <- rowSums(cpm(mat) > 1) >= 3
    mat <- DGEList(mat[keep, ])
    mat <- calcNormFactors(mat)
  } else if (data_type == 'Nonlesional_mRNA') {
    mat <- mRNA[, grepl('Nonlesional', colnames(mRNA))]
    keep <- rowSums(cpm(mat) > 1) >= 3
    mat <- DGEList(mat[keep, ])
    mat <- calcNormFactors(mat)
  } else if (data_type == 'Blood_miRNA') {
    mat <- miRNA
    keep <- rowSums(cpm(mat) > 1) >= 3
    mat <- DGEList(mat[keep, ])
    mat <- calcNormFactors(mat)
  } else if (data_type == 'Blood_Proteomics') {
    mat <- log2(prot)
  }
  
  # Filter, median centre data
  top <- round(0.5 * nrow(mat))
  if (supervised) {                              # By limma p-value
    des <- model.matrix(~ y)
    if (grepl('RNA', data_type)) {
      mat <- voom(mat, des)
    } 
    fit <- eBayes(lmFit(mat, des))
    hits <- topTable(fit, coef = 'y', number = top, sort.by = 'p')
    hits <- rownames(hits)
    mat <- mat[hits, ]
    if (grepl('RNA', data_type)) {
      mat <- mat$E
    }
    mat <- sweep(mat, 1, apply(mat, 1, median))
    dm <- dist.matrix(t(mat), method = 'euclidean')
  } else {                                       # By leading fold change
    if (grepl('RNA', data_type)) {
      mat <- cpm(mat, log = TRUE, prior.count = 1)
    }
    mat <- sweep(mat, 1, apply(mat, 1, median))
    dm <- matrix(nrow = ncol(mat), ncol = ncol(mat))
    for (i in 2:ncol(mat)) {
      for (j in 1:(i - 1)) {
        top_idx <- nrow(mat) - top + 1
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
  colnames(out)[colnames(out) == 'new'] <- data_type
  return(out)
  
}

# Execute in parallel
combo <- function(x, y) inner_join(x, y, by = 'Subject')
data_types <- c('Blood_Proteomics', 'Blood_miRNA',
                'Blood_mRNA', 'Lesional_mRNA', 'Nonlesional_mRNA')
foreach(d = data_types, .combine = combo) %dopar% 
  clusters(d, supervised = TRUE) %>%
  fwrite('./Results/SupervisedClusters.csv')
foreach(d = data_types, .combine = combo) %dopar% 
  clusters(d, supervised = FALSE) %>%
  fwrite('./Results/UnsupervisedClusters.csv')


