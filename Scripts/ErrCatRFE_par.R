### PREDICTIVE POWER BY DATA TYPE
### Categorical Response

# Load libraries, register cores 
library(data.table)
library(plyr)
library(randomForest)
library(caret)
library(dplyr)
library(doMC)
registerDoMC(cores = 20)

# Set seeds
set.seed(123)
tr_seeds <- vector(mode = 'list', length = 11)
for (i in seq_len(10)) {
  tr_seeds[[i]] <- sample.int(1000, 3)
}
tr_seeds[[11]] <- sample.int(1000, 1)
rfe_seeds <- vector(mode = 'list', length = 11)
for (i in seq_len(10)) {
  rfe_seeds[[i]] <- sample.int(1000, 20)
}
rfe_seeds[[11]] <- sample.int(1000, 1)

# Import data
mRNA <- readRDS('./Data/mRNA_RawCounts.rds')
mRNA <- mRNA[, grepl('wk00', colnames(mRNA))]
miRNA <- readRDS('./Data/miRNA_RawCounts.rds')
miRNA <- miRNA[, grepl('wk00', colnames(miRNA))]
prot <- readRDS('./Data/Prot_RawDat.rds')
prot <- prot[, grepl('wk00', colnames(prot))]

# Define response
clin <- fread('./Data/Clin_Baseline.csv') 
y <- as.factor(clin$PASI_75)

# Helper functions
trCtrl <- trainControl(method = 'cv', summaryFunction = mnLogLoss, 
                       classProbs = TRUE, seeds = tr_seeds)
my_rfFuncs <- rfFuncs
my_rfFuncs$summary <- mnLogLoss
rfeCtrl <- rfeControl(functions = my_rfFuncs, rerank = TRUE, 
                      method = 'cv', seeds = rfe_seeds)
subsets <- function(x) {
  out <- 10 + ((x - 10) / 400) * seq_len(19)^2
  return(round(out))
}
fill <- function(data_type, x) {
  out <- rep(NA, 1000)
  success <- FALSE
  while (!success) {
    if (data_type == 'Clinical') {
      fit <- train(x, y, method = 'rf', trControl = trCtrl, 
                   tuneGrid = data.frame(.mtry = 4:6), metric = 'logLoss')
    } else {
      fit <- rfe(x, y, sizes = subsets(ncol(x)), rfeControl = rfeCtrl, 
                 metric = 'logLoss', maximize = FALSE)
      cat(paste('Following RFE, the optimal random forest model for', data_type,
                'contains', fit$bestSubset, 'probes.'))
    }
    NA_idx <- which(is.na(out))
    if (length(NA_idx) > 0) {
      i <- min(NA_idx)
      if (is(fit, 'rfe')) {
        res <- fit$resample %>% filter(Variables == fit$bestSubset)
      } else {
        res <- fit$resample
      }
      out[i:(i + 9)] <- res$logLoss
    } else {
      success <- TRUE
    }
  }
  return(out)
}

# Error by data type
loss <- function(data_type) {
  
  if (data_type == 'Clinical') {
    
    x <- clin %>% select(-Subject, -DeltaPASI, -PASI_75)
    x <- model.matrix(~., data = x)
    x <- x[, -1]
    out <- data_frame(Clinical = fill(data_type, x))
    return(out)
    
  } else if (data_type == 'Blood_mRNA') {
    
    mat <- mRNA[, grepl('Blood', colnames(mRNA))]
    keep <- rowSums(cpm(mat) > 1) >= 3
    mat <- DGEList(mat[keep, ])
    mat <- calcNormFactors(mat)
    mat <- cpm(mat, log = TRUE, prior.count = 1)
    x <- t(mat)
    out <- data_frame(Blood_mRNA = fill(data_type, x))
    return(out)
    
  } else if (data_type == 'Lesional_mRNA') {
    
    mat <- mRNA[, grepl('Lesional', colnames(mRNA))]
    keep <- rowSums(cpm(mat) > 1) >= 3
    mat <- DGEList(mat[keep, ])
    mat <- calcNormFactors(mat)
    mat <- cpm(mat, log = TRUE, prior.count = 1)
    x <- t(mat)
    out <- data_frame(Lesional_mRNA = fill(data_type, x))
    return(out)
    
  } else if (data_type == 'Nonlesional_mRNA') {
    
    mat <- mRNA[, grepl('Nonlesional', colnames(mRNA))]
    keep <- rowSums(cpm(mat) > 1) >= 3
    mat <- DGEList(mat[keep, ])
    mat <- calcNormFactors(mat)
    mat <- cpm(mat, log = TRUE, prior.count = 1)
    x <- t(mat)
    out <- data_frame(Nonlesional_mRNA = fill(data_type, x))
    return(out)
    
  } else if (data_type == 'Blood_miRNA') {
    
    mat <- miRNA
    keep <- rowSums(cpm(mat) > 1) >= 3
    mat <- DGEList(mat[keep, ])
    mat <- calcNormFactors(mat)
    mat <- cpm(mat, log = TRUE, prior.count = 1)
    x <- t(mat)
    out <- data_frame(Blood_miRNA = fill(data_type, x))
    return(out)
    
  } else if (data_type == 'Blood_Proteomics') {
    
    mat <- log2(prot + 1)
    x <- t(mat)
    out <- data_frame(Blood_Proteomics = fill(data_type, x))
    return(out)
    
  }
  
}

# Execute in parallel
data_types <- c('Clinical', 'Blood_Proteomics', 'Blood_miRNA',
                'Blood_mRNA', 'Lesional_mRNA', 'Nonlesional_mRNA')
out <- foreach(d = data_types, .combine = cbind) %dopar% loss(d)
fwrite(out, './Results/ErrCatRFE.csv')


