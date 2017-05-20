### ERROR BY DATA TYPE
### Continuous Response

# Load libraries, register cores
library(data.table)
library(plyr)
library(randomForest)
library(caret)
library(edgeR)
library(dplyr)
library(doMC)
registerDoMC(4)

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

# Helper functions
trCtrl <- trainControl(method = 'cv', seeds = tr_seeds)
my_rfFuncs <- rfFuncs
my_rfFuncs$fit <- function(x, y, first, last, ...) {
  randomForest(x, y, ntree = 1000, importance = TRUE, ...)
}
my_rfFuncs$rank <- function(object, x, y) {
  imp <- importance(object, type = 1, scale = FALSE)
  df <- data_frame(Overall = imp[, 1], 
                   var = rownames(imp))
  if (any(!complete.cases(df))) {
    df$Overall[which(is.na(df$Overall))] <- median(df$Overall, na.rm = TRUE)
  }
  df %>% arrange(desc(Overall))
}
rfeCtrl <- rfeControl(functions = my_rfFuncs, rerank = TRUE, 
                      method = 'cv', seeds = rfe_seeds)
subsets <- function(x) {
  round(10 + ((x - 10) / 400) * seq_len(19)^2)
}
fill <- function(data_type, x) {
    if (data_type == 'Clinical') {
      fit <- train(x, y, method = 'rf', tuneGrid = data.frame(.mtry = 8:10),
                   trControl = trCtrl)
    } else {
      fit <- rfe(x, y, sizes = subsets(ncol(x)), rfeControl = rfeCtrl)
    }
    if (is(fit, 'rfe')) {
      res <- fit$resample %>% filter(Variables == fit$bestSubset)
      loss <- res$RMSE
      tune <- res$Variables
    } else {
      loss <- fit$resample$RMSE 
      tune <- fit$bestTune$mtry
    }
  out <- data_frame(Loss = loss, Tune = tune)
  colnames(out) <- paste(data_type, colnames(out), sep = '_')
  return(out)
}

# Error by data type
loss <- function(data_type) {
  
  if (data_type == 'Clinical') {
    
    x <- clin %>% select(-Subject, -DeltaPASI, -PASI_75)
    x <- model.matrix(~., data = x)
    x <- x[, -1]
    return(fill(data_type, x))
    
  } else if (data_type == 'Blood_mRNA') {
    
    mat <- mRNA[, grepl('Blood', colnames(mRNA))]
    keep <- rowSums(cpm(mat) > 1) >= 3
    mat <- DGEList(mat[keep, ])
    mat <- calcNormFactors(mat)
    mat <- cpm(mat, log = TRUE, prior.count = 1)
    x <- t(mat)
    return(fill(data_type, x))
    
  } else if (data_type == 'Lesional_mRNA') {
    
    mat <- mRNA[, grepl('Lesional', colnames(mRNA))]
    keep <- rowSums(cpm(mat) > 1) >= 3
    mat <- DGEList(mat[keep, ])
    mat <- calcNormFactors(mat)
    mat <- cpm(mat, log = TRUE, prior.count = 1)
    x <- t(mat)
    return(fill(data_type, x))
    
  } else if (data_type == 'Nonlesional_mRNA') {
    
    mat <- mRNA[, grepl('Nonlesional', colnames(mRNA))]
    keep <- rowSums(cpm(mat) > 1) >= 3
    mat <- DGEList(mat[keep, ])
    mat <- calcNormFactors(mat)
    mat <- cpm(mat, log = TRUE, prior.count = 1)
    x <- t(mat)
    return(fill(data_type, x))
    
  } else if (data_type == 'Blood_miRNA') {
    
    mat <- miRNA
    keep <- rowSums(cpm(mat) > 1) >= 3
    mat <- DGEList(mat[keep, ])
    mat <- calcNormFactors(mat)
    mat <- cpm(mat, log = TRUE, prior.count = 1)
    x <- t(mat)
    return(fill(data_type, x))
    
  } else if (data_type == 'Blood_Proteomics') {
    
    mat <- log2(prot)
    x <- t(mat)
    return(fill(data_type, x))
    
  }
  
}

# Execute in parallel
data_types <- c('Clinical', 'Blood_Proteomics', 'Blood_miRNA',
                'Blood_mRNA', 'Lesional_mRNA', 'Nonlesional_mRNA')
out <- foreach(d = data_types, .combine = cbind) %dopar% loss(d)
fwrite(out, './Results/ErrCont_RFE.csv')


