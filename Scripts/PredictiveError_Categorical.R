### ERROR BY DATA TYPE
### Categorical Response

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
y <- as.factor(clin$PASI_75)

# Helper functions
cross_entropy <- function(data, lev = NULL, model = NULL) {
  data <- data[complete.cases(data), ]
  pred <- data[, lev[2]]
  eps <- 1e-15
  pred <- pmin(pmax(pred, eps), 1 - eps)
  obs <- ifelse(data$obs == lev[2], 1, 0)
  ce <- -(sum(obs * log(pred) + (1 - obs) * log(1 - pred))) / length(obs)
  c(cross_entropy = ce)
}
trCtrl <- trainControl(method = 'cv', summaryFunction = cross_entropy, 
                       classProbs = TRUE, seeds = tr_seeds)
my_rfFuncs <- rfFuncs
my_rfFuncs$summary <- cross_entropy
my_rfFuncs$fit <- function(x, y, first, last, ...) {
  randomForest(x, y, ntree = 1000, importance = TRUE, ...)
}
my_rfFuncs$rank <- function(object, x, y) {
  imp <- importance(object, type = 1, scale = FALSE)
  data_frame(Overall = imp[, 1], 
             var = rownames(imp)) %>%
    arrange(desc(Overall))
}
rfeCtrl <- rfeControl(functions = my_rfFuncs, method = 'cv', seeds = rfe_seeds)
subsets <- function(x) {
  tmp <- round(10 + ((x - 10) / 50^2.5) * seq_len(49)^2.5)
  out <- tmp[c(seq_len(30), seq(32, 48, 2))]
  return(out)
}
fill <- function(data_type, x) {
  if (data_type == 'Clinical') {
    fit <- train(x, y, method = 'rf', trControl = trCtrl, 
                 tuneGrid = data.frame(.mtry = 4:6), 
                 metric = 'cross_entropy', maximize = FALSE)
  } else {
    fit <- rfe(x, y, sizes = subsets(ncol(x)), rfeControl = rfeCtrl, 
               metric = 'cross_entropy', maximize = FALSE)
  }
  loss <- fit$resample$cross_entropy
  if (is(fit, 'rfe')) {
    tune <- fit$bestSubset
  } else {
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
fwrite(out, './Results/ErrCat_RFE.csv')


### ADD: TX_3TISSUE, BLOOD_3PLATFORM, AND ALL
### May have to use a different variable ranking procedure for this


