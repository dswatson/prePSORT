### PREDICTIVE POWER BY DATA TYPE
### Categorical Response

# Load libraries, register cores, set seed
library(data.table)
library(limma)
library(edgeR)
library(caret)
library(dplyr)
library(doMC)
registerDoMC(cores = 10)
set.seed(123)

# Define response
clin <- fread('./Data/Clin_Baseline.csv') 
y <- as.factor(clin$PASI_75)

# Helper functions
trCtrl <- trainControl(method = 'cv', summaryFunction = mnLogLoss, 
                       classProbs = TRUE)
fill <- function(data_type) {
  out <- rep(NA, 1000)
  success <- FALSE
  while (!success) {
    if (data_type == 'Clinical') {
      fit <- train(x, y, method = 'rf', trControl = trCtrl, 
                      tuneGrid = data.frame(.mtry = 4:6), metric = 'logLoss')
    } else if (grepl('mRNA', data_type)) {
      fit <- train(x, y, method = 'pam', trControl = trCtrl,
                   tuneLength = 30, metric = 'logLoss')
    } else {
      fit <- sbf(x, y, sbfControl = sbfCtrl)
    }
    NA_idx <- which(is.na(out))
    if (length(NA_idx) > 0) {
      i <- min(NA_idx)
      out[i:(i + 9)] <- fit$resample$logLoss
    } else {
      success <- TRUE
    }
  }
  return(out)
}

# Clinical
x <- clin %>% select(-Subject, -DeltaPASI, -PASI_75)
x <- model.matrix(~., data = x)
x <- x[, -1]
err_cat <- data_frame(Clinical = fill('Clinical'))

# Blood_mRNA
mRNA <- readRDS('./Data/mat_mRNA.rds')
mat <- mRNA[, grepl('Blood.wk00', colnames(mRNA))]
x <- t(mat)
err_cat$Blood_mRNA <- fill('Blood_mRNA')

# Lesional_mRNA
mat <- mRNA[, grepl('Lesional.wk00', colnames(mRNA))]
x <- t(mat)
err_cat$Lesional_mRNA <- fill('Lesional_mRNA')

# Nonlesional_mRNA
mat <- mRNA[, grepl('Nonlesional.wk00', colnames(mRNA))]
x <- t(mat)
err_cat$Nonlesional_mRNA <- fill('Nonlesional_mRNA')

# Blood_miRNA
filter_fns <- list(
  summary = mnLogLoss,
  fit = function(x, y) {
    x <- miRNA_mat[colnames(x), rownames(x)]
    default_mtry <- floor(sqrt(nrow(x)))
    tunes <- data.frame(.mtry = c(round(0.9 * default_mtry), 
                                  default_mtry, 
                                  round(1.1 * default_mtry)))
    train(t(x), y, method = 'rf', trControl = trCtrl, tuneGrid = tunes, 
          metric = 'logLoss')
  },
  pred = function(object, x) {
    cbind(data.frame(pred = predict(object, x)), 
          predict(object, x, type = 'prob'))
  },
  score = function(x, y) {  
    des <- model.matrix(~ y)
    v <- voom(t(x), des)
    mod <- eBayes(lmFit(v, des), robust = TRUE)
    top <- topTable(mod, coef = 2, number = Inf)
    p <- top$P.Value
    names(p) <- rownames(top)
    return(p)
  },
  filter = function(score, x, y) {
    score <= 0.05
  }
)
sbfCtrl <- sbfControl(method = 'cv', functions = filter_fns, multivariate = TRUE)
miRNA <- read.csv('./Data/miRNA.csv', row.names = 1) %>% as.matrix()
keep <- rowSums(cpm(miRNA) > 1) >= 3
miRNA <- DGEList(miRNA[keep, ])
miRNA <- calcNormFactors(miRNA)
mat <- miRNA[, grepl('wk00', colnames(miRNA))]
miRNA_mat <- readRDS('./Data/mat_miRNA.rds')
miRNA_mat <- miRNA_mat[, grepl('wk00', colnames(miRNA_mat))]
x <- t(mat$counts) / mat$samples$norm.factors
err_cat$Blood_miRNA <- fill('Blood_miRNA')

# Blood_Proteomics
filter_fns$fit <- function(x, y) {
  default_mtry <- floor(sqrt(ncol(x)))
  tunes <- data.frame(.mtry = c(round(0.9 * default_mtry), 
                                default_mtry,
                                round(1.1 * default_mtry)))
  train(x, y, method = 'rf', trControl = trCtrl, tuneGrid = tunes,
        metric = 'logLoss')
}
filter_fns$score <- function(x, y) {
  mod <- eBayes(lmFit(t(x), model.matrix(~ y)), robust = TRUE)
  top <- topTable(mod, coef = 2, number = Inf)
  p <- top$P.Value
  names(p) <- rownames(top)
  return(p)
}
sbfCtrl <- sbfControl(method = 'cv', functions = filter_fns, multivariate = TRUE)
prot <- readRDS('./Data/mat_prot.rds')
x <- t(prot[, grepl('wk00', colnames(prot))])
err_cat$Blood_Proteomics <- fill('Blood_Proteomics')

# Export
fwrite(err_cat, './Results/err_cat.csv')


