### PREDICTIVE POWER BY DATA TYPE
### Continuous Response

# Load libraries, register cores, set seed
library(data.table)
library(tximport)
library(limma)
library(edgeR)
library(caret)
library(dplyr)
library(doMC)
registerDoMC(cores = 10)
set.seed(123)

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
trCtrl <- trainControl(method = 'cv')
fill <- function(data_type) {
  out <- rep(NA, 1000)
  success <- FALSE
  while (!success) {
    if (data_type == 'Clinical') {
      fit <- train(x, y, method = 'rf', trControl = trCtrl, 
                   tuneGrid = data.frame(.mtry = 4:6))
    } else {
      fit <- sbf(x, y, sbfControl = sbfCtrl)
    }
    NA_idx <- which(is.na(out))
    if (length(NA_idx) > 0) {
      i <- min(NA_idx)
      out[i:(i + 9)] <- fit$resample$RMSE
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
err_cont <- data_frame(Clinical = fill('Clinical'))

#mRNA
filter_fns <- list(
  summary = defaultSummary,
  fit = function(x, y) {
    x <- mRNA_mat[colnames(x), rownames(x)]
    default_mtry <- floor(sqrt(nrow(x)))
    tunes <- data.frame(.mtry = c(round(0.9 * default_mtry), 
                                  default_mtry, 
                                  round(1.1 * default_mtry)))
    train(t(x), y, method = 'rf', trControl = trCtrl, tuneGrid = tunes)
  },
  pred = function(object, x) {
    predict(object, x)
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
    score <= 0.01
  }
)
sbfCtrl <- sbfControl(method = 'cv', functions = filter_fns, multivariate = TRUE)
clin <- fread('./Data/Clinical.csv')
t2g <- fread('./Data/Hs79.t2g.csv')
files <- file.path('./Data/RawCounts', clin$Sample, 'abundance.tsv')
txi <- tximport(files, type = 'kallisto', tx2gene = t2g, reader = fread, 
                countsFromAbundance = 'lengthScaledTPM')
mRNA <- txi$counts
colnames(mRNA) <- clin$Sample
mRNA_mat <- readRDS('./Data/mat_mRNA.rds')

# Blood_mRNA
mat <- mRNA[, grepl('Blood.wk00', colnames(mRNA))]
keep <- rowSums(cpm(mat) > 1) >= 3
mat <- DGEList(mat[keep, ])
mat <- calcNormFactors(mat)
x <- t(mat$counts) / mat$samples$norm.factors
err_cont$Blood_mRNA <- fill('Blood_mRNA')

# Lesional_mRNA
mat <- mRNA[, grepl('Lesional.wk00', colnames(mRNA))]
keep <- rowSums(cpm(mat) > 1) >= 3
mat <- DGEList(mat[keep, ])
mat <- calcNormFactors(mat)
x <- t(mat$counts) / mat$samples$norm.factors
err_cont$Lesional_mRNA <- fill('Lesional_mRNA')

# Nonlesional_mRNA
mat <- mRNA[, grepl('Nonlesional.wk00', colnames(mRNA))]
keep <- rowSums(cpm(mat) > 1) >= 3
mat <- DGEList(mat[keep, ])
mat <- calcNormFactors(mat)
x <- t(mat$counts) / mat$samples$norm.factors
err_cont$Nonlesional_mRNA <- fill('Nonlesional_mRNA')

# Blood_miRNA
filter_fns$fit <- function(x, y) {
  x <- miRNA_mat[colnames(x), rownames(x)]
  default_mtry <- floor(sqrt(nrow(x)))
  tunes <- data.frame(.mtry = c(round(0.9 * default_mtry), 
                                default_mtry, 
                                round(1.1 * default_mtry)))
  train(t(x), y, method = 'rf', trControl = trCtrl, tuneGrid = tunes)
}
filter_fns$filter <- function(score, x, y) {
  score <= 0.05
}
sbfCtrl <- sbfControl(method = 'cv', functions = filter_fns, multivariate = TRUE)
miRNA <- read.csv('./Data/miRNA.csv', row.names = 1) %>% as.matrix()
keep <- rowSums(cpm(miRNA) > 1) >= 3
miRNA <- DGEList(miRNA[keep, ])
miRNA <- calcNormFactors(miRNA)
mat <- miRNA[, grepl('wk00', colnames(miRNA))]
miRNA_mat <- readRDS('./Data/mat_miRNA.rds')
miRNA_mat <- miRNA_mat[, grepl('wk00', colnames(miRNA_mat))]
x <- t(mat$counts) / mat$samples$norm.factors
err_cont$Blood_miRNA <- fill('Blood_miRNA')

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
err_cont$Blood_Proteomics <- fill('Blood_Proteomics')

# Export
fwrite(err_cont, './Results/err_cont.csv')


