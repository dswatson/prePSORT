### CLUSTER CONCORDANCE PLOTS ###

# Load libraries
library(data.table)
library(infotheo)
library(corrplot)
library(dplyr)

# Import clinical stratifier
clin <- fread('./Data/Clinical.csv') %>%
  distinct(Subject, PASI_75) 

### Supervised Clusters ###

# Load cluster assignments
df <- read.csv('./Results/Clusters/Supervised.csv', row.names = 1) %>%
  mutate(PASI_75 = clin$PASI_75)

# Create mutual information matrix
mat <- matrix(nrow = ncol(df), ncol = ncol(df), 
              dimnames = list(colnames(df), colnames(df)))
for (i in 2:nrow(mat)) {
  for (j in 1:(i - 1)) {
    mat[i, j] <- mutinformation(df[[i]], df[[j]])
  }
}
mat <- natstobits(mat)
mat[is.na(mat)] <- max(mat, na.rm = TRUE)

# Plot
corrplot(mat, method = 'shade', is.corr = FALSE, type = 'lower', diag = FALSE,  
         title = 'Cluster Concordance: Supervised', mar = c(0, 1, 1, 0),
         addgrid.col = 'black', addCoef.col = 'black',
         tl.cex = 0.65, tl.col = 'black', tl.srt = 45, cl.lim = c(0, 1))

### Unsupervised Clusters ###

# Load cluster assignments
df <- read.csv('./Results/Clusters/Unsupervised.csv', row.names = 1) %>%
  mutate(PASI_75 = clin$PASI_75)

# Create mutual information matrix
mat <- matrix(nrow = ncol(df), ncol = ncol(df), 
              dimnames = list(colnames(df), colnames(df)))
for (i in 2:nrow(mat)) {
  for (j in 1:(i - 1)) {
    mat[i, j] <- mutinformation(df[[i]], df[[j]])
  }
}
mat <- natstobits(mat)
mat[is.na(mat)] <- max(mat, na.rm = TRUE)

# Plot
corrplot(mat, method = 'shade', is.corr = FALSE, type = 'lower', diag = FALSE,  
         title = 'Cluster Concordance: Unsupervised', mar = c(0, 1, 1, 0),
         addgrid.col = 'black', addCoef.col = 'black',
         tl.cex = 0.65, tl.col = 'black', tl.srt = 45, cl.lim = c(0, 1))


