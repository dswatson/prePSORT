# Load libraries
library(corrplot)

# Import cluster assignments
df <- read.csv('./Results/ConsensusClusters/clusters.csv', row.names = 1)

# Create p-value matrix
mat <- matrix(nrow = ncol(df), ncol = ncol(df), 
              dimnames = list(colnames(df), colnames(df)))
dimnames(mat) <- lapply(dimnames(mat), function(i) gsub('.Response', '', i))
for (i in 2:nrow(mat)) {
  for (j in 1:(i - 1)) {
    mat[i, j] <- fisher.test(df[[i]], df[[j]])$p.value
  }
}

# Transform, fill NAs
mat <- -log10(mat)
mat[is.na(mat)] <- ceiling(max(mat, na.rm = TRUE))

# Plot
corrplot(mat, is.corr = FALSE, type = 'lower', diag = FALSE, method = 'shade', 
         title = 'Cluster Concordance', addgrid.col = 'black', mar = c(0, 0, 1, 0), 
         tl.cex = 0.75, tl.col = 'black', tl.srt = 45, cl.lim = c(0, 3))





### Using mutual information ###

# Load libraries
library(infotheo)
library(corrplot)

# Import cluster assignments
clin <- fread('./Data/Clinical.csv') %>%
  select(Subject, PASI_75) %>%
  unique()
df <- read.csv('./Results/ConsensusClusters/clusters_km.csv', row.names = 1) %>%
  mutate(PASI_75 = clin$PASI_75)

# Create p-value matrix
mat <- matrix(nrow = ncol(df), ncol = ncol(df), 
              dimnames = list(colnames(df), colnames(df)))
dimnames(mat) <- lapply(dimnames(mat), function(i) gsub('.Response', '', i))
for (i in 2:nrow(mat)) {
  for (j in 1:(i - 1)) {
    mat[i, j] <- mutinformation(df[[i]], df[[j]])
  }
}
mat <- natstobits(mat)

# Fill NAs
#mat[is.na(mat)] <- ceiling(max(mat, na.rm = TRUE))
mat[is.na(mat)] <- max(mat, na.rm = TRUE)

# Plot
corrplot(mat, is.corr = FALSE, type = 'lower', diag = FALSE, method = 'shade', 
         title = 'Cluster Concordance', addgrid.col = 'black', mar = c(0, 0, 1, 0), 
         tl.cex = 0.75, tl.col = 'black', tl.srt = 45, cl.lim = c(0, 1))




### TO DO:
# Rerun ConsensusClustering with maxK = 4
# Baseline numbers for all platforms and PASI75
# Mutual info with scale up to max possible MI given n and k


