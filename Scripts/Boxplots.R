# Load libraries
library(tidyverse)
library(ggsci)

# Import data
df_cat <- read_csv('./Results/ErrCat_RFE.csv')
df_cat <- df_cat[, grepl('Loss', colnames(df_cat))]
colnames(df_cat) <- gsub('_Loss', '', colnames(df_cat))
df_cont <- read_csv('./Results/ErrCont_RFE.csv')
df_cont <- df_cont[, grepl('Loss', colnames(df_cont))]
colnames(df_cont) <- gsub('_Loss', '', colnames(df_cont))
brks <- c('Clinical', 'Blood_Proteomics', 'Blood_miRNA',
          'Blood_mRNA', 'Lesional_mRNA', 'Nonlesional_mRNA')

# Plotting function
boxes <- function(dat, resp) {
  
  # Melt
  df <- gather(dat, Data, Loss) %>%
    mutate(Data = factor(Data, levels = brks))
  
  # Plot
  if (resp == 'Categorical') {
    ylab <- 'Cross Entropy'
  } else {
    ylab <- 'RMSE'
  }
  ggplot(df, aes(Data, Loss, fill = Data)) +
    geom_boxplot() + 
    labs(title = paste('Predictive Error:\n', resp, 'Response'),
         y = ylab) + 
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(breaks = brks,
                      values = c(pal_npg()(6)[c(2, 1, 3:6)]))
  
  # Export
  ggsave(paste0('./Results/Figures/PredictiveError_', resp, '.pdf'))
  
}

# Run
boxes(df_cat, 'Categorical')
boxes(df_cont, 'Continuous')



