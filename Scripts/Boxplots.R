# Load libraries
library(tidyverse)
library(ggsci)

# Import data
df_cat <- read_csv('./Results/err_cat.csv')
df_cont <- read_csv('./Results/err_cont.csv')
brks <- c('Clinical', 'Blood_Proteomics', 'Blood_miRNA',
          'Blood_mRNA', 'Lesional_mRNA', 'Nonlesional_mRNA')

### CATEGORICAL ###
df <- gather(df_cat, Data, Loss) %>%
  mutate(Data = factor(Data, levels = brks))

# Build plot
ggplot(df, aes(Data, Loss, fill = Data)) +
  geom_boxplot() + 
  labs(title = 'Test Error:\nCategorical Response',
       y = 'Cross Entropy') + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(breaks = brks,
                     values = c(pal_npg()(6)[c(2, 1, 3:6)]))

# These distributions have some extreme outliers...let's winsorise
winsorise <- function(x, multiple = 3) {
  y <- x - median(x)
  lim <- mad(y, center = 0) * multiple
  y[y > lim] <- lim
  y[y < -lim] <- -lim
  y <- y + median(x)
  return(y)
}
df <- df %>% 
  group_by(Data) %>%
  mutate(Loss = winsorise(Loss)) %>%
  ungroup()
ggplot(df, aes(Data, Loss, fill = Data)) +
  geom_boxplot() + 
  labs(title = 'Test Error:\nCategorical Response',
       y = 'Cross Entropy') + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(breaks = brks,
                    values = c(pal_npg()(6)[c(2, 1, 3:6)]))
# Much better!


### CONTINUOUS ###
df <- gather(df_cont, Data, Loss) %>%
  mutate(Data = factor(Data, levels = brks))

# Build plot
ggplot(df, aes(Data, Loss, fill = Data)) +
  geom_boxplot() + 
  labs(title = 'Test Error:\nContinuous Response',
       y = 'RMSE') + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(breaks = brks,
                    values = c(pal_npg()(6)[c(2, 1, 3:6)]))



