# Load libraries
library(data.table)
library(dplyr)

# Combo function
combo <- function(tissue) {
  
  # Load baseline
  wk00 <- fread(paste0('./Results/Response/', tissue, '.wk00.Response.txt')) %>%
            rename(p.value_wk00 = p.value,
                   logFC_wk00   = logFC) %>%
            select(EnsemblID, GeneSymbol, logFC_wk00, p.value_wk00)
  
  # Load deltas
  delta01 <- fread(paste0('./Results/Time/', tissue, '.Delta01.txt')) %>%
    rename(p.value_delta01 = p.value,
           logFC_delta01   = logFC) %>%
    select(EnsemblID, logFC_delta01, p.value_delta01)
  delta12 <- fread(paste0('./Results/Time/', tissue, '.Delta12.txt')) %>%
    rename(p.value_delta12 = p.value,
           logFC_delta12   = logFC) %>%
    select(EnsemblID, logFC_delta12, p.value_delta12)
  
  # Combos
  mix01 <- inner_join(wk00, delta01, by = 'EnsemblID') %>%
    rowwise() %>%
    mutate(p.value_combo = pchisq(-2 * sum(log(c(p.value_delta01, 
                                                 p.value_wk00))),
                                  df = 4, lower = FALSE)) 
  mix01 %>%
    mutate(q.value = p.adjust(p.value_combo, method = 'BH')) %>%
    arrange(p.value_combo) %>%
    fwrite(paste0('./Results/Combos/',
                  tissue, '.Combo_Baseline.Delta01.txt'), sep = '\t')
  mix12 <- inner_join(wk00, delta12, by = 'EnsemblID') %>%
    rowwise() %>%
    mutate(p.value_combo = pchisq(-2 * sum(log(c(p.value_delta12, 
                                                 p.value_wk00))),
                                  df = 4, lower = FALSE))
  mix12 %>%
    mutate(q.value = p.adjust(p.value_combo, method = 'BH')) %>%
    arrange(p.value_combo) %>%
    fwrite(paste0('./Results/Combos/', 
                  tissue, '.Combo_Baseline.Delta12.txt'), sep = '\t')
  
}

# Run
for (tissue in c('Blood', 'Lesional', 'Nonlesional')) combo(tissue)


