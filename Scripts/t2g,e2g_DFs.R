##########################################################################################
### This script generates the t2g and e2g data frames that are respectively saved as #####
### Ensembl.Hs79.Tx.csv.zip and Ensembl.Hs79.GeneSymbols.csv.zip in the Data directory ###
##########################################################################################

# Load libraries
library(AnnotationHub)
library(ensembldb)

# Load genome build GRCh38, v.79
ah <- AnnotationHub()
hs79 <- query(ah, c('ensembl', 'gtf', '79', 'homo sapiens'))
gr <- hs79[[names(hs79)]]

# Create ensembledb
db <- ensDbFromGRanges(gr, organism='Homo_sapiens', version=79, genomeVersion='GRCh38')
edb <- EnsDb(db)

# Create, export t2g df
t2g <- transcripts(edb, columns=c('tx_id', 'gene_id'), return.type='data.frame')
write.csv(t2g, 'Ensembl.Hs79.Tx.csv', row.names=FALSE)

# Create, export e2g df
e2g <- transcripts(edb, columns=c('gene_id', 'gene_name'), return.type='data.frame') %>% unique(.)
write.csv(e2g, 'Ensembl.Hs79.GeneSymbols.csv', row.names=FALSE)