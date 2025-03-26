library(xCell)
library(dplyr)
library(data.table)
library(stringr)

# Load RNA matrix
rna_t = read.csv('mRNA_RSEM_UQ_log2_Tumor.cct', sep = '\t', row.names = 1)
colnames(rna_t) = paste0(colnames(rna_t), '_T')

# Load Metadata
meta = readRDS('results/figure3/meta.rds')
meta = meta[meta$cts_id %in% colnames(rna_t),]

# Kept only genes with zero counts less than 50%
rna_t_subset = rna_t[rowSums(rna_t==0)/ncol(rna_t) < 0.5,]

# xCell Deconvolution
results = xCellAnalysis(rna_t_subset)

results%>%apply(.,1, function(x) x < 1) %>% unlist %>% all

# filtered out immune cell types with zero readout in > 80% of samples
apply(results, 1, function(x){
sum(x<=1e-10)/length(x)
}) %>% sort

apply(results, 1, function(x){
sum(x <=1e-10)/length(x) <= 0.8
}) -> row_filter
row_filter%>%table


results_filtered = results[row_filter,]
saveRDS(results_filtered, 'results/figure6/xcell.rds')