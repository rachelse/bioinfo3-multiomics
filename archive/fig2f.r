# pacman::p_load(WGCNA)
# # allowWGCNAThreads()

# # Read in patient list
# patients <- read.table("105.list", header=FALSE)
# # replace "-" to "."
# patients$V1 <- gsub("-", ".", patients$V1)

# # Read in the gene list
# genes <- read.table("genes_exist.uniq", header=FALSE)

# # Read in cnv data
# cnv <- read.table("SCNA_log2_gene_level.cct", header=TRUE, sep="\t", row.names=1)
# # get columns of patients
# cnv <- cnv[, patients$V1]
# # sort rows by genes
# cnv <- cnv[match(genes$V1, rownames(cnv)), ]
# print(dim(cnv))

# # Read in mRNA data
# mRNA <- read.table("mRNA_RSEM_UQ_log2_Tumor.cct", header=TRUE)
# # get columns of patients
# mRNA <- mRNA[, patients$V1]
# # sort rows of genes
# mRNA <- mRNA[match(genes$V1, rownames(mRNA)), ]
# print(dim(mRNA))

# # Read in proteomics data
# proteomics <- read.table("proteomics_gene_level_MD_abundance_tumor.cct", header=TRUE)
# # get columns of patients
# proteomics <- proteomics[, patients$V1]
# # sort rows of genes
# proteomics <- proteomics[match(genes$V1, rownames(proteomics)), ]
# print(dim(proteomics))

# # assert the rownames of cnv and mRNA are the same
# stopifnot(all.equal(rownames(cnv), rownames(mRNA)))
# stopifnot(all.equal(rownames(cnv), rownames(proteomics)))

# mrna.vs.cna <- corAndPvalue(mRNA, cnv, method="spearman")
# pome.vs.cna <- corAndPvalue(proteomics, cnv, method="spearman")

# # save the results
# write.csv(mrna.vs.cna, "mrna-vs-cna.csv")
# write.csv(pome.vs.cna, "pome-vs-cna.csv")
pvalue <- 0.05

# Read in the matrix (tsv format)
mrna.vs.cna <- read.csv ('matrix_cnv_mrna_small.tsv', header=TRUE, row.names=1, sep="\t")
pome.vs.cna <- read.csv ('matrix_cnv_prot_small.tsv', header=TRUE, row.names=1, sep="\t")

# For each col, FDR correction by Benjamini-Hochberg
cat ('Adjusting p-values ...\n')
for (i in 1:nrow(mrna.vs.cna)) mrna.vs.cna[i,] <- p.adjust (mrna.vs.cna[i,], method='BH')
for (i in 1:nrow(pome.vs.cna)) pome.vs.cna[i,] <- p.adjust (pome.vs.cna[i,], method='BH')

# check if the p-values are less than 0.05
cna_mrna <- mrna.vs.cna < pvalue
cna_pome <- pome.vs.cna < pvalue

# save the results
write.table (cna_mrna, 'cna_mrna_small.tsv', sep="\t", quote=FALSE)
write.table (cna_pome, 'cna_pome_small.tsv', sep="\t", quote=FALSE)

# From the table, get the indices where the diagonal value is 1
mrna_diag <- diag (cna_mrna)
pome_diag <- diag (cna_pome)

# Save the results
# row names are row names of cna_mrna
write.table (mrna_diag, 'mrna_diag_small.tsv', sep="\t", quote=FALSE, col.names=FALSE)
write.table (pome_diag, 'pome_diag_small.tsv', sep="\t", quote=FALSE, col.names=FALSE)