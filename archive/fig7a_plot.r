# import complexheatmap library
library(ComplexHeatmap)
library(circlize)

# Read in membership
membership <- read.table("member3.tsv", header=TRUE, row.names=1)

member2 <- membership[membership$Name == 1,]
member1 <- membership[membership$Name == 2,]

# sort each row by the membership column in reverse order
member2 <- member2[order(-member2$membership1),]
member1 <- member1[order(-member1$membership2),]

row_names2 <- rownames(member2)
row_names1 <- rownames(member1)

row_names_all <- c(row_names2, row_names1)

# concatenate the two membership tables
membership <- rbind(member2, member1)

# Read in mmc1.csv
mmc1 <- read.csv("mmc1.csv", header=TRUE, row.names=1)
# change "-" in the row names to "."
rownames(mmc1) <- gsub("-", ".", rownames(mmc1))
# filter out the rows and order the rows in the order of row_names_all
mmc1 <- mmc1[row_names_all,]
# Get columns named 'Bailey', 'Collisson', and 'Moffitt'
mmc1 <- mmc1[,c("Bailey", "Collisson", "Moffitt")]

# Read in clinical table
clinical <- read.csv("clinical_table_140.tsv", header=TRUE, row.names=1, sep="\t")
# change "-" in the row names to "."
rownames(clinical) <- gsub("-", ".", rownames(clinical))
# filter out the rows and order the rows in the order of row_names_all
clinical <- clinical[row_names_all,]
# Get 'vital_status' column
clinical <- clinical[,c("vital_status")]

# Read in matrices
cnv <- read.table("fig7a_scna.txt", header=TRUE, row.names=1)
# filter out the columns and order the columns in the order of row_names_all
cnv <- cnv[,row_names_all]
# # Normalize each row
# cnv <- t(apply(cnv, 1, function(x) (x - mean(x)) / sd(x)))
print(dim(cnv))
# row_clustering_cnv <- hclust(dist(cnv[,1:52]))
# cnv <- cnv[row_clustering_cnv$order,]
# Get the number of negative values in each row from column 1 to 52
cnv_neg <- apply(cnv[,1:52], 1, function(x) sum(x < 0))
# sort the rows by the number of negative values
cnv <- cnv[order(cnv_neg),]
# Get the top 5 and bottom 5 rows
cnv <- rbind(cnv[1:5,], cnv[(nrow(cnv)-4):nrow(cnv),])
prot <- read.table("fig7a_proteomics.txt", header=TRUE, row.names=1)
# filter out the columns and order the columns in the order of row_names_all
prot <- prot[,row_names_all]
# # Normalize each row
# prot <- t(apply(prot, 1, function(x) (x - mean(x)) / sd(x)))
print(dim(prot))
# row_clustering_prot <- hclust(dist(prot[,1:52]))
# prot <- prot[row_clustering_prot$order,]
# Get the number of negative values in each row from column 1 to 52
prot_neg <- apply(prot[,1:52], 1, function(x) sum(x < 0))
# sort the rows by the number of negative values
prot <- prot[order(prot_neg),]
# Get the top 10 and bottom 10 rows
prot <- rbind(prot[1:10,], prot[(nrow(prot)-9):nrow(prot),])
rna <- read.table("fig7a_mRNA.txt", header=TRUE, row.names=1)
# filter out the columns and order the columns in the order of row_names_all
rna <- rna[,row_names_all]
# # Normalize each row
# rna <- t(apply(rna, 1, function(x) (x - mean(x)) / sd(x)))
print(dim(rna))
# row_clustering_rna <- hclust(dist(rna[,1:52]))
# rna <- rna[row_clustering_rna$order,]
# Get the number of negative values in each row from column 1 to 52
rna_neg <- apply(rna[,1:52], 1, function(x) sum(x < 0))
# sort the rows by the number of negative values
rna <- rna[order(rna_neg),]
# Get the top 100 and bottom 100 rows
rna <- rbind(rna[1:100,], rna[(nrow(rna)-99):nrow(rna),])
glyco <- read.table("fig7a_glycosylation.txt", header=TRUE, row.names=1)
# filter out the columns and order the columns in the order of row_names_all
glyco <- glyco[,row_names_all]
# # Normalize each row
# glyco <- t(apply(glyco, 1, function(x) (x - mean(x)) / sd(x)))
print(dim(glyco))
# row_clustering_glyco <- hclust(dist(glyco[,1:52]))
# glyco <- glyco[row_clustering_glyco$order,]
# Get the number of negative values in each row from column 1 to 52
glyco_neg <- apply(glyco[,1:52], 1, function(x) sum(x < 0))
# sort the rows by the number of negative values
glyco <- glyco[order(glyco_neg),]
# Get the top 15 and bottom 15 rows
glyco <- rbind(glyco[1:15,], glyco[(nrow(glyco)-14):nrow(glyco),])
phospho <- read.table("fig7a_phosphoproteomics.txt", header=TRUE, row.names=1)
# filter out the columns and order the columns in the order of row_names_all
phospho <- phospho[,row_names_all]
# # Normalize each row
# phospho <- t(apply(phospho, 1, function(x) (x - mean(x)) / sd(x)))
print(dim(phospho))
# row_clustering_phospho <- hclust(dist(phospho[,1:52]))
# phospho <- phospho[row_clustering_phospho$order,]
# Get the number of negative values in each row from column 1 to 52
phospho_neg <- apply(phospho[,1:52], 1, function(x) sum(x < 0))
# sort the rows by the number of negative values
phospho <- phospho[order(phospho_neg),]
# Get the top 15 and bottom 15 rows
phospho <- rbind(phospho[1:15,], phospho[(nrow(phospho)-14):nrow(phospho),])

# vector size of 105, filled with 1
group <- rep(1, 105)
# from index 53 to 105, fill with 2
group[53:105] <- 2

# make new column in membership named "membership_plot"
membership$membership_plot <- membership$membership1
# from index 53 to 105, fill with membership2
membership$membership_plot[53:105] <- membership$membership2[53:105]

# Start to draw the heatmap
# heatmap of group. 1 to blue, 2 to orange
# heatmap of membership ranging from 0 to 1. Gradient from white to blue
# heatmap of vital_status, brown for "Deceased", green for "Living"
# heatpmap of collisson, purple for "quasimesenchymal", orange for "classical", blue for "exocrine-like"
# heatmap of Bailey, blue for "squamous", orange for "immunogenic", sky blue for "pancreatic progenitor", brown for "ADEX"
# heatmap of Moffitt, red for "BASAL-LIKE", orange for "CLASSICAL"
hmap_clu <- HeatmapAnnotation(df = data.frame(group = group, 
                                membership = membership$membership_plot,
                                vital_status = clinical,
                                collisson = mmc1[,"Collisson"],
                                Bailey = mmc1[,"Bailey"],
                                Moffitt = mmc1[,"Moffitt"]),
                              col = list(group = c("1" = "blue", "2" = "orange"),
                                        membership = circlize::colorRamp2(c(0.4, 1), c("white", "blue")),
                                        vital_status = c("Deceased" = "brown", "Living" = "green"),
                                        collisson = c("quasimesenchymal" = "purple", "classical" = "orange", "exocrine-like" = "blue"),
                                        Bailey = c("squamous" = "blue", "immunogenic" = "orange", "pancreatic progenitor" = "sky blue", "ADEX" = "brown"),
                                        Moffitt = c("BASAL-LIKE" = "red", "CLASSICAL" = "orange")),
                              which = "column",
                              na_col = "gray",
                              show_legend = TRUE)

# # heatmap of cnv, proteomics, rna, glyco, phospho
# heatmap_cnv <- Heatmap(cnv, name = "CNV", col = circlize::colorRamp2(c(-5, 0, 5), c("blue", "white", "red")), show_row_names = FALSE, show_column_names = FALSE)
# heatmap_prot <- Heatmap(prot, name = "Proteomics", col = circlize::colorRamp2(c(-5, 0, 5), c("blue", "white", "red")), show_row_names = FALSE, show_column_names = FALSE)
# heatmap_rna <- Heatmap(rna, name = "RNA", col = circlize::colorRamp2(c(-5, 0, 5), c("blue", "white", "red")), show_row_names = FALSE, show_column_names = FALSE)
# heatmap_glyco <- Heatmap(glyco, name = "Glycosylation", col = circlize::colorRamp2(c(-5, 0, 5), c("blue", "white", "red")), show_row_names = FALSE, show_column_names = FALSE)
# heatmap_phospho <- Heatmap(phospho, name = "Phosphoproteomics", col = circlize::colorRamp2(c(-5, 0, 5), c("blue", "white", "red")), show_row_names = FALSE, show_column_names = FALSE)
# save figure
pdf("fig7a_heatmap_upper.pdf", width=10, height=20)
# combine the heatmaps
# heatmap_all = heatmap_cnv %v% heatmap_prot %v% heatmap_rna %v% heatmap_glyco %v% heatmap_phospho
# draw(heatmap_all, ht_gap = unit(1, "cm"))
draw(hmap_clu)
dev.off()