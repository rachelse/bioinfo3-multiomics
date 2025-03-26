library(NMF)
library(doParallel)

# Read in list of patients
patients <- read.table("105.list", header=FALSE)
# substitute "-" to "." in the patient names
patients$V1 <- gsub("-", ".", patients$V1)
# print(patients[1:5,])

# Read in SCNA data
scna <- read.table("SCNA_log2_gene_level.cct", header=TRUE, row.names=1, sep="\t")
# only keep the columns with the  patients in the list
scna <- scna[,patients$V1]
# sort the column in the order of patients$V1
scna <- scna[,order(colnames(scna))]
# If there are any non-numeric values or NA, remove the row
scna <- scna[apply(scna, 1, function(x) all(is.finite(x))),]
# subtract by the median
scna_median <- apply(scna, 1, median)
scna <- scna - scna_median
# # check if there is any negative value
# print(dim(scna))
# print("Read in SCNA data")
# number of rows
scna_rows <- nrow(scna)

# Read in proteomics data, nonnegative
proteomics <- read.table("proteomics_gene_level_MD_abundance_tumor.cct", header=TRUE, row.names=1)
# only keep the columns with the  patients in the list
proteomics <- proteomics[,patients$V1]
# sort the column in the order of patients$V1
proteomics <- proteomics[,order(colnames(proteomics))]
# Filter out rows that only contain non-infinite numeric values and no NA
proteomics <- proteomics[apply(proteomics, 1, function(x) all(is.finite(x))),]
# divide by the median of each row
prot_median <- apply(proteomics, 1, median)
proteomics <- proteomics / prot_median
# log2
proteomics <- log2(proteomics)
# Check if there is any negative value
# print(dim(proteomics))
# print("Read in proteomics data")
# number of rows
prot_rows <- nrow(proteomics)

# Read in mRNA data, nonnegative
mRNA <- read.table("mRNA_RSEM_UQ_log2_Tumor.cct", header=TRUE, row.names=1)
# only keep the columns with the  patients in the list
mRNA <- mRNA[,patients$V1]
# sort the column in the order of patients$V1
mRNA <- mRNA[,order(colnames(mRNA))]
# If there are any non-numeric values or NA, remove the row
mRNA <- mRNA[apply(mRNA, 1, function(x) all(is.finite(x))),]
# Also remove the rows if there is any 0
mRNA <- mRNA[rowSums(mRNA == 0) == 0,]
# subtract with median of each row
mrna_median <- apply(mRNA, 1, median)
mRNA <- mRNA - mrna_median
# check if there is any negative value
# print(dim(mRNA))
# # print(mRNA[1:5,1:5])
# print("Read in mRNA data")
# number of rows
mRNA_rows <- nrow(mRNA)

# Read in glycosylation data
glycosylation <- read.table("N-glycoproteomics_Site_level_ratio_tumor.cct", header=TRUE, row.names=1)
# only keep the columns with the  patients in the list
glycosylation <- glycosylation[,patients$V1]
# sort the column in the order of patients$V1
glycosylation <- glycosylation[,order(colnames(glycosylation))]
# If there are any non-numeric values or NA, remove the row
glycosylation <- glycosylation[apply(glycosylation, 1, function(x) all(is.finite(x))),]
# # divide by the median of each row
# glyc_median <- apply(glycosylation, 1, median)
# glycosylation <- glycosylation - glyc_median
# Check if there is any negative value
# print(dim(glycosylation))
# print("Read in glycosylation data")
# number of rows
glyc_rows <- nrow(glycosylation)

# Read in phosphoproteomics data, nonnegative
phosphoproteomics <- read.table("phosphoproteomics_site_level_MD_abundance_tumor.cct", header=TRUE, row.names=1)
# phosphoproteomics <- read.table("phosphoproteomics_MultiSite_level_MD_abundance_tumor.cct", header=TRUE, row.names=1, sep="\t")
# only keep the columns with the  patients in the list
phosphoproteomics <- phosphoproteomics[,patients$V1]
# sort the column in the order of patients$V1
phosphoproteomics <- phosphoproteomics[,order(colnames(phosphoproteomics))]
# If there are any non-numeric values or NA, remove the row
phosphoproteomics <- phosphoproteomics[apply(phosphoproteomics, 1, function(x) all(is.finite(x))),]
# divide by the median of each row
phos_median <- apply(phosphoproteomics, 1, median)
phosphoproteomics <- phosphoproteomics / phos_median
# log2
phosphoproteomics <- log2(phosphoproteomics)
# Check if there is any negative value
# print(dim(phosphoproteomics))
# print("Read in phosphoproteomics data")
# number of rows
phos_rows <- nrow(phosphoproteomics)

# If there is any warning, see warnings
if (length(warnings()) > 0) {
  print(warnings())
}

# Combine all data above columnwise
data <- rbind(scna, mRNA, proteomics, glycosylation, phosphoproteomics)
print(scna_rows + mRNA_rows + prot_rows + glyc_rows + phos_rows)
print(dim(data))

# Remove bottom 5% row of the data with the lowest std
feat_sd <- apply(data, 1, sd)
# print("Feature standard deviation")
# print(dim(feat_sd))
sd_cutoff <- quantile(feat_sd, 0.05)
idx.keep <- which(feat_sd > sd_cutoff)
data <- data[idx.keep,]
# print(dim(data))
# print("Combined all data")

# Recalculate the number of rows based on the idx.keep
scna_rows_new <- sum(idx.keep <= scna_rows)
mRNA_rows_new <- sum(idx.keep > scna_rows & idx.keep <= scna_rows + mRNA_rows)
prot_rows_new <- sum(idx.keep > scna_rows + mRNA_rows & idx.keep <= scna_rows + mRNA_rows + prot_rows)
glyc_rows_new <- sum(idx.keep > scna_rows + mRNA_rows + prot_rows & idx.keep <= scna_rows + mRNA_rows + prot_rows + glyc_rows)
phos_rows_new <- sum(idx.keep > scna_rows + mRNA_rows + prot_rows + glyc_rows)
# check if the sum of rows is the same as the number of rows in data
if (scna_rows_new + mRNA_rows_new + prot_rows_new + glyc_rows_new + phos_rows_new != nrow(data)) {
  print("Error: number of rows does not match")
  print(scna_rows_new + mRNA_rows_new + prot_rows_new + glyc_rows_new + phos_rows_new)
  print(nrow(data))
  quit()
}

# save idx.keep
write.table(idx.keep, "idx.keep.txt", sep="\t", quote=FALSE)

# z-score for each row
# data <- t(apply(data, 1, function(x) (x - mean(x)) / sd(x)))
# column wise
data <- apply(data, 2, function(x) (x - mean(x)) / sd(x))
# save data
write.table(data, "data.txt", sep="\t", quote=FALSE)
# save data split by type
write.table(data[1:scna_rows_new,], "fig7a_scna.txt", sep="\t", quote=FALSE)
write.table(data[(scna_rows_new + 1):(scna_rows_new + mRNA_rows_new),], "fig7a_mRNA.txt", sep="\t", quote=FALSE)
write.table(data[(scna_rows_new + mRNA_rows_new + 1):(scna_rows_new + mRNA_rows_new + prot_rows_new),], "fig7a_proteomics.txt", sep="\t", quote=FALSE)
write.table(data[(scna_rows_new + mRNA_rows_new + prot_rows_new + 1):(scna_rows_new + mRNA_rows_new + prot_rows_new + glyc_rows_new),], "fig7a_glycosylation.txt", sep="\t", quote=FALSE)
write.table(data[(scna_rows_new + mRNA_rows_new + prot_rows_new + glyc_rows_new + 1):(scna_rows_new + mRNA_rows_new + prot_rows_new + glyc_rows_new + phos_rows_new),], "fig7a_phosphoproteomics.txt", sep="\t", quote=FALSE)
data <- t(data)

# Create one data matrix with all negative numbers zeroed and all the rest of the data retained
positive_only_matrix <- data
for (i in 1:ncol(positive_only_matrix)) {
  positive_only_matrix[, i] <- ifelse(positive_only_matrix[, i] < 0, 0, positive_only_matrix[, i])
}

# Step 2: Create another data matrix with all positive numbers zeroed and signs of all negative numbers removed
negative_only_matrix <- data
for (i in 1:ncol(negative_only_matrix)) {
  negative_only_matrix[, i] <- ifelse(negative_only_matrix[, i] > 0, 0, abs(negative_only_matrix[, i]))
}

# Step 3: Concatenate both matrices horizontally
combined_matrix <- cbind(positive_only_matrix, negative_only_matrix)
# transpose
combined_matrix <- t(combined_matrix)

# Remove all the rows with 0
combined_matrix <- combined_matrix[rowSums(combined_matrix) != 0,]
print(dim(combined_matrix))
print("Non-negative matrix")

# NMF
# Run NMF, 500 iterations
res <- nmf(as.matrix(combined_matrix), 2, .options=list(verbose = TRUE), maxIter=500)

W <- res@fit@W # basis matrix (samples x components)
H <- res@fit@H # coefficient matrix (components x features)

# Save the results
write.table(W, "W.txt", sep="\t", quote=FALSE)
write.table(H, "H.txt", sep="\t", quote=FALSE)
