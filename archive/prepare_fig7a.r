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
# Create "Description" column as the first column
scna <- cbind(Description=rownames(scna), scna)
print(dim(scna))
# save the scna data without first row
write.table(scna, "panoply_scna.gct", sep="\t", quote=FALSE)

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
# Create "Description" column as the first column
proteomics <- cbind(Description=rownames(proteomics), proteomics)
print(dim(proteomics))
# save the proteomics data
write.table(proteomics, "panoply_proteomics.gct", sep="\t", quote=FALSE)

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
# Create "Description" column as the first column
mRNA <- cbind(Description=rownames(mRNA), mRNA)
print(dim(mRNA))
# save the mRNA data
write.table(mRNA, "panoply_mRNA.gct", sep="\t", quote=FALSE)

# Read in glycosylation data
glycosylation <- read.table("N-glycoproteomics_Site_level_ratio_tumor.cct", header=TRUE, row.names=1)
# only keep the columns with the  patients in the list
glycosylation <- glycosylation[,patients$V1]
# sort the column in the order of patients$V1
glycosylation <- glycosylation[,order(colnames(glycosylation))]
# If there are any non-numeric values or NA, remove the row
glycosylation <- glycosylation[apply(glycosylation, 1, function(x) all(is.finite(x))),]
# Create "Description" column as the first column
glycosylation <- cbind(Description=rownames(glycosylation), glycosylation)
print(dim(glycosylation))
# save the glycosylation data
write.table(glycosylation, "panoply_glycosylation.gct", sep="\t", quote=FALSE)

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
# Create "Description" column as the first column
phosphoproteomics <- cbind(Description=rownames(phosphoproteomics), phosphoproteomics)
print(dim(phosphoproteomics))
# save the phosphoproteomics data
write.table(phosphoproteomics, "panoply_phosphoproteomics.gct", sep="\t", quote=FALSE)