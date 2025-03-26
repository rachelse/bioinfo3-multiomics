# Read in the patient list
patient_list <- read.table("bioinfo3/kras_vaf_longlived.list", header = FALSE)

# Read in the clinical data
clinical <- read.table("bioinfo3/clinical_table_survival.tsv", header = TRUE, sep = "\t")
# Only get the rows that are in the patient list
clinical <- clinical[clinical$case_id %in% patient_list$V1, ]
# Reset the row index
rownames(clinical) <- 1:nrow(clinical)
# Exchange "-" with "." in the case_id
clinical$case_id <- gsub("-", ".", clinical$case_id)
# From vital_status, change Deceased to 1 and Living to 0
clinical$vital_status <- as.numeric(clinical$vital_status == "Deceased")

# Replace "-" with "."
patient_list$V1 <- gsub("-", ".", patient_list$V1)

# Read in the data
proteomics <- read.table("bioinfo3/N-glycoproteomics_Site_level_ratio_tumor_cox.cct", header=TRUE)
# Combine the row names and the "Gene" column with ":"
rownames(proteomics) <- paste(rownames(proteomics), proteomics$Gene, sep = ":")
# Only get the columns that are in the patient list
proteomics <- proteomics[, colnames(proteomics) %in% patient_list$V1]

# Get row ends with "N2H9F0S0G0:APOD"
proteomics <- proteomics[grep("N2H9F0S0G0:APOD", rownames(proteomics)), ]
# Get the second row
proteomics <- proteomics[2, ]

# Transpose
proteomics <- t(proteomics)
# Change column name to "N2H9F0S0G0_APOD"
colnames(proteomics) <- "N2H9F0S0G0_APOD"
# Combine the clinical data with the proteomics data
proteomics <- cbind(proteomics, clinical[match(rownames(proteomics), clinical$case_id), ])
# Remove if there is any NA
proteomics <- proteomics[complete.cases(proteomics), ]
# Get the median of "N2H9F0S0G0_APOD", and divide the data into two groups
proteomics$N2H9F0S0G0_APOD <- as.numeric(proteomics$N2H9F0S0G0_APOD > median(proteomics$N2H9F0S0G0_APOD))
print(dim(proteomics))

library(survival)
library(ggplot2)
library(ggfortify)
library(survminer)
# Plot Kaplan-Meier curve
fit <- survfit(Surv(follow_up_days, vital_status) ~ N2H9F0S0G0_APOD, data = proteomics)
summary_fit <- summary(fit)
pval <- surv_pvalue(fit)
print(pval)

# Plot the Kaplan-Meier curve with p-value, white background, without CI
autoplot(fit, conf.int = FLASE) + theme_minimal() + theme(legend.position = "none") + ggtitle(paste("p-value: ", pval$pval, sep = "")) + 
theme(panel.background = element_rect(fill = "white")) + theme(axis.line = element_line(colour = "black")) + theme(axis.text = element_text(size = 12)) + theme(axis.title = element_text(size = 14))

# Save the plot
ggsave("fig3f_1.png", width = 6, height = 6)