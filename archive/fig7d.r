# Read in the C1_classical.list and C2_classical.list
c1_list <- read.table("C1_classical.list", header = FALSE)
c2_list <- read.table("C2_classical.list", header = FALSE)

# concatenate the two lists
patient_list <- rbind(c1_list, c2_list)

# Read in the clinical data
clinical <- read.table("clinical_table_survival.tsv", header = TRUE, sep = "\t")
# From case_id, replace "-" with "."
clinical$case_id <- gsub("-", ".", clinical$case_id)
# Only get the rows that are in the patient list
clinical <- clinical[clinical$case_id %in% patient_list$V1, ]
# Reset the row index
rownames(clinical) <- 1:nrow(clinical)
# From vital_status, change Deceased to 1 and Living to 0
clinical$vital_status <- as.numeric(clinical$vital_status == "Deceased")

# two clusters: one in c1_list, the other not in c1_list
clinical$c1_list <- as.numeric(clinical$case_id %in% c1_list$V1)

library(survival)
library(ggplot2)
library(ggfortify)
library(survminer)
# Plot Kaplan-Meier curve
fit <- survfit(Surv(follow_up_days, vital_status) ~ c1_list, data = clinical)
summary_fit <- summary(fit)
pval <- surv_pvalue(fit)
print(pval)

# Plot the Kaplan-Meier curve with p-value, white background, no standard error
autoplot(fit, conf.int = FLASE) + theme_minimal() + theme(legend.position = "none") + ggtitle(paste("p-value: ", pval$pval, sep = "")) + 
theme(panel.background = element_rect(fill = "white")) + theme(axis.line = element_line(colour = "black")) + theme(axis.text = element_text(size = 12)) + theme(axis.title = element_text(size = 14))

# Save the plot
ggsave("fig7d.png", width = 6, height = 6)