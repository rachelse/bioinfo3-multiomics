# Get C1.list and C2.list
c1_list <- read.table("C1.list", header = FALSE)
c2_list <- read.table("C2.list", header = FALSE)

# combine c1_list and c2_list's row_names
row_names_all <- c(c1_list$V1, c2_list$V1)

# Read in mmc1.csv
mmc1 <- read.csv("mmc1.csv", header=TRUE, row.names=1)
# change "-" in the row names to "."
rownames(mmc1) <- gsub("-", ".", rownames(mmc1))
# filter out the rows and order the rows in the order of row_names_all
mmc1 <- mmc1[row_names_all,]
# Get columns named 'Bailey', 'Collisson', and 'Moffitt'
mmc1 <- mmc1[,c("Bailey", "Collisson", "Moffitt")]
print(dim(mmc1))

# Read in clinical table
clinical <- read.csv("clinical_table_140.tsv", header=TRUE, row.names=1, sep="\t")
# change "-" in the row names to "."
rownames(clinical) <- gsub("-", ".", rownames(clinical))
# filter out the rows and order the rows in the order of row_names_all
clinical <- clinical[row_names_all,]
# # Get 'vital_status' column
# clinical <- clinical[,c("vital_status")]

# make the contingency table based on C2 and collisson=="quasimesenchymal" and do the fisher's exact test
# make 2X2 contingency table
contingency_table <- matrix(0, nrow=2, ncol=2)
# iterate through row_names_all
for (name in row_names_all) {
  # if the row is in c2_list
  if (name %in% c2_list$V1) {
    # if the row is in collisson and the value is "quasimesenchymal"
    if (!is.na(mmc1[name, "Collisson"]) && mmc1[name, "Collisson"] == "quasimesenchymal") {
        # increase the value of the first row and the first column
        contingency_table[1, 1] <- contingency_table[1, 1] + 1
    } else {
        # increase the value of the first row and the second column
        contingency_table[1, 2] <- contingency_table[1, 2] + 1
    }
  } else {
    # if the row is in collisson and the value is "quasimesenchymal"
    if (!is.na(mmc1[name, "Collisson"]) && mmc1[name, "Collisson"] == "quasimesenchymal") {
        # increase the value of the second row and the first column
        contingency_table[2, 1] <- contingency_table[2, 1] + 1
    } else {
        # increase the value of the second row and the second column
        contingency_table[2, 2] <- contingency_table[2, 2] + 1
    }
  }
}
# do fisher's exact test
fisher_result <- fisher.test(contingency_table)
# get the p-value
pval_quasi <- fisher_result$p.value
log_pval_quasi <- -log10(pval_quasi)

# do same with the C2 and squamous
contingency_table <- matrix(0, nrow=2, ncol=2)
# iterate through row_names_all
for (name in row_names_all) {
  # if the row is in c2_list
  if (name %in% c2_list$V1) {
    # if the row is in collisson and the value is "squamous"
    if (!is.na(mmc1[name, "Bailey"]) && mmc1[name, "Bailey"] == "squamous") {
        # increase the value of the first row and the first column
        contingency_table[1, 1] <- contingency_table[1, 1] + 1
    } else {
        # increase the value of the first row and the second column
        contingency_table[1, 2] <- contingency_table[1, 2] + 1
    }
  } else {
    # if the row is in collisson and the value is "squamous"
    if (!is.na(mmc1[name, "Bailey"]) && mmc1[name, "Bailey"] == "squamous") {
        # increase the value of the second row and the first column
        contingency_table[2, 1] <- contingency_table[2, 1] + 1
    } else {
        # increase the value of the second row and the second column
        contingency_table[2, 2] <- contingency_table[2, 2] + 1
    }
  }
}
# do fisher's exact test
fisher_result <- fisher.test(contingency_table)
# get the p-value
pval_squamous <- fisher_result$p.value
log_pval_squamous <- -log10(pval_squamous)

# do same with the C1 and moffitt classical
contingency_table <- matrix(0, nrow=2, ncol=2)
# iterate through row_names_all
for (name in row_names_all) {
  # if the row is in c1_list
  if (name %in% c1_list$V1) {
    # if the row is in moffitt and the value is "CLASSICAL"
    if (!is.na(mmc1[name, "Moffitt"]) && mmc1[name, "Moffitt"] == "CLASSICAL") {
        # increase the value of the first row and the first column
        contingency_table[1, 1] <- contingency_table[1, 1] + 1
    } else {
        # increase the value of the first row and the second column
        contingency_table[1, 2] <- contingency_table[1, 2] + 1
    }
  } else {
    # if the row is in moffitt and the value is "CLASSICAL"
    if (!is.na(mmc1[name, "Moffitt"]) && mmc1[name, "Moffitt"] == "CLASSICAL") {
        # increase the value of the second row and the first column
        contingency_table[2, 1] <- contingency_table[2, 1] + 1
    } else {
        # increase the value of the second row and the second column
        contingency_table[2, 2] <- contingency_table[2, 2] + 1
    }
  }
}
# do fisher's exact test
fisher_result <- fisher.test(contingency_table)
# get the p-value
pval_classical <- fisher_result$p.value
log_pval_classical <- -log10(pval_classical)

# do same with the C2 and moffitt basal-like
contingency_table <- matrix(0, nrow=2, ncol=2)
# iterate through row_names_all
for (name in row_names_all) {
  # if the row is in c2_list
  if (name %in% c2_list$V1) {
    # if the row is in moffitt and the value is "BASAL"
    if (!is.na(mmc1[name, "Moffitt"]) && mmc1[name, "Moffitt"] == "BASAL-LIKE") {
        # increase the value of the first row and the first column
        contingency_table[1, 1] <- contingency_table[1, 1] + 1
    } else {
        # increase the value of the first row and the second column
        contingency_table[1, 2] <- contingency_table[1, 2] + 1
    }
  } else {
    # if the row is in moffitt and the value is "BASAL"
    if (!is.na(mmc1[name, "Moffitt"]) && mmc1[name, "Moffitt"] == "BASAL-LIKE") {
        # increase the value of the second row and the first column
        contingency_table[2, 1] <- contingency_table[2, 1] + 1
    } else {
        # increase the value of the second row and the second column
        contingency_table[2, 2] <- contingency_table[2, 2] + 1
    }
  }
}
# do fisher's exact test
fisher_result <- fisher.test(contingency_table)
# get the p-value
pval_basal <- fisher_result$p.value
log_pval_basal <- -log10(pval_basal)

# do same with the C1 and Bailey pancreatic progenitor
contingency_table <- matrix(0, nrow=2, ncol=2)
# iterate through row_names_all
for (name in row_names_all) {
  # if the row is in c1_list
  if (name %in% c1_list$V1) {
    # if the row is in bailey and the value is "pancreatic progenitor"
    if (!is.na(mmc1[name, "Bailey"]) && mmc1[name, "Bailey"] == "pancreatic progenitor") {
        # increase the value of the first row and the first column
        contingency_table[1, 1] <- contingency_table[1, 1] + 1
    } else {
        # increase the value of the first row and the second column
        contingency_table[1, 2] <- contingency_table[1, 2] + 1
    }
  } else {
    # if the row is in bailey and the value is "pancreatic progenitor"
    if (!is.na(mmc1[name, "Bailey"]) && mmc1[name, "Bailey"] == "pancreatic progenitor") {
        # increase the value of the second row and the first column
        contingency_table[2, 1] <- contingency_table[2, 1] + 1
    } else {
        # increase the value of the second row and the second column
        contingency_table[2, 2] <- contingency_table[2, 2] + 1
    }
  }
}
# do fisher's exact test
fisher_result <- fisher.test(contingency_table)
# get the p-value
pval_progenitor <- fisher_result$p.value
log_pval_progenitor <- -log10(pval_progenitor)

# do same with the C1 and Living
contingency_table <- matrix(0, nrow=2, ncol=2)
# iterate through row_names_all
for (name in row_names_all) {
  # if the row is in c1_list
  if (name %in% c1_list$V1) {
    # if the row is in clinical and the value is "Dead"
    if (!is.na(clinical[name, "vital_status"]) && clinical[name, "vital_status"] == "Living") {
        # increase the value of the first row and the first column
        contingency_table[1, 1] <- contingency_table[1, 1] + 1
    } else {
        # increase the value of the first row and the second column
        contingency_table[1, 2] <- contingency_table[1, 2] + 1
    }
  } else {
    # if the row is in clinical and the value is "Dead"
    if (!is.na(clinical[name, "vital_status"]) && clinical[name, "vital_status"] == "Living") {
        # increase the value of the second row and the first column
        contingency_table[2, 1] <- contingency_table[2, 1] + 1
    } else {
        # increase the value of the second row and the second column
        contingency_table[2, 2] <- contingency_table[2, 2] + 1
    }
  }
}
# do fisher's exact test
fisher_result <- fisher.test(contingency_table)
# get the p-value
pval_death <- fisher_result$p.value
log_pval_death <- -log10(pval_death)

# do same with the C2 and Deceased
contingency_table <- matrix(0, nrow=2, ncol=2)
# iterate through row_names_all
for (name in row_names_all) {
  # if the row is in c2_list
  if (name %in% c2_list$V1) {
    # if the row is in clinical and the value is "Dead"
    if (!is.na(clinical[name, "vital_status"]) && clinical[name, "vital_status"] == "Deceased") {
        # increase the value of the first row and the first column
        contingency_table[1, 1] <- contingency_table[1, 1] + 1
    } else {
        # increase the value of the first row and the second column
        contingency_table[1, 2] <- contingency_table[1, 2] + 1
    }
  } else {
    # if the row is in clinical and the value is "Dead"
    if (!is.na(clinical[name, "vital_status"]) && clinical[name, "vital_status"] == "Deceased") {
        # increase the value of the second row and the first column
        contingency_table[2, 1] <- contingency_table[2, 1] + 1
    } else {
        # increase the value of the second row and the second column
        contingency_table[2, 2] <- contingency_table[2, 2] + 1
    }
  }
}
# do fisher's exact test
fisher_result <- fisher.test(contingency_table)
# get the p-value
pval_death_c2 <- fisher_result$p.value
log_pval_death_c2 <- -log10(pval_death_c2)

# do same with the C1 and Collison exocrine-like
contingency_table <- matrix(0, nrow=2, ncol=2)
# iterate through row_names_all
for (name in row_names_all) {
  # if the row is in c1_list
  if (name %in% c1_list$V1) {
    # if the row is in collisson and the value is "exocrine-like"
    if (!is.na(mmc1[name, "Collisson"]) && mmc1[name, "Collisson"] == "exocrine-like") {
        # increase the value of the first row and the first column
        contingency_table[1, 1] <- contingency_table[1, 1] + 1
    } else {
        # increase the value of the first row and the second column
        contingency_table[1, 2] <- contingency_table[1, 2] + 1
    }
  } else {
    # if the row is in collisson and the value is "exocrine-like"
    if (!is.na(mmc1[name, "Collisson"]) && mmc1[name, "Collisson"] == "exocrine-like") {
        # increase the value of the second row and the first column
        contingency_table[2, 1] <- contingency_table[2, 1] + 1
    } else {
        # increase the value of the second row and the second column
        contingency_table[2, 2] <- contingency_table[2, 2] + 1
    }
  }
}
# do fisher's exact test
fisher_result <- fisher.test(contingency_table)
# get the p-value
pval_exocrine <- fisher_result$p.value
log_pval_exocrine <- -log10(pval_exocrine)

# do same with the C1 and Bailey ADEX
contingency_table <- matrix(0, nrow=2, ncol=2)
# iterate through row_names_all
for (name in row_names_all) {
  # if the row is in c1_list
  if (name %in% c1_list$V1) {
    # if the row is in bailey and the value is "ADEX"
    if (!is.na(mmc1[name, "Bailey"]) && mmc1[name, "Bailey"] == "ADEX") {
        # increase the value of the first row and the first column
        contingency_table[1, 1] <- contingency_table[1, 1] + 1
    } else {
        # increase the value of the first row and the second column
        contingency_table[1, 2] <- contingency_table[1, 2] + 1
    }
  } else {
    # if the row is in bailey and the value is "ADEX"
    if (!is.na(mmc1[name, "Bailey"]) && mmc1[name, "Bailey"] == "ADEX") {
        # increase the value of the second row and the first column
        contingency_table[2, 1] <- contingency_table[2, 1] + 1
    } else {
        # increase the value of the second row and the second column
        contingency_table[2, 2] <- contingency_table[2, 2] + 1
    }
  }
}
# do fisher's exact test
fisher_result <- fisher.test(contingency_table)
# get the p-value
pval_adex <- fisher_result$p.value
log_pval_adex <- -log10(pval_adex)

# do same with the C2 and Bailey immunogenic
contingency_table <- matrix(0, nrow=2, ncol=2)
# iterate through row_names_all
for (name in row_names_all) {
  # if the row is in c2_list
  if (name %in% c2_list$V1) {
    # if the row is in bailey and the value is "immunogenic"
    if (!is.na(mmc1[name, "Bailey"]) && mmc1[name, "Bailey"] == "immunogenic") {
        # increase the value of the first row and the first column
        contingency_table[1, 1] <- contingency_table[1, 1] + 1
    } else {
        # increase the value of the first row and the second column
        contingency_table[1, 2] <- contingency_table[1, 2] + 1
    }
  } else {
    # if the row is in bailey and the value is "immunogenic"
    if (!is.na(mmc1[name, "Bailey"]) && mmc1[name, "Bailey"] == "immunogenic") {
        # increase the value of the second row and the first column
        contingency_table[2, 1] <- contingency_table[2, 1] + 1
    } else {
        # increase the value of the second row and the second column
        contingency_table[2, 2] <- contingency_table[2, 2] + 1
    }
  }
}
# do fisher's exact test
fisher_result <- fisher.test(contingency_table)
# get the p-value
pval_immunogenic <- fisher_result$p.value
log_pval_immunogenic <- -log10(pval_immunogenic)

# get all the log p-values into dataframe
df <- data.frame(Group = rep(c("C1", "C2"), 10),
                Feature = rep(c("Collison: quasimesenchymal", "Bailey: squamous", "Moffitt: CLASSICAL", "Moffitt: BASAL-LIKE", "Bailey: pancreatic progenitor", "Living", "Deceased", "Collison: exocrine-like", "Bailey: ADEX", "Bailey: immunogenic"), each=2),
                LogPValue = c(NA, log_pval_quasi, NA, log_pval_squamous, -log_pval_classical, NA, NA, log_pval_basal, -log_pval_progenitor, NA, -log_pval_death, NA, NA, log_pval_death_c2, -log_pval_exocrine, NA, -log_pval_adex, NA, NA, log_pval_immunogenic),
                PValue = c(NA, pval_quasi, NA, pval_squamous, pval_classical, NA, NA, pval_basal, pval_progenitor, NA, pval_death, NA, NA, pval_death_c2, pval_exocrine, NA, pval_adex, NA, NA, pval_immunogenic))

df$Feature <- factor(df$Feature, levels = rev(unique(df$Feature)))

library(ggplot2)
# Plot, y ticks is Feature
ggplot(df, aes(x = LogPValue, y = Feature, size = PValue)) +
  geom_point(aes(color = Group), alpha = 0.7) +
  scale_size_continuous(trans = 'log10', range = c(2, 10), 
                        breaks = rev(c(1, 0.01, 0.0001, 1e-06, 1e-08)),
                        labels = rev(c("1", "0.01", "1e-04", "1e-06", "1e-08"))) +
  scale_x_continuous(position = "bottom", limits = c(-10,10)) +
  labs(size = "p-value", x = "+/- log10(p-value)", y = "") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  geom_vline(xintercept = log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed")

# Save the plot
ggsave("fig7c.png", width = 8, height = 6, dpi = 300)