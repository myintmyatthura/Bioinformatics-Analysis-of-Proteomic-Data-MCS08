library(limma)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

data <- read.csv("ProteinquantitationHE&CI.csv", row.names = 1)

# filter immunoglobulin
data <- data[!grepl("immunoglobulin", data$ProteinDescriptions, ignore.case = TRUE), ]

# separate the data based on groups for data cleaning 
CI_cols <- grep("^CI", colnames(data))
HE_cols <- grep("^HE", colnames(data))
CI_data <- data[, CI_cols]
HE_data <- data[, HE_cols]

# impute data using mean of the group
impute_mean <- function(group_data) {
  imputed <- apply(group_data, 1, function(impute){
    impute[is.na(impute)] <- mean(impute, na.rm = TRUE)
    return(impute)
  })
  return((t(imputed)))
}

# apply the function to the groups
CI_data_imputed <- impute_mean(CI_data)
HE_data_imputed <- impute_mean(HE_data)

# combine and remove row with NaN values
# data contains only numeric values
data_preprocessed <- na.omit(cbind(CI_data_imputed, HE_data_imputed))

# convert into matrix format
data_matrix <- as.matrix(data_preprocessed)

# normalization 
log_data_matrix <- log2(data_matrix + 1)
normalized_matrix <- normalizeQuantiles(log_data_matrix)

# create experimental group and design matrix model
group <- factor(c(rep("CI", 4), rep("HE", 4)))
design <- model.matrix(~0 + group) 

# define contrast with the 2 columns
contrast.matrix <- makeContrasts(groupCI - groupHE, levels = design)

# fit the model
fit <- lmFit(normalized_matrix, design)
# apply contrast
fitlc <- contrasts.fit(fit, contrast.matrix)
fitlc <- eBayes(fitlc)

# extract result
result <- topTable(fitlc, adjust="BH", number = Inf)

# filter based on log fold change and adjusted p value
filtered_results <- result[abs(result$logFC) > 1 & result$adj.P.Val < 0.05, ]

# Create a significance column for classification
result$Significance <- ifelse(result$adj.P.Val < 0.05 & abs(result$logFC) > 1, 
                              ifelse(result$logFC > 1, "Upregulated", "Downregulated"), 
                              "Not Significant")

# Select only significant proteins for labeling
significant_proteins <- subset(result, Significance != "Not Significant")

# Add protein names (assuming row names contain the protein names)
significant_proteins$Protein <- rownames(significant_proteins)

# Create the volcano plot
volcano_plot <- ggplot(result, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
  theme_minimal() +
  labs(title = "Volcano Plot",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-Value") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  theme(legend.position = "right") +
  
  # Add labels only for significant proteins
  geom_text(data = significant_proteins, aes(label = Protein), 
            vjust = -0.5, hjust = 0.5, size = 3, check_overlap = TRUE)

# Display the plot
print(volcano_plot)

# Filter the normalized matrix to only include significant proteins
filtered_expression_matrix <- normalized_matrix[significant_proteins$Protein, ]

pheatmap(
  filtered_expression_matrix, 
  scale = "row",
  clustering_distance_rows = "euclidean", 
  clustering_method = "ward.D2",
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "logFC Heatmap with Clustering"
)





