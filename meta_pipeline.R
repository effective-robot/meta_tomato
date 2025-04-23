# Load necessary libraries
library(sva)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(metaRNASeq)

### Step 1: Load Data ###
counts <- read.csv("counts_matrix.csv", row.names = 1)
metadata <- read.csv("metadata.csv")

# ✅ Create output directory if it doesn't exist
output_dir <- "meta_tomato"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# ✅ Select which stresses to keep (Modify this list before running)
stresses_to_keep <- c("clvb", "rs", "eb", "cf", "cf2",
                      "lb", "lb2", "nem", "nem2", "pvx", "tswv", "tswv2")  # Modify this list as needed

# ✅ Filter metadata and counts for selected stresses
metadata_filtered <- metadata[metadata$Stress %in% stresses_to_keep, ]
counts_filtered <- counts[, metadata_filtered$Sample]

# Ensure row names match sample names
rownames(metadata_filtered) <- metadata_filtered$Sample
all(rownames(metadata_filtered) == colnames(counts_filtered))  # Should return TRUE

### Step 2: Perform ComBat-Seq (Batch Effect Removal) ###
# Convert counts to matrix format
counts_matrix <- as.matrix(counts_filtered)

# Apply ComBat-Seq to remove batch effects
corrected_counts <- ComBat_seq(counts_matrix, batch = metadata_filtered$Batch)

# Ensure corrected counts are non-negative
corrected_counts[corrected_counts < 0] <- 0

# Save batch-corrected counts
write.csv(corrected_counts, file.path(output_dir, "batch_corrected_counts.csv"), row.names = TRUE)

# ✅ Find and remove genes (rows) with zero variance
zero_variance_genes <- rownames(corrected_counts)[apply(corrected_counts, 1, var) == 0]
cat("Genes with zero variance:", length(zero_variance_genes), "\n")
corrected_counts <- corrected_counts[apply(corrected_counts, 1, var) > 0, ]

### Step 3: PCA on Batch-Corrected Counts ###
pca_res <- prcomp(t(corrected_counts), scale. = TRUE)
metadata_filtered$PC1 <- pca_res$x[,1]
metadata_filtered$PC2 <- pca_res$x[,2]

# Save PCA plot
pca_plot <- ggplot(metadata_filtered, aes(x = PC1, y = PC2, color = Stress)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA of Batch-Corrected Counts", x = "PC1", y = "PC2")
ggsave(file.path(output_dir, "PCA_plot.png"), pca_plot, width = 7, height = 5)

### Step 4: Run DEA on Selected Stresses ###
dea_results <- list()

for (stress in stresses_to_keep) {
  if (stress %in% metadata_filtered$Stress) {
    cat("Running DEA for:", stress, "\n")
    
    # Filter metadata and counts for this stress
    metadata_stress <- metadata_filtered[metadata_filtered$Stress == stress, ]
    counts_stress <- corrected_counts[, metadata_stress$Sample]
    
    # Create DESeq dataset
    dds_stress <- DESeqDataSetFromMatrix(countData = counts_stress,
                                         colData = metadata_stress,
                                         design = ~ Condition)
    
    # Run DEA
    dds_stress <- DESeq(dds_stress)
    results_stress <- results(dds_stress)
    
    # Save results
    dea_results[[stress]] <- results_stress$pvalue
    write.csv(as.data.frame(results_stress), file.path(output_dir, paste0("DEA_results_", stress, ".csv")), row.names = TRUE)
  }
}

### Step 5: Perform Meta-Analysis ###
meta_res <- fishercomb(dea_results)

meta_res_df <- data.frame(
  Gene = rownames(corrected_counts),
  rawpval = meta_res$rawpval,
  adjpval = meta_res$adjpval
)

# Filter significant genes
significant_genes <- meta_res_df[meta_res_df$adjpval < 0.05, ]

if (nrow(significant_genes) > 0) {
  write.csv(significant_genes, file.path(output_dir, "meta_analysis_results.csv"), row.names = FALSE)
  cat("Meta-analysis results saved.\n")
} else {
  cat("No significant genes found.\n")
}

write.table(meta_res_df$Gene[meta_res_df$adjpval < 0.05], file.path(output_dir, "genes.csv"), row.names = FALSE, col.names = FALSE, quote = FALSE)

cat("Pipeline completed successfully! Results saved in:", output_dir, "\n")


