# ============================================================
# E. coli RNA-Seq Differential Expression Analysis
# Water vs Ciprofloxacin
# ============================================================


# ============================================================
# 1. Load Required Packages
# ============================================================

library(DESeq2)
library(pheatmap)
library(ggplot2)


# ============================================================
# 2. Define Input and Output Directories
# ============================================================

count_dir <- "results/counts"
plot_dir <- "results/plots"

if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}


# ============================================================
# 3. Read Gene Count Files
# ============================================================

# Read the first Water sample
sample1 <- read.table(
  file.path(count_dir, "SRR22578513_gene_counts.txt"),
  header = TRUE,
  sep = "\t",
  comment.char = "#"
)

# Read the second Water sample
sample2 <- read.table(
  file.path(count_dir, "SRR22578514_gene_counts.txt"),
  header = TRUE,
  sep = "\t",
  comment.char = "#"
)

# Read the first Ciprofloxacin sample
sample3 <- read.table(
  file.path(count_dir, "SRR22578533_gene_counts.txt"),
  header = TRUE,
  sep = "\t",
  comment.char = "#"
)

# Read the second Ciprofloxacin sample
sample4 <- read.table(
  file.path(count_dir, "SRR22578536_gene_counts.txt"),
  header = TRUE,
  sep = "\t",
  comment.char = "#"
)


# ============================================================
# 4. Construct Gene Count Matrix
# ============================================================

countData <- data.frame(
  row.names = sample1$Geneid,
  Water_1 = sample1[, 7],
  Water_2 = sample2[, 7],
  Cipro_1 = sample3[, 7],
  Cipro_2 = sample4[, 7]
)

# Check the count matrix
head(countData)

# Check dimensions
dim(countData)


# ============================================================
# 5. Define Experimental Conditions
# ============================================================

colData <- data.frame(
  condition = factor(
    c(
      "Water",
      "Water",
      "Ciprofloxacin",
      "Ciprofloxacin"
    ),
    levels = c("Ciprofloxacin", "Water")
  )
)

# Assign sample names to metadata
rownames(colData) <- colnames(countData)

# Display sample information
colData


# ============================================================
# 6. Create DESeq2 Dataset
# ============================================================

dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData = colData,
  design = ~ condition
)


# ============================================================
# 7. Filter Low-Count Genes
# ============================================================

dds <- dds[
  rowSums(counts(dds)) >= 10,
]

# Check dimensions after filtering
dim(dds)


# ============================================================
# 8. Run DESeq2 Differential Expression Analysis
# ============================================================

dds <- DESeq(dds)


# ============================================================
# 9. Extract Water vs Ciprofloxacin Results
# ============================================================

res <- results(
  dds,
  contrast = c(
    "condition",
    "Water",
    "Ciprofloxacin"
  )
)

# Display DESeq2 summary
summary(res)


# ============================================================
# 10. Count Significant Genes
# ============================================================

significant_gene_count <- sum(
  res$padj < 0.05,
  na.rm = TRUE
)

cat(
  "Number of significant genes:",
  significant_gene_count,
  "\n"
)


# ============================================================
# 11. Create Significant Gene Table
# ============================================================

sig <- subset(
  res,
  !is.na(padj) & padj < 0.05
)

# Sort genes by adjusted p-value
sig <- sig[
  order(sig$padj),
]

# Convert to data frame
sig_df <- as.data.frame(sig)

# Add gene IDs
sig_df$Gene <- rownames(sig_df)

# Reorder columns
sig_df <- sig_df[
  ,
  c(
    "Gene",
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj"
  )
]

# Save significant genes
write.csv(
  sig_df,
  "Significant_Genes_Water_vs_Ciprofloxacin.csv",
  row.names = FALSE
)


# ============================================================
# 12. Separate Genes by Direction of Expression Change
# ============================================================

# Positive log2 fold change:
# Higher expression in Water relative to Ciprofloxacin

water_up <- subset(
  sig_df,
  log2FoldChange > 0
)

# Negative log2 fold change:
# Higher expression in Ciprofloxacin relative to Water

cipro_up <- subset(
  sig_df,
  log2FoldChange < 0
)

# Display counts
cat(
  "Genes higher in Water:",
  nrow(water_up),
  "\n"
)

cat(
  "Genes higher in Ciprofloxacin:",
  nrow(cipro_up),
  "\n"
)


# Save directional gene lists
write.csv(
  water_up,
  "Water_Higher_Expression_Genes.csv",
  row.names = FALSE
)

write.csv(
  cipro_up,
  "Ciprofloxacin_Higher_Expression_Genes.csv",
  row.names = FALSE
)


# ============================================================
# 13. MA Plot
# ============================================================

png(
  file.path(
    plot_dir,
    "MA_Plot_Water_vs_Ciprofloxacin.png"
  ),
  width = 1000,
  height = 800,
  res = 120
)

plotMA(
  res,
  ylim = c(-10, 10),
  alpha = 0.05,
  main = "MA Plot: Water vs Ciprofloxacin"
)

dev.off()


# Display MA plot in R
plotMA(
  res,
  ylim = c(-10, 10),
  alpha = 0.05,
  main = "MA Plot: Water vs Ciprofloxacin"
)


# ============================================================
# 14. Prepare Data for Volcano Plot
# ============================================================

volcano_data <- as.data.frame(res)

# Add gene IDs
volcano_data$Gene <- rownames(volcano_data)

# Remove genes with missing p-values
volcano_data <- volcano_data[
  !is.na(volcano_data$pvalue),
]

# Define significance
volcano_data$Significance <- "Not Significant"

volcano_data$Significance[
  !is.na(volcano_data$padj) &
    volcano_data$padj < 0.05
] <- "Significant"


# ============================================================
# 15. Volcano Plot
# ============================================================

png(
  file.path(
    plot_dir,
    "Volcano_Plot_Water_vs_Ciprofloxacin.png"
  ),
  width = 1000,
  height = 800,
  res = 120
)

plot(
  volcano_data$log2FoldChange,
  -log10(volcano_data$pvalue),
  pch = 16,
  xlab = "Log2 Fold Change",
  ylab = "-Log10 P-value",
  main = "Volcano Plot: Water vs Ciprofloxacin"
)

# Add log2 fold-change reference line
abline(
  v = 0,
  lty = 2
)

# Add p-value reference line
abline(
  h = -log10(0.05),
  lty = 2
)

# Highlight significant genes
points(
  volcano_data$log2FoldChange[
    volcano_data$Significance == "Significant"
  ],
  -log10(
    volcano_data$pvalue[
      volcano_data$Significance == "Significant"
    ]
  ),
  pch = 16
)

dev.off()


# Display Volcano plot in R
plot(
  volcano_data$log2FoldChange,
  -log10(volcano_data$pvalue),
  pch = 16,
  xlab = "Log2 Fold Change",
  ylab = "-Log10 P-value",
  main = "Volcano Plot: Water vs Ciprofloxacin"
)

abline(
  v = 0,
  lty = 2
)

abline(
  h = -log10(0.05),
  lty = 2
)

points(
  volcano_data$log2FoldChange[
    volcano_data$Significance == "Significant"
  ],
  -log10(
    volcano_data$pvalue[
      volcano_data$Significance == "Significant"
    ]
  ),
  pch = 16
)


# ============================================================
# 16. Variance Stabilizing Transformation
# ============================================================

vsd <- vst(
  dds,
  blind = FALSE
)


# ============================================================
# 17. PCA Analysis
# ============================================================

pca_plot <- plotPCA(
  vsd,
  intgroup = "condition"
) +
  ggtitle(
    "PCA Plot: Water vs Ciprofloxacin"
  )

# Display PCA plot
pca_plot

# Save PCA plot
ggsave(
  filename = file.path(
    plot_dir,
    "PCA_Plot_Water_vs_Ciprofloxacin.png"
  ),
  plot = pca_plot,
  width = 8,
  height = 6,
  dpi = 300
)


# ============================================================
# 18. Select Top 50 Significant Genes
# ============================================================

# Select the top 50 genes based on adjusted p-value
top_genes <- rownames(sig)[
  1:min(50, nrow(sig))
]


# ============================================================
# 19. Prepare Data for Heatmap
# ============================================================

# Extract normalized counts
normalized_counts <- counts(
  dds,
  normalized = TRUE
)

# Select top significant genes
heatmap_data <- normalized_counts[
  top_genes,
]

# Apply log2 transformation
heatmap_data_log <- log2(
  heatmap_data + 1
)


# ============================================================
# 20. Create Sample Annotation
# ============================================================

annotation_col <- data.frame(
  Condition = colData(dds)$condition
)

# Assign sample names
rownames(annotation_col) <- colnames(
  heatmap_data_log
)


# ============================================================
# 21. Heatmap
# ============================================================

png(
  file.path(
    plot_dir,
    "Heatmap_Top50_DEGs.png"
  ),
  width = 1200,
  height = 1000,
  res = 120
)

pheatmap(
  heatmap_data_log,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  main = "Top 50 Differentially Expressed Genes"
)

dev.off()


# Display Heatmap in R
pheatmap(
  heatmap_data_log,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  main = "Top 50 Differentially Expressed Genes"
)


# ============================================================
# 22. Final Analysis Summary
# ============================================================

cat("\n")
cat("============================================\n")
cat("RNA-Seq Differential Expression Summary\n")
cat("============================================\n")
cat(
  "Total genes after filtering:",
  nrow(res),
  "\n"
)
cat(
  "Significant genes (padj < 0.05):",
  significant_gene_count,
  "\n"
)
cat(
  "Genes higher in Water:",
  nrow(water_up),
  "\n"
)
cat(
  "Genes higher in Ciprofloxacin:",
  nrow(cipro_up),
  "\n"
)
cat("============================================\n")
