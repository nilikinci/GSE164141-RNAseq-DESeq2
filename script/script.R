# ============================================================
# RNA-seq Differential Expression Analysis
# GSE164141
# ============================================================

# 1. Setup and Packages

library(DESeq2)
library(apeglm)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)

# 2. Data Import

count_file <- "../data/GSE164141_BT549_RawCounts_Matrix.csv.gz"

stopifnot(file.exists(count_file))

count_data <- read.csv(
  gzfile(count_file),
  check.names = FALSE)

colnames(count_data)[1] <- "gene"

stopifnot(!anyNA(count_data$gene))
stopifnot(!anyDuplicated(count_data$gene))

# 3. Metadata Preparation

metadata <- data.frame(sample = c(
  "BT549_DOXminus_1", "BT549_DOXminus_2", "BT549_DOXminus_3",
  "BT549_DOXplus_1", "BT549_DOXplus_2", "BT549_DOXplus_3"),
  condition = c(
    "DOXminus", "DOXminus", "DOXminus", "DOXplus", "DOXplus", "DOXplus"))

metadata$condition <- factor(
  metadata$condition,
  levels = c("DOXminus", "DOXplus"))

metadata <- metadata |>
  column_to_rownames(var = "sample")

# 4. Count Matrix Preparation

deseq_counts <- count_data |> dplyr::select(gene,
                                            BT549_DOXminus_1, BT549_DOXminus_2, BT549_DOXminus_3, BT549_DOXplus_1, BT549_DOXplus_2, BT549_DOXplus_3)

deseq_counts <- deseq_counts |> column_to_rownames(var = "gene")

stopifnot(
  identical(colnames(deseq_counts), rownames(metadata)))

stopifnot(
  all(!is.na(as.matrix(deseq_counts))),
  all(as.matrix(deseq_counts) >= 0),
  all(as.matrix(deseq_counts) == round(as.matrix(deseq_counts))))

# 5. DESeq2 Dataset Construction

dds <- DESeqDataSetFromMatrix(
  countData = deseq_counts,
  colData = metadata,
  design = ~ condition)

dds <- DESeq(dds)

stopifnot(
  "condition_DOXplus_vs_DOXminus" %in% resultsNames(dds)) 

results_data <- results(
  dds,
  contrast = c("condition", "DOXplus", "DOXminus"),
  alpha = 0.05)

results_df <- as.data.frame(results_data) |>
  rownames_to_column(var = "gene") |>
  mutate(
    fold_change = 2^log2FoldChange)

# 6. LFC Shrinkage

shrunk_results <- lfcShrink(
  dds,
  coef = "condition_DOXplus_vs_DOXminus",
  type = "apeglm")

shrunk_results_df <- as.data.frame(shrunk_results) |>
  rownames_to_column(var = "gene") |>
  rename(
    shrunk_log2FoldChange = log2FoldChange)

comparison_df <- results_df |>
  left_join(
    shrunk_results_df |>
      dplyr::select(gene, shrunk_log2FoldChange),
    by = "gene")

stopifnot(nrow(comparison_df) == nrow(results_df))

# 7. Significant DEG Identification and Classification (Raw log2FC)

significant_genes <- results_df |> filter(padj < 0.05, abs(log2FoldChange) > 1)

significant_genes <- significant_genes |> mutate(regulation = case_when(
  log2FoldChange > 1 ~ "upregulated", log2FoldChange < -1 ~ "downregulated") )

regulation_counts <- significant_genes |> group_by(regulation) |> summarise(n = n(), .groups = "drop")

# 8. DEG Classification Using Shrunken log2FC

significant_shrunk_genes <- shrunk_results_df |> filter(padj < 0.05, abs(shrunk_log2FoldChange) > 1)

significant_shrunk_genes <- significant_shrunk_genes |> mutate(regulation = case_when(
  shrunk_log2FoldChange > 1 ~ "upregulated", shrunk_log2FoldChange < -1 ~ "downregulated") )

regulation_shrunk_counts <- significant_shrunk_genes |> group_by(regulation) |> 
  summarise(n = n(), .groups = "drop")

# 9. Top Differentially Expressed Genes (Raw log2FC)

top10_up <- significant_genes |> filter(regulation == "upregulated") |> 
  slice_max(order_by = log2FoldChange, n = 10, with_ties = FALSE)

top10_down <- significant_genes |> filter(regulation == "downregulated") |> 
  slice_min(order_by = log2FoldChange, n = 10, with_ties = FALSE)

# 10. Top Differentially Expressed Genes Using Shrunken log2FC

top10_up_shrunk <- significant_shrunk_genes |> filter(regulation == "upregulated") |> 
  slice_max(order_by = shrunk_log2FoldChange, n = 10, with_ties = FALSE)

top10_down_shrunk <- significant_shrunk_genes |> filter(regulation == "downregulated") |> 
  slice_min(order_by = shrunk_log2FoldChange, n = 10, with_ties = FALSE)

# 11. Quality Control and Exploratory Analysis

#PCA

vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)

percent_var <- round(100 * attr(pca_data, "percentVar") )  

pca_plot <- ggplot(pca_data,
                   aes(
                     x = PC1,
                     y = PC2,
                     color = condition) ) + geom_point(aes(shape = condition) ) +   
  theme_classic() + 
  geom_text_repel(aes(label = name), vjust = -0.8, show.legend = FALSE ) +  
  labs(
    x = paste0(
      "PC1 (", percent_var[1], "% variance)" ),
    y = paste0(
      "PC2 (", percent_var[2], "% variance)" ) )

ggsave(filename = "../figures/PCA_all_conditions.png", plot = pca_plot, width = 8, height = 6, dpi = 300)

# Sample-to-sample distance

sample_dist <- dist(t(assay(vsd)))

sample_dist_matrix <- as.matrix(sample_dist)

heatmap_annotation <- data.frame(
  Group = factor(
    metadata$condition,
    levels = c("DOXminus", "DOXplus"),
    labels = c("DOX-", "DOX+")),
  row.names = rownames(metadata))

pheatmap(sample_dist_matrix, annotation_col = heatmap_annotation, 
         filename = "../figures/sample_distance_heatmap.png", width = 8, height = 7)

# Library size

library_sizes <- colSums(assay(dds))

# Expression distribution

vsd_long <- as.data.frame(assay(vsd)) |> rownames_to_column(var = "gene") |>
  pivot_longer(
    cols = -gene,
    names_to = "sample",
    values_to = "expression")

expression_distribution_plot <- ggplot(vsd_long, aes(
  x = sample, y = expression) ) +
  geom_boxplot(outlier.shape = NA) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1) )

ggsave(filename =
         "../figures/expression_distribution_boxplot.png", plot =
         expression_distribution_plot, width = 8, height = 6, dpi = 300)

# Top variable genes

gene_variance <- apply(assay(vsd), 1, var)

top10_idx <- order(gene_variance, decreasing = TRUE)[1:10]   

top10_matrix <- assay(vsd)[top10_idx, ]

pheatmap(top10_matrix, scale = "row", annotation_col = heatmap_annotation, 
         filename = "../figures/top10_variable_genes_heatmap_v2.png",
         width = 9, height = 9, fontsize_row = 8, fontsize_col = 9, angle_col = 90)

# 12. Volcano Plot

# Raw log2FC

volcano_df <- results_df |>
  filter(
    !is.na(padj),
    !is.na(log2FoldChange)) |>
  mutate(
    interpretation = case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "upregulated",
      padj < 0.05 & log2FoldChange < -1 ~ "downregulated",
      TRUE ~ "not DEG"),
    neg_log10_padj = -log10(pmax(padj, .Machine$double.xmin)))

volcano_plot <- ggplot(
  volcano_df,
  aes(
    x = log2FoldChange,
    y = neg_log10_padj,
    color = interpretation)) +
  geom_point(
    alpha = 0.7,
    size = 1.5) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed") + geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed") +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.12))) + labs(
      title = "Volcano Plot - DOXplus vs DOXminus (Raw log2FC)",
      subtitle = "DEGs: padj < 0.05 and |log2FC| > 1",
      x = "log2 Fold Change",
      y = "-log10 adjusted p-value (capped for plotting)",
      color = "Interpretation") +
  theme_classic()

ggsave(
  "../figures/volcano_plot.png",
  plot = volcano_plot,
  width = 12,
  height = 8,
  dpi = 300)

# Shrunken log2FC

shrunk_volcano_df <- shrunk_results_df |>
  filter(
    !is.na(padj),
    !is.na(shrunk_log2FoldChange)) |>
  mutate(
    interpretation = case_when(
      padj < 0.05 & shrunk_log2FoldChange > 1 ~ "upregulated",
      padj < 0.05 & shrunk_log2FoldChange < -1 ~ "downregulated",
      TRUE ~ "not DEG"),
    neg_log10_padj = -log10(pmax(padj, .Machine$double.xmin)))

shrunk_volcano_plot <- ggplot(
  shrunk_volcano_df,
  aes(
    x = shrunk_log2FoldChange,
    y = neg_log10_padj,
    color = interpretation)) +
  geom_point(
    alpha = 0.7,
    size = 1.5) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed") + geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed") +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.12))) + labs(
      title = "Volcano Plot - DOXplus vs DOXminus (Shrunken log2FC)",
      subtitle = "DEGs: padj < 0.05 and |Shrunken log2FC| > 1",
      x = "Shrunken log2 Fold Change",
      y = "-log10 adjusted p-value (capped for plotting)",
      color = "Interpretation") +
  theme_classic()

ggsave(
  "../figures/shrunk_volcano_plot.png",
  plot = shrunk_volcano_plot,
  width = 12,
  height = 8,
  dpi = 300)

# 13. Raw vs Shrunken log2FC Comparison

lfc_comparison_plot <- comparison_df |>
  filter(
    !is.na(log2FoldChange),
    !is.na(shrunk_log2FoldChange)) |>
  ggplot(
    aes(
      x = log2FoldChange,
      y = shrunk_log2FoldChange)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed") +
  coord_equal() +
  labs(
    title = "Raw vs Shrunken log2 Fold Change",
    x = "Raw log2 Fold Change",
    y = "Shrunken log2 Fold Change") +
  theme_classic()

ggsave(
  "../figures/raw_vs_shrunk_log2FC.png",
  plot = lfc_comparison_plot,
  width = 8,
  height = 8,
  dpi = 300)

# 14. Export Results

write.csv(
  results_df,
  "../results/DESeq2_all_results.csv",
  row.names = FALSE)

write.csv(
  comparison_df,
  "../results/DESeq2_shrunk_vs_raw_results.csv",
  row.names = FALSE)

write.csv(
  shrunk_results_df,
  "../results/DESeq2_shrunk_results.csv",
  row.names = FALSE)

write.csv(
  significant_shrunk_genes,
  "../results/shrunk_significant_DEGs.csv",
  row.names = FALSE)

write.csv(
  significant_genes,
  "../results/significant_DEGs.csv",
  row.names = FALSE)

write.csv(
  regulation_counts,
  "../results/DEG_counts.csv",
  row.names = FALSE)

write.csv(
  regulation_shrunk_counts,
  "../results/DEG_counts_shrunk.csv",
  row.names = FALSE)

write.csv(
  top10_up,
  "../results/top10_upregulated_genes.csv",
  row.names = FALSE)

write.csv(
  top10_up_shrunk,
  "../results/top10_shrunk_upregulated_genes.csv",
  row.names = FALSE)

write.csv(
  top10_down,
  "../results/top10_downregulated_genes.csv",
  row.names = FALSE)

write.csv(
  top10_down_shrunk,
  "../results/top10_shrunk_downregulated_genes.csv",
  row.names = FALSE)

capture.output(
  sessionInfo(),
  file = "../results/sessionInfo.txt")