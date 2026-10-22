
# ============================================================
# GSE54723: Mouse Liver Microarray Analysis
# Comparison: Palatinose vs Sucrose
# ============================================================

required_packages <- c(
  "GEOquery", "Biobase", "limma",
  "AnnotationDbi", "illuminaMousev2.db",
  "ggplot2", "pheatmap"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Install missing packages: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(GEOquery)
library(Biobase)
library(limma)
library(AnnotationDbi)
library(illuminaMousev2.db)
library(ggplot2)
library(pheatmap)

# Use the project root as the working directory.
project_dir <- getwd()

fig_dir <- file.path(project_dir, "results", "figures")
results_dir <- file.path(project_dir, "results", "tables")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. Download GEO data and assign sample groups
# ------------------------------------------------------------

gse <- getGEO(
  "GSE54723",
  GSEMatrix = TRUE,
  getGPL = FALSE
)

if (length(gse) > 1) {
  idx <- grep("GPL6885", names(gse))
  if (length(idx) != 1) {
    stop("Could not uniquely identify platform GPL6885.")
  }
} else {
  idx <- 1
}

eset <- gse[[idx]]
expr_matrix <- exprs(eset)
sample_info <- pData(eset)

sample_info$group <- factor(
  sub(
    "^treatment: ",
    "",
    sample_info$characteristics_ch1
  ),
  levels = c("sucrose", "palatinose")
)

if (!identical(rownames(sample_info), colnames(expr_matrix))) {
  stop("Sample metadata and expression columns do not match.")
}

if (anyNA(sample_info$group) ||
    !identical(as.integer(table(sample_info$group)), c(7L, 7L))) {
  stop("Expected 7 sucrose and 7 palatinose samples.")
}

if (any(!is.finite(expr_matrix)) || any(expr_matrix <= 0)) {
  stop("Expression values need review before log transformation.")
}

# ------------------------------------------------------------
# 2. Log transformation
# ------------------------------------------------------------
# GEO2R's automatic heuristic applies log2 to this dataset.
# This matches GEO2R's procedure but does not independently
# validate the complete original preprocessing history.

log_expr <- log2(expr_matrix)

# Save the sample information used in the analysis.
write.csv(
  data.frame(
    GEO_ID = rownames(sample_info),
    Title = sample_info$title,
    Group = sample_info$group
  ),
  file.path(results_dir, "sample_metadata.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 3. Quality-control plots
# ------------------------------------------------------------

# 3a. Expression distributions before and after log2.
png(
  file.path(fig_dir, "01_expression_distribution.png"),
  width = 1200, height = 800, res = 120
)

par(mar = c(8, 4, 3, 1))
boxplot(
  expr_matrix,
  las = 2,
  outline = FALSE,
  main = "GSE54723: Original expression values",
  ylab = "Expression value",
  xlab = "GEO sample"
)
dev.off()

png(
  file.path(fig_dir, "02_log2_expression_distribution.png"),
  width = 1200, height = 800, res = 120
)

par(mar = c(8, 4, 3, 1))
boxplot(
  log_expr,
  las = 2,
  outline = FALSE,
  main = "GSE54723: log2-transformed expression",
  ylab = "log2(expression)",
  xlab = "GEO sample"
)
dev.off()

# 3b. PCA: each point represents one sample.
pca <- prcomp(t(log_expr), center = TRUE, scale. = FALSE)

pca_variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)

pca_data <- data.frame(
  GEO_ID = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Group = sample_info[rownames(pca$x), "group"]
)

p_pca <- ggplot(
  pca_data,
  aes(x = PC1, y = PC2, color = Group, label = GEO_ID)
) +
  geom_point(size = 3) +
  labs(
    title = "PCA of GSE54723 liver samples",
    x = sprintf("PC1 (%.2f%%)", pca_variance[1]),
    y = sprintf("PC2 (%.2f%%)", pca_variance[2]),
    color = "Diet group"
  ) +
  theme_classic()

ggsave(
  file.path(fig_dir, "03_pca_samples.png"),
  p_pca, width = 8, height = 6, dpi = 300
)

# 3c. Sample-to-sample Euclidean distance heatmap.
sample_dist <- dist(t(log_expr))
distance_matrix <- as.matrix(sample_dist)

annotation_col <- data.frame(
  Group = sample_info[colnames(log_expr), "group"]
)
rownames(annotation_col) <- colnames(log_expr)

pheatmap(
  distance_matrix,
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  main = "Sample-to-sample Euclidean distances",
  filename = file.path(fig_dir, "04_sample_distance_heatmap.png"),
  width = 9,
  height = 8
)

# ------------------------------------------------------------
# 4. Differential expression using limma
# ------------------------------------------------------------
# Coefficient grouppalatinose means palatinose relative
# to sucrose. Negative logFC indicates lower expression
# in palatinose; positive logFC indicates higher expression.

design <- model.matrix(~ group, data = sample_info)

fit <- lmFit(log_expr, design)
fit <- eBayes(fit)

de_results <- topTable(
  fit,
  coef = "grouppalatinose",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

de_results$ProbeID <- rownames(de_results)

# ------------------------------------------------------------
# 5. Probe annotation — retain one row per probe
# ------------------------------------------------------------

probe_annotation <- AnnotationDbi::select(
  illuminaMousev2.db,
  keys = unique(de_results$ProbeID),
  keytype = "PROBEID",
  columns = c("SYMBOL", "GENENAME", "ENTREZID")
)

# Remove exact duplicate annotation records.
probe_annotation <- unique(probe_annotation)

# Keep one annotation record per probe for a one-row-per-probe
# output table. Ambiguous mappings require further review.
probe_annotation <- probe_annotation[
  !duplicated(probe_annotation$PROBEID),
]

annotated_results <- merge(
  de_results,
  probe_annotation,
  by.x = "ProbeID",
  by.y = "PROBEID",
  all.x = TRUE,
  sort = FALSE
)

# Restore the original differential-expression ranking.
result_order <- match(
  annotated_results$ProbeID,
  de_results$ProbeID
)
annotated_results <- annotated_results[order(result_order), ]

# Validate the output.
stopifnot(
  nrow(annotated_results) == nrow(de_results),
  !anyDuplicated(annotated_results$ProbeID)
)

write.csv(
  de_results,
  file.path(results_dir, "differential_expression_all_probes.csv"),
  row.names = FALSE
)

write.csv(
  annotated_results,
  file.path(results_dir, "differential_expression_annotated.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 6. Volcano plot
# ------------------------------------------------------------

volcano_data <- de_results
volcano_data$neg_log10_p <- -log10(
  pmax(volcano_data$P.Value, .Machine$double.xmin)
)
volcano_data$FDR_significant <- (
  volcano_data$adj.P.Val < 0.05
)

p_volcano <- ggplot(
  volcano_data,
  aes(x = logFC, y = neg_log10_p, color = FDR_significant)
) +
  geom_point(alpha = 0.65, size = 1.4) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dotted"
  ) +
  scale_color_manual(
    values = c("FALSE" = "grey55", "TRUE" = "red"),
    labels = c("FALSE" = "Not FDR-significant",
               "TRUE" = "FDR-significant"),
    name = "Result"
  ) +
  labs(
    title = "Differential expression: Palatinose vs Sucrose",
    subtitle = "Horizontal line shows raw p = 0.05",
    x = "log2 fold change",
    y = "-log10(raw p-value)"
  ) +
  theme_classic()

ggsave(
  file.path(fig_dir, "05_volcano_plot.png"),
  p_volcano, width = 8, height = 6, dpi = 300
)

# ------------------------------------------------------------
# 7. Summary table, including SOCS2
# ------------------------------------------------------------

socs2_rows <- annotated_results[
  !is.na(annotated_results$SYMBOL) &
    annotated_results$SYMBOL == "Socs2",
]

if (nrow(socs2_rows) > 0) {
  socs2_row <- socs2_rows[1, ]
  socs2_probe <- socs2_row$ProbeID
  socs2_logFC <- socs2_row$logFC
  socs2_p <- socs2_row$P.Value
  socs2_fdr <- socs2_row$adj.P.Val
} else {
  socs2_probe <- NA_character_
  socs2_logFC <- NA_real_
  socs2_p <- NA_real_
  socs2_fdr <- NA_real_
}

summary_table <- data.frame(
  Dataset = "GSE54723",
  Comparison = "Palatinose vs Sucrose",
  Samples = ncol(log_expr),
  Sucrose_samples = sum(sample_info$group == "sucrose"),
  Palatinose_samples = sum(sample_info$group == "palatinose"),
  Probes_tested = nrow(de_results),
  FDR_significant = sum(de_results$adj.P.Val < 0.05, na.rm = TRUE),
  SOCS2_ProbeID = socs2_probe,
  SOCS2_logFC = socs2_logFC,
  SOCS2_raw_p = socs2_p,
  SOCS2_adjusted_p = socs2_fdr
)

write.csv(
  summary_table,
  file.path(results_dir, "analysis_summary.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 8. Reproducibility information
# ------------------------------------------------------------

writeLines(
  capture.output(sessionInfo()),
  file.path(project_dir, "sessionInfo.txt")
)

message("Analysis completed.")
message("Review preprocessing assumptions before biological interpretation.")
message(
  "FDR-significant probes: ",
  sum(de_results$adj.P.Val < 0.05, na.rm = TRUE)
)

