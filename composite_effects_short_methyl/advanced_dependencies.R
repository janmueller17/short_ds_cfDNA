library(readr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(ggridges)
library(ggthemes)

# Define constants
base_path = "X:/db05/Jan/scfDNA_intro_paper/analysis/for_github/composite_effects_short_methyl/"
# base_path = "/base/path/for/analysis/"

pseudocount = 0.01

# Load data
full_df = read_csv(paste0(base_path, "/unified_count_matrices_with_pseudocount_covariates.csv")) %>%
  column_to_rownames(var = "Geneid")

# Define subsets of features
feature_patterns = c("GC_pct_", "mean_vals", "rel_std_vals", "_S42_", "_S44_", "_S45_", "_S46_")
feature_dfs = map(feature_patterns, ~ full_df %>% select(contains(.x)))

# GC content adjustment function
adjust_for_GC = function(df_in, gc_in, val_in, type, scaling) {
  work_df = df_in
  lowessFit = lowess(x = work_df[[gc_in]], y = work_df[[val_in]], iter = 10) %>%
    as.data.frame()
  colnames(lowessFit) = c("GC", "fit_val")
  lowessFit$fit_val[lowessFit$fit_val < 0] = min(lowessFit$fit_val[lowessFit$fit_val > 0])

  if (type == "OoE") {
    OoE = work_df[[val_in]] / lowessFit$fit_val
    OoE[OoE < 0] = min(OoE[OoE > 0])
    OoE = log2(OoE)
    df_in = df_in %>% mutate(!!paste0(val_in, "_GC_adj_OoE") := OoE)
  } else if (type == "residual") {
    resid_val = work_df[[val_in]] - lowessFit$fit_val
    if (scaling == "log") {
      resid_val = log10(resid_val)
    } else if (scaling == "sqrt") {
      resid_val = sqrt(resid_val)
    }
    df_in = df_in %>% mutate(!!paste0(val_in, "_GC_adj_resid") := resid_val)
  }

  return(df_in)
}

# Apply GC content adjustment
mean_vals = full_df %>% select(ends_with("mean_vals"))
covariates = full_df %>% select(starts_with("GC_pct"))
mean_and_gc_df = merge(mean_vals, covariates, by = 0, all = TRUE)
row.names(mean_and_gc_df) = mean_and_gc_df$Row.names
mean_and_gc_df$Row.names = NULL

# List of GC-content and mean value pairs to adjust
adjustment_list = list(
  c("GC_pct_Gencode_gene_bodies", "RNA_Gencode_gene_bodies_mean_vals"),
  c("GC_pct_Gencode_5UTR", "RNA_Gencode_5UTR_mean_vals"),
  c("GC_pct_Gencode_3UTR", "RNA_Gencode_3UTR_mean_vals"),
  c("GC_pct_Gencode_gene_bodies", "short_cf_Gencode_gene_bodies_mean_vals"),
  c("GC_pct_Gencode_5UTR", "short_cf_Gencode_5UTR_mean_vals"),
  c("GC_pct_Gencode_3UTR", "short_cf_Gencode_3UTR_mean_vals"),
  c("GC_pct_epd_TSS_up_1000", "short_cf_epd_TSS_up_1000_mean_vals"),
  c("GC_pct_epd_TSS_up_2to5", "short_cf_epd_TSS_up_2to5_mean_vals"),
  c("GC_pct_epd_TSS_down_1000", "short_cf_epd_TSS_down_1000_mean_vals"),
  c("GC_pct_epd_core_promoters", "short_cf_epd_core_promoters_mean_vals"),
  c("GC_pct_CpG_islands_matched_EPD_TSS", "short_cf_CpG_islands_matched_EPD_TSS_mean_vals"),
  c("GC_pct_Gencode_gene_bodies", "MBD_Gencode_gene_bodies_mean_vals"),
  c("GC_pct_Gencode_5UTR", "MBD_Gencode_5UTR_mean_vals"),
  c("GC_pct_Gencode_3UTR", "MBD_Gencode_3UTR_mean_vals"),
  c("GC_pct_epd_TSS_up_1000", "MBD_epd_TSS_up_1000_mean_vals"),
  c("GC_pct_epd_TSS_up_2to5", "MBD_epd_TSS_up_2to5_mean_vals"),
  c("GC_pct_epd_TSS_down_1000", "MBD_epd_TSS_down_1000_mean_vals"),
  c("GC_pct_epd_core_promoters", "MBD_epd_core_promoters_mean_vals"),
  c("GC_pct_CpG_islands_matched_EPD_TSS", "MBD_CpG_islands_matched_EPD_TSS_mean_vals")
)

# Adjust all pairs
for (pair in adjustment_list) {
  mean_and_gc_df = adjust_for_GC(df_in = mean_and_gc_df, gc_in = pair[1], val_in = pair[2], type = "OoE", scaling = "none")
}

# Exclude problematic genes
exclude_genes = c("ENSG00000002745", "ENSG00000002726")
adj_mean_df = mean_and_gc_df[!(row.names(mean_and_gc_df) %in% exclude_genes), ]

# Heatmap and clustering preparation
heatmap_mat = as.matrix(cbind(
  adj_mean_df$MBD_CpG_islands_matched_EPD_TSS_mean_vals_GC_adj_OoE,
  adj_mean_df$short_cf_epd_core_promoters_mean_vals_GC_adj_OoE,
  adj_mean_df$RNA_Gencode_gene_bodies_mean_vals
))

row.names(heatmap_mat) = row.names(adj_mean_df)
colnames(heatmap_mat) = c("MBD_CpGi_adj", "short_cf_core_prom_adj", "RNA_gene_bodies")

# Plot histograms for the distributions
hist(heatmap_mat[, 1], breaks = 100, freq = TRUE)
hist(heatmap_mat[, 2], breaks = 100, freq = TRUE)
hist(heatmap_mat[, 3], breaks = 100, freq = TRUE)

# Log transform RNA gene bodies
heatmap_mat[, 3] = log10(heatmap_mat[, 3])
heatmap_mat = heatmap_mat[complete.cases(heatmap_mat), ]

# Prepare annotation matrix
anno_mat = heatmap_mat[, 3]
heatmap_mat = heatmap_mat[, -3]

# Define color functions
col_fun = colorRamp2(c(quantile(heatmap_mat, c(0.02)), quantile(heatmap_mat, c(0.50)), quantile(heatmap_mat, c(0.98))), c("blue", "gray95", "red"))
col_fun_anno = colorRamp2(c(-2, 0, 4), c("green", "gray95", "purple4"))

# Heatmap plotting
set.seed(123)  # Ensure reproducibility
HM = Heatmap(heatmap_mat,
              cluster_columns = FALSE,
              cluster_rows = TRUE,
              show_row_names = FALSE,
              column_names_rot = 45,
              show_row_dend = TRUE,
              row_split = 10,
              border = TRUE,
              col = col_fun,
              name = "GC adjust",
              right_annotation = rowAnnotation("RNA(logTPM)" = anno_mat, col = list("RNA(logTPM)" = col_fun_anno), border = TRUE),
              clustering_distance_rows = "euclidean",
              clustering_method_rows = "ward.D")

pdf(paste0(base_path, "/heatmap_MBD_epd_CpGi_adj_short_cf_epd_core_prom_adj_RNA-anno.pdf"), width = 4, height = 6)
draw(HM)
dev.off()

# Extract cluster assignment of genes
cluster_assign = data.frame(GeneID = character(), Cluster = character())
for (i in seq_along(row_order(HM))) {
  cluster_genes = row.names(heatmap_mat[row_order(HM)[[i]], ])
  cluster_assign = rbind(cluster_assign, data.frame(GeneID = cluster_genes, Cluster = paste("cluster", i, sep = "")))
}

# Mean annotation values per cluster
anno_df = data.frame(anno_mat)
anno_df$GeneID = rownames(anno_df)
cluster_and_anno = inner_join(cluster_assign, anno_df, by = "GeneID")
write.table(cluster_and_anno, paste0(base_path, "/heatmap_MBD_epd_CpGi_adj_short_cf_epd_core_prom_adj_RNA-anno_cluster_assignment.csv"), quote = FALSE, sep = ",", row.names = FALSE)


# Plot annotation distributions per cluster
cluster_count = table(cluster_and_anno$Cluster)
cluster_and_anno = cluster_and_anno[cluster_and_anno$Cluster %in% names(cluster_count[cluster_count >= 200]), ] #remove cluster with less than x genes

pdf(paste0(base_path, "/heatmap_annotation_distribution_per_cluster.pdf"), width = 4, height = 4)
ggplot(cluster_and_anno, aes(x = anno_mat, y = fct_reorder(.f = Cluster, .x = anno_mat, .fun = median), fill = after_stat(x))) +
   geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01, quantile_lines = TRUE, quantiles = 2) +
   scale_fill_gradient2(low = "green", mid = "gray95", high = "purple4", midpoint = 0) +
   labs(x = 'Gene expression [log10(TPM+0.01)]', y = "Cluster") +
   lims(x = c(-2,4)) +
   theme_clean() +
   theme(
     legend.position="none",
     panel.spacing = unit(0.1, "lines"),
     strip.text.x = element_text(size = 8)
     )
dev.off()


