# ===============================================================================
# Single-Cell RNA-seq Analysis: Human Melanoma Dataset
# Author: Jason Summers
# Date: 6/27/2025
# Description: Comprehensive preprocessing, quality control, and clustering analysis 
#              of 10k melanoma cells using Seurat workflow
# ===============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(patchwork)
  library(ggplot2)
  library(cowplot)
})

# Set reproducible seed
set.seed(12)

# Define constants and directories
DATA_DIR <- "data/"
RESULTS_DIR <- "results/"
FIGURES_DIR <- file.path(RESULTS_DIR, "figures/")

# Create output directories if they don't exist
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)

# ===============================================================================
# DATA LOADING AND INITIAL SETUP
# ===============================================================================

cat("Loading melanoma dataset...\n")
melanoma_data <- Read10X_h5(file.path(DATA_DIR, "sc5p_v2_hs_melanoma_10k_filtered_feature_bc_matrix.h5"))

# Create Seurat object with initial filtering
melanoma <- CreateSeuratObject(
  counts = melanoma_data,
  project = "Melanoma_10X",
  min.cells = 3,      # Filter genes expressed in < 3 cells
  min.features = 200  # Filter cells with < 200 genes
)

cat(sprintf("Initial dataset: %d cells x %d genes\n", 
            ncol(melanoma), nrow(melanoma)))

# ===============================================================================
# QUALITY CONTROL METRICS
# ===============================================================================

cat("Calculating QC metrics...\n")

# Calculate mitochondrial gene percentage
melanoma[["percent.mt"]] <- PercentageFeatureSet(melanoma, pattern = "^MT-")

# Add cell complexity (log10 genes per UMI)
melanoma[["log10GenesPerUMI"]] <- log10(melanoma[["nFeature_RNA"]]) / log10(melanoma[["nCount_RNA"]])

# Visualize QC metrics before filtering
qc_metrics_before <- c("nFeature_RNA", "nCount_RNA", "percent.mt")

vln_plot_before_qc <- VlnPlot(
  melanoma,
  features = qc_metrics_before,
  ncol = 3,
  pt.size = 0.1
) + 
  plot_annotation(title = "QC Metrics Before Filtering",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 16)))

# Save QC plots before filtering
ggsave(
  file.path(FIGURES_DIR, "01_qc_metrics_before_filtering.png"),
  vln_plot_before_qc,
  width = 12,
  height = 10,
  dpi = 300,
  bg = "white"
)

# Create scatter plots to identify outliers
scatter_plot <- FeatureScatter(melanoma, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  FeatureScatter(melanoma, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

ggsave(
  file.path(FIGURES_DIR, "01_feature_scatter_before_filtering.png"),
  scatter_plot,
  width = 12,
  height = 6,
  dpi = 300,
  bg = "white"
)

# ===============================================================================
# CELL FILTERING BASED ON QC METRICS
# ===============================================================================

cat("Applying quality control filters...\n")

# Define filtering thresholds
min_features <- 200
max_features <- 6000
max_mt_percent <- 10

# Apply filters
melanoma_filtered <- subset(
  melanoma,
  subset = nFeature_RNA > min_features & 
    nFeature_RNA < max_features & 
    percent.mt < max_mt_percent
)

cat(sprintf("After filtering: %d cells x %d genes\n", 
            ncol(melanoma_filtered), nrow(melanoma_filtered)))
cat(sprintf("Removed %d cells (%.1f%%)\n", 
            ncol(melanoma) - ncol(melanoma_filtered),
            100 * (ncol(melanoma) - ncol(melanoma_filtered)) / ncol(melanoma)))

# Visualize QC metrics after filtering
vln_plot_after_qc <- VlnPlot(
  melanoma_filtered,
  features = qc_metrics_before,
  ncol = 3,
  pt.size = 0.1
) + 
  plot_annotation(title = "QC Metrics After Filtering",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 16)))

ggsave(
  file.path(FIGURES_DIR, "02_qc_metrics_after_filtering.png"),
  vln_plot_after_qc,
  width = 12,
  height = 10,
  dpi = 300,
  bg = "white"
)

# Update main object
melanoma <- melanoma_filtered
rm(melanoma_filtered)

# ===============================================================================
# NORMALIZATION AND FEATURE SELECTION
# ===============================================================================

cat("Normalizing data...\n")

# Log-normalize the data
melanoma <- NormalizeData(
  melanoma,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

cat("Finding variable features...\n")

# Identify highly variable features
melanoma <- FindVariableFeatures(
  melanoma,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

# Get top 10 most variable genes for labeling
top10_variable <- head(VariableFeatures(melanoma), 10)

# Plot variable features
variable_features_plot <- VariableFeaturePlot(melanoma) +
  ggtitle("Highly Variable Genes (Top 2000)") +
  theme(plot.title = element_text(hjust = 0.5))

variable_features_labeled <- LabelPoints(
  plot = variable_features_plot,
  points = top10_variable,
  repel = TRUE,
  max.overlaps = 20
)

ggsave(
  file.path(FIGURES_DIR, "03_variable_features.png"),
  variable_features_labeled,
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)

# ===============================================================================
# DATA SCALING AND PCA
# ===============================================================================

cat("Scaling data...\n")

# Scale data for PCA
melanoma <- ScaleData(
  melanoma,
  features = rownames(melanoma),
  verbose = FALSE
)

cat("Running PCA...\n")

# Perform PCA
melanoma <- RunPCA(
  melanoma,
  features = VariableFeatures(melanoma),
  verbose = FALSE
)

# Visualize PCA results
pca_dim_plot <- DimPlot(melanoma, reduction = "pca") +
  ggtitle("PCA: PC1 vs PC2") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  file.path(FIGURES_DIR, "04_pca_plot.png"),
  pca_dim_plot,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)

# Create PCA heatmap
png(
  file.path(FIGURES_DIR, "04_pca_heatmap.png"),
  width = 14,
  height = 12,
  units = "in",
  res = 300
)
DimHeatmap(
  melanoma,
  dims = 1:12,
  cells = 500,
  balanced = TRUE,
  reduction = "pca",
  ncol = 3
)
dev.off()

# Determine dimensionality with elbow plot
elbow_plot <- ElbowPlot(melanoma, ndims = 30) +
  ggtitle("PCA Elbow Plot") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = 15, linetype = "dashed", color = "red", alpha = 0.7)

ggsave(
  file.path(FIGURES_DIR, "04_elbow_plot.png"),
  elbow_plot,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)

# ===============================================================================
# CLUSTERING ANALYSIS
# ===============================================================================

cat("Performing clustering analysis...\n")

# Build SNN graph and find clusters
n_pcs <- 15  # Based on elbow plot

melanoma <- FindNeighbors(
  melanoma,
  dims = 1:n_pcs,
  verbose = FALSE
)

# Perform clustering at resolution 0.5
melanoma <- FindClusters(
  melanoma,
  resolution = 0.5,
  random.seed = 12,
  verbose = FALSE
)

cat(sprintf("Found %d clusters at resolution 0.5\n", 
            length(unique(Idents(melanoma)))))

# ===============================================================================
# DIMENSIONALITY REDUCTION - UMAP
# ===============================================================================

cat("Running UMAP...\n")

# Run UMAP
melanoma <- RunUMAP(
  melanoma,
  dims = 1:n_pcs,
  seed.use = 12,
  verbose = FALSE
)

# Create UMAP plot
umap_clusters <- DimPlot(
  melanoma,
  reduction = "umap",
  label = TRUE,
  label.size = 6,
  pt.size = 0.5
) + 
  ggtitle("Cell Clusters (UMAP)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(override.aes = list(size = 3)))

ggsave(
  file.path(FIGURES_DIR, "05_umap_clusters.png"),
  umap_clusters,
  width = 12,
  height = 10,
  dpi = 300,
  bg = "white"
)

# ===============================================================================
# SAVE PROCESSED DATA
# ===============================================================================

cat("Saving processed Seurat object...\n")

# Save the processed Seurat object
saveRDS(melanoma, file = file.path(RESULTS_DIR, "melanoma_processed.rds"))

# Save session info
writeLines(capture.output(sessionInfo()), 
           file.path(RESULTS_DIR, "session_info.txt"))

# Create summary statistics
summary_stats <- data.frame(
  Metric = c("Initial cells", "Initial genes", "Cells after QC", "Genes after filtering",
             "Highly variable genes", "PCs used", "Final clusters"),
  Value = c(ncol(Read10X_h5(file.path(DATA_DIR, "sc5p_v2_hs_melanoma_10k_filtered_feature_bc_matrix.h5"))),
            nrow(Read10X_h5(file.path(DATA_DIR, "sc5p_v2_hs_melanoma_10k_filtered_feature_bc_matrix.h5"))),
            ncol(melanoma),
            nrow(melanoma),
            length(VariableFeatures(melanoma)),
            n_pcs,
            length(unique(Idents(melanoma))))
)

write.csv(summary_stats, 
          file.path(RESULTS_DIR, "analysis_summary.csv"), 
          row.names = FALSE)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Results saved to:", RESULTS_DIR, "\n")
cat("Figures saved to:", FIGURES_DIR, "\n")
cat("Processed Seurat object:", file.path(RESULTS_DIR, "melanoma_processed.rds"), "\n")