# ===============================================================================
# Single-Cell RNA-seq Analysis: Human Breast Cancer Dataset
# Author: Jason Summers
# Date: 6/27/2025
# Description: Comprehensive preprocessing, quality control, and clustering analysis 
#              of 10k breast cancer cells using Seurat workflow
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

cat("Loading breast cancer dataset...\n")
breast_cancer_data <- Read10X(data.dir = file.path(DATA_DIR, "raw_feature_bc_matrix"))

# Create Seurat object with initial filtering
breast_cancer <- CreateSeuratObject(
  counts = breast_cancer_data,
  project = "breast_cancer_10X",
  min.cells = 3,      # Filter genes expressed in < 3 cells
  min.features = 200  # Filter cells with < 200 genes
)

cat(sprintf("Initial dataset: %d cells x %d genes\n", 
            ncol(breast_cancer), nrow(breast_cancer)))

# ===============================================================================
# QUALITY CONTROL METRICS
# ===============================================================================

cat("Calculating QC metrics...\n")

# Calculate mitochondrial gene percentage
breast_cancer[["percent.mt"]] <- PercentageFeatureSet(breast_cancer, pattern = "^MT-")

# Visualize QC metrics before filtering
qc_metrics_before <- c("nFeature_RNA", "nCount_RNA", "percent.mt")

vln_plot_before_qc <- VlnPlot(
  breast_cancer,
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
scatter_plot <- FeatureScatter(breast_cancer, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  FeatureScatter(breast_cancer, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

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
max_mt_percent <- 25

# Apply filters
breast_cancer_filtered <- subset(
  breast_cancer,
  subset = nFeature_RNA > min_features & 
    nFeature_RNA < max_features & 
    percent.mt < max_mt_percent
)

cat(sprintf("After filtering: %d cells x %d genes\n", 
            ncol(breast_cancer_filtered), nrow(breast_cancer_filtered)))
cat(sprintf("Removed %d cells (%.1f%%)\n", 
            ncol(breast_cancer) - ncol(breast_cancer_filtered),
            100 * (ncol(breast_cancer) - ncol(breast_cancer_filtered)) / ncol(breast_cancer)))

# Visualize QC metrics after filtering
vln_plot_after_qc <- VlnPlot(
  breast_cancer_filtered,
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
breast_cancer <- breast_cancer_filtered
rm(breast_cancer_filtered)

# ===============================================================================
# NORMALIZATION AND FEATURE SELECTION
# ===============================================================================

cat("Normalizing data...\n")

# Log-normalize the data
breast_cancer <- NormalizeData(
  breast_cancer,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

cat("Finding variable features...\n")

# Identify highly variable features
breast_cancer <- FindVariableFeatures(
  breast_cancer,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

# Get top 10 most variable genes for labeling
top10_variable <- head(VariableFeatures(breast_cancer), 10)

# Plot variable features
variable_features_plot <- VariableFeaturePlot(breast_cancer) +
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
breast_cancer <- ScaleData(
  breast_cancer,
  features = rownames(breast_cancer),
  verbose = FALSE
)

cat("Running PCA...\n")

# Perform PCA
breast_cancer <- RunPCA(
  breast_cancer,
  features = VariableFeatures(breast_cancer),
  verbose = FALSE
)

# Visualize PCA results
pca_dim_plot <- DimPlot(breast_cancer, reduction = "pca") +
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
  breast_cancer,
  dims = 1:12,
  cells = 500,
  balanced = TRUE,
  reduction = "pca",
  ncol = 3
)
dev.off()

# Determine dimensionality with elbow plot
elbow_plot <- ElbowPlot(breast_cancer, ndims = 30) +
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

breast_cancer <- FindNeighbors(
  breast_cancer,
  dims = 1:n_pcs,
  verbose = FALSE
)

# Perform clustering at resolution 0.5
breast_cancer <- FindClusters(
  breast_cancer,
  resolution = 0.5,
  random.seed = 12,
  verbose = FALSE
)

cat(sprintf("Found %d clusters at resolution 0.5\n", 
            length(unique(Idents(breast_cancer)))))

# ===============================================================================
# DIMENSIONALITY REDUCTION - UMAP
# ===============================================================================

cat("Running UMAP...\n")

# Run UMAP
breast_cancer <- RunUMAP(
  breast_cancer,
  dims = 1:n_pcs,
  seed.use = 12,
  verbose = FALSE
)

# Create UMAP plot
umap_clusters <- DimPlot(
  breast_cancer,
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
saveRDS(breast_cancer, file = file.path(RESULTS_DIR, "breast_cancer_processed.rds"))

# Save session info
writeLines(capture.output(sessionInfo()), 
           file.path(RESULTS_DIR, "session_info.txt"))

# Create summary statistics
summary_stats <- data.frame(
  Metric = c("Initial cells", "Initial genes", "Cells after QC", "Genes after filtering",
             "Highly variable genes", "PCs used", "Final clusters"),
  Value = c(ncol(Read10X(data.dir = "data/raw_feature_bc_matrix/")),
            nrow(Read10X(data.dir = "data/raw_feature_bc_matrix/")),
            ncol(breast_cancer),
            nrow(breast_cancer),
            length(VariableFeatures(breast_cancer)),
            n_pcs,
            length(unique(Idents(breast_cancer))))
)

write.csv(summary_stats, 
          file.path(RESULTS_DIR, "analysis_summary.csv"), 
          row.names = FALSE)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Results saved to:", RESULTS_DIR, "\n")
cat("Figures saved to:", FIGURES_DIR, "\n")
cat("Processed Seurat object:", file.path(RESULTS_DIR, "breast_cancer_processed.rds"), "\n")
