# Load required libraries
library(Seurat)
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(dplyr)

# Set reproducible seed
set.seed(12)

# Define constants and directories
DATA_DIR <- "data/"
RESULTS_DIR <- "results/"
FIGURES_DIR <- file.path(RESULTS_DIR, "figures/")

# Load the processed Seurat object
breast_cancer <- readRDS(file.path(RESULTS_DIR, "breast_cancer_processed.rds"))

# Convert Seurat object to SingleCellExperiment
sce <- as.SingleCellExperiment(breast_cancer)

# Set Human Primary Cell Atlas as reference dataset
ref <- celldex::BlueprintEncodeData()

# Extract cluster information
clusters <- breast_cancer$seurat_clusters

# Run SingleR to annotate clusters
cluster_annotations <- SingleR(
  test = sce,
  ref = ref,
  labels = ref$label.main,
  clusters = clusters
)

# Create a data frame of cluster annotations
cluster_assignments <- data.frame(
  Cluster = rownames(cluster_annotations),
  CellType = cluster_annotations$labels,
  stringsAsFactors = FALSE
)

# Write results
write.csv(
  cluster_assignments,
  file = file.path(RESULTS_DIR, "cluster_annotations.csv"),
  row.names = FALSE
)

# Print summary
cat("Cluster annotations completed!\n")
print(cluster_assignments)