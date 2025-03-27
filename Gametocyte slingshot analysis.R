BiocManager::install("slingshot")
BiocManager::install("DelayedMatrixStats")
aBiocManager::install("tradeSeq", force = TRUE)
detach("package:tradeSeq", unload = TRUE)
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(tradeSeq)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(SingleCellExperiment)
library(Matrix)
library(Cairo)
library(slingshot)
library(DelayedMatrixStats)

setwd("/storage/work/skt5723/Single Cell Gametocyte stuff")
Sys.setenv(DISPLAY=":0")  # Prevent X11 errors
options(device = "pdf")  
Dogga_data <- readRDS("/storage/group/mul27/default/scRNAseq_share/v3lab.ref.sce.RDS")
assay_data <- as(counts(Dogga_data), "dgCMatrix")
assay_data <- as(logcounts(Dogga_data), "CsparseMatrix")
meta_data <- as.data.frame(colData(Dogga_data))
Dogga_data_seurat <- CreateSeuratObject(
  counts = assay_data,  
  meta.data = meta_data
)

DG_gametocytes <- subset(Dogga_data_seurat, subset = stageLR %in% c("gametocyte (developing)", "gametocyte (male)", "gametocyte (female)"))
DG_gametocytes_NF54 <- subset(DG_gametocytes, subset = ident %in% c("NF54"))
DG_gametocytes_NF54 <- NormalizeData(DG_gametocytes_NF54)
DG_gametocytes_NF54 <- FindVariableFeatures(DG_gametocytes_NF54, nfeatures = 1000)
DG_gametocytes_NF54 <- ScaleData(DG_gametocytes_NF54)
DG_gametocytes_NF54 <- RunPCA(DG_gametocytes_NF54, npcs = 50)
CairoPNG("elbow_plot.png")
ElbowPlot(DG_gametocytes_NF54)
dev.off() 
DG_gametocytes_NF54 <- FindNeighbors(DG_gametocytes_NF54, dims = 1:10)
DG_gametocytes_NF54 <- FindClusters(DG_gametocytes_NF54, resolution = 0.25)  # Adjust resolution for finer or coarser clusters
DG_gametocytes_NF54 <- RunUMAP(
  DG_gametocytes_NF54,
  dims = 1:4,                # Use the first 4 PCs
  n.components = 2,          # 2D UMAP
  min.dist = 0.3,            # Minimum distance between embedded points
  n.neighbors = 50,          # Number of neighbors
  spread = 1,                # Effective scale of embedded points
  local.connectivity = 1,    # Local connectivity
  a = 70,                    # Parameter for embedding
  b = 1                      # Parameter for embedding
)
dg_umap_plot <- DimPlot(DG_gametocytes_NF54, reduction = "umap", group.by = "stageHR") +
  ggplot2::ggtitle("UMAP with Automatic Clusters")
dg_umap_plot_strain <- DimPlot(DG_gametocytes, reduction = "umap", group.by = "ident") +
  ggplot2::ggtitle("UMAP with Automatic Clusters")
# Save the plot to a file
ggplot2::ggsave("dg_umap_plot.png", plot = dg_umap_plot, width = 8, height = 6, dpi = 300)

# **Slingshot -----
# Convert Seurat object to SingleCellExperiment
DG_sce <- as.SingleCellExperiment(DG_gametocytes_NF54)
start_cluster <- c("committed")
# Run Slingshot using PCA embeddings
DG_sce <- slingshot(DG_sce, reducedDim = "PCA", clusterLabels = DG_sce$stageHR, start.clus = start_cluster)
saveRDS(DG_sce, file = "DG_sce.RDS")


DG_sce <- readRDS("/storage/work/skt5723/Single Cell Gametocyte stuff/DG_sce.RDS")

# Extract UMAP embeddings
umap_embeddings <- reducedDim(DG_sce, "UMAP")

# Create a data frame for UMAP embeddings and stageHR
umap_data <- data.frame(
  UMAP1 = umap_embeddings[, 1],
  UMAP2 = umap_embeddings[, 2],
  Stage = DG_sce$stageHR
)


# Create a ggplot object
p <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = Stage)) +
  geom_point(size = 1.5, alpha = 0.8) +
  theme_minimal() +
  labs(title = "Slingshot Trajectory on UMAP", x = "UMAP1", y = "UMAP2")


# Add Slingshot trajectories to the plot
slingshot_data <- SlingshotDataSet(DG_sce)
for (curve in slingshot_data@curves) {
  # Extract curve points in PCA space
  curve_points_pca <- curve$s[curve$ord, ]
  
  # Map PCA points to UMAP space using nearest neighbors in PCA space
  pca_embeddings <- reducedDim(DG_sce, "PCA")  # Get PCA embeddings
  nn_indices <- FNN::get.knnx(pca_embeddings, curve_points_pca, k = 1)$nn.index  # Find nearest neighbors
  curve_points_umap <- umap_embeddings[nn_indices, ]  # Use UMAP embeddings of nearest neighbors
  
  # Convert to data frame and set column names
  curve_points_umap <- as.data.frame(curve_points_umap)
  colnames(curve_points_umap) <- c("UMAP1", "UMAP2")
  
  # Add trajectory to the plot
  p <- p + geom_path(data = curve_points_umap, aes(x = UMAP1, y = UMAP2), color = "black", linewidth = 1)
}

# Save the UMAP plot with trajectories
ggsave("slingshot_trajectory_umap.png", plot = p, width = 8, height = 6, dpi = 300)
lineages <- slingshot::slingLineages(DG_sce)
print(lineages)

# Extract pseudotime values
pseudotime <- slingshot::slingPseudotime(DG_sce)

# Prepare data for tradeSeq
counts <- counts(DG_sce)  # Raw counts

cell_weights <- slingshot::slingCurveWeights(DG_sce)  # Cell weights

# Evaluate optimal number of knots
set.seed(123)  # For reproducibility

# Evaluate knots
DG_knots <- evaluateK(
  counts = counts,
  pseudotime = pseudotime,
  cellWeights = cell_weights,
  k = 3:10,  # Test 3 to 10 knots
  nGenes = 200  # Subsample 200 genes for evaluation
)
print(knots)  # Choose the optimal number of knots (e.g., 8)

# Transpose counts matrix and convert to matrix
counts <- t(counts)
counts <- as.matrix(counts)

# Extract male and female lineage pseudotime
pseudotime_male <- pseudotime[, 1, drop = FALSE]
pseudotime_female <- pseudotime[, 2, drop = FALSE]

# Extract male and female lineage cell weights
cell_weights_male <- cell_weights[, 1, drop = FALSE]
cell_weights_female <- cell_weights[, 2, drop = FALSE]

# Replace NA values in pseudotime and cell_weights with 0
pseudotime_male[is.na(pseudotime_male)] <- 0
pseudotime_female[is.na(pseudotime_female)] <- 0
cell_weights_male[is.na(cell_weights_male)] <- 0
cell_weights_female[is.na(cell_weights_female)] <- 0

# Create logical indices for cells with non-zero weights
male_cells <- cell_weights_male > 0
female_cells <- cell_weights_female > 0

# Subset counts for male and female lineages
counts_male <- counts[male_cells, , drop = FALSE]
counts_male <- t(counts_male)
counts_female <- counts[female_cells, , drop = FALSE]
counts_female <- t(counts_female)

# Subset pseudotime and cell_weights for male and female lineages
pseudotime_male <- pseudotime_male[male_cells, , drop = FALSE]
pseudotime_female <- pseudotime_female[female_cells, , drop = FALSE]
cell_weights_male <- cell_weights_male[male_cells, , drop = FALSE]
cell_weights_female <- cell_weights_female[female_cells, , drop = FALSE]

# Check dimensions for male
cat("Male counts:", dim(counts_male), "\n")
cat("Male pseudotime:", dim(pseudotime_male), "\n")
cat("Male cell weights:", dim(cell_weights_male), "\n")

# Check dimensions for female
cat("Female counts:", dim(counts_female), "\n")
cat("Female pseudotime:", dim(pseudotime_female), "\n")
cat("Female cell weights:", dim(cell_weights_female), "\n")

# Check for rows with zero counts
zero_sum_rows <- rowSums(counts_male) == 0
sum(zero_sum_rows)

# Function to check input dimensions
check_input_dimensions <- function(counts, pseudotime, cellWeights) {
  if (!is.null(dim(pseudotime)) & !is.null(dim(cellWeights))) {
    if (!identical(nrow(pseudotime), ncol(counts))) {
      stop("❌ ERROR: 'pseudotime' and count matrix must have equal number of cells. ",
           "Found: ", nrow(pseudotime), " vs ", ncol(counts))
    }
    if (!identical(nrow(cellWeights), ncol(counts))) {
      stop("❌ ERROR: 'cellWeights' and count matrix must have equal number of cells. ",
           "Found: ", nrow(cellWeights), " vs ", ncol(counts))
    }
  } else {
    warning("⚠️ Warning: 'pseudotime' or 'cellWeights' might be NULL. Ensure they are correctly initialized.")
  }
  message("✅ Input dimensions check passed! Ready for fitGAM().")
}

# Example Usage
check_input_dimensions(counts_male, pseudotime_male, cell_weights_male)

# Adjust small cell weights
cell_weights_male[cell_weights_male <= 0] <- 1e-6

# Fit GAM models
gam_fit_male <- fitGAM(
  counts = counts_male,
  pseudotime = as.vector(pseudotime_male),
  cellWeights = cell_weights_male,
  nknots = 8
)

# Test for DEGs in male lineage
degs_male <- startVsEndTest(
  models = gam_fit_male,
  global = TRUE,
  lineages = 1,
  l2fc = 1
)

# Filter significant DEGs
degs_male <- degs_male[degs_male$waldStat > 60 & abs(degs_male$log2fc) > 1, ]

# Save DEGs
write.csv(degs_male, file = "degs_male_lineage.csv", row.names = TRUE)

# Sort and extract top genes
degs_male_sorted <- degs_male[order(abs(degs_male$waldStat), decreasing = TRUE), ]
top_genes_male <- rownames(degs_male_sorted)[1:1000]

# Heatmap setup
heatmap_data_male <- gam_fit_male[top_genes_male, ]
heatmap_data_male <- heatmap_data_male[, order(pseudotime_male[, 1])]
heatmap_data_male_matrix <- assay(heatmap_data_male, "counts")

# Create and save heatmap
heatmap_male <- Heatmap(
  heatmap_data_male_matrix,
  name = "Expression",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  column_title = "Cells ordered by pseudotime",
  row_title = "Top Genes"
)

CairoPNG("heatmap_male.png", width = 800, height = 800)
draw(heatmap_male)
dev.off()

ggsave("heatmap_male_lineage.png", plot = heatmap_male, width = 8, height = 6, dpi = 300)
