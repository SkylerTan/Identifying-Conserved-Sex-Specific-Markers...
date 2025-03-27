install.packages("ggVennDiagram")
BiocManager::install("org.Pf.plasmo.db")
BiocManager::install("clusterProfiler")
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
library(ggVennDiagram)
library(org.Pf.plasmo.db)
library(clusterProfiler)


Dogga_data <- readRDS("/storage/group/mul27/default/scRNAseq_share/v3lab.ref.sce.RDS")
assay_data <- as(counts(Dogga_data), "dgCMatrix")
assay_data <- as(logcounts(Dogga_data), "CsparseMatrix")
meta_data <- as.data.frame(colData(Dogga_data))
Dogga_data_seurat <- CreateSeuratObject(
  counts = assay_data,  
  meta.data = meta_data
)
table(Dogga_data_seurat@meta.data$stageHR)

DG_gametocytes_c_l <- subset(Dogga_data_seurat, subset = stageHR %in% c("committed", "early stalk", "late stalk"))
DG_gametocytes_c_l_NF54 <- subset(DG_gametocytes_c_l, subset = ident %in% c("NF54"))
DG_gametocytes_c_l_NF54 <- NormalizeData(DG_gametocytes_c_l_NF54)
DG_gametocytes_c_l_NF54 <- FindVariableFeatures(DG_gametocytes_c_l_NF54, nfeatures = 1000)
DG_gametocytes_c_l_NF54 <- ScaleData(DG_gametocytes_c_l_NF54)
DG_gametocytes_c_l_NF54 <- RunPCA(DG_gametocytes_c_l_NF54, npcs = 50)
CairoPNG("elbow_plot.png")
ElbowPlot(DG_gametocytes_c_l_NF54)
dev.off() 
DG_gametocytes_c_l_NF54 <- FindNeighbors(DG_gametocytes_c_l_NF54, dims = 1:10)
DG_gametocytes_c_l_NF54 <- FindClusters(DG_gametocytes_c_l_NF54, resolution = 0.25)  # Adjust resolution for finer or coarser clusters
DG_gametocytes_c_l_NF54 <- RunUMAP(
  DG_gametocytes_c_l_NF54,
  dims = 1:4,                # Use the first 4 PCs
  n.components = 2,          # 2D UMAP
  min.dist = 0.3,            # Minimum distance between embedded points
  n.neighbors = 50,          # Number of neighbors
  spread = 1,                # Effective scale of embedded points
  local.connectivity = 1,    # Local connectivity
  a = 70,                    # Parameter for embedding
  b = 1                      # Parameter for embedding
)
dg_umap_plot <- DimPlot(DG_gametocytes_c_l_NF54, reduction = "umap", group.by = "stageHR") +
  ggplot2::ggtitle("UMAP with Automatic Clusters")
# Save the plot to a file
ggplot2::ggsave("dg_umap_plot.png", plot = dg_umap_plot, width = 8, height = 6, dpi = 300)

# **Slingshot -----
# Convert Seurat object to SingleCellExperiment
DG_c_l_sce <- as.SingleCellExperiment(DG_gametocytes_c_l_NF54)
start_cluster <- c("committed")
# Run Slingshot using PCA embeddings
DG_c_l_sce <- slingshot(DG_c_l_sce, reducedDim = "PCA", clusterLabels = DG_c_l_sce$stageHR, start.clus = start_cluster)
saveRDS(DG_c_l_sce, file = "DG_c_l_sce.RDS")

# Early to late
DG_gametocytes_e_l <- subset(Dogga_data_seurat, subset = stageHR %in% c("early stalk", "late stalk"))
DG_gametocytes_e_l_NF54 <- subset(DG_gametocytes_e_l, subset = ident %in% c("NF54"))
DG_gametocytes_e_l_NF54 <- NormalizeData(DG_gametocytes_e_l_NF54)
DG_gametocytes_e_l_NF54 <- FindVariableFeatures(DG_gametocytes_e_l_NF54, nfeatures = 1000)
DG_gametocytes_e_l_NF54 <- ScaleData(DG_gametocytes_e_l_NF54)
DG_gametocytes_e_l_NF54 <- RunPCA(DG_gametocytes_e_l_NF54, npcs = 50)
CairoPNG("elbow_plot.png")
ElbowPlot(DG_gametocytes_e_l_NF54)
dev.off() 
DG_gametocytes_e_l_NF54 <- FindNeighbors(DG_gametocytes_e_l_NF54, dims = 1:10)
DG_gametocytes_e_l_NF54 <- FindClusters(DG_gametocytes_e_l_NF54, resolution = 0.25)  # Adjust resolution for finer or coarser clusters
DG_gametocytes_e_l_NF54 <- RunUMAP(
  DG_gametocytes_e_l_NF54,
  dims = 1:4,                # Use the first 4 PCs
  n.components = 2,          # 2D UMAP
  min.dist = 0.3,            # Minimum distance between embedded points
  n.neighbors = 50,          # Number of neighbors
  spread = 1,                # Effective scale of embedded points
  local.connectivity = 1,    # Local connectivity
  a = 70,                    # Parameter for embedding
  b = 1                      # Parameter for embedding
)
dg_umap_plot <- DimPlot(DG_gametocytes_e_l_NF54, reduction = "umap", group.by = "stageHR") +
  ggplot2::ggtitle("UMAP with Automatic Clusters")
# Save the plot to a file
ggplot2::ggsave("dg_umap_plot.png", plot = dg_umap_plot, width = 8, height = 6, dpi = 300)

# **Slingshot -----
# Convert Seurat object to SingleCellExperiment
DG_e_l_sce <- as.SingleCellExperiment(DG_gametocytes_e_l_NF54)
start_cluster <- c("early stalk")
# Run Slingshot using PCA embeddings
DG_e_l_sce <- slingshot(DG_e_l_sce, reducedDim = "PCA", clusterLabels = DG_e_l_sce$stageHR, start.clus = start_cluster)
saveRDS(DG_e_l_sce, file = "DG_e_l_sce.RDS")

DG_c_l_sce <- readRDS("/storage/work/skt5723/Single Cell Gametocyte stuff/DG_c_l_sce.RDS")
DG_e_l_sce <- readRDS("/storage/work/skt5723/Single Cell Gametocyte stuff/DG_e_l_sce.RDS")
#HERE
# Extract UMAP embeddings
DG_c_l_sce_umap_embeddings <- reducedDim(DG_c_l_sce, "UMAP")
DG_e_l_sce_umap_embeddings <- reducedDim(DG_e_l_sce, "UMAP")

# Create a data frame for UMAP embeddings and stageHR
DG_c_l_umap_data <- data.frame(
  UMAP1 = DG_c_l_sce_umap_embeddings[, 1],
  UMAP2 = DG_c_l_sce_umap_embeddings[, 2],
  Stage = DG_c_l_sce$stageHR
)
DG_e_l_umap_data <- data.frame(
  UMAP1 = DG_e_l_sce_umap_embeddings[, 1],
  UMAP2 = DG_e_l_sce_umap_embeddings[, 2],
  Stage = DG_e_l_sce$stageHR
)


# Create a ggplot object
p1 <- ggplot(DG_c_l_umap_data, aes(x = UMAP1, y = UMAP2, color = Stage)) +
  geom_point(size = 1.5, alpha = 0.8) +
  theme_minimal() +
  labs(title = "Slingshot Trajectory on UMAP", x = "UMAP1", y = "UMAP2")
p2 <- ggplot(DG_e_l_umap_data, aes(x = UMAP1, y = UMAP2, color = Stage)) +
  geom_point(size = 1.5, alpha = 0.8) +
  theme_minimal() +
  labs(title = "Slingshot Trajectory on UMAP", x = "UMAP1", y = "UMAP2")


# Add Slingshot trajectories to the plot
DG_c_l_slingshot_data <- SlingshotDataSet(DG_c_l_sce)
for (curve in DG_c_l_slingshot_data@curves) {
  # Extract curve points in PCA space
  curve_points_pca <- curve$s[curve$ord, ]
  
  # Map PCA points to UMAP space using nearest neighbors in PCA space
  pca_embeddings <- reducedDim(DG_c_l_sce, "PCA")  # Get PCA embeddings
  nn_indices <- FNN::get.knnx(pca_embeddings, curve_points_pca, k = 1)$nn.index  # Find nearest neighbors
  curve_points_umap <- DG_c_l_sce_umap_embeddings[nn_indices, ]  # Use UMAP embeddings of nearest neighbors
  
  # Convert to data frame and set column names
  curve_points_umap <- as.data.frame(curve_points_umap)
  colnames(curve_points_umap) <- c("UMAP1", "UMAP2")
  
  # Add trajectory to the plot
  p1 <- p1 + geom_path(data = curve_points_umap, aes(x = UMAP1, y = UMAP2), color = "black", linewidth = 1)
}

# Save the UMAP plot with trajectories
ggsave("DG_c_l_sce_slingshot_trajectory_umap.png", plot = p1, width = 8, height = 6, dpi = 300)
DG_c_l_sce_lineages <- slingshot::slingLineages(DG_c_l_sce)
print(DG_c_l_sce_lineages)

# Extract pseudotime values
DG_c_l_sce_pseudotime <- slingshot::slingPseudotime(DG_c_l_sce)

# Prepare data for tradeSeq
DG_c_l_sce_counts <- counts(DG_c_l_sce)  # Raw counts

DG_c_l_sce_cell_weights <- slingshot::slingCurveWeights(DG_c_l_sce)  # Cell weights

# Evaluate optimal number of knots
set.seed(123)  # For reproducibility

# Evaluate knots
DG_c_l_sce_knots <- evaluateK(
  counts = DG_c_l_sce_counts,
  pseudotime = DG_c_l_sce_pseudotime,
  cellWeights = DG_c_l_sce_cell_weights,
  k = 3:10,  # Test 3 to 10 knots
  nGenes = 200  # Subsample 200 genes for evaluation
)
print(DG_c_l_sce_knots)  # Choose the optimal number of knots (e.g., 8)

#EL PLOTS
# Add Slingshot trajectories to the plot
DG_e_l_slingshot_data <- SlingshotDataSet(DG_e_l_sce)
for (curve in DG_e_l_slingshot_data@curves) {
  # Extract curve points in PCA space
  curve_points_pca <- curve$s[curve$ord, ]
  
  # Map PCA points to UMAP space using nearest neighbors in PCA space
  pca_embeddings <- reducedDim(DG_e_l_sce, "PCA")  # Get PCA embeddings
  nn_indices <- FNN::get.knnx(pca_embeddings, curve_points_pca, k = 1)$nn.index  # Find nearest neighbors
  curve_points_umap <- DG_e_l_sce_umap_embeddings[nn_indices, ]  # Use UMAP embeddings of nearest neighbors
  
  # Convert to data frame and set column names
  curve_points_umap <- as.data.frame(curve_points_umap)
  colnames(curve_points_umap) <- c("UMAP1", "UMAP2")
  
  # Add trajectory to the plot
  p2 <- p2 + geom_path(data = curve_points_umap, aes(x = UMAP1, y = UMAP2), color = "black", linewidth = 1)
}

# Save the UMAP plot with trajectories
ggsave("DG_e_l_sce_slingshot_trajectory_umap.png", plot = p2, width = 8, height = 6, dpi = 300)
DG_e_l_sce_lineages <- slingshot::slingLineages(DG_e_l_sce)
print(DG_e_l_sce_lineages)

# Extract pseudotime values
DG_e_l_sce_pseudotime <- slingshot::slingPseudotime(DG_e_l_sce)

# Prepare data for tradeSeq
DG_e_l_sce_counts <- counts(DG_e_l_sce)  # Raw counts

DG_e_l_sce_cell_weights <- slingshot::slingCurveWeights(DG_e_l_sce)  # Cell weights

# Evaluate optimal number of knots
set.seed(123)  # For reproducibility

# Evaluate knots
DG_e_l_sce_knots <- evaluateK(
  counts = DG_e_l_sce_counts,
  pseudotime = DG_e_l_sce_pseudotime,
  cellWeights = DG_e_l_sce_cell_weights,
  k = 3:10,  # Test 3 to 10 knots
  nGenes = 200  # Subsample 200 genes for evaluation
)
print(DG_e_l_sce_knots)  # Choose the optimal number of knots (e.g., 8)

# Transpose counts matrix and convert to matrix
DG_c_l_sce_counts <- as.matrix(DG_c_l_sce_counts)
DG_e_l_sce_counts <- as.matrix(DG_e_l_sce_counts)

DG_c_l_sce_pseudotime <- DG_c_l_sce_pseudotime[, 1, drop = FALSE]
DG_e_l_sce_pseudotime <- DG_e_l_sce_pseudotime[, 1, drop = FALSE]

DG_c_l_sce_cell_weights <- DG_c_l_sce_cell_weights[, 1, drop = FALSE]
DG_e_l_sce_cell_weights <- DG_e_l_sce_cell_weights[, 1, drop = FALSE]

DG_c_l_sce_pseudotime[is.na(DG_c_l_sce_pseudotime)] <- 0
DG_c_l_sce_cell_weights[is.na(DG_c_l_sce_cell_weights)] <- 0
DG_e_l_sce_pseudotime[is.na(DG_e_l_sce_pseudotime)] <- 0
DG_e_l_sce_cell_weights[is.na(DG_e_l_sce_cell_weights)] <- 0

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

check_input_dimensions(DG_c_l_sce_counts, DG_c_l_sce_pseudotime, DG_c_l_sce_cell_weights)
check_input_dimensions(DG_e_l_sce_counts, DG_e_l_sce_pseudotime, DG_e_l_sce_cell_weights)

DG_c_l_sce_gam_fit <- fitGAM(
  counts = DG_c_l_sce_counts,
  pseudotime = as.vector(DG_c_l_sce_pseudotime),
  cellWeights = DG_c_l_sce_cell_weights,
  nknots = 8
)
DG_e_l_sce_gam_fit <- fitGAM(
  counts = DG_e_l_sce_counts,
  pseudotime = as.vector(DG_e_l_sce_pseudotime),
  cellWeights = DG_e_l_sce_cell_weights,
  nknots = 8
)

# Save each GAM fit as an RDS file
saveRDS(DG_c_l_sce_gam_fit, file = "DG_c_l_sce_gam_fit.rds")
saveRDS(DG_e_l_sce_gam_fit, file = "DG_e_l_sce_gam_fit.rds")

# To load them back later:
DG_c_l_sce_gam_fit <- readRDS("DG_c_l_sce_gam_fit.rds")
DG_e_l_sce_gam_fit <- readRDS("DG_e_l_sce_gam_fit.rds")


degs_DG_c_l_sce <- startVsEndTest(
  models = DG_c_l_sce_gam_fit,
  global = TRUE,
  lineages = 1,
  l2fc = 1
)
degs_DG_e_l_sce <- startVsEndTest(
  models = DG_e_l_sce_gam_fit,
  global = TRUE,
  lineages = 1,
  l2fc = 1
)

write.csv(degs_DG_e_l_sce, file = "degs_DG_e_l_sce.csv", row.names = TRUE)

de <- degs_DG_c_l_sce
de$logfc <- de$logFClineage1
de$wald <- de$waldStat
de$pvalue <- de$pvalue_lineage1
de$diffexpressed <- "NO"
de$diffexpressed[de$logfc > 1 & de$wald > 20] <- "UP"
de$diffexpressed[de$logfc < -1 & de$wald > 20] <- "DOWN"
def <- de[which(de$diffexpressed != "NO"),]
def$geneid2 <- rownames(def)
committed_latestalk <- def


de <- degs_DG_e_l_sce
de$logfc <- de$logFClineage1
de$wald <- de$waldStat
de$pvalue <- de$pvalue_lineage1
de$diffexpressed <- "NO"
de$diffexpressed[de$logfc > 1 & de$wald > 5] <- "UP"
de$diffexpressed[de$logfc < -1 & de$wald > 5] <- "DOWN"
def <- de[which(de$diffexpressed != "NO"),]
def$geneid2 <- rownames(def)
earlystalk_latestalk <- def

fig_s5_B <- ggVennDiagram(list(committed_latestalk$geneid2,earlystalk_latestalk$geneid2), label_alpha = 0, label = "count", label_size = 7,
                          category.names = c("committed_latestalk","earlystalk_latestalk")
) +   ggplot2::scale_fill_gradient(low="white",high = "thistle")
ggsave("fig_s5_B.png", plot = fig_s5_B, width = 8, height = 6, dpi = 300)



stalk_both_60 <- data.frame(geneid2 = union(committed_latestalk$geneid2,earlystalk_latestalk$geneid2),order = c(1:length(union(committed_latestalk$geneid2,earlystalk_latestalk$geneid2))))
stalk_both_60 <- left_join(stalk_both_60,pfdb56, by = "geneid2")


### fig S5 D

geneList = keys(org.Pf.plasmo.db)
genelist <- stalk_both_60
genelist$geneid2 <- gsub("-", "_", genelist$geneid2)

# module3
enr <- enrichGO(
  genelist$geneid2,
  universe = geneList,
  OrgDb = 'org.Pf.plasmo.db',
  pvalueCutoff = 0.05,
  keyType = "SYMBOL",  
  ont = "BP", 
  qvalueCutoff = 0.2
)

fig_s5_D <- clusterProfiler::dotplot(enr)
go_res <- enr@result

ggsave("fig_s5_D.png", fig_s5_D ,
       width=7, height=2.7, dpi=300, bg = "transparent")
