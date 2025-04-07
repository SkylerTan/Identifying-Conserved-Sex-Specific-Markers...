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
table(Dogga_data_seurat@meta.data$stageHR)
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
table(DG_gametocytes_NF54@meta.data$stageHR)
DG_gametocytes_NF54 <- RunUMAP(
  DG_gametocytes_NF54,
  dims = 1:4,               
  n.components = 2,         
  min.dist = 0.3,          
  n.neighbors = 50,         
  spread = 1,               
  local.connectivity = 1,  
  a = 70,                   
  b = 1                     
)
stage_colors <- c(
  "early trophozoite" = "#fa8f87",
  "late trophozoite" = "#d46d66",  
  "early schizont" = "#7B87B6",
  "late schizont" = "#3D517E",
  "early ring" = "#C3E796",
  "late ring" = "#5AC870",
  "committed" = "#E8C9CF",
  "early stalk" = "#CB97A2",
  "late stalk" = "#AB747F", 
  "branching" = "#6C424C",
  "early female" = "#88d1b3",  
  "early male" = "#FFEDA1",  
  "late female" = "#6cad98",
  "late male" = "#FFAE2C"
)
dg_umap_plot <- DimPlot(DG_gametocytes_NF54, reduction = "umap", group.by = "stageHR") +
  scale_color_manual(values = stage_colors) + ggplot2::ggtitle("UMAP with Automatic Clusters")
dg_umap_plot_strain <- DimPlot(DG_gametocytes, reduction = "umap", group.by = "ident") +
  ggplot2::ggtitle("UMAP with Automatic Clusters")
dg_umap_plot
ggplot2::ggsave("dg_umap_plot.png", plot = dg_umap_plot, width = 8, height = 6, dpi = 300)

# slingshot
DG_sce <- as.SingleCellExperiment(DG_gametocytes_NF54)
start_cluster <- c("committed")

DG_sce <- slingshot(DG_sce, reducedDim = "PCA", clusterLabels = DG_sce$stageHR, start.clus = start_cluster)
saveRDS(DG_sce, file = "DG_sce.RDS")


DG_sce <- readRDS("/storage/work/skt5723/Single Cell Gametocyte stuff/DG_sce.RDS")

umap_embeddings <- reducedDim(DG_sce, "UMAP")

umap_data <- data.frame(
  UMAP1 = umap_embeddings[, 1],
  UMAP2 = umap_embeddings[, 2],
  Stage = DG_sce$stageHR
)

p <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = Stage)) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_manual(values = stage_colors) +
  theme_minimal() +
  labs(title = "Slingshot Trajectory on UMAP", x = "UMAP1", y = "UMAP2")

slingshot_data <- SlingshotDataSet(DG_sce)
for (curve in slingshot_data@curves) {
  curve_points_pca <- curve$s[curve$ord, ]
  
  pca_embeddings <- reducedDim(DG_sce, "PCA")  
  nn_indices <- FNN::get.knnx(pca_embeddings, curve_points_pca, k = 1)$nn.index  
  umap_embeddings <- reducedDim(DG_sce, "UMAP")  
  curve_points_umap <- umap_embeddings[nn_indices, , drop = FALSE]  

  curve_points_umap <- as.data.frame(curve_points_umap)
  colnames(curve_points_umap) <- c("UMAP1", "UMAP2")
  curve_points_umap$index <- 1:nrow(curve_points_umap)

  smoothed_umap1 <- loess(UMAP1 ~ index, data = curve_points_umap, span = 0.4)
  smoothed_umap2 <- loess(UMAP2 ~ index, data = curve_points_umap, span = 0.4)
  curve_points_umap$UMAP1 <- predict(smoothed_umap1)
  curve_points_umap$UMAP2 <- predict(smoothed_umap2)

  p <- p + geom_path(data = curve_points_umap, aes(x = UMAP1, y = UMAP2), color = "black", linewidth = 1)
}
p

ggsave("slingshot_trajectory_umap.png", plot = p, width = 8, height = 6, dpi = 300)
lineages <- slingshot::slingLineages(DG_sce)
print(lineages)
pseudotime <- slingPseudotime(slingshot_data)  # assuming your slingshot object is called 'sds'
umap_data$pseudotime_lineage1 <- pseudotime[,1]  # first lineage
umap_data$pseudotime_lineage2 <- pseudotime[,2]  # second lineage (if present)
pt_all <- c(umap_data$pseudotime_lineage1, umap_data$pseudotime_lineage2)
pt_range <- range(pt_all, na.rm = TRUE)
L1 <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = pseudotime_lineage1)) +
  geom_point(size = 1.5) +
  scale_color_viridis_c(limits = pt_range) +
  theme_minimal() +
  labs(title = "Pseudotime - Lineage 1")
L1
L2 <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = pseudotime_lineage2)) +
  geom_point(size = 1.5) +
  scale_color_viridis_c(limits = pt_range) +
  theme_minimal() +
  labs(title = "Pseudotime - Lineage 2")
L2
umap_long <- umap_data %>%
  pivot_longer(
    cols = starts_with("pseudotime_lineage"),
    names_to = "Lineage",
    values_to = "Pseudotime"
  )
ggplot(umap_long, aes(x = Stage, y = Pseudotime, fill = Lineage)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Pseudotime distribution across stages for both lineages")

ggsave("DG_female_lineage_umap.png", plot = L1, width = 6, height = 5)
ggsave("DG_male_lineage_umap.png", plot = L2, width = 6, height = 5)
# Extract pseudotime values
pseudotime <- slingshot::slingPseudotime(DG_sce)

counts <- counts(DG_sce) 

cell_weights <- slingshot::slingCurveWeights(DG_sce) 

set.seed(123) 
DG_knots <- evaluateK(
  counts = counts,
  pseudotime = pseudotime,
  cellWeights = cell_weights,
  k = 3:10, 
  nGenes = 200 
)
print(knots)

counts <- t(counts)
counts <- as.matrix(counts)

pseudotime_male <- pseudotime[, 1, drop = FALSE]
pseudotime_female <- pseudotime[, 2, drop = FALSE]

cell_weights_male <- cell_weights[, 1, drop = FALSE]
cell_weights_female <- cell_weights[, 2, drop = FALSE]

pseudotime_male[is.na(pseudotime_male)] <- 0
pseudotime_female[is.na(pseudotime_female)] <- 0
cell_weights_male[is.na(cell_weights_male)] <- 0
cell_weights_female[is.na(cell_weights_female)] <- 0

male_cells <- cell_weights_male > 0
female_cells <- cell_weights_female > 0

counts_male <- counts[male_cells, , drop = FALSE]
counts_male <- t(counts_male)
counts_female <- counts[female_cells, , drop = FALSE]
counts_female <- t(counts_female)

pseudotime_male <- pseudotime_male[male_cells, , drop = FALSE]
pseudotime_female <- pseudotime_female[female_cells, , drop = FALSE]
cell_weights_male <- cell_weights_male[male_cells, , drop = FALSE]
cell_weights_female <- cell_weights_female[female_cells, , drop = FALSE]

cat("Male counts:", dim(counts_male), "\n")
cat("Male pseudotime:", dim(pseudotime_male), "\n")
cat("Male cell weights:", dim(cell_weights_male), "\n")

cat("Female counts:", dim(counts_female), "\n")
cat("Female pseudotime:", dim(pseudotime_female), "\n")
cat("Female cell weights:", dim(cell_weights_female), "\n")

zero_sum_rows <- rowSums(counts_male) == 0
sum(zero_sum_rows)

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

check_input_dimensions(counts_male, pseudotime_male, cell_weights_male)

cell_weights_male[cell_weights_male <= 0] <- 1e-6

gam_fit_male <- fitGAM(
  counts = counts_male,
  pseudotime = as.vector(pseudotime_male),
  cellWeights = cell_weights_male,
  nknots = 8
)

degs_male <- startVsEndTest(
  models = gam_fit_male,
  global = TRUE,
  lineages = 1,
  l2fc = 1
)

degs_male <- degs_male[degs_male$waldStat > 60 & abs(degs_male$log2fc) > 1, ]

write.csv(degs_male, file = "degs_male_lineage.csv", row.names = TRUE)

degs_male_sorted <- degs_male[order(abs(degs_male$waldStat), decreasing = TRUE), ]
top_genes_male <- rownames(degs_male_sorted)[1:1000]

heatmap_data_male <- gam_fit_male[top_genes_male, ]
heatmap_data_male <- heatmap_data_male[, order(pseudotime_male[, 1])]
heatmap_data_male_matrix <- assay(heatmap_data_male, "counts")

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
