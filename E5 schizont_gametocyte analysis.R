BiocManager::install("GenomicFeatures", type = "source")
BiocManager::install("txdbmaker")
BiocManager::install("ComplexHeatmap")
install.packages("reshape2")
library(ComplexHeatmap)
library(tradeSeq)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(SingleCellExperiment)
library(Matrix)
library(slingshot)
library(DelayedMatrixStats)
library(ggVennDiagram)
library(org.Pf.plasmo.db)
library(clusterProfiler)
library(dplyr)
library(txdbmaker)
library(reshape2)
setwd("/storage/work/skt5723/Single Cell Gametocyte stuff")
Dogga_data <- readRDS("/storage/group/mul27/default/scRNAseq_share/v3lab.ref.sce.RDS")
assay_data <- as(counts(Dogga_data), "dgCMatrix")
assay_data <- as(logcounts(Dogga_data), "CsparseMatrix")
meta_data <- as.data.frame(colData(Dogga_data))
Dogga_data_seurat <- CreateSeuratObject(
  counts = assay_data,  
  meta.data = meta_data
)
Dogga_data_seurat[["percent.mt"]] <- PercentageFeatureSet(Dogga_data_seurat, pattern = "MIT")
max(Dogga_data_seurat@meta.data$percent.mt)
head(Dogga_data_seurat@meta.data, 5)

plot1 <- FeatureScatter(Dogga_data_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Dogga_data_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1
plot2
VlnPlot(Dogga_data_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
VlnPlot(Dogga_data_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.2)

Dogga_data_seurat <- subset(Dogga_data_seurat, subset = nCount_RNA > 1000 & nCount_RNA < 10000 & nFeature_RNA > 300 & nFeature_RNA < 3000 & percent.mt < 2)
Dogga_data_seurat <- NormalizeData(Dogga_data_seurat)
Dogga_data_seurat <- FindVariableFeatures(Dogga_data_seurat, nfeatures = 1000)
Dogga_data_seurat <- ScaleData(Dogga_data_seurat)
Dogga_data_seurat <- RunPCA(Dogga_data_seurat, npcs = 50)
#ElbowPlot(Dogga_data_seurat)
Dogga_data_seurat <- FindNeighbors(Dogga_data_seurat, dims = 1:20)
Dogga_data_seurat <- FindClusters(Dogga_data_seurat, resolution = 0.2)
Dogga_data_seurat <- FindNeighbors(Dogga_data_seurat, dims = 1:20)
Dogga_data_seurat <- FindClusters(Dogga_data_seurat, resolution = 0.25)
DG_schizonts <- subset(Dogga_data_seurat,
  subset = stageLR %in% c("early schizont", "late schizont"))
DG_rings <- subset(Dogga_data_seurat, subset = stageLR %in% c("early ring", "late ring"))
table(Dogga_data_seurat@meta.data$stageLR)
# Calculate the mean expression for each of the genes
#Schizonts
mean_expr_1222600 <- mean(GetAssayData(DG_schizonts, layer = "counts")["PF3D7-1222600", ][GetAssayData(DG_schizonts, layer = "counts")["PF3D7-1222600", ] > 0])
mean_expr_0935400 <- mean(GetAssayData(DG_schizonts, layer = "counts")["PF3D7-0935400", ][GetAssayData(DG_schizonts, layer = "counts")["PF3D7-0935400", ] > 0])
mean_expr_1328800 <- mean(GetAssayData(DG_schizonts, layer = "counts")["PF3D7-1328800", ][GetAssayData(DG_schizonts, layer = "counts")["PF3D7-1328800", ] > 0])
#Rings
mean_expr_1222600_rings <- mean(GetAssayData(DG_rings, layer = "counts")["PF3D7-1222600", ][GetAssayData(DG_rings, layer = "counts")["PF3D7-1222600", ] > 0])
mean_expr_0406500 <- mean(GetAssayData(DG_rings, layer = "counts")["PF3D7-0406500", ][GetAssayData(DG_rings, layer = "counts")["PF3D7-0406500", ] > 0])
mean_expr_0936500 <- mean(GetAssayData(DG_rings, layer = "counts")["PF3D7-0936500", ][GetAssayData(DG_rings, layer = "counts")["PF3D7-0936500", ] > 0])

DG_committed <- subset(
  Dogga_data_seurat,
  subset = stageLR %in% c("gametocyte (developing)", "gametocyte (male)", "gametocyte (female)") |
  (stageLR %in% c("early schizont", "late schizont") &
      (GetAssayData(Dogga_data_seurat, layer = "counts")["PF3D7-1222600", ] > mean_expr_1222600 |
        GetAssayData(Dogga_data_seurat, layer = "counts")["PF3D7-0935400", ] > mean_expr_0935400 |
        GetAssayData(Dogga_data_seurat, layer = "counts")["PF3D7-1328800", ] > mean_expr_1328800)) |
  (stageLR %in% c("early ring", "late ring") &
      (GetAssayData(Dogga_data_seurat, layer = "counts")["PF3D7-1222600", ] > mean_expr_1222600_rings |
        GetAssayData(Dogga_data_seurat, layer = "counts")["PF3D7-0406500", ] > mean_expr_0406500 |
        GetAssayData(Dogga_data_seurat, layer = "counts")["PF3D7-0936500", ] > mean_expr_0936500)) 
)
table(DG_committed@meta.data$stageLR)
Idents(DG_committed) <- DG_committed$stageLR

# Get counts matrix
counts <- GetAssayData(DG_committed, layer = "counts")

# Define the genes
schizont_genes <- c("PF3D7-1222600", "PF3D7-0935400", "PF3D7-1328800")
ring_genes <- c("PF3D7-1222600", "PF3D7-0406500", "PF3D7-0936500")

# Calculate mean expression externally (do not put into Seurat metadata)
mean_expr_schizont <- colMeans(counts[schizont_genes, , drop = FALSE])
mean_expr_ring <- colMeans(counts[ring_genes, , drop = FALSE])

# Create a dataframe for plotting
plot_df <- data.frame(
  cell = colnames(DG_committed),
  stage = DG_committed$stageLR,
  mean_expr = ifelse(DG_committed$stageLR %in% c("early ring", "late ring"),
                     mean_expr_ring,
              ifelse(DG_committed$stageLR %in% c("early schizont", "late schizont"),
                     mean_expr_schizont, NA)) 
)

# Keep only rings and schizonts
plot_df <- plot_df[plot_df$stage %in% c("early ring", "late ring", "early schizont", "late schizont"), ]

# Plot
stage_plot <- ggplot(plot_df, aes(x = stage, y = mean_expr, fill = stage)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Mean expression of ring and schizont genes",
       x = "Stage",
        y = "Mean Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, NA)


# Display the plot
stage_plot

# Save the plot
ggsave("stage_boxplot.png", plot = stage_plot, width = 6, height = 8, dpi = 300)

DG_committed <- subset(DG_committed, subset = ident %in% c("NF54"))
DG_committed <- NormalizeData(DG_committed)
DG_committed <- FindVariableFeatures(DG_committed, nfeatures = 1000)
DG_committed <- ScaleData(DG_committed)
DG_committed <- RunPCA(DG_committed, npcs = 50)
#ElbowPlot(DG_committed)
DG_committed <- FindNeighbors(DG_committed, dims = 1:20)
DG_committed <- FindClusters(DG_committed, resolution = 0.25)  # Adjust resolution for finer or coarser clusters
DG_committed <- RunUMAP(
  DG_committed,
  dims = 1:20,
  n.components = 2,
  min.dist = 0.1,
  n.neighbors = 50,
  spread = 1,
  local.connectivity = 1,
  a = 50,
  b = 2
)
stage_colors <- c(
  "early trophozoite" = "#FFEDA1",
  "late trophozoite" = "#FFAE2C",  
  "early schizont" = "#b69e00",
  "late schizont" = "#7aad00",
  "early ring" = "#f8746b",
  "late ring" = "#df8c00",
  "committed" = "#467063",
  "early stalk" = "#00c08b",
  "late stalk" = "#23872a", 
  "branching" = "#00b3f0",
  "early female" = "#619bff",  
  "early male" = "#fa87ed",  
  "late female" = "#c67bff",
  "late male" = "#ff64b0"
)

dg_umap_plot <- DimPlot(DG_committed, reduction = "umap", group.by = "stageHR") +
  scale_color_manual(values = stage_colors) + ggplot2::ggtitle("UMAP with Automatic Clusters")
dg_umap_plot
ggplot2::ggsave("dg_comit_umap_plot.png", plot = dg_umap_plot, width = 8, height = 6, dpi = 300)

# **Slingshot -----
# Convert Seurat object to SingleCellExperiment
DG_committed_sce <- as.SingleCellExperiment(DG_committed, assay = "RNA")
start_cluster <- c("early schizont")
# Run Slingshot using PCA embeddings
DG_committed_sce <- slingshot(DG_committed_sce, reducedDim = "PCA", clusterLabels = DG_committed_sce$stageHR, start.clus = start_cluster)
saveRDS(DG_committed_sce, file = "DG_committed_sce.RDS")
DG_committed_sce <- readRDS("DG_committed_sce.RDS")
pfdb56 <- read.csv("/storage/work/skt5723/Single Cell Gametocyte stuff/GenesByTaxon_GeneModelDump.csv")
umap_embeddings <- reducedDim(DG_committed_sce, "UMAP")

# Create a data frame for UMAP embeddings and stageHR
umap_data <- data.frame(
  UMAP1 = umap_embeddings[, 1],
  UMAP2 = umap_embeddings[, 2],
  Stage = DG_committed_sce$stageHR
)

#Slingshot trajectory analysis
p <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = Stage)) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_manual(values = stage_colors) +  
  theme_minimal() +
  labs(title = "Slingshot Trajectory on UMAP", x = "UMAP1", y = "UMAP2")


# Add Slingshot trajectories to the plot
slingshot_data <- SlingshotDataSet(DG_committed_sce)
for (curve in slingshot_data@curves) {
  # Extract curve points in PCA space
  curve_points_pca <- curve$s[curve$ord, ]
  
  # Map PCA points to UMAP space using nearest neighbors in PCA space
  pca_embeddings <- reducedDim(DG_committed_sce, "PCA")  # Get PCA embeddings
  nn_indices <- FNN::get.knnx(pca_embeddings, curve_points_pca, k = 1)$nn.index  # Find nearest neighbors
  curve_points_umap <- umap_embeddings[nn_indices, ]  # Use UMAP embeddings of nearest neighbors
  
  # Convert to data frame and set column names
  curve_points_umap <- as.data.frame(curve_points_umap)
  colnames(curve_points_umap) <- c("UMAP1", "UMAP2")
  
  # Smooth the curve points using loess
  curve_points_umap$index <- 1:nrow(curve_points_umap)
  smoothed_umap1 <- loess(UMAP1 ~ index, data = curve_points_umap, span = 0.4)
  smoothed_umap2 <- loess(UMAP2 ~ index, data = curve_points_umap, span = 0.4)
  curve_points_umap$UMAP1 <- predict(smoothed_umap1)
  curve_points_umap$UMAP2 <- predict(smoothed_umap2)
  
  # Add smoother trajectory to the plot
  p <- p + geom_path(data = curve_points_umap, aes(x = UMAP1, y = UMAP2), color = "black", linewidth = 1)
}
p
# Save the UMAP plot with trajectories
ggsave("schizont slingshot_trajectory_umap.png", plot = p, width = 8, height = 6, dpi = 300)
lineages <- slingshot::slingLineages(DG_committed_sce)
print(lineages)
rm(p)

# Subset into sections
DG_r_s_sce <- DG_committed_sce[, colData(DG_committed_sce)$stageHR %in% c("early ring", "late ring", "early schizont", "late schizont", "committed")]
DG_er_ls_sce <- DG_committed_sce[, colData(DG_committed_sce)$stageHR %in% c("early ring", "late schizont", "late ring", "committed")]
DG_e_l_sce <- DG_committed_sce[, colData(DG_committed_sce)$stageHR %in% c("early schizont", "late schizont")]
# Preparing pseudotime for fitGAM
# Remove cells where pseudotime or cell weight is NA or 0

# DG_r_s_sce
valid_cells_r_s <- !is.na(DG_r_s_sce_pseudotime) & DG_r_s_sce_pseudotime != 0 & 
  !is.na(DG_r_s_sce_cell_weights) & DG_r_s_sce_cell_weights != 0
DG_r_s_sce_pseudotime <- DG_r_s_sce_pseudotime[valid_cells_r_s]
DG_r_s_sce_cell_weights <- DG_r_s_sce_cell_weights[valid_cells_r_s]
DG_r_s_sce_counts <- DG_r_s_sce_counts[, valid_cells_r_s]

# DG_e_l_sce
valid_cells_e_l <- !is.na(DG_e_l_sce_pseudotime) & DG_e_l_sce_pseudotime != 0 & 
  !is.na(DG_e_l_sce_cell_weights) & DG_e_l_sce_cell_weights != 0
DG_e_l_sce_pseudotime <- DG_e_l_sce_pseudotime[valid_cells_e_l]
DG_e_l_sce_cell_weights <- DG_e_l_sce_cell_weights[valid_cells_e_l]
DG_e_l_sce_counts <- DG_e_l_sce_counts[, valid_cells_e_l]

# DG_er_ls_sce
valid_cells_er_ls <- !is.na(DG_er_ls_sce_pseudotime) & DG_er_ls_sce_pseudotime != 0 & 
  !is.na(DG_er_ls_sce_cell_weights) & DG_er_ls_sce_cell_weights != 0
DG_er_ls_sce_pseudotime <- DG_er_ls_sce_pseudotime[valid_cells_er_ls]
DG_er_ls_sce_cell_weights <- DG_er_ls_sce_cell_weights[valid_cells_er_ls]
DG_er_ls_sce_counts <- DG_er_ls_sce_counts[, valid_cells_er_ls]


DG_r_s_knots <- evaluateK(
  counts = DG_r_s_sce_counts,
  pseudotime = DG_r_s_sce_pseudotime,
  cellWeights = DG_r_s_sce_cell_weights,
  k = 3:10,  # Test 3 to 10 knots
  nGenes = 200  # Subsample 200 genes for evaluation
)

print(DG_r_s_knots)
DG_e_l_knots <- evaluateK(
  counts = DG_e_l_sce_counts,
  pseudotime = DG_e_l_sce_pseudotime,
  cellWeights = DG_e_l_sce_cell_weights,
  k = 3:10,  # Test 3 to 10 knots
  nGenes = 200  # Subsample 200 genes for evaluation
)

print(DG_e_l_knots)
DG_er_ls_knots <- evaluateK(
  counts = DG_er_ls_sce_counts,
  pseudotime = DG_er_ls_sce_pseudotime,
  cellWeights = DG_er_ls_sce_cell_weights,
  k = 3:10,  # Test 3 to 10 knots
  nGenes = 200  # Subsample 200 genes for evaluation
)
print(DG_er_ls_knots)

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

check_input_dimensions(DG_r_s_sce_counts, DG_r_s_sce_pseudotime, DG_r_s_sce_cell_weights)
check_input_dimensions(DG_e_l_sce_counts, DG_e_l_sce_pseudotime, DG_e_l_sce_cell_weights)
check_input_dimensions(DG_er_ls_sce_counts, DG_er_ls_sce_pseudotime, DG_er_ls_sce_cell_weights)

DG_r_s_sce_gam_fit <- fitGAM(
  counts = DG_r_s_sce_counts,
  pseudotime = as.vector(DG_r_s_sce_pseudotime),
  cellWeights = DG_r_s_sce_cell_weights,
  nknots = 9
)

DG_e_l_sce_gam_fit <- fitGAM(
  counts = DG_e_l_sce_counts,
  pseudotime = as.vector(DG_e_l_sce_pseudotime),
  cellWeights = DG_e_l_sce_cell_weights,
  nknots = 8
)

DG_er_ls_sce_gam_fit <- fitGAM(
  counts = DG_er_ls_sce_counts,
  pseudotime = as.vector(DG_er_ls_sce_pseudotime),
  cellWeights = DG_er_ls_sce_cell_weights,
  nknots = 9
)

saveRDS(DG_r_s_sce_gam_fit, file = "DG_r_s_sce_gam_fit.rds")
saveRDS(DG_e_l_sce_gam_fit, file = "DG_e_l_sce_gam_fit.rds")
saveRDS(DG_er_ls_sce_gam_fit, file = "DG_r_s_sce_gam_fit.rds")

DG_r_s_sce_gam_fit <- readRDS("DG_r_s_sce_gam_fit.rds")
DG_e_l_sce_gam_fit <- readRDS("DG_e_l_sce_gam_fit.rds")
DG_er_ls_sce_gam_fit <- readRDS("DG_r_s_sce_gam_fit.rds")

degs_DG_r_s_sce <- startVsEndTest(
  models = DG_r_s_sce_gam_fit,
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

degs_DG_er_ls_sce <- startVsEndTest(
  models = DG_r_s_sce_gam_fit,
  global = TRUE,
  lineages = 1,
  l2fc = 1
)

write.csv(degs_DG_r_s_sce, file = "degs_DG_r_s_sce.csv", row.names = TRUE)
write.csv(degs_DG_e_l_sce, file = "degs_DG_e_l_sce.csv", row.names = TRUE)
write.csv(degs_DG_er_ls_sce, file = "degs_DG_er_ls_sce.csv", row.names = TRUE)

de <- degs_DG_r_s_sce
de$logfc <- de$logFClineage1
de$wald <- de$waldStat
de$pvalue <- de$pvalue_lineage1
de$diffexpressed <- "NO"
de$diffexpressed[de$logfc > 1 & de$wald > 30] <- "UP"
de$diffexpressed[de$logfc < -1 & de$wald > 30] <- "DOWN"
def <- de[which(de$diffexpressed != "NO"),]
def$geneid2 <- rownames(def)
earlyschizont_committed <- def

de <- degs_DG_e_l_sce
de$logfc <- de$logFClineage1
de$wald <- de$waldStat
de$pvalue <- de$pvalue_lineage1
de$diffexpressed <- "NO"
de$diffexpressed[de$logfc > 1 & de$wald > 30] <- "UP"
de$diffexpressed[de$logfc < -1 & de$wald > 30] <- "DOWN"
def <- de[which(de$diffexpressed != "NO"),]
def$geneid2 <- rownames(def)
earlyschizont_lateschizont <- def

de <- degs_DG_er_ls_sce
de$logfc <- de$logFClineage1
de$wald <- de$waldStat
de$pvalue <- de$pvalue_lineage1
de$diffexpressed <- "NO"
de$diffexpressed[de$logfc > 1 & de$wald > 30] <- "UP"
de$diffexpressed[de$logfc < -1 & de$wald > 30] <- "DOWN"
def <- de[which(de$diffexpressed != "NO"),]
def$geneid2 <- rownames(def)
lateschizont_committed <- def

fig_s5_B <- ggVennDiagram(list(earlyschizont_committed$geneid2,earlyschizont_lateschizont$geneid2, lateschizont_committed$geneid2), label_alpha = 0, label = "count", label_size = 7,
                          category.names = c("earlyschizont_committed","earlyschizont_lateschizont", "lateschizont_committed")
) +  ggplot2::scale_fill_gradient(low="white",high = "thistle")
print(fig_s5_B)
ggsave("fig_s5_B.png", plot = fig_s5_B, width = 8, height = 6, dpi = 300)
# Find genes present in all three sets
common_genes <- Reduce(intersect, list(
  earlyschizont_committed$geneid2,
  earlyschizont_lateschizont$geneid2,
  lateschizont_committed$geneid2
))

# Print the common genes
print(common_genes)


all_genes <- data.frame(
  geneid2 = union(
  earlyschizont_committed$geneid2,
  union(earlyschizont_lateschizont$geneid2, lateschizont_committed$geneid2)
  )
  order = 1:length(
  union(
    earlyschizont_committed$geneid2,
    union(earlyschizont_lateschizont$geneid2, lateschizont_committed$geneid2)
  )
  )
)

all_genes <- left_join(
  all_genes,
  pfdb56,
  by = c("geneid2" = "Gene.ID")
)


### fig S5 D

geneList = keys(org.Pf.plasmo.db)
genelist <- all_genes
genelist$geneid2 <- gsub("-", "_", genelist$geneid2)

# module3
enr <- enrichGO(
  genelist$geneid2,
  universe = geneList,
  OrgDb = 'org.Pf.plasmo.db',
  pvalueCutoff = 0.5,
  keyType = "SYMBOL",  
  ont = "BP", 
  qvalueCutoff = 0.5
)

fig_s5_D <- clusterProfiler::dotplot(enr)
print(fig_s5_D)
go_res <- enr@result

ggsave("fig_s5_D.png", fig_s5_D , width=7, height=2.7, dpi=300, bg = "transparent")

#DG_committed_sce
lineages <- slingshot::slingLineages(DG_committed_sce)
print(lineages)
# Rerun of script but with committed trophs in smootherplots
DG_schizonts <- subset(Dogga_data_seurat, subset = stageLR %in% c("early schizont", "late schizont"))
DG_rings <- subset(Dogga_data_seurat, subset = stageLR %in% c("early ring", "late ring"))
DG_trophs <- subset(Dogga_data_seurat, subset = stageLR %in% c("late trophozoite"))
# Calculate the mean expression for each of the genes
#Schizonts
mean_expr_1222600 <- mean(GetAssayData(DG_schizonts, layer = "counts")["PF3D7-1222600", ][GetAssayData(DG_schizonts, layer = "counts")["PF3D7-1222600", ] > 0])
mean_expr_0935400 <- mean(GetAssayData(DG_schizonts, layer = "counts")["PF3D7-0935400", ][GetAssayData(DG_schizonts, layer = "counts")["PF3D7-0935400", ] > 0])
mean_expr_1328800 <- mean(GetAssayData(DG_schizonts, layer = "counts")["PF3D7-1328800", ][GetAssayData(DG_schizonts, layer = "counts")["PF3D7-1328800", ] > 0])
#Rings
mean_expr_1222600_rings <- mean(GetAssayData(DG_rings, layer = "counts")["PF3D7-1222600", ][GetAssayData(DG_rings, layer = "counts")["PF3D7-1222600", ] > 0])
mean_expr_0406500 <- mean(GetAssayData(DG_rings, layer = "counts")["PF3D7-0406500", ][GetAssayData(DG_rings, layer = "counts")["PF3D7-0406500", ] > 0])
mean_expr_0936500 <- mean(GetAssayData(DG_rings, layer = "counts")["PF3D7-0936500", ][GetAssayData(DG_rings, layer = "counts")["PF3D7-0936500", ] > 0])
#Trophs
mean_expr_1222600_trophs <- mean(GetAssayData(DG_trophs, layer = "counts")["PF3D7-1222600", ][GetAssayData(DG_trophs, layer = "counts")["PF3D7-1222600", ] > 0])
mean_expr_0935400_trophs <- mean(GetAssayData(DG_trophs, layer = "counts")["PF3D7-0935400", ][GetAssayData(DG_trophs, layer = "counts")["PF3D7-0935400", ] > 0])

DG_committed <- subset(
  Dogga_data_seurat,
  subset = stageLR %in% c("gametocyte (developing)", "gametocyte (male)", "gametocyte (female)") |
    (stageLR %in% c("early schizont", "late schizont") &
      (GetAssayData(Dogga_data_seurat, layer = "counts")["PF3D7-1222600", ] > mean_expr_1222600 |
        GetAssayData(Dogga_data_seurat, layer = "counts")["PF3D7-0935400", ] > mean_expr_0935400 |
        GetAssayData(Dogga_data_seurat, layer = "counts")["PF3D7-1328800", ] > mean_expr_1328800)) 
  |
  (stageLR %in% c("early ring", "late ring") &
      (GetAssayData(Dogga_data_seurat, layer = "counts")["PF3D7-1222600", ] > mean_expr_1222600_rings |
        GetAssayData(Dogga_data_seurat, layer = "counts")["PF3D7-0406500", ] > mean_expr_0406500 |
        GetAssayData(Dogga_data_seurat, layer = "counts")["PF3D7-0936500", ] > mean_expr_0936500)) 
  |
  (stageLR %in% c("late trophozoite") &
      (GetAssayData(Dogga_data_seurat, layer = "counts")["PF3D7-1222600", ] > mean_expr_1222600_trophs |
        GetAssayData(Dogga_data_seurat, layer = "counts")["PF3D7-0935400", ] > mean_expr_0935400_trophs) 
))
table(DG_committed@meta.data$stageLR)
table(Dogga_data_seurat@meta.data$stageLR)
Idents(DG_committed) <- DG_committed$stageLR

# Get counts matrix
counts <- GetAssayData(DG_committed, layer = "counts")

# Define the genes
schizont_genes <- c("PF3D7-1222600", "PF3D7-0935400", "PF3D7-1328800")
ring_genes <- c("PF3D7-1222600", "PF3D7-0406500", "PF3D7-0936500")
troph_genes <- c("PF3D7-1222600", "PF3D7-0935400")

mean_expr_schizont <- colMeans(counts[schizont_genes, , drop = FALSE])
mean_expr_ring <- colMeans(counts[ring_genes, , drop = FALSE])
mean_expr_troph <- colMeans(counts[troph_genes, , drop = FALSE])


plot_df <- data.frame(
  cell = colnames(DG_committed),
  stage = DG_committed$stageLR,
  mean_expr = ifelse(DG_committed$stageLR %in% c("early ring", "late ring"),
                     mean_expr_ring,
              ifelse(DG_committed$stageLR %in% c("early schizont", "late schizont"),
                     mean_expr_schizont,
              ifelse(DG_committed$stageLR %in% c("late trophozoite"),
                      mean_expr_troph,NA )))
)

# Keep only rings and schizonts
plot_df <- plot_df[plot_df$stage %in% c("early ring", "late ring", "early schizont", "late schizont", "late trophozoite"), ]

# Plot
stage_plot <- ggplot(plot_df, aes(x = stage, y = mean_expr, fill = stage)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Mean expression of ring and schizont genes",
    x = "Stage",
    y = "Mean Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, NA)


# Display the plot
stage_plot

# Save the plot
ggsave("stage_boxplot_w_trops.png", plot = stage_plot, width = 6, height = 8, dpi = 300)

DG_committed <- subset(DG_committed, subset = ident %in% c("NF54"))
DG_committed <- NormalizeData(DG_committed)
DG_committed <- FindVariableFeatures(DG_committed, nfeatures = 1000)
DG_committed <- ScaleData(DG_committed)
DG_committed <- RunPCA(DG_committed, npcs = 50)
ElbowPlot(DG_committed)
DG_committed <- FindNeighbors(DG_committed, dims = 1:19)
DG_committed <- FindClusters(DG_committed, resolution = 0.25)  # Adjust resolution for finer or coarser clusters
DG_committed <- RunUMAP(
  DG_committed,
  dims = 1:19,
  n.components = 2,
  min.dist = 0.1,
  n.neighbors = 50,
  spread = 1,
  local.connectivity = 1,
  a = 50,
  b = 2
)

dg_umap_plot <- DimPlot(DG_committed, reduction = "umap", group.by = "stageHR") +
  scale_color_manual(values = stage_colors) + ggplot2::ggtitle("UMAP with Automatic Clusters")
dg_umap_plot
ggplot2::ggsave("dg_troph_comit_umap_plot.png", plot = dg_umap_plot, width = 8, height = 6, dpi = 300)

# **Slingshot -----
# Convert Seurat object to SingleCellExperiment
DG_committed_sce <- as.SingleCellExperiment(DG_committed, assay = "RNA")
start_cluster <- c("early schizont") #use late trophozoite if you want to reattempt the trophozoite slingshot analysis
# Run Slingshot using PCA embeddings
DG_committed_sce <- slingshot(DG_committed_sce, reducedDim = "PCA", clusterLabels = DG_committed_sce$stageHR, start.clus = start_cluster)
saveRDS(DG_committed_sce, file = "DG_committed_sce.RDS")
pfdb56 <- read.csv("/storage/work/skt5723/Single Cell Gametocyte stuff/GenesByTaxon_GeneModelDump.csv")
umap_embeddings <- reducedDim(DG_committed_sce, "UMAP")

# Create a data frame for UMAP embeddings and stageHR
umap_data <- data.frame(
  UMAP1 = umap_embeddings[, 1],
  UMAP2 = umap_embeddings[, 2],
  Stage = DG_committed_sce$stageHR
)

#Slingshot trajectory analysis
p <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = Stage)) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_manual(values = stage_colors) +  
  theme_minimal() +
  labs(title = "Slingshot Trajectory on UMAP", x = "UMAP1", y = "UMAP2")


# Add Slingshot trajectories to the plot
slingshot_data <- SlingshotDataSet(DG_committed_sce)
for (curve in slingshot_data@curves) {
  # Extract curve points in PCA space
  curve_points_pca <- curve$s[curve$ord, ]
  
  # Map PCA points to UMAP space using nearest neighbors in PCA space
  pca_embeddings <- reducedDim(DG_committed_sce, "PCA")  # Get PCA embeddings
  nn_indices <- FNN::get.knnx(pca_embeddings, curve_points_pca, k = 1)$nn.index  # Find nearest neighbors
  curve_points_umap <- umap_embeddings[nn_indices, ]  # Use UMAP embeddings of nearest neighbors
  
  # Convert to data frame and set column names
  curve_points_umap <- as.data.frame(curve_points_umap)
  colnames(curve_points_umap) <- c("UMAP1", "UMAP2")
  
  # Smooth the curve points using loess
  curve_points_umap$index <- 1:nrow(curve_points_umap)
  smoothed_umap1 <- loess(UMAP1 ~ index, data = curve_points_umap, span = 0.4)
  smoothed_umap2 <- loess(UMAP2 ~ index, data = curve_points_umap, span = 0.4)
  curve_points_umap$UMAP1 <- predict(smoothed_umap1)
  curve_points_umap$UMAP2 <- predict(smoothed_umap2)
  
  # Add smoother trajectory to the plot
  p <- p + geom_path(data = curve_points_umap, aes(x = UMAP1, y = UMAP2), color = "black", linewidth = 1)
}
p
# Save the UMAP plot with trajectories
ggsave("schizont_slingshot_trajectory_umap.png", plot = p, width = 8, height = 6, dpi = 300)
lineages <- slingshot::slingLineages(DG_committed_sce)
print(lineages)
rm(p)

# Extract pseudotime values
DG_committed_pseudotime <- slingshot::slingPseudotime(DG_committed_sce)
DG_committed_pseudotime[is.na(DG_committed_pseudotime)] <- 1e-6

# Prepare data for tradeSeq
DG_committed_counts <- counts(DG_committed_sce)  # Raw counts

DG_committed_cell_weights <- slingshot::slingCurveWeights(DG_committed_sce)  # Cell weights

check_input_dimensions(DG_committed_counts, DG_committed_pseudotime, DG_committed_cell_weights)

DG_sex_sce_gam_fit <- fitGAM(
  counts = DG_committed_counts,
  pseudotime = DG_committed_pseudotime,
  cellWeights = DG_committed_cell_weights,
  nknots = 8
)
saveRDS(DG_sex_sce_gam_fit, file = "DG_committed_sex_sce_gam_fit.rds")
DG_sex_sce_gam_fit <- readRDS("/storage/work/skt5723/Single Cell Gametocyte stuff/DG_committed_sex_sce_gam_fit.rds")

pseudotime_values <- as.data.frame(slingshot::slingPseudotime(DG_committed_sce))

# Get cell annotations
cell_annotations <- colData(DG_committed_sce)$stageHR

# Combine into a single data frame
pseudotime_annotations_df <- data.frame(
  Cell = colnames(DG_committed_sce),
  Annotation = cell_annotations,
  Pseudotime = pseudotime_values[, 1] # Assuming you’re checking the first lineage
)

# Filter the data for each annotation and extract pseudotime ranges
stalk_pseudotime <- pseudotime_annotations_df$Pseudotime[pseudotime_annotations_df$Annotation == "early stalk"]
branching_pseudotime <- pseudotime_annotations_df$Pseudotime[pseudotime_annotations_df$Annotation == "branching"]
Fem_pseudotime <- pseudotime_annotations_df$Pseudotime[pseudotime_annotations_df$Annotation == "early female"]
Mal_pseudotime <- pseudotime_annotations_df$Pseudotime[pseudotime_annotations_df$Annotation == "early male"]
Ring_pseudotime <- pseudotime_annotations_df$Pseudotime[pseudotime_annotations_df$Annotation == "early ring"]
Schizont_pseudotime <- pseudotime_annotations_df$Pseudotime[pseudotime_annotations_df$Annotation == "early schizont"]
# Output basic stats to see ranges
summary(Schizont_pseudotime)
summary(Ring_pseudotime)
summary(stalk_pseudotime)
summary(branching_pseudotime)
summary(Fem_pseudotime)
summary(Mal_pseudotime)

upregulated_earlyschizont_committed <- earlyschizont_committed$geneid2[earlyschizont_committed$diffexpressed == "UP"]
upregulated_earlyschizont_lateschizont <- earlyschizont_lateschizont$geneid2[earlyschizont_lateschizont$diffexpressed == "UP"]
upregulated_lateschizont_committed <- lateschizont_committed$geneid2[lateschizont_committed$diffexpressed == "UP"]
#after filtering them out of R via comparison to malaria cell atlas
genesused <- c("PF3D7-1328800", "PF3D7-1222600", "PF3D7-0935400", "PF3D7-1033200", "PF3D7-1335000", "PF3D7-1227600",  "PF3D7-0820900", 
               "PF3D7-1222400", "PF3D7-1134100", "PF3D7-0113800", "PF3D7-1460500", "PF3D7-1433800",  "PF3D7-0207700", "PF3D7-0624400", 
               "PF3D7-0813200", "PF3D7-1019300", "PF3D7-0219800", "PF3D7-0316700", "PF3D7-0216800", "PF3D7-1466200", "PF3D7-0830800", 
               "PF3D7-0416300", "PF3D7-1451800", "PF3D7-1232200", "PF3D7-1465100", "PF3D7-0615500", "PF3D7-0815500", "PF3D7-1410700",
               "PF3D7-1433900", "PF3D7-1111900", "PF3D7-1431500", "PF3D7-1013000", "PF3D7-0914800", "PF3D7-0731600", "PF3D7-0828300",
               "PF3D7-1334800", "PF3D7-0413700", "PF3D7-1035000", "PF3D7-0606800", "PF3D7-1145600", "PF3D7-1011400", "PF3D7-1104300",
               "PF3D7-0902800", "PF3D7-1328000", "PF3D7-1248600", "PF3D7-0829500", "PF3D7-0912400", "PF3D7-0715900", "PF3D7-1323300",
               "PF3D7-0301200", "PF3D7-0811100", "PF3D7-0807500", "PF3D7-1305900", "PF3D7-1330400.1", "PF3D7-0911600", "PF3D7-1432900", 
               "PF3D7-0512500", "PF3D7-0114000", "PF3D7-0936100", "PF3D7-1469600", "PF3D7-1437500", "PF3D7-1315800", "PF3D7-1207200",
               "PF3D7-0307500", "PF3D7-1420600", "PF3D7-1237000", "PF3D7-0306500", "PF3D7-0915700", "PF3D7-1121900", "PF3D7-0715300",
               "PF3D7-0515100", "PF3D7-1301900", "PF3D7-1024500", "PF3D7-1029800", "PF3D7-0618500", "PF3D7-1333900", "PF3D7-0826200",
               "PF3D7-1306200", "PF3D7-0804400", "PF3D7-1467600", "PF3D7-0213000", "PF3D7-1415400", "PF3D7-1353800", "PF3D7-1015100",
               "PF3D7-0810000", "PF3D7-0623000", "PF3D7-0704200", "PF3D7-0317400", "PF3D7-1417700", "PF3D7-0207900", "PF3D7-1439000",
               "PF3D7-1106200", "PF3D7-0808000", "PF3D7-1021300", "PF3D7-1205700", "PF3D7-1350400", "PF3D7-1457800", "PF3D7-0605500",
               "PF3D7-0704900" )

genesused <- gsub("_","-",genesused)

fig_s9_AB <- list()
if(j %in% rownames(counts(DG_sex_sce_gam_fit))) {
  print(paste("Gene", j, "is found in the counts data"))
} else {
  print(paste("Gene", j, "is not found in the counts data"))
}

for(j in genesused){
  p <- plotSmoothers(DG_sex_sce_gam_fit, counts(DG_sex_sce_gam_fit), gene = j, xlab = paste(j)) + geom_vline(xintercept=0.74, linetype="dashed", color = "red", linewidth=1) +
   geom_vline(xintercept=46, linetype="dashed", color = "red", linewidth=1) +
   geom_vline(xintercept=70, linetype="dashed", color = "red", linewidth=1) +
   geom_vline(xintercept=75, linetype="dashed", color = "red", linewidth=1) + theme_classic(base_size = 25) + theme(legend.position="none")
  
  plot_path <- paste0("/storage/work/skt5723/Single Cell Gametocyte stuff/new_smoother_plots_2/", j, ".png")
  ggsave(plot_path, plot = p, width = 7, height = 5, dpi = 600, bg = "transparent")
  
  # Print the file path to check if it's correct
  print(paste("Saved plot to:", plot_path))
  
  fig_s9_AB[[j]] <- p
}

genesused2 <- c("PF3D7-1439000", "PF3D7-1232200", "PF3D7-1403200", "PF3D7-1415400")

fig_s9_AB <- list()
if(j %in% rownames(counts(DG_sex_sce_gam_fit))) {
  print(paste("Gene", j, "is found in the counts data"))
} else {
  print(paste("Gene", j, "is not found in the counts data"))
}

for (j in genesused2) {
  p <- plotSmoothers(
    DG_sex_sce_gam_fit, 
    counts(DG_sex_sce_gam_fit), 
    gene = j, 
    xlab = paste(j),
    curvesCols = curvesCols,
    sample = 0  # <- disables plotting individual data points
  ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
    geom_vline(xintercept = 46, linetype = "dashed", color = "red", linewidth = 1) +
    theme_classic(base_size = 25) +
    theme(legend.position = "none")
  
  plot_path <- paste0("/storage/work/skt5723/Single Cell Gametocyte stuff/Schizont_smoother_plots/", j, ".png")
  ggsave(plot_path, plot = p, width = 7, height = 5, dpi = 600, bg = "transparent")
  
  print(paste("Saved plot to:", plot_path))
  
  fig_s9_AB[[j]] <- p
}
