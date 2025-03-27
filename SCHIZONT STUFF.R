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
library(Cairo)
library(slingshot)
library(DelayedMatrixStats)
library(ggVennDiagram)
library(org.Pf.plasmo.db)
library(clusterProfiler)
library(dplyr)
library(txdbmaker)
setwd("/storage/work/skt5723/Single Cell Gametocyte stuff")
Dogga_data <- readRDS("/storage/group/mul27/default/scRNAseq_share/v3lab.ref.sce.RDS")
assay_data <- as(counts(Dogga_data), "dgCMatrix")
assay_data <- as(logcounts(Dogga_data), "CsparseMatrix")
meta_data <- as.data.frame(colData(Dogga_data))
Dogga_data_seurat <- CreateSeuratObject(
  counts = assay_data,  
  meta.data = meta_data
)
DG_schizonts <- subset(Dogga_data_seurat,
  subset = stageLR %in% c("early schizont", "late schizont"))
DG_rings <- subset(Dogga_data_seurat,
                       subset = stageLR %in% c("early ring", "late ring"))
table(Dogga_data_seurat@meta.data$stageLR)
table(DG_schizonts@meta.data$stageLR)
# Calculate the mean expression for each of the genes
#Schizonts
mean_expr_1222600 <- mean(GetAssayData(DG_schizonts, layer = "counts")["PF3D7-1222600", ][GetAssayData(DG_schizonts, layer = "counts")["PF3D7-1222600", ] > 0])
mean_expr_0935400 <- mean(GetAssayData(DG_schizonts, layer = "counts")["PF3D7-0935400", ][GetAssayData(DG_schizonts, layer = "counts")["PF3D7-0935400", ] > 0])
mean_expr_1328800 <- mean(GetAssayData(DG_schizonts, layer = "counts")["PF3D7-1328800", ][GetAssayData(DG_schizonts, layer = "counts")["PF3D7-1328800", ] > 0])
#Rings
mean_expr_1222600_rings <- mean(GetAssayData(DG_rings, layer = "counts")["PF3D7-1222600", ][GetAssayData(DG_rings, layer = "counts")["PF3D7-0406500", ] > 0])
mean_expr_0406500 <- mean(GetAssayData(DG_rings, layer = "counts")["PF3D7-0406500", ][GetAssayData(DG_rings, layer = "counts")["PF3D7-0406500", ] > 0])
mean_expr_0936500 <- mean(GetAssayData(DG_rings, layer = "counts")["PF3D7-0936500", ][GetAssayData(DG_rings, layer = "counts")["PF3D7-0936500", ] > 0])
overall_mean <- mean(c(mean_expr_1222600, mean_expr_0935400, mean_expr_1328800))
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

# Subset the metadata for only early and late schizonts
schizont_genes <- c("PF3D7-1222600", "PF3D7-0935400", "PF3D7-1328800")
ring_genes <- c("PF3D7-1222600", "PF3D7-0406500", "PF3D7-0936500")

# Combine the two gene sets and subset expression data
all_genes <- unique(c(schizont_genes, ring_genes))
expr_data <- GetAssayData(DG_committed, layer = "counts")[all_genes, ]

# Convert data to long format for ggplot
expr_long <- melt(as.data.frame(t(expr_data)))
colnames(expr_long) <- c("Cell", "Gene", "Expression")
expr_long$Stage <- DG_committed$stageLR[match(expr_long$Cell, rownames(DG_committed@meta.data))]
filtered_expr <- expr_long[expr_long$Stage %in% c("early schizont", "late schizont", "early ring", "late ring"), ]
gene_stage_plot <- ggplot(filtered_expr, aes(x = Stage, y = Expression, fill = Gene)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Expression of selected genes in schizont and ring stages",
       x = "Stage",
       y = "Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ Gene, scales = "free_y")  # One facet per gene
print(gene_stage_plot)
ggsave("gene_stage_boxplot.png", plot = gene_stage_plot, width = 8, height = 6, dpi = 300)

DG_committed <- subset(DG_committed, subset = ident %in% c("NF54"))
DG_committed <- NormalizeData(DG_committed)
DG_committed <- FindVariableFeatures(DG_committed, nfeatures = 1000)
DG_committed <- ScaleData(DG_committed)
DG_committed <- RunPCA(DG_committed, npcs = 50)
ElbowPlot(DG_committed)
DG_committed <- FindNeighbors(DG_committed, dims = 1:20)
DG_committed <- FindClusters(DG_committed, resolution = 0.25)  # Adjust resolution for finer or coarser clusters
DG_committed <- RunUMAP(
  DG_committed,
  dims = 1:20,                # Use the first 4 PCs
  n.components = 2,          # 2D UMAP
  min.dist = 0.1,            # Minimum distance between embedded points
  n.neighbors = 50,          # Number of neighbors
  spread = 1,                # Effective scale of embedded points
  local.connectivity = 1,    # Local connectivity
  a = 50,                    # Parameter for embedding
  b = 2                      # Parameter for embedding
)
dg_umap_plot <- DimPlot(DG_committed, reduction = "umap", group.by = "stageHR") +
  ggplot2::ggtitle("UMAP with Automatic Clusters")
dg_umap_plot
ggplot2::ggsave("dg_comit_umap_plot.png", plot = dg_umap_plot, width = 8, height = 6, dpi = 300)

# **Slingshot -----
# Convert Seurat object to SingleCellExperiment
DG_committed_sce <- as.SingleCellExperiment(DG_committed, assay = "RNA")
start_cluster <- c("early schizont")
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
rm(p)

# Save the UMAP plot with trajectories
ggsave("schizont slingshot_trajectory_umap.png", plot = p, width = 8, height = 6, dpi = 300)
lineages <- slingshot::slingLineages(DG_committed_sce)
print(lineages)

# Extract pseudotime values
pseudotime <- slingshot::slingPseudotime(DG_committed_sce)
pseudotime[is.na(pseudotime)] <- 0

# Prepare data for tradeSeq
counts <- counts(DG_committed_sce)  # Raw counts

cell_weights <- slingshot::slingCurveWeights(DG_committed_sce)  # Cell weights

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
print(DG_knots)

# Subset into sections
DG_r_s_sce <- DG_committed_sce[, colData(DG_committed_sce)$stageHR %in% c("early ring", "late ring", "early schizont", "late schizont", "committed")]
DG_er_ls_sce <- DG_committed_sce[, colData(DG_committed_sce)$stageHR %in% c("early ring", "late schizont", "late ring", "committed")]
DG_e_l_sce <- DG_committed_sce[, colData(DG_committed_sce)$stageHR %in% c("early schizont", "late schizont")]
# Preparing pseudotime for fitGAM
DG_r_s_sce_pseudotime <- slingshot::slingPseudotime(DG_r_s_sce)
DG_r_s_sce_pseudotime[is.na(DG_r_s_sce_pseudotime)] <- 0
DG_e_l_sce_pseudotime <- slingshot::slingPseudotime(DG_e_l_sce)
DG_e_l_sce_pseudotime[is.na(DG_e_l_sce_pseudotime)] <- 0
DG_er_ls_sce_pseudotime <- slingshot::slingPseudotime(DG_er_ls_sce)
DG_er_ls_sce_pseudotime[is.na(DG_er_ls_sce_pseudotime)] <- 0
# Preparing counts for fitGAM
DG_r_s_sce_counts <- counts(DG_r_s_sce)
DG_e_l_sce_counts <- counts(DG_e_l_sce)
DG_er_ls_sce_counts <- counts(DG_er_ls_sce)
DG_r_s_sce_counts <- as.matrix(DG_r_s_sce_counts)
DG_e_l_sce_counts <- as.matrix(DG_e_l_sce_counts)
DG_er_ls_sce_counts <- as.matrix(DG_er_ls_sce_counts)
# Preparing cell weights for fitGAM
DG_r_s_sce_cell_weights <- slingshot::slingCurveWeights(DG_r_s_sce) 
DG_e_l_sce_cell_weights <- slingshot::slingCurveWeights(DG_e_l_sce) 
DG_er_ls_sce_cell_weights <- slingshot::slingCurveWeights(DG_er_ls_sce) 

DG_r_s_sce_pseudotime <- DG_r_s_sce_pseudotime[, 1, drop = FALSE]
DG_e_l_sce_pseudotime <- DG_e_l_sce_pseudotime[, 1, drop = FALSE]
DG_er_ls_sce_pseudotime <- DG_er_ls_sce_pseudotime[, 1, drop = FALSE]

DG_r_s_sce_cell_weights <- DG_r_s_sce_cell_weights[, 1, drop = FALSE]
DG_r_s_sce_cell_weights[is.na(DG_r_s_sce_cell_weights)] <- 0
DG_e_l_sce_cell_weights <- DG_e_l_sce_cell_weights[, 1, drop = FALSE]
DG_e_l_sce_cell_weights[is.na(DG_e_l_sce_cell_weights)] <- 0
DG_er_ls_sce_cell_weights <- DG_er_ls_sce_cell_weights[, 1, drop = FALSE]
DG_r_s_sce_cell_weights[is.na(DG_r_s_sce_cell_weights)] <- 0

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
  nknots = 10
)

DG_e_l_sce_gam_fit <- fitGAM(
  counts = DG_e_l_sce_counts,
  pseudotime = as.vector(DG_e_l_sce_pseudotime),
  cellWeights = DG_e_l_sce_cell_weights,
  nknots = 10
)

DG_er_ls_sce_gam_fit <- fitGAM(
  counts = DG_er_ls_sce_counts,
  pseudotime = as.vector(DG_er_ls_sce_pseudotime),
  cellWeights = DG_er_ls_sce_cell_weights,
  nknots = 10
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
) +   ggplot2::scale_fill_gradient(low="white",high = "thistle")
print(fig_s5_B)
ggsave("fig_s5_B.png", plot = fig_s5_B, width = 8, height = 6, dpi = 300)

all_genes <- data.frame(
  geneid2 = union(
    earlyschizont_committed$geneid2,
    union(earlyschizont_lateschizont$geneid2, lateschizont_committed$geneid2)
  ),
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

ggsave("fig_s5_D.png", fig_s5_D ,
       width=7, height=2.7, dpi=300, bg = "transparent")

#DG_committed_sce
lineages <- slingshot::slingLineages(DG_committed_sce)
print(lineages)

# Extract pseudotime values
DG_committed_pseudotime <- slingshot::slingPseudotime(DG_committed_sce)
DG_committed_pseudotime[is.na(DG_committed_pseudotime)] <- 0

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



# Get the upregulated genes from each dataset
increasing_genes <- intersect(
  earlyschizont_committed$geneid2[earlyschizont_committed$diffexpressed == "UP"],
  earlyschizont_lateschizont$geneid2[earlyschizont_lateschizont$diffexpressed == "UP"]
)

upregulated_earlyschizont_committed <- earlyschizont_committed$geneid2[earlyschizont_committed$diffexpressed == "UP"]
upregulated_earlyschizont_lateschizont <- earlyschizont_lateschizont$geneid2[earlyschizont_lateschizont$diffexpressed == "UP"]
upregulated_lateschizont_committed <- lateschizont_committed$geneid2[lateschizont_committed$diffexpressed == "UP"]

# Find the intersection of these upregulated genes
common_upregulated_genes <- intersect(upregulated_earlyschizont_committed, upregulated_lateschizont_committed)

# Extract the data for these genes from each dataset
common_upregulated_genes_data <- degs_DG_r_s_sce[degs_DG_r_s_sce$geneid2 %in% upregulated_earlyschizont_lateschizont, ]

genesused <- c("PF3D7_1466800", "PF3D7_1146800", "PF3D7_0110000", "PF3D7_1369300", "PF3D7_1403200", "PF3D7_0114000", "PF3D7_1016900", "PF3D7_1416400", "PF3D7_0416100", "PF3D7_0406200", "PF3D7_1368800", "PF3D7_1312700", "PF3D7_1410000", "PF3D7_0513000", "PF3D7_0513700", "PF3D7_0421900", "PF3D7_1470600", "PF3D7_0314700", "PF3D7_0406500", "PF3D7_1429300", "PF3D7_1423400", "PF3D7_0513600", "PF3D7_1245000", "PF3D7_0517500", "PF3D7_0904200", "PF3D7_0813500", "PF3D7_1369400")
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
  
  plot_path <- paste0("/storage/work/skt5723/Single Cell Gametocyte stuff/other_smoother_plots/", j, ".png")
  ggsave(plot_path, plot = p, width = 7, height = 5, dpi = 600, bg = "transparent")
  
  # Print the file path to check if it's correct
  print(paste("Saved plot to:", plot_path))
  
  fig_s9_AB[[j]] <- p
}
