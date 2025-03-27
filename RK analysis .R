install.packages("Cairo")
install.packages("seriation")
install.packages("RColorBrewer")
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
library(seriation)
library(RColorBrewer)

setwd("/storage/work/skt5723/Single Cell Gametocyte stuff")
Sys.setenv(DISPLAY=":0")  # Prevent X11 errors
options(device = "pdf")  
Rajat_data <- readRDS("/storage/work/skt5723/RKsc/PfE5_SingleR.rds")

Rajat_data_gametocytes <- subset(Rajat_data, subset = SingleR_labels %in% c("branching", "committed", "early female", "early male", "early stalk", "late female", "late male", "late stalk"))
Rajat_data_gametocytes <- subset(Rajat_data_gametocytes, subset = my_ann %in% c("Early and Late Female", "Branching", "Stalk", "Late Male"))
Rajat_data_gametocytes[["percent.mt"]] <- PercentageFeatureSet(Rajat_data_gametocytes, pattern = "MIT")
max(Rajat_data_gametocytes@meta.data$percent.mt)
head(Rajat_data_gametocytes@meta.data, 5)

plot1 <- FeatureScatter(Rajat_data_gametocytes, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Rajat_data_gametocytes, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1
plot2
VlnPlot(Rajat_data_gametocytes, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
VlnPlot(Rajat_data_gametocytes, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.2)

Rajat_data_gametocytes <- subset(Rajat_data_gametocytes, subset = nCount_RNA > 1000 & nCount_RNA < 10000 & nFeature_RNA > 300 & nFeature_RNA < 3000 & percent.mt < 2)
table(Rajat_data_gametocytes@meta.data$SingleR_labels)
Rajat_data_gametocytes <- NormalizeData(Rajat_data_gametocytes)
Rajat_data_gametocytes <- FindVariableFeatures(Rajat_data_gametocytes, nfeatures = 1000)
Rajat_data_gametocytes <- ScaleData(Rajat_data_gametocytes)
Rajat_data_gametocytes <- RunPCA(Rajat_data_gametocytes, npcs = 50)
RK_ElbowPlot <- ElbowPlot(Rajat_data_gametocytes)
print(RK_ElbowPlot)
ggplot2::ggsave("rk_elbow_plot.png", plot = RK_ElbowPlot, width = 8, height = 6, dpi = 300)
Rajat_data_gametocytes <- FindNeighbors(Rajat_data_gametocytes, dims = 1:20)
Rajat_data_gametocytes <- FindClusters(Rajat_data_gametocytes, resolution = 0.25)  # Adjust resolution for finer or coarser clusters
Rajat_data_gametocytes <- RunUMAP(
  Rajat_data_gametocytes,
  dims = 1:4,              
  n.components = 2,       
  min.dist = 0.3,          
  n.neighbors = 50,        
  spread = 0.5,              
  local.connectivity = 2,  
  a = 100,                  
  b = 1                       
)
rk_umap_plot_myann <- DimPlot(Rajat_data_gametocytes, reduction = "umap", group.by = "my_ann") +
  ggplot2::ggtitle("RK UMAP myann")
rk_umap_plot_labels <- DimPlot(Rajat_data_gametocytes, reduction = "umap", group.by = "SingleR_labels") +
  ggplot2::ggtitle("RK UMAP SRL")
rk_umap_plot_labels
# Save the plot to a file
ggplot2::ggsave("rk_gam_umap_myann.png", plot = rk_umap_plot_myann, width = 8, height = 6, dpi = 300)
ggplot2::ggsave("rk_gam_umap_SRL.png", plot = rk_umap_plot_labels, width = 8, height = 6, dpi = 300)
# **Slingshot -----
# Convert Seurat object to SingleCellExperiment
RK_sce <- as.SingleCellExperiment(Rajat_data_gametocytes)
start_cluster <- c("committed")
# Run Slingshot using PCA embeddings
RK_sce <- slingshot(RK_sce, reducedDim = "PCA", clusterLabels = RK_sce$SingleR_labels, start.clus = start_cluster)
saveRDS(RK_sce, file = "RK_sce.RDS")

RK_sce <- readRDS("/storage/work/skt5723/Single Cell Gametocyte stuff/RK_sce.RDS")
pfdb56 <- read.csv("/storage/work/skt5723/Single Cell Gametocyte stuff/GenesByTaxon_GeneModelDump.csv")
umap_embeddings <- reducedDim(RK_sce, "UMAP")

# Create a data frame for UMAP embeddings and stageHR
umap_data <- data.frame(
  UMAP1 = umap_embeddings[, 1],
  UMAP2 = umap_embeddings[, 2],
  Stage = RK_sce$my_ann
)

#Slingshot trajectory analysis
p <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = Stage)) +
  geom_point(size = 1.5, alpha = 0.8) +
  theme_minimal() +
  labs(title = "Slingshot Trajectory on UMAP", x = "UMAP1", y = "UMAP2")


# Add Slingshot trajectories to the plot
slingshot_data <- SlingshotDataSet(RK_sce)
for (curve in slingshot_data@curves) {
  # Extract curve points in PCA space
  curve_points_pca <- curve$s[curve$ord, ]
  
  # Map PCA points to UMAP space using nearest neighbors in PCA space
  pca_embeddings <- reducedDim(RK_sce, "PCA")  # Get PCA embeddings
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

# Save the UMAP plot with trajectories
ggsave("RK_slingshot_trajectory_umap.png", plot = p, width = 8, height = 6, dpi = 300)
lineages <- slingshot::slingLineages(RK_sce)
print(lineages)

RK_c_l_sce <- RK_sce[, colData(RK_sce)$SingleR_labels %in% c("committed", "early stalk", "late stalk")]
RK_e_l_sce <- RK_sce[, colData(RK_sce)$SingleR_labels %in% c("early stalk", "late stalk")]
# Extract pseudotime values
RK_c_l_pseudotime <- slingshot::slingPseudotime(RK_c_l_sce)
RK_c_l_pseudotime[is.na(RK_c_l_pseudotime)] <- 0
RK_e_l_pseudotime <- slingshot::slingPseudotime(RK_e_l_sce)
RK_e_l_pseudotime[is.na(RK_e_l_pseudotime)] <- 0
# Prepare data for tradeSeq
RK_c_l_counts <- counts(RK_c_l_sce)  # Raw counts
RK_e_l_counts <- counts(RK_e_l_sce)
RK_c_l_cell_weights <- slingshot::slingCurveWeights(RK_c_l_sce)  # Cell weights
RK_e_l_cell_weights <- slingshot::slingCurveWeights(RK_e_l_sce)
# Evaluate optimal number of knots
set.seed(123)  # For reproducibility
# Evaluate knots
RK_c_l_knots <- evaluateK(
  counts = RK_c_l_counts,
  pseudotime = RK_c_l_pseudotime,
  cellWeights = RK_c_l_cell_weights,
  k = 3:10,  # Test 3 to 10 knots
  nGenes = 200  # Subsample 200 genes for evaluation
)
#7
RK_e_l_knots <- evaluateK(
  counts = RK_e_l_counts,
  pseudotime = RK_e_l_pseudotime,
  cellWeights = RK_e_l_cell_weights,
  k = 3:10,  # Test 3 to 10 knots
  nGenes = 200  # Subsample 200 genes for evaluation
)

# Extract pseudotime values
RK_c_l_sce_pseudotime <- slingshot::slingPseudotime(RK_c_l_sce)
RK_c_l_sce_pseudotime[is.na(RK_c_l_sce_pseudotime)] <- 0
RK_e_l_sce_pseudotime <- slingshot::slingPseudotime(RK_e_l_sce)
RK_e_l_sce_pseudotime[is.na(RK_e_l_sce_pseudotime)] <- 0
# Prepare data for tradeSeq
RK_c_l_sce_counts <- counts(RK_c_l_sce)
RK_e_l_sce_counts <- counts(RK_e_l_sce)

RK_c_l_sce_cell_weights <- slingshot::slingCurveWeights(RK_c_l_sce) 
RK_e_l_sce_cell_weights <- slingshot::slingCurveWeights(RK_e_l_sce) 
RK_c_l_sce_cell_weights <- RK_c_l_sce_cell_weights[, 1, drop = FALSE]
RK_c_l_sce_cell_weights[RK_c_l_sce_cell_weights <= 0] <- 1e-6 
RK_e_l_sce_cell_weights <- RK_e_l_sce_cell_weights[, 1, drop = FALSE]
RK_e_l_sce_cell_weights[RK_e_l_sce_cell_weights <= 0] <- 1e-6 


RK_c_l_sce_counts <- as.matrix(RK_c_l_sce_counts)
RK_e_l_sce_counts <- as.matrix(RK_e_l_sce_counts)

RK_c_l_sce_pseudotime <- RK_c_l_sce_pseudotime[, 1, drop = FALSE]
RK_c_l_sce_pseudotime[is.na(RK_c_l_sce_pseudotime)] <- 0
RK_e_l_sce_pseudotime <- RK_e_l_sce_pseudotime[, 1, drop = FALSE]
RK_e_l_sce_pseudotime[is.na(RK_e_l_sce_pseudotime)] <- 0



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

check_input_dimensions(RK_c_l_sce_counts, RK_c_l_sce_pseudotime, RK_c_l_sce_cell_weights)
check_input_dimensions(RK_e_l_sce_counts, RK_e_l_sce_pseudotime, RK_e_l_sce_cell_weights)

RK_c_l_sce_gam_fit <- fitGAM(
  counts = RK_c_l_sce_counts,
  pseudotime = RK_c_l_sce_pseudotime,
  cellWeights = RK_c_l_sce_cell_weights,
  nknots = 7
)
RK_e_l_sce_gam_fit <- fitGAM(
  counts = RK_e_l_sce_counts,
  pseudotime = as.vector(RK_e_l_sce_pseudotime),
  cellWeights = RK_e_l_sce_cell_weights,
  nknots = 7
)

# Save each GAM fit as an RDS file
saveRDS(RK_c_l_sce_gam_fit, file = "RK_c_l_sce_gam_fit.rds")
saveRDS(RK_e_l_sce_gam_fit, file = "RK_e_l_sce_gam_fit.rds")

# To load them back later:
RK_c_l_sce_gam_fit <- readRDS("RK_c_l_sce_gam_fit.rds")
RK_e_l_sce_gam_fit <- readRDS("RK_e_l_sce_gam_fit.rds")

degs_RK_c_l_sce1 <- startVsEndTest(
  models = RK_c_l_sce_gam_fit,
  global = TRUE,
  lineages = 1,
  l2fc = 1
)
degs_RK_e_l_sce1 <- startVsEndTest(
  models = RK_e_l_sce_gam_fit,
  global = TRUE,
  lineages = 1,
  l2fc = 1
)

de <- degs_RK_c_l_sce1
de$logfc <- de$logFClineage1
de$wald <- de$waldStat
de$pvalue <- de$pvalue_lineage1
de$diffexpressed <- "NO"
de$diffexpressed[de$logfc > 2 & de$wald > 30] <- "UP"
de$diffexpressed[de$logfc < -2 & de$wald > 30] <- "DOWN"
def <- de[which(de$diffexpressed != "NO"),]
def$geneid2 <- rownames(def)
committed_latestalk <- def


de <- degs_RK_e_l_sce1
de$logfc <- de$logFClineage1
de$wald <- de$waldStat
de$pvalue <- de$pvalue_lineage1
de$diffexpressed <- "NO"
de$diffexpressed[de$logfc > 2 & de$wald > 30] <- "UP"
de$diffexpressed[de$logfc < -2 & de$wald > 30] <- "DOWN"
def <- de[which(de$diffexpressed != "NO"),]
def$geneid2 <- rownames(def)
earlystalk_latestalk <- def

rk_fig_s5_B <- ggVennDiagram(list(committed_latestalk$geneid2,earlystalk_latestalk$geneid2), label_alpha = 0, label = "count", label_size = 7,
                          category.names = c("committed_latestalk","earlystalk_latestalk")
) +   ggplot2::scale_fill_gradient(low="white",high = "thistle")
rk_fig_s5_B
ggsave("rk_fig_s5_B.png", plot = rk_fig_s5_B, width = 8, height = 6, dpi = 300)

stalk_both_60 <- data.frame(geneid2 = union(committed_latestalk$geneid2,earlystalk_latestalk$geneid2),order = c(1:length(union(committed_latestalk$geneid2,earlystalk_latestalk$geneid2))))
stalk_both_60 <- left_join(
  stalk_both_60,
  pfdb56,
  by = c("geneid2" = "Gene.ID")
)


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

rk_fig_s5_D <- clusterProfiler::dotplot(enr)
go_res <- enr@result
rk_fig_s5_D
ggsave("rk_fig_s5_D.png", rk_fig_s5_D ,
       width=7, height=2.7, dpi=300, bg = "transparent")
genelist <- stalk_both_60

md <- Rajat_data_gametocytes@meta.data
n <- 100
p <- 100

RK_sex_sce_gam_fit@colData$stagenew <- md$SingleR_labels
RK_sce <- subset(RK_sex_sce_gam_fit, , stagenew %in% c('committed','early stalk','late stalk'))

glist <- genelist

yhatSmoothx <- predictSmooth(RK_sce, gene = glist@geneid2, nPoints = n, tidy = TRUE)
yhatSmooth <- predictSmooth(RK_sce, gene = glist@geneid2, nPoints = n, tidy = FALSE)
x <- t(scale(t(yhatSmooth[, 1:n])))
o <- seriate(dist(x), method = "R2E")
glist <- glist[get_order(o),]
x <- t(scale(t(yhatSmooth[, 1:p])))[glist$geneid2,]

yhatSmoothx <- yhatSmoothx %>% mutate(Group =
                                        case_when(lineage == '1' & time <= 0.74 ~ "committed", 
                                                  lineage == '1' & time > 0.74 & time <= 1.165 ~ "early stalk",
                                                  lineage == '1' & time > 1.165 ~ "late stalk",
                                                  lineage == '2' & time <= 0.74 ~ "committed", 
                                                  lineage == '2' & time > 0.74 & time <= 1.165 ~ "early stalk",
                                                  lineage == '2' & time > 1.165 ~ "late stalk")
)

annotation_col = data.frame(
  Time = yhatSmoothx$time[1:p],
  Stages = yhatSmoothx$Group[1:p]
)

ann_colors = list(
  Time = colors <- colorRampPalette(brewer.pal(11,'Blues')[-6])(n),
  Stages = colors_v3all)

# annotation_row = data.frame(
# GeneClass = factor(c(rep("S1",179)))
# )
# rownames(annotation_row) = c(1:nrow(genelist))

annotation_row = data.frame(
  GeneClass = factor(c(rep("S1",25),rep("S2",8),rep("S3",146)))
)
rownames(annotation_row) = glist$geneid2

fig_s5_hm <- ComplexHeatmap::pheatmap(x,
                                      cluster_cols = FALSE,
                                      cluster_rows = FALSE,
                                      show_rownames = TRUE,
                                      show_colnames = FALSE,
                                      labels_row = glist$genename2,
                                      annotation_col = annotation_col,
                                      fontsize_row = 6)
fig_s5_hm

RK_sex_sce <- RK_sce
lineages <- slingshot::slingLineages(RK_sex_sce)
print(lineages)

# Extract pseudotime values
RK_sex_pseudotime <- slingshot::slingPseudotime(RK_sex_sce)
RK_sex_pseudotime[is.na(RK_sex_pseudotime)] <- 0

# Prepare data for tradeSeq
RK_sex_counts <- counts(RK_sex_sce)  # Raw counts

RK_sex_cell_weights <- slingshot::slingCurveWeights(RK_sex_sce)  # Cell weights

check_input_dimensions(RK_sex_counts, RK_sex_pseudotime, RK_sex_cell_weights)
# View your cell annotations
table(colData(RK_sex_sce)$my_ann)

# View pseudotime values
head(slingshot::slingPseudotime(RK_sex_sce))
# Extract pseudotime for each cell (if multiple lineages, it gives a matrix)
pseudotime_values <- as.data.frame(slingshot::slingPseudotime(RK_sex_sce))

# Get cell annotations
cell_annotations <- colData(RK_sex_sce)$my_ann

# Combine into a single data frame
pseudotime_annotations_df <- data.frame(
  Cell = colnames(RK_sex_sce),
  Annotation = cell_annotations,
  Pseudotime = pseudotime_values[, 1] # Assuming you’re checking the first lineage
)

# Filter the data for each annotation and extract pseudotime ranges
stalk_pseudotime <- pseudotime_annotations_df$Pseudotime[pseudotime_annotations_df$Annotation == "Stalk"]
branching_pseudotime <- pseudotime_annotations_df$Pseudotime[pseudotime_annotations_df$Annotation == "Branching"]
Fem_pseudotime <- pseudotime_annotations_df$Pseudotime[pseudotime_annotations_df$Annotation == "Early and Late Female"]
Mal_pseudotime <- pseudotime_annotations_df$Pseudotime[pseudotime_annotations_df$Annotation == "Late Male"]
# Output basic stats to see ranges
summary(stalk_pseudotime)
summary(branching_pseudotime)
summary(Fem_pseudotime)
summary(Mal_pseudotime)

RK_sex_knots <- evaluateK(
  counts = RK_sex_counts,
  pseudotime = RK_sex_pseudotime,
  cellWeights = RK_sex_cell_weights,
  k = 3:10,  # Test 3 to 10 knots
  nGenes = 200  # Subsample 200 genes for evaluation
)

ggplot(pseudotime_annotations_df, aes(x = Pseudotime, color = Annotation)) +
  geom_density() +
  theme_minimal() +
  labs(title = "Annotations Along Pseudotime")

RK_sex_sce_gam_fit <- fitGAM(
  counts = RK_sex_counts,
  pseudotime = RK_sex_pseudotime,
  cellWeights = RK_sex_cell_weights,
  nknots = 7
)

saveRDS(RK_sex_sce_gam_fit, file = "RK_sex_sce_gam_fit.rds")
RK_sex_sce_gam_fit <- readRDS("RK_sex_sce_gam_fit.rds")

degs_RK_fsex <- startVsEndTest(
  models = RK_sex_sce_gam_fit,
  global = TRUE,
  lineages = 1,
  l2fc = 1
)
degs_RK_msex <- startVsEndTest(
  models = RK_sex_sce_gam_fit,
  global = TRUE,
  lineages = 2,
  l2fc = 1
)

lineages <- slingshot::slingLineages(RK_sex_sce_gam_fit)
print(lineages)
de <- degs_RK_msex
de$logfc <- de$logFClineage2
de$wald <- de$waldStat
de$pvalue <- de$pvalue_lineage2
de$diffexpressed <- "NO"
de$diffexpressed[de$logfc > 1 & de$wald > 30] <- "UP"
de$diffexpressed[de$logfc < -1 & de$wald > 30] <- "DOWN"
def <- de[which(de$diffexpressed != "NO"),]
def$geneid2 <- rownames(def)
committed_male <- def
table(committed_male$diffexpressed)


de <- degs_RK_fsex
de$logfc <- de$logFClineage1
de$wald <- de$waldStat
de$pvalue <- de$pvalue_lineage1
de$diffexpressed <- "NO"
de$diffexpressed[de$logfc > 1 & de$wald > 30] <- "UP"
de$diffexpressed[de$logfc < -1 & de$wald > 30] <- "DOWN"
def <- de[which(de$diffexpressed != "NO"),]
def$geneid2 <- rownames(def)
committed_female <- def
table(committed_female$diffexpressed)

up_male <- (
  committed_male$geneid2[committed_male$diffexpressed == "UP"]
)

up_female <- (
  committed_female$geneid2[committed_female$diffexpressed == "UP"]
)

cm_up_male <- setdiff(
  committed_male$geneid2[committed_male$diffexpressed == "UP"], 
  committed_female$geneid2[committed_female$diffexpressed == "UP"]
)
cf_up_female <- setdiff(
  committed_female$geneid2[committed_female$diffexpressed == "UP"], 
  committed_male$geneid2[committed_male$diffexpressed == "UP"]
)

committed_latestalkf <- intersect(
  committed_latestalk$geneid2, cf_up_female
)
committed_latestalkm <- intersect(
  committed_latestalk$geneid2, cm_up_male
)
earlystalk_latestalkf <- intersect(
  committed_latestalk$geneid2, cf_up_female
)
earlystalk_latestalkm <- intersect(
  committed_latestalk$geneid2, cm_up_male
)
genesused <- c("PF3D7-0110000", "PF3D7-1146800", 
              "PF3D7-1416400", "PF3D7-1403200","PF3D7-1466800", "PF3D7-0422300", "PF3D7-0905300", "PF3D7-1122900", "PF3D7−1202300", "PF3D7-1227500")

fig_s9_AB <- list()
for(j in genesused){
  p <- plotSmoothers(RK_sex_sce_gam_fit, counts(RK_sex_sce_gam_fit), gene = j, xlab = paste(j)) + geom_vline(xintercept=0.74, linetype="dashed", color = "red", linewidth=1) +
    geom_vline(xintercept=25, linetype="dashed", color = "red", linewidth=1) +
    geom_vline(xintercept=28, linetype="dashed", color = "red", linewidth=1) +
    geom_vline(xintercept=45, linetype="dashed", color = "red", linewidth=1) + theme_classic(base_size = 25) + theme(legend.position="none")
  
  plot_path <- paste0("/storage/work/skt5723/Single Cell Gametocyte stuff/other_smoother_plots/", j, ".png")
  ggsave(plot_path, plot = p, width = 7, height = 5, dpi = 600, bg = "transparent")
  
  # Print the file path to check if it's correct
  print(paste("Saved plot to:", plot_path))
  
  fig_s9_AB[[j]] <- p
}
sc_annotated_data <- readRDS("/storage/work/skt5723/RKsc/PfE5_SingleR.rds")
sc_annotated_data[["percent.mt"]] <- PercentageFeatureSet(sc_annotated_data, pattern = "MIT")
max(sc_annotated_data@meta.data$percent.mt)
head(sc_annotated_data@meta.data, 5)

plot1 <- FeatureScatter(sc_annotated_data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sc_annotated_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1
plot2
VlnPlot(sc_annotated_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
VlnPlot(sc_annotated_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.2)
Rajat_data_gametocytes <- subset(Rajat_data_gametocytes, subset = nCount_RNA > 1000 & nCount_RNA < 10000 & nFeature_RNA > 300 & nFeature_RNA < 3000 & percent.mt < 2)
sc_annotated_data <- NormalizeData(sc_annotated_data)
sc_annotated_data <- FindVariableFeatures(sc_annotated_data, nfeatures = 1000)
sc_annotated_data <- ScaleData(sc_annotated_data)
sc_annotated_data <- RunPCA(sc_annotated_data, npcs = 50)
elbow_plot <- ElbowPlot(sc_annotated_data)
ggsave("elbow_plot.png", plot = elbow_plot, width = 6, height = 4, dpi = 300)
sc_annotated_data <- FindNeighbors(sc_annotated_data, dims = 1:20)
sc_annotated_data <- FindClusters(sc_annotated_data, resolution = 0.2)  # Adjust resolution for finer or coarser clusters
stage_colors <- c(
  "branching" = "#1f77b4",  
  "committed" = "#c0c2ff",  
  "early female" = "#ff7878",  
  "early male" = "#7dff7b",  
  "early ring" = "#feffc0",
  "early schizont" = "#fed8a3",
  "early stalk" = "#6b70fc",
  "early trophozoite" = "#ffe075",
  "late female" = "#ff0000",
  "late male" = "#006100",
  "late ring" = "#f7fa00",
  "late schizont" = "#fa9000",
  "late stalk" = "#131bff",
  "late trophozoite" = "#f9c100" 
)
FeaturePlot(sc_annotated_data, features = "PF3D7-1205700", cols = c("lightgrey", "blue"))
sc_annotated_data <- RunUMAP(
  sc_annotated_data,
  dims = 1:20,                # Use the first 4 PCs
  n.components = 2,          # 2D UMAP
  min.dist = 0.3,            # Minimum distance between embedded points
  n.neighbors = 50,          # Number of neighbors
  spread = 1,                # Effective scale of embedded points
  local.connectivity = 1,    # Local connectivity
  a = 70,                    # Parameter for embedding
  b = 1                      # Parameter for embedding
)

dg_umap_plot <- DimPlot(sc_annotated_data, reduction = "umap", group.by = "SingleR_labels") +
  scale_color_manual(values = stage_colors) +
  ggplot2::ggtitle("UMAP with Automatic Clusters")
dg_umap_plot
ggplot2::ggsave("RK_whole_umap_plot.png", plot = dg_umap_plot, width = 8, height = 6, dpi = 300)
print(lineages)
male_genes <- cm_up_male
female_genes <- cf_up_female
late_stalk_cells <- subset(sc_annotated_data, subset = my_ann == "Stalk")

femgenes <- c("PF3D7-1457100", "PF3D7-1357600", 
              "PF3D7-1316300", "PF3D7-1310700","PF3D7-1250100", "PF3D7-0920300", "PF3D7-0914400")
malgenes <- c("PF3D7-1205700", "PF3D7-1122900", "PF3D7-0905300")

Idents(sc_annotated_data) <- "SingleR_labels"
lsm <- as.data.frame(rowSums(FetchData(late_stalk_cells, vars = malgenes, layer = "data")))
colnames(lsm) <- "lsm"
lsm$rank <- c(1:nrow(lsm))
mcells <- lsm[order(lsm$lsm, decreasing = TRUE),]
mcells$rank <- c(1:nrow(mcells))
mcellsx <- mcells[which(mcells$lsm >= mean(lsm$lsm)),] %>% rownames()
fcellsx <- tail(rownames(mcells), length(mcellsx))

lsf <- as.data.frame(rowSums(FetchData(late_stalk_cells, vars = femgenes, layer = "data")))
colnames(lsf) <- "lsf"
lsf$rank <- c(1:nrow(lsf))
fcells <- lsf[order(lsf$lsf, decreasing = TRUE),]
fcells$rank <- c(1:nrow(fcells))
fcellsy <- fcells[which(fcells$lsf >= mean(lsf$lsf)),] %>% rownames()
mcellsy <- tail(rownames(fcells), length(fcellsy))

intersect(mcellsx, fcellsy) %>% length()
mcellsz <- setdiff(mcellsx,fcellsy)
fcellsz <- setdiff(fcellsy,mcellsx)
late_stalk_cells@meta.data$cellid <- rownames(late_stalk_cells@meta.data)

combined_df <- data.frame(
  mcellsz = c(mcellsz, rep(NA, max(0, length(fcellsz) - length(mcellsz)))),
  fcellsz = c(fcellsz, rep(NA, max(0, length(mcellsz) - length(fcellsz))))
)
write.csv(combined_df, "mcellsz_fcellsz.csv", row.names = FALSE)

late_stalk_cells@meta.data$fm <- "other"
late_stalk_cells@meta.data[which(late_stalk_cells@meta.data$SingleR_labels == "late stalk"),]$fm <- "other late stalk"
late_stalk_cells@meta.data$cellids <- rownames(late_stalk_cells@meta.data)
late_stalk_cells@meta.data[which(late_stalk_cells@meta.data$cellid %in% mcellsz),]$fm <- "male-like"
late_stalk_cells@meta.data[which(late_stalk_cells@meta.data$cellid %in% fcellsz),]$fm <- "female-like"
View(late_stalk_cells@meta.data)
table(late_stalk_cells$fm)

sc_annotated_data@meta.data$cellid <- rownames(sc_annotated_data@meta.data)
sc_annotated_data@meta.data$fm <- "other"
sc_annotated_data@meta.data[which(sc_annotated_data@meta.data$SingleR_labels == "late stalk"),]$fm <- "other late stak"
sc_annotated_data@meta.data[which(sc_annotated_data@meta.data$SingleR_labels == "late stalk" & sc_annotated_data@meta.data$cellid %in% mcellsx),]$fm <- "male-like"
sc_annotated_data@meta.data[which(sc_annotated_data@meta.data$SingleR_labels == "late stalk" & sc_annotated_data@meta.data$cellid %in% fcellsx),]$fm <- "female-like"
table(sc_annotated_data$fm)

saveRDS(late_stalk_cells, file = "RK_classified_late_stalk_cells.rds")
saveRDS(sc_annotated_data, file = "Classified_RK.rds")
write.csv(late_stalk_cells@meta.data, "late_stalk_cells.csv", row.names = FALSE)

late_stalk_cells <- readRDS("/storage/work/skt5723/Single Cell Gametocyte stuff/RK_classified_late_stalk_cells.rds")
sc_annotated_data <- readRDS("/storage/work/skt5723/Single Cell Gametocyte stuff/Classified_RK.rds")
Idents(late_stalk_cells) <- late_stalk_cells$fm
late_stalk_cells <- NormalizeData(late_stalk_cells)

defm <- FindMarkers(late_stalk_cells,ident.1 = mcellsz, ident.2 = fcellsz, only.pos = FALSE, test.use = "MAST")
defm$geneid2 <- rownames(defm)

defm <- defm %>% mutate(genename2 =
                          case_when(genename == 'N/A' ~ geneid2,
                                    TRUE ~ genename))

defm$updown <- "other" 
defm[which(defm$avg_log2FC > 0.4 & defm$p_val_adj < 0.05),]$updown <- "up"
defm[which(defm$avg_log2FC < -0.4 & defm$p_val_adj < 0.05),]$updown <- "down"
table(defm$updown)
write.csv(defm, "/storage/work/skt5723/Single Cell Gametocyte stuff/RK_defm_results.csv", row.names = TRUE)
updown <- table(defm$updown)
p <- ggplot(data=defm, aes(x=avg_log2FC, y=-log10(p_val_adj), col=updown)) + geom_point(size=2.5) + theme_minimal()
# Add lines as before...
fig_s9_C1 <- p + geom_vline(xintercept=c(-0.4, 0.4), col="red") +
  geom_hline(yintercept=1.30103, col="red") + 
  annotate("label", x =c(-1.5,0.8), y = 25, label = c(paste0(updown[[1]]," genes"), paste0(updown[[3]]," genes")), col=c("red","steelblue"), size = 6) +
  annotate("text", x =c(-1), y = 17, label = c(paste0("log2FC > 0.4, padj < 0.05")), size = 6) +
  theme(legend.position="none") + 
  theme(axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 16))
fig_s9_C1
ggplot2::ggsave("rk_fig_s9_c1.png", plot = fig_s9_C1, width = 8, height = 6, dpi = 300)
rm(p)

data <- cbind(sc_annotated_data@meta.data[,c("orig.ident","fm")], sc_annotated_data@reductions$umap@cell.embeddings)
data$label <- "select"
data[which(data$fm != "female-like"),]$label <- "other"
data$order <- ifelse(data$label=="other", 1, 2)
p <- ggplot(data) + geom_point(aes(x = umap_1, y = umap_2, color = label),size=3, shape=16) + 
  geom_point(data = subset(data, label == "select"), aes(x = umap_1, y = umap_2, color = label),size=3, shape=16) + 
  theme_classic() +
  ggtitle("Female Like Late Stalk") +
  theme(legend.title = element_text(colour="Black", size=12, face="bold"),
        legend.text = element_text(colour="Black", size=16))  +   scale_color_manual(values = c("other"="grey","select"="#a3767f")) + theme(legend.position="none") 
print(p)
ggplot2::ggsave("rk_fig_s9_c2.png", plot = p, width = 8, height = 6, dpi = 300)
rm(p)

data <- cbind(sc_annotated_data@meta.data[,c("orig.ident","fm")], sc_annotated_data@reductions$umap@cell.embeddings)
data$label <- "select"
data[which(data$fm != "male-like"),]$label <- "other"
data$order <- ifelse(data$label=="other", 1, 2)
p <- ggplot(data) + geom_point(aes(x = umap_1, y = umap_2, color = label),size=3, shape=16) + 
  geom_point(data = subset(data, label == "select"), aes(x = umap_1, y = umap_2, color = label),size=3, shape=16) + 
  theme_classic() +
  ggtitle("Male Like Late Stalk") +
  theme(legend.title = element_text(colour="Black", size=12, face="bold"),
        legend.text = element_text(colour="Black", size=16))  +   scale_color_manual(values = c("other"="grey","select"="#a3767f")) + theme(legend.position="none")
print(p)
ggplot2::ggsave("rk_fig_s9_c3.png", plot = p, width = 8, height = 6, dpi = 300)
rm(p)

scm <- as.data.frame(rowSums(FetchData(sc_annotated_data, vars = malgenes, layer = "data")))
colnames(scm) <- "scm"
scm$rank <- c(1:nrow(scm))
mcells <- scm[order(scm$scm, decreasing = TRUE),]
mcells$rank <- c(1:nrow(mcells))
mcellsx <- mcells[which(mcells$scm >= mean(scm$scm)),] %>% rownames()
fcellsx <- tail(rownames(mcells), length(mcellsx))

scf <- as.data.frame(rowSums(FetchData(sc_annotated_data, vars = femgenes, layer = "data")))
colnames(scf) <- "scf"
scf$rank <- c(1:nrow(scf))
fcells <- scf[order(scf$scf, decreasing = TRUE),]
fcells$rank <- c(1:nrow(fcells))
fcellsy <- fcells[which(fcells$scf >= mean(scf$scf)),] %>% rownames()
mcellsy <- tail(rownames(fcells), length(fcellsy))

intersect(mcellsx, fcellsy) %>% length()
mcellsz <- setdiff(mcellsx,fcellsy)
fcellsz <- setdiff(fcellsy,mcellsx)
sc_annotated_data@meta.data$cellid <- rownames(sc_annotated_data@meta.data)

seur <- RK_sex_sce
seur$newfm <- seur$fm
seur$newfm <- gsub("other late stak","other",seur$newfm)

sel <- defm[which(defm$updown == "up"),]
gene.set <- sel$geneid2
scale_color_manual(values = c("trophozoite"="#FEEEAA","gametocyte (developing)"="thistle",
                              "gametocyte (male)"="mediumpurple","gametocyte (female)"="",
                              "not assigned"="grey","ring"="#78C679","schizont"="#85B1D3"))


data <- cbind(sc_annotated_data@meta.data[,c("cellid","fm")], sc_annotated_data@reductions$umap@cell.embeddings)
data$label <- "other"
data[which(data$fm == "female-like"),]$label <- "female-like"
data[which(data$fm == "male-like"),]$label <- "male-like"
data$order <- ifelse(data$label=="other", 1, 2)
sc_annotated_data <- RunPCA(sc_annotated_data, npcs = 50)
RK_ElbowPlot <- ElbowPlot(sc_annotated_data)
print(RK_ElbowPlot)
sc_annotated_data <- FindNeighbors(sc_annotated_data, dims = 1:20)
sc_annotated_data <- FindClusters(sc_annotated_data, resolution = 0.25)  # Adjust resolution for finer or coarser clusters
sc_annotated_data <- RunUMAP(
  sc_annotated_data,
  dims = 1:4,              
  n.components = 2,       
  min.dist = 0.3,          
  n.neighbors = 50,        
  spread = 0.5,              
  local.connectivity = 2,  
  a = 50,                  
  b = 1                       
)
rkg_umap_plot_myann <- DimPlot(Rajat_data_gametocytes, reduction = "umap", group.by = "my_ann") +
  ggplot2::ggtitle("RK UMAP myann")
print(rkg_umap_plot_myann)
a3 <- ggplot(data) + geom_point(aes(x = umap_1, y = umap_2, color = label),size=3, shape=16) + 
  geom_point(data = subset(data, label == "select"), aes(x = umap_1, y = umap_2, color = label),size=3, shape=16) + 
  theme_classic() +
  theme(legend.title = element_text(colour="Black", size=12, face="bold"),
        legend.text = element_text(colour="Black", size=16))  +   scale_color_manual(values = c("other"="grey","female-like"="#551A8B","male-like"="mediumpurple")) + theme(legend.position="none") 

print(a3)

FeaturePlot(sc_annotated_data, features = "PF3D7-0728100") 

data <- cbind(sc_annotated_data@meta.data[,c("orig.ident","fm")], sc_annotated_data@reductions$umap@cell.embeddings)
data$label <- "select"
data[which(data$fm != "female-like"),]$label <- "other"
data$order <- ifelse(data$label=="other", 1, 2)
a4 <- ggplot(data) + geom_point(aes(x = umap_1, y = umap_2, color = label),size=3, shape=16) + 
  geom_point(data = subset(data, label == "select"), aes(x = umap_1, y = umap_2, color = label),size=3, shape=16) + 
  theme_classic() +
  theme(legend.title = element_text(colour="Black", size=12, face="bold"),
        legend.text = element_text(colour="Black", size=16))  +   scale_color_manual(values = c("other"="grey","select"="#551A8B")) + theme(legend.position="none") 

print(a4)

data <- cbind(v3sex@meta.data[,c("orig.ident","fm")], v3sex@reductions$umap@cell.embeddings)
data$label <- "select"
data[which(data$fm != "male-like"),]$label <- "other"
data$order <- ifelse(data$label=="other", 1, 2)
a5 <- ggplot(data) + geom_point(aes(x = UMAP_2, y = UMAP_3, color = label),size=3, shape=16) + 
  geom_point(data = subset(data, label == "select"), aes(x = UMAP_2, y = UMAP_3, color = label),size=3, shape=16) + 
  theme_classic() +
  theme(legend.title = element_text(colour="Black", size=12, face="bold"),
        legend.text = element_text(colour="Black", size=16))  +   scale_color_manual(values = c("other"="grey","select"="mediumpurple")) + theme(legend.position="none")

a5 <- blankit(a5)
a5

