install.packages("htmlwidgets")
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
library(plotly)
library(htmlwidgets)
setwd("/storage/work/skt5723/Single Cell Gametocyte stuff")
Dd2_data <- readRDS("/storage/work/skt5723/Single Cell Gametocyte stuff/Plasmodium_Dd2.rds")
table(Dd2_data@meta.data$my_ann)
table(Dd2_data@meta.data$StageHR_v3lab)
Dd2_data_gametocytes <- subset(Dd2_data, subset = StageHR_v3lab %in% c("committed", "early stalk", "late stalk", "branching", "early female", "late female", "early male", "late male"))
table(Dd2_data_gametocytes@meta.data$my_ann)
table(Dd2_data_gametocytes@meta.data$StageHR_v3lab)
Dd2_data_gametocytes <- NormalizeData(Dd2_data_gametocytes)
Dd2_data_gametocytes <- FindVariableFeatures(Dd2_data_gametocytes, nfeatures = 1000)
Dd2_data_gametocytes <- ScaleData(Dd2_data_gametocytes)
Dd2_data_gametocytes <- RunPCA(Dd2_data_gametocytes, npcs = 50)
Dd2_ElbowPlot <- ElbowPlot(Dd2_data_gametocytes)
print(Dd2_ElbowPlot)
Dd2_data_gametocytes <- FindNeighbors(Dd2_data_gametocytes, dims = 1:20)
Dd2_data_gametocytes <- FindClusters(Dd2_data_gametocytes, resolution = 0.25) 
Dd2_data_gametocytes <- RunUMAP(
  Dd2_data_gametocytes,
  dims = 1:20,              
  n.components = 3,       
  min.dist = 0.01,          
  n.neighbors = 50,        
  spread = 0.2,              
  local.connectivity = 2,  
  a = 1000,                  
  b = 1                       
)
umap_3d <- Embeddings(Dd2_data_gametocytes, "umap")

umap_plot <- plot_ly(
  x = ~umap_3d[,1],
  y = ~umap_3d[,2],
  z = ~umap_3d[,3],
  color = ~Idents(Dd2_data_gametocytes),
  colors = RColorBrewer::brewer.pal(n = length(unique(Idents(Dd2_data_gametocytes))), "Set1"),
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 3, opacity = 0.8)
) %>%
  layout(scene = list(
    xaxis = list(title = "UMAP 1"),
    yaxis = list(title = "UMAP 2"),
    zaxis = list(title = "UMAP 3")
  ))

saveWidget(umap_plot, "3D_UMAP_plot.html")
Dd2_umap_plot_myann <- DimPlot(Dd2_data_gametocytes, reduction = "umap", group.by = "my_ann") +
  ggplot2::ggtitle("Dd2 UMAP myann")
Dd2_umap_plot_labels <- DimPlot(Dd2_data_gametocytes, reduction = "umap", group.by = "SingleR_labels") +
  ggplot2::ggtitle("Dd2 UMAP SRL")
print(Dd2_umap_plot_myann)
print(Dd2_umap_plot_labels)
ggplot2::ggsave("Dd2_gam_umap_myann.png", plot = Dd2_umap_plot_myann, width = 8, height = 6, dpi = 300)
ggplot2::ggsave("Dd2_gam_umap_SRL.png", plot = Dd2_umap_plot_labels, width = 8, height = 6, dpi = 300)
# **Slingshot -----
Dd2_sce <- as.SingleCellExperiment(Dd2_data_gametocytes)
start_cluster <- c("committed")
Dd2_sce <- slingshot(Dd2_sce, reducedDim = "PCA", clusterLabels = Dd2_sce$StageHR_v3lab, start.clus = start_cluster)
saveRDS(Dd2_sce, file = "Dd2_sce.RDS")

Dd2_sce <- readRDS("/storage/work/skt5723/Single Cell Gametocyte stuff/Dd2_sce.RDS")
pfdb56 <- read.csv("/storage/work/skt5723/Single Cell Gametocyte stuff/GenesByTaxon_GeneModelDump.csv")
umap_embeddings <- reducedDim(Dd2_sce, "UMAP")

umap_data <- data.frame(
  UMAP1 = umap_embeddings[, 1],
  UMAP2 = umap_embeddings[, 2],
  Stage = Dd2_sce$my_ann
)

p <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = Stage)) +
  geom_point(size = 1.5, alpha = 0.8) +
  theme_minimal() +
  labs(title = "Slingshot Trajectory on UMAP", x = "UMAP1", y = "UMAP2")

slingshot_data <- SlingshotDataSet(Dd2_sce)
for (curve in slingshot_data@curves) {
  curve_points_pca <- curve$s[curve$ord, ]
  
  pca_embeddings <- reducedDim(Dd2_sce, "PCA")  
  nn_indices <- FNN::get.knnx(pca_embeddings, curve_points_pca, k = 1)$nn.index
  curve_points_umap <- umap_embeddings[nn_indices, ]  
  
  curve_points_umap <- as.data.frame(curve_points_umap)
  colnames(curve_points_umap) <- c("UMAP1", "UMAP2")
  
  p <- p + geom_path(data = curve_points_umap, aes(x = UMAP1, y = UMAP2), color = "black", linewidth = 1)
}

ggsave("Dd2_slingshot_trajectory_umap.png", plot = p, width = 8, height = 6, dpi = 300)
lineages <- slingshot::slingLineages(Dd2_sce)
print(lineages)

pseudotime <- slingshot::slingPseudotime(Dd2_sce)
pseudotime[is.na(pseudotime)] <- 0

counts <- counts(Dd2_sce) 

cell_weights <- slingshot::slingCurveWeights(Dd2_sce) 

set.seed(123) 

Dd2_knots <- evaluateK(
  counts = counts,
  pseudotime = pseudotime,
  cellWeights = cell_weights,
  k = 3:10, 
  nGenes = 200  
)
print(Dd2_knots)
#9
Dd2_c_l_sce <- Dd2_sce[, colData(Dd2_sce)$SingleR_labels %in% c("committed", "early stalk", "late stalk")]
Dd2_e_l_sce <- Dd2_sce[, colData(Dd2_sce)$SingleR_labels %in% c("early stalk", "late stalk")]

Dd2_c_l_sce_pseudotime <- slingshot::slingPseudotime(Dd2_c_l_sce)
Dd2_c_l_sce_pseudotime[is.na(Dd2_c_l_sce_pseudotime)] <- 0
Dd2_e_l_sce_pseudotime <- slingshot::slingPseudotime(Dd2_e_l_sce)
Dd2_e_l_sce_pseudotime[is.na(Dd2_e_l_sce_pseudotime)] <- 0

Dd2_c_l_sce_counts <- counts(Dd2_c_l_sce)
Dd2_e_l_sce_counts <- counts(Dd2_e_l_sce)

Dd2_c_l_sce_cell_weights <- slingshot::slingCurveWeights(Dd2_c_l_sce) 
Dd2_e_l_sce_cell_weights <- slingshot::slingCurveWeights(Dd2_e_l_sce) 
Dd2_c_l_sce_cell_weights <- Dd2_c_l_sce_cell_weights[, 1, drop = FALSE]
Dd2_c_l_sce_cell_weights[Dd2_c_l_sce_cell_weights <= 0] <- 1e-6 
Dd2_e_l_sce_cell_weights <- Dd2_e_l_sce_cell_weights[, 1, drop = FALSE]
Dd2_e_l_sce_cell_weights[Dd2_e_l_sce_cell_weights <= 0] <- 1e-6 


Dd2_c_l_sce_counts <- as.matrix(Dd2_c_l_sce_counts)
Dd2_e_l_sce_counts <- as.matrix(Dd2_e_l_sce_counts)

Dd2_c_l_sce_pseudotime <- Dd2_c_l_sce_pseudotime[, 1, drop = FALSE]
Dd2_c_l_sce_pseudotime[is.na(Dd2_c_l_sce_pseudotime)] <- 0
Dd2_e_l_sce_pseudotime <- Dd2_e_l_sce_pseudotime[, 1, drop = FALSE]
Dd2_e_l_sce_pseudotime[is.na(Dd2_e_l_sce_pseudotime)] <- 0



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
  message("✅ ")
}

check_input_dimensions(Dd2_c_l_sce_counts, Dd2_c_l_sce_pseudotime, Dd2_c_l_sce_cell_weights)
check_input_dimensions(Dd2_e_l_sce_counts, Dd2_e_l_sce_pseudotime, Dd2_e_l_sce_cell_weights)

Dd2_c_l_sce_gam_fit <- fitGAM(
  counts = Dd2_c_l_sce_counts,
  pseudotime = Dd2_c_l_sce_pseudotime,
  cellWeights = Dd2_c_l_sce_cell_weights,
  nknots = 10
)
Dd2_e_l_sce_gam_fit <- fitGAM(
  counts = Dd2_e_l_sce_counts,
  pseudotime = as.vector(Dd2_e_l_sce_pseudotime),
  cellWeights = Dd2_e_l_sce_cell_weights,
  nknots = 9
)

saveRDS(Dd2_c_l_sce_gam_fit, file = "Dd2_c_l_sce_gam_fit.rds")
saveRDS(Dd2_e_l_sce_gam_fit, file = "Dd2_e_l_sce_gam_fit.rds")

Dd2_c_l_sce_gam_fit <- readRDS("Dd2_c_l_sce_gam_fit.rds")
Dd2_e_l_sce_gam_fit <- readRDS("Dd2_e_l_sce_gam_fit.rds")

degs_Dd2_c_l_sce1 <- startVsEndTest(
  models = Dd2_c_l_sce_gam_fit,
  global = TRUE,
  lineages = 1,
  l2fc = 1
)
degs_Dd2_e_l_sce1 <- startVsEndTest(
  models = Dd2_e_l_sce_gam_fit,
  global = TRUE,
  lineages = 1,
  l2fc = 1
)

de <- degs_Dd2_c_l_sce1
de$logfc <- de$logFClineage1
de$wald <- de$waldStat
de$pvalue <- de$pvalue_lineage1
de$diffexpressed <- "NO"
de$diffexpressed[de$logfc > 1 & de$wald > 30] <- "UP"
de$diffexpressed[de$logfc < -1 & de$wald > 30] <- "DOWN"
def <- de[which(de$diffexpressed != "NO"),]
def$geneid2 <- rownames(def)
committed_latestalk <- def


de <- degs_Dd2_e_l_sce1
de$logfc <- de$logFClineage1
de$wald <- de$waldStat
de$pvalue <- de$pvalue_lineage1
de$diffexpressed <- "NO"
de$diffexpressed[de$logfc > 1 & de$wald > 30] <- "UP"
de$diffexpressed[de$logfc < -1 & de$wald > 30] <- "DOWN"
def <- de[which(de$diffexpressed != "NO"),]
def$geneid2 <- rownames(def)
earlystalk_latestalk <- def

Dd2_fig_s5_B <- ggVennDiagram(list(committed_latestalk$geneid2,earlystalk_latestalk$geneid2), label_alpha = 0, label = "count", label_size = 7,
                             category.names = c("committed_latestalk","earlystalk_latestalk")
) +   ggplot2::scale_fill_gradient(low="white",high = "thistle")
print(Dd2_fig_s5_B)
ggsave("Dd2_fig_s5_B.png", plot = Dd2_fig_s5_B, width = 8, height = 6, dpi = 300)



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
  pvalueCutoff = 0.5,
  keyType = "SYMBOL",  
  ont = "BP", 
  qvalueCutoff = 0.5
)

Dd2_fig_s5_D <- clusterProfiler::dotplot(enr)
go_res <- enr@result
ggsave("Dd2_fig_s5_D.png", Dd2_fig_s5_D ,
       width=7, height=2.7, dpi=300, bg = "transparent")

Dd2_sex_sce <- Dd2_sce
lineages <- slingshot::slingLineages(Dd2_sex_sce)
print(lineages)

Dd2_sex_pseudotime <- slingshot::slingPseudotime(Dd2_sex_sce)
Dd2_sex_pseudotime[is.na(Dd2_sex_pseudotime)] <- 0
Dd2_sex_counts <- counts(Dd2_sex_sce)  
Dd2_sex_cell_weights <- slingshot::slingCurveWeights(Dd2_sex_sce) 
check_input_dimensions(Dd2_sex_counts, Dd2_sex_pseudotime, Dd2_sex_cell_weights)
head(slingshot::slingPseudotime(Dd2_sex_sce))
pseudotime_values <- as.data.frame(slingshot::slingPseudotime(Dd2_sex_sce))

cell_annotations <- colData(Dd2_sex_sce)$StageHR_v3lab
table(colData(Dd2_sex_sce)$StageHR_v3lab)
pseudotime_annotations_df <- data.frame(
  Cell = colnames(Dd2_sex_sce),
  Annotation = cell_annotations,
  Pseudotime = pseudotime_values[, 1]
)

early_stalk_pseudotime <- pseudotime_annotations_df$Pseudotime[pseudotime_annotations_df$Annotation == "early stalk"]
late_stalk_pseudotime <- pseudotime_annotations_df$Pseudotime[pseudotime_annotations_df$Annotation == "late stalk"]
branching_pseudotime <- pseudotime_annotations_df$Pseudotime[pseudotime_annotations_df$Annotation == "branching"]
eFem_pseudotime <- pseudotime_annotations_df$Pseudotime[pseudotime_annotations_df$Annotation == "early female"]
lFem_pseudotime <- pseudotime_annotations_df$Pseudotime[pseudotime_annotations_df$Annotation == "late female"]
eMal_pseudotime <- pseudotime_annotations_df$Pseudotime[pseudotime_annotations_df$Annotation == "early male"]
lMal_pseudotime <- pseudotime_annotations_df$Pseudotime[pseudotime_annotations_df$Annotation == "late male"]
summary(early_stalk_pseudotime)
summary(late_stalk_pseudotime)
summary(branching_pseudotime)
summary(eFem_pseudotime)
summary(lFem_pseudotime)
summary(eMal_pseudotime)
summary(lMal_pseudotime)


ggplot(pseudotime_annotations_df, aes(x = Pseudotime, color = Annotation)) +
  geom_density() +
  theme_minimal() +
  labs(title = "Annotations Along Pseudotime")

Dd2_sex_sce_gam_fit <- fitGAM(
  counts = Dd2_sex_counts,
  pseudotime = Dd2_sex_pseudotime,
  cellWeights = Dd2_sex_cell_weights,
  nknots = 10
)

saveRDS(Dd2_sex_sce_gam_fit, file = "Dd2_sex_sce_gam_fit.rds")
Dd2_sex_sce_gam_fit <- readRDS("Dd2_sex_sce_gam_fit.rds")

degs_Dd2_fsex <- startVsEndTest(
  models = Dd2_sex_sce_gam_fit,
  global = TRUE,
  lineages = 1,
  l2fc = 1
)
degs_Dd2_msex <- startVsEndTest(
  models = Dd2_sex_sce_gam_fit,
  global = TRUE,
  lineages = 2,
  l2fc = 1
)

lineages <- slingshot::slingLineages(Dd2_sex_sce)
print(lineages)
de <- degs_Dd2_msex
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


de <- degs_Dd2_fsex
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

committed_stalkf <- intersect(
  stalk_both_60$geneid2, cf_up_female
)
committed_stalkm <- intersect(
  stalk_both_60$geneid2, cm_up_male
)
genesused <- committed_stalkm
fig_s9_AB <- list()

for(j in genesused){
  p <- plotSmoothers(Dd2_sex_sce_gam_fit, counts(Dd2_sex_sce_gam_fit), gene = j, xlab = paste(j)) + geom_vline(xintercept=0.74, linetype="dashed", color = "red", linewidth=1) +
    geom_vline(xintercept=35, linetype="dashed", color = "red", linewidth=1) +
    geom_vline(xintercept=40, linetype="dashed", color = "red", linewidth=1) +
    geom_vline(xintercept=45, linetype="dashed", color = "red", linewidth=1) + theme_classic(base_size = 25) + theme(legend.position="none")
  
  plot_path <- paste0("/storage/work/skt5723/Single Cell Gametocyte stuff/Dd2_smoother_plots/", j, ".png")
  ggsave(plot_path, plot = p, width = 7, height = 5, dpi = 600, bg = "transparent")
  
  print(paste("Saved plot to:", plot_path))
  
  fig_s9_AB[[j]] <- p
}

sc_annotated_data <- readRDS("/storage/work/skt5723/Single Cell Gametocyte stuff/Plasmodium_Dd2.rds")
sc_annotated_data <- NormalizeData(sc_annotated_data)
sc_annotated_data <- FindVariableFeatures(sc_annotated_data, nfeatures = 1000)
sc_annotated_data <- ScaleData(sc_annotated_data)
sc_annotated_data <- RunPCA(sc_annotated_data, npcs = 50)
elbow_plot <- ElbowPlot(sc_annotated_data)
ggsave("elbow_plot.png", plot = elbow_plot, width = 6, height = 4, dpi = 300)
sc_annotated_data <- FindNeighbors(sc_annotated_data, dims = 1:20)
sc_annotated_data <- FindClusters(sc_annotated_data, resolution = 0.2) 
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
sc_annotated_data <- RunUMAP(
  sc_annotated_data,
  dims = 1:4,              
  n.components = 2,        
  min.dist = 0.3,           
  n.neighbors = 50,         
  spread = 1,              
  local.connectivity = 1,  
  a = 70,                   
  b = 2                  
)

dg_umap_plot <- DimPlot(sc_annotated_data, reduction = "umap", group.by = "SingleR_labels") +
  scale_color_manual(values = stage_colors) +
  ggplot2::ggtitle("UMAP with Automatic Clusters")
print(dg_umap_plot)
ggplot2::ggsave("Dd2_whole_umap_plot.png", plot = dg_umap_plot, width = 8, height = 6, dpi = 300)

male_genes <- cm_up_male
female_genes <- cf_up_female
table(sc_annotated_data@meta.data$my_ann)
late_stalk_cells <- subset(sc_annotated_data, subset = my_ann == "Stalk and Branching")

female_genes <- intersect(female_genes, rownames(late_stalk_cells))
male_genes <- intersect(male_genes, rownames(late_stalk_cells))

Idents(sc_annotated_data) <- "SingleR_labels"
lsm <- as.data.frame(rowSums(FetchData(late_stalk_cells, vars = male_genes, layer = "data")))
colnames(lsm) <- "lsm"
lsm$rank <- c(1:nrow(lsm))
mcells <- lsm[order(lsm$lsm, decreasing = TRUE),]
mcells$rank <- c(1:nrow(mcells))
mcellsx <- mcells[which(mcells$lsm >= mean(lsm$lsm)),] %>% rownames()
fcellsx <- tail(rownames(mcells), length(mcellsx))

lsf <- as.data.frame(rowSums(FetchData(late_stalk_cells, vars = female_genes, layer = "data")))
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

saveRDS(late_stalk_cells, file = "Dd2_classified_late_stalk_cells.rds")
saveRDS(sc_annotated_data, file = "Classified_Dd2.rds")
write.csv(late_stalk_cells@meta.data, "late_stalk_cells.csv", row.names = FALSE)

late_stalk_cells <- readRDS("/storage/work/skt5723/Single Cell Gametocyte stuff/Dd2_classified_late_stalk_cells.rds")
sc_annotated_data <- readRDS("/storage/work/skt5723/Single Cell Gametocyte stuff/Classified_Dd2.rds")
Idents(late_stalk_cells) <- late_stalk_cells$fm
late_stalk_cells <- NormalizeData(late_stalk_cells)

defm <- FindMarkers(late_stalk_cells,ident.1 = mcellsz, ident.2 = fcellsz, only.pos = FALSE, test.use = "MAST")
defm$geneid2 <- rownames(defm)

defm <- defm %>% mutate(genename2 =
                          case_when(genename == 'N/A' ~ geneid2,
                                    TRUE ~ genename))

defm$updown <- "other" 
defm[which(defm$avg_log2FC > 2 & defm$p_val_adj < 0.05),]$updown <- "up"
defm[which(defm$avg_log2FC < -2 & defm$p_val_adj < 0.05),]$updown <- "down"
table(defm$updown)
write.csv(defm, "/storage/work/skt5723/Single Cell Gametocyte stuff/Dd2_defm_results.csv", row.names = TRUE)
updown <- table(defm$updown)
p <- ggplot(data=defm, aes(x=avg_log2FC, y=-log10(p_val_adj), col=updown)) + geom_point(size=2.5) + theme_minimal()
fig_s9_C1 <- p + geom_vline(xintercept=c(-0.4, 0.4), col="red") +
  geom_hline(yintercept=1.30103, col="red") + 
  annotate("label", x =c(-1.5,0.8), y = 25, label = c(paste0(updown[[1]]," genes"), paste0(updown[[3]]," genes")), col=c("red","steelblue"), size = 6) +
  annotate("text", x =c(-1), y = 17, label = c(paste0("log2FC > 0.4, padj < 0.05")), size = 6) +
  theme(legend.position="none") + 
  theme(axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 16))
ggplot2::ggsave("Dd2_fig_s9_c1.png", plot = fig_s9_C1, width = 8, height = 6, dpi = 300)
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
ggplot2::ggsave("Dd2_fig_s9_c2.png", plot = p, width = 8, height = 6, dpi = 300)
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

ggplot2::ggsave("Dd2_fig_s9_c3.png", plot = p, width = 8, height = 6, dpi = 300)
rm(p)

scm <- as.data.frame(rowSums(FetchData(sc_annotated_data, vars = male_genes, layer = "data")))
colnames(scm) <- "scm"
scm$rank <- c(1:nrow(scm))
mcells <- scm[order(scm$scm, decreasing = TRUE),]
mcells$rank <- c(1:nrow(mcells))
mcellsx <- mcells[which(mcells$scm >= mean(scm$scm)),] %>% rownames()
fcellsx <- tail(rownames(mcells), length(mcellsx))

scf <- as.data.frame(rowSums(FetchData(sc_annotated_data, vars = female_genes, layer = "data")))
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

seur <- Dd2_sex_sce
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
Dd2_ElbowPlot <- ElbowPlot(sc_annotated_data)
print(Dd2_ElbowPlot)
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
rkg_umap_plot_myann <- DimPlot(Dd2_data_gametocytes, reduction = "umap", group.by = "my_ann") +
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
