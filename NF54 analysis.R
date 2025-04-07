install.packages("Cairo")
install.packages("seriation")
install.packages("RColorBrewer")
install.packages("tidyverse")
BiocManager::install("scmap")
install.packages("reticulate")
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
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
library(tidyverse)
library(DoubletFinder)
library(scmap)
setwd("/storage/work/skt5723/Single Cell Gametocyte stuff")
Rajat_data <- readRDS("/storage/work/skt5723/RKsc/PfE5_SingleR.rds")
Rajat_data <- NormalizeData(Rajat_data)
Rajat_data <- FindVariableFeatures(Rajat_data, nfeatures = 1000)
Rajat_data <- ScaleData(Rajat_data)
Rajat_data <- RunPCA(Rajat_data, npcs = 50)
nExp <- round(0.075 * ncol(Rajat_data))  # Adjust 0.075 based on expected doublet rate
RK_ElbowPlot <- ElbowPlot(Rajat_data)
print(RK_ElbowPlot)
Rajat_data <- FindNeighbors(Rajat_data, dims = 1:20)
Rajat_data <- FindClusters(Rajat_data, resolution = 0.25)  
Rajat_data <- RunUMAP(
  Rajat_data,
  dims = 1:20,              
  n.components = 2,       
  min.dist = 0.3,          
  n.neighbors = 50,        
  spread = 0.5,              
  local.connectivity = 2,  
  a = 100,                  
  b = 1                       
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

rk_umap_plot_myann <- DimPlot(Rajat_data, reduction = "umap", group.by = "my_ann") + 
  ggplot2::ggtitle("RK UMAP myann")
rk_umap_plot_myann
# pK identification (no ground-truth)
sweep.list <- paramSweep(Rajat_data, PCs = 1:20)
sweep.stats <- summarizeSweep(sweep.list)
bcmvn <- find.pK(sweep.stats)

optimal.pk <- sweep.stats %>%
  dplyr::filter(BCreal == max(BCreal)) %>%
  dplyr::select(pK)
optimal.pk <- as.numeric(as.character(optimal.pk[[1]]))
print(optimal.pk)

summary(sweep.stats$BCreal)

# Define the output directory for saving the doublet detection plot
output_dir_doublet <- "Doublet_finder"
dir.create(output_dir_doublet, showWarnings = FALSE)  # Create the directory if it doesn't exist

# Create the Bimodality Coefficient vs. pK plot
doublet_plot <- ggplot(sweep.stats, aes(x = as.numeric(as.character(pK)), y = BCreal)) +
  geom_line() +
  labs(x = "pK", y = "BCreal", title = "Bimodality Coefficient vs. pK")

# Save the plot in the "Doublet detection and removal" folder
ggsave(file.path(output_dir_doublet, "Bimodality_Coefficient_vs_pK.png"), plot = doublet_plot, width = 8, height = 6, dpi = 300)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic doublet proportion estimate
annotations <- Rajat_data@meta.data$seurat_clusters # use clusters as user-defined cell types
homotypic.prop <- modelHomotypic(annotations) # estimate proportion of homotypic doublets

# Set 24% as the doublet formation rate
multiplet_rates_10x <- data.frame('Multiplet_rate'= c(0.004, 0.008, 0.0160, 0.023, 0.031, 0.039, 0.046, 0.054, 0.061, 0.069, 0.076, 0.16, 0.24),
                                  'Loaded_cells' = c(800, 1600, 3200, 4800, 6400, 8000, 9600, 11200, 12800, 14400, 16000, 40000, 60000),
                                  'Recovered_cells' = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 20000, 30000))

print(multiplet_rates_10x)

#**Change doublet_rate based on cell loading/recovery**
doublet_rate <- 0.09

# Compute the number of expected doublets
nExp.poi <- round(doublet_rate * nrow(Rajat_data@meta.data)) 

# Adjust for homotypic doublet proportion
nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop)) 

# Print results for confirmation
print(paste("Estimated number of doublets (nExp.poi):", nExp.poi))
print(paste("Adjusted number of doublets (nExp.poi.adj):", nExp.poi.adj))

Rajat_data <- doubletFinder(seu = Rajat_data, 
                                        PCs = 1:20, 
                                        pK = optimal.pk,
                                        nExp = nExp.poi.adj)
umap_plot_doublet_finder <- DimPlot(Rajat_data, reduction = "umap", group.by = "DF.classifications_0.25_0.005_1456", 
                                    pt.size = 0.2, label = TRUE) + 
  labs(title = "UMAP with Single and Doublets") + 
  scale_color_manual(values = c("Singlet" = "#00BA38", "Doublet" = "#F8766D"))
umap_plot_doublet_finder
table(Rajat_data$DF.classifications_0.25_0.005_1456)

Doublet_score_cluster <- ggplot(Rajat_data@meta.data, aes(x = seurat_clusters, y = pANN_0.25_0.005_1456, fill = DF.classifications_0.25_0.005_1456)) +
  geom_violin(trim = FALSE) +
  labs(title = "Doublet Scores by Cluster",
       x = "Seurat Cluster",
       y = "Doublet Score (pANN)") +
  theme_minimal() +
  scale_fill_manual(values = c("Singlet" = "#32CD32", "Doublet" = "#F8766D")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

Doublet_score_cluster
# Save the plot as a .png file
ggsave(file.path(output_dir_doublet, "Doublet_score_seurat_cluster.png"), plot = Doublet_score_cluster, width = 15, height = 10, dpi = 300)

# Save Rajat_data object as an RDS file
saveRDS(Rajat_data, file = "Rajat_data_singlet.rds")

# Subset for single cells
Rajat_data <- readRDS("/storage/work/skt5723/Single Cell Gametocyte stuff/Rajat_data_singlet.rds")
Rajat_data_singlet <- subset(Rajat_data, subset = DF.classifications_0.25_0.005_1456 == "Singlet")

rk_umap_plot_labels <- DimPlot(Rajat_data_singlet, reduction = "umap", group.by = "SingleR_labels") + 
  scale_color_manual(values = stage_colors) + 
  ggplot2::ggtitle("RK UMAP singlet")
rk_umap_plot_labels
ggplot2::ggsave("rk_whole_umap_singlet.png", plot = rk_umap_plot_labels, width = 8, height = 6, dpi = 300)

Rajat_data_gametocytes <- subset(Rajat_data_singlet, subset = SingleR_labels %in% c("branching", "committed", "early female", "early male", "early stalk", "late female", "late male", "late stalk"))
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
Rajat_data_gametocytes <- FindNeighbors(Rajat_data_gametocytes, dims = 1:15)
Rajat_data_gametocytes <- FindClusters(Rajat_data_gametocytes, resolution = 0.25)  

Rajat_data_gametocytes <- RunUMAP(
  Rajat_data_gametocytes,
  dims = 1:6,              
  n.components = 2,       
  min.dist = 0.1,          
  n.neighbors = 20,        
  spread = 0.1,              
  local.connectivity = 2,  
  a = 100,                  
  b = 1                       
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

rk_umap_plot_myann <- DimPlot(Rajat_data_gametocytes, reduction = "umap", group.by = "my_ann") + 
  ggplot2::ggtitle("RK UMAP myann")
rk_umap_plot_labels <- DimPlot(Rajat_data_gametocytes, reduction = "umap", group.by = "SingleR_labels") + 
  scale_color_manual(values = stage_colors) + 
  ggplot2::ggtitle("RK UMAP SRL")
rk_umap_plot_labels
FeaturePlot(Rajat_data_gametocytes, features = "PF3D7-1466800") 
ggplot2::ggsave("rk_gam_umap_SRL.png", plot = rk_umap_plot_labels, width = 8, height = 6, dpi = 300)

v3lab.sce <- readRDS("/storage/group/mul27/default/scRNAseq_share/v3lab.ref.sce.RDS") 
md.v3lab.seur <- readRDS("/storage/group/mul27/default/scRNAseq_share/v3lab.ref.seur.metadata.RDS") 

seur <- Rajat_data_gametocytes
#Convert to SCE
sce <- as.SingleCellExperiment(seur)
logcounts(sce) <- log2(SingleCellExperiment::counts(sce) + 1)
rowData(sce)$feature_symbol <- rownames(sce)
sce <- sce[!duplicated(rownames(sce)), ]

##Projection
sce_mapped <- scmapCell(sce,list(
  v3lab = metadata(v3lab.sce)$scmap_cell_index
))

topsim_v3lab <- sce_mapped$v3lab$similarities[1, ]

# Assigning stage labels only when the top 3 most closest reference cells to a query cell all have a similar stage label and a cosine similarity greater than, say 0.4 for example
##Assigning stage labels using thresholds
sce_mapped_lbls <- scmapCell2Cluster(sce_mapped, list(
  as.character(SummarizedExperiment::colData(v3lab.sce)$stageHR)),w = 3, threshold = 0.4)

# Obtaining the reference cell that is closest to each query cell in neighbourhood space regardless of cosine similarity. For those where the top 10 most closest cells constitute several stages, we also obtain the second most similar stage.
sc_cells_v3lab <- sce_mapped$v3lab$cells %>% 
  as.data.frame()  %>%
  mutate(across(everything(), ~as.character(SummarizedExperiment::colData(v3lab.sce)$stageHR)[.])) %>% 
  dplyr::select(where(~n_distinct(.) > 1))

##Get stages for cell that have similar stage assignment by scmap
sc_cells_hom_v3lab <- sce_mapped$v3lab$cells %>% 
  as.data.frame()  %>%
  mutate(across(everything(), ~as.character(SummarizedExperiment::colData(v3lab.sce)$stageHR)[.])) %>% 
  dplyr::select(where(~n_distinct(.) == 1)) %>%
  slice_head(., n = 1)

##Get similaraties for cells that have similar stage assignment by scmap and row bind
scmap_cs_v3lab <- sce_mapped$v3lab$similarities %>% 
  as.data.frame()  %>%
  dplyr::select(colnames(sc_cells_hom_v3lab)) %>%
  slice_head(., n = 1) %>%
  mutate(across(everything(), ~as.character(.))) %>%
  bind_rows(.,sc_cells_hom_v3lab) %>% 
  t() %>% 
  as.data.frame() %>% 
  setNames(c('sim_v3lab','stage_v3lab'))

### Then
sc_top_cell_id_v3lab <- sce_mapped$v3lab$cells %>%
  .[1,] %>% 
  data.frame(.,fix.empty.names = TRUE) %>%
  mutate(across(`.`, ~colnames(v3lab.sce)[.])) %>% 
  setNames('v3lab_ref_cells')

# Adding annotations to MSC query datasets
stg_refnd_lvls_v3lab = c("early ring","late ring","early trophozoite","late trophozoite","early schizont","late schizont", "committed", "early stalk", "late stalk", "branching",
                         "early female","late female","early male",
                         "late male","unassigned")

##Adding annotations to query cells
seur_anotd <- sce_mapped_lbls$scmap_cluster_labs %>% 
  data.frame()  %>% 
  mutate(cell_bc = colnames(sce_mapped$v3lab$cells))%>%
  column_to_rownames("cell_bc") %>%
  dplyr::rename("Stage_scmap_v3lab" = v3lab) %>% 
  mutate(across(Stage_scmap_v3lab, ~droplevels(factor(., levels = stg_refnd_lvls_v3lab)))) %>%
  AddMetaData(seur, .) %>%
  AddMetaData(., scmap_cs_v3lab)%>%
  AddMetaData(., sc_top_cell_id_v3lab)

df_v3lab <- sce_mapped$v3lab$cells %>% 
  as.data.frame()  %>%
  mutate(across(everything(), ~as.character(SummarizedExperiment::colData(v3lab.sce)$stageHR)[.])) 
seur_anotd$v3lab_lab1 <- df_v3lab[1,] %>% as.character()
seur_anotd$v3lab_lab2 <- df_v3lab[2,] %>% as.character()
seur_anotd$v3lab_lab3 <- df_v3lab[3,] %>% as.character()

dfsim_v3lab <- sce_mapped$v3lab$cells %>% 
  as.data.frame()   
seur_anotd$v3lab_sim1 <- dfsim_v3lab[1,] %>% as.numeric()
seur_anotd$v3lab_sim2 <- dfsim_v3lab[2,] %>% as.numeric()
seur_anotd$v3lab_sim3 <- dfsim_v3lab[3,] %>% as.numeric()

qmd <- seur_anotd@meta.data
qmd <- dplyr::left_join(qmd, md.v3lab.seur, by = c("v3lab_ref_cells" = "BCrev2"))
seur_anotd@meta.data[,c("StageHR_v3lab",
                        "umap_v3lab_1","umap_v3lab_2","umap_v3lab_3")] <- qmd[,c("stageHR","UMAP_1","UMAP_2","UMAP_3")]
seur_anotd$topsim_v3lab <- topsim_v3lab

if(ncol(seur_anotd)<100){np <- 10; nn = 5} else{np <- 30; nn = 30}

seur_anotd <- seur_anotd %>% Seurat::NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 500, verbose = FALSE) %>% ScaleData(., features = rownames(.), verbose = FALSE) %>% RunPCA(., npcs = np, features = VariableFeatures(object = .), verbose = FALSE)

seur_anotd <- seur_anotd %>% RunUMAP(dims = 1:5, n.neighbors = nn, n.components = 3L, verbose=F, seed.use = 42) %>% FindNeighbors(dims = 1:20, verbose = FALSE) %>% FindClusters(resolution = 0.2)
#seur_anotd <- seur_anotd %>% RunUMAP(dims = 1:20, n.neighbors = nn, n.components = 3L, verbose=F, seed.use = 42) %>% FindNeighbors(dims = 1:20, verbose = FALSE) %>% FindClusters(resolution = 0.5)

p1_MCA <- UMAPPlot(seur_anotd, dims = c(2,3), label = TRUE) + theme(legend.position = "none")
p2_MCA <- UMAPPlot(seur_anotd, dims = c(2,3), group.by = "StageHR_v3lab", cols = stage_colors) + theme(legend.position = "none")
p3_MCA <- ggplot(md.v3lab.seur, aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(color='grey80') +
  geom_point(aes(x=umap_v3lab_2, y=umap_v3lab_3, colour=StageHR_v3lab), seur_anotd@meta.data) +
  theme_classic() + scale_color_manual(values = stage_colors) + guides(colour = guide_legend(override.aes = list(size=7))) + theme(legend.position = "right")

print(p1_MCA)
print(p2_MCA)
print(p3_MCA)

ggplot2::ggsave("MCA_UMAP_dims2_3_custom_colors_StageHR_v3lab.png", plot = p3_MCA, width = 8, height = 6, dpi = 300)

# Save full version of mapped object
saveRDS(seur_anotd, file = file.path(output_dir_RDS, 'W2_obj_MCA_SingleR_my_ann.rds'))

# Remove specific columns from the meta.data of seur_anotd
colnames(seur_anotd@meta.data)
seur_anotd@meta.data <- seur_anotd@meta.data[, !(colnames(seur_anotd@meta.data) %in% c("Stage_scmap_v3lab", "sim_v3lab", "stage_v3lab", "v3lab_ref_cells", "v3lab_lab1", "v3lab_lab2", "v3lab_lab3", "v3lab_sim1", "v3lab_sim2", "v3lab_sim3", "RNA_snn_res.0.2"))]

# Save reduced version of mapped object
saveRDS(seur_anotd, file = file.path(output_dir_RDS, 'W2_obj_MCA_SingleR_my_ann_clean.rds'))

#slingshot
RK_sce <- as.SingleCellExperiment(Rajat_data_singlet)
start_cluster <- c("committed")
RK_sce <- slingshot(RK_sce, reducedDim = "PCA", clusterLabels = RK_sce$SingleR_labels, start.clus = start_cluster)
saveRDS(RK_sce, file = "RK_sce3.RDS")
RK_sce <- readRDS("/storage/work/skt5723/Single Cell Gametocyte stuff/RK_sce3.RDS")
pfdb56 <- read.csv("/storage/work/skt5723/Single Cell Gametocyte stuff/GenesByTaxon_GeneModelDump.csv")
umap_embeddings <- reducedDim(RK_sce, "UMAP")

umap_data <- data.frame(
  UMAP1 = umap_embeddings[, 1],
  UMAP2 = umap_embeddings[, 2],
  Stage = RK_sce$SingleR_labels
)

#slingshot trajectory analysis
p <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = Stage)) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_manual(values = stage_colors) +  
  theme_minimal() +
  labs(title = "Slingshot Trajectory on UMAP", x = "UMAP1", y = "UMAP2")

p

slingshot_data <- SlingshotDataSet(RK_sce)
for (curve in slingshot_data@curves) {
  curve_points_pca <- curve$s[curve$ord, ]

  pca_embeddings <- reducedDim(RK_sce, "PCA") 
  nn_indices <- FNN::get.knnx(pca_embeddings, curve_points_pca, k = 1)$nn.index  
  curve_points_umap <- umap_embeddings[nn_indices, ] 
  
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

ggsave("RK_slingshot_trajectory_umap.png", plot = p, width = 8, height = 6, dpi = 300)
lineages <- slingshot::slingLineages(RK_sce)
print(lineages)
pseudotime <- slingPseudotime(slingshot_data) 
umap_data$pseudotime_lineage1 <- pseudotime[,1]  
umap_data$pseudotime_lineage2 <- pseudotime[,2] 
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
# Compute pseudotime
pseudotime <- slingPseudotime(RK_sce)

# Extract cells belonging to Lineage 1 (L1) and Lineage 2 (L2)
L1_cells <- colnames(RK_sce)[!is.na(pseudotime[,1])]  # Cells with non-NA pseudotime for L1
L2_cells <- colnames(RK_sce)[!is.na(pseudotime[,2])]  # Cells with non-NA pseudotime for L2

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

ggsave("RK_female_lineage_umap.png", plot = L1, width = 6, height = 5)
ggsave("RK_male_lineage_umap.png", plot = L2, width = 6, height = 5)

RK_c_l_sce <- RK_sce[, colData(RK_sce)$SingleR_labels %in% c("committed", "early stalk", "late stalk")]
RK_e_l_sce <- RK_sce[, colData(RK_sce)$SingleR_labels %in% c("early stalk", "late stalk")]

valid_cells_c_l <- rowSums(is.na(RK_c_l_pseudotime) | RK_c_l_pseudotime == 0) == 0 &
  rowSums(is.na(RK_c_l_cell_weights) | RK_c_l_cell_weights == 0) == 0

RK_c_l_pseudotime <- RK_c_l_pseudotime[valid_cells_c_l, , drop = FALSE]
RK_c_l_counts <- RK_c_l_counts[, valid_cells_c_l]
RK_c_l_cell_weights <- RK_c_l_cell_weights[valid_cells_c_l, , drop = FALSE]

valid_cells_e_l <- rowSums(is.na(RK_e_l_pseudotime) | RK_e_l_pseudotime == 0) == 0 &
  rowSums(is.na(RK_e_l_cell_weights) | RK_e_l_cell_weights == 0) == 0

RK_e_l_pseudotime <- RK_e_l_pseudotime[valid_cells_e_l, , drop = FALSE]
RK_e_l_counts <- RK_e_l_counts[, valid_cells_e_l]
RK_e_l_cell_weights <- RK_e_l_cell_weights[valid_cells_e_l, , drop = FALSE]

set.seed(123)  

RK_c_l_knots <- evaluateK(
  counts = RK_c_l_counts,
  pseudotime = RK_c_l_pseudotime,
  cellWeights = RK_c_l_cell_weights,
  k = 3:10,  
  nGenes = 200  
)
#8
RK_e_l_knots <- evaluateK(
  counts = RK_e_l_counts,
  pseudotime = RK_e_l_pseudotime,
  cellWeights = RK_e_l_cell_weights,
  k = 3:10,  
  nGenes = 200  
)
#9

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
    warning("Warning: 'pseudotime' or 'cellWeights' might be NULL.")
  }
  message("✅ ")
}

check_input_dimensions(RK_c_l_counts, RK_c_l_pseudotime, RK_c_l_cell_weights)
check_input_dimensions(RK_e_l_counts, RK_e_l_pseudotime, RK_e_l_cell_weights)

RK_c_l_sce_gam_fit <- fitGAM(
  counts = RK_c_l_counts,
  pseudotime = RK_c_l_pseudotime,
  cellWeights = RK_c_l_cell_weights,
  nknots = 8
)

RK_e_l_sce_gam_fit <- fitGAM(
  counts = RK_e_l_counts,
  pseudotime = RK_e_l_pseudotime,
  cellWeights = RK_e_l_cell_weights,
  nknots = 9
)

saveRDS(RK_c_l_sce_gam_fit, file = "RK_c_l_sce_gam_fit.rds")
saveRDS(RK_e_l_sce_gam_fit, file = "RK_e_l_sce_gam_fit.rds")

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
de$diffexpressed[de$logfc > 1 & de$wald > 15] <- "UP"
de$diffexpressed[de$logfc < -1 & de$wald > 15] <- "DOWN"
def <- de[which(de$diffexpressed != "NO"),]
def$geneid2 <- rownames(def)
committed_latestalk <- def


de <- degs_RK_e_l_sce1
de$logfc <- de$logFClineage1
de$wald <- de$waldStat
de$pvalue <- de$pvalue_lineage1
de$diffexpressed <- "NO"
de$diffexpressed[de$logfc > 1 & de$wald > 15] <- "UP"
de$diffexpressed[de$logfc < -1 & de$wald > 15] <- "DOWN"
def <- de[which(de$diffexpressed != "NO"),]
def$geneid2 <- rownames(def)
earlystalk_latestalk1 <- def

rk_fig_s5_B <- ggVennDiagram(list(committed_latestalk$geneid2,earlystalk_latestalk1$geneid2), label_alpha = 0, label = "count", label_size = 7,
                          category.names = c("committed_latestalk","earlystalk_latestalk")
) +   ggplot2::scale_fill_gradient(low="white",high = "thistle")
rk_fig_s5_B
ggsave("rk_fig_s5_B.png", plot = rk_fig_s5_B, width = 8, height = 6, dpi = 300)

stalk_both_60 <- data.frame(geneid2 = union(committed_latestalk$geneid2,earlystalk_latestalk1$geneid2),order = c(1:length(union(committed_latestalk$geneid2,earlystalk_latestalk1$geneid2))))
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
RK_sex_pseudotime <- slingshot::slingPseudotime(RK_sex_sce)
RK_sex_pseudotime[is.na(RK_sex_pseudotime)] <- 0
RK_sex_counts <- counts(RK_sex_sce)
RK_sex_cell_weights <- slingshot::slingCurveWeights(RK_sex_sce)  # Cell weights
lineages <- slingshot::slingLineages(RK_sex_sce)
print(lineages)
pseudotime_values <- as.data.frame(slingshot::slingPseudotime(RK_sex_sce))

cell_annotations <- colData(RK_sex_sce)$SingleR_labels

pseudotime_annotations_df <- data.frame(
  Cell = colnames(RK_sex_sce),
  Annotation = cell_annotations,
  Pseudotime = pseudotime_values[, 1] 
)
committed_pseudotime <- pseudotime_annotations_df$Pseudotime[pseudotime_annotations_df$Annotation == "committed"]
early_stalk_pseudotime <- pseudotime_annotations_df$Pseudotime[pseudotime_annotations_df$Annotation == "early stalk"]
late_stalk_pseudotime <- pseudotime_annotations_df$Pseudotime[pseudotime_annotations_df$Annotation == "late stalk"]
branching_pseudotime <- pseudotime_annotations_df$Pseudotime[pseudotime_annotations_df$Annotation == "branching"]
eFem_pseudotime <- pseudotime_annotations_df$Pseudotime[pseudotime_annotations_df$Annotation == "early female"]
Fem_pseudotime <- pseudotime_annotations_df$Pseudotime[pseudotime_annotations_df$Annotation == "late female"]
eMal_pseudotime <- pseudotime_annotations_df$Pseudotime[pseudotime_annotations_df$Annotation == "early male"]
Mal_pseudotime <- pseudotime_annotations_df$Pseudotime[pseudotime_annotations_df$Annotation == "late male"]
summary(committed_pseudotime)
summary(early_stalk_pseudotime)
summary(late_stalk_pseudotime)
summary(branching_pseudotime)
summary(eFem_pseudotime)
summary(Fem_pseudotime)
summary(eMal_pseudotime)
summary(Mal_pseudotime)

check_input_dimensions(RK_sex_counts, RK_sex_pseudotime, RK_sex_cell_weights)
table(colData(RK_sex_sce)$my_ann)

head(slingshot::slingPseudotime(RK_sex_sce))
pseudotime_values <- as.data.frame(slingshot::slingPseudotime(RK_sex_sce))

cell_annotations <- colData(RK_sex_sce)$my_ann


RK_sex_knots <- evaluateK(
  counts = RK_sex_counts,
  pseudotime = RK_sex_pseudotime,
  cellWeights = RK_sex_cell_weights,
  k = 3:10,  
  nGenes = 200  
)
print(RK_sex_knots)
ggplot(pseudotime_annotations_df, aes(x = Pseudotime, color = Annotation)) +
  geom_density() +
  theme_minimal() +
  labs(title = "Annotations Along Pseudotime")

RK_sex_sce_gam_fit <- fitGAM(
  counts = RK_sex_counts,
  pseudotime = RK_sex_pseudotime,
  cellWeights = RK_sex_cell_weights,
  nknots = 8
)

saveRDS(RK_sex_sce_gam_fit, file = "RK_sex_sce_gam_fit.rds")
RK_sex_sce_gam_fit <- readRDS("RK_sex_sce_gam_fit.rds")

degs_RK_fsex <- startVsEndTest(
  models = RK_sex_sce_gam_fit,
  global = TRUE,
  lineages = 1,
  l2fc = 1
)

genesused <- (stalk_both_60$geneid2)

fig_s9_AB <- list()
for(j in genesused){
  p <- plotSmoothers(RK_sex_sce_gam_fit, counts(RK_sex_sce_gam_fit), gene = j, xlab = paste(j)) + geom_vline(xintercept=0.74, linetype="dashed", color = "red", linewidth=1) +
    geom_vline(xintercept=25, linetype="dashed", color = "red", linewidth=1) +
    geom_vline(xintercept=45, linetype="dashed", color = "red", linewidth=1) + theme_classic(base_size = 25) + theme(legend.position="none")
  
  plot_path <- paste0("/storage/work/skt5723/Single Cell Gametocyte stuff/RK_other_smoother_plots/", j, ".png")
  ggsave(plot_path, plot = p, width = 7, height = 5, dpi = 600, bg = "transparent")
  
  print(paste("Saved plot to:", plot_path))
  
  fig_s9_AB[[j]] <- p
}
genesused <- c("PF3D7-1227500", "PF3D7-0930200", "PF3D7-0110000")
fig_s9_AB <- list()
curvesCols <- c("purple", "red")
for (j in genesused) {
  p <- plotSmoothers(
    RK_sex_sce_gam_fit, 
    counts(RK_sex_sce_gam_fit), 
    gene = j, 
    xlab = paste(j),
    curvesCols = curvesCols,
    sample = 0  # <- disables plotting individual data points
  ) +
    geom_vline(xintercept = 19, linetype = "dashed", color = "red", linewidth = 1) +
    geom_vline(xintercept = 48, linetype = "dashed", color = "red", linewidth = 1) +
    theme_classic(base_size = 25) +
    theme(legend.position = "none")
  
  plot_path <- paste0("/storage/work/skt5723/Single Cell Gametocyte stuff/RK_other_smoother_plots/", j, ".png")
  ggsave(plot_path, plot = p, width = 7, height = 5, dpi = 600, bg = "transparent")
  
  print(paste("Saved plot to:", plot_path))
  
  fig_s9_AB[[j]] <- p
}
fig_s9_AB

# Rajat
plot_gene_expression_diff <- function(gene_name, pseudotime1, pseudotime2, expression_data, lineage1, lineage2) {
  # Extract expression data for the gene
  gene_expression_1 <- expression_data[gene_name, ]
  gene_expression_2 <- expression_data[gene_name, ]
  # Create a data frame for Lineage 1
  df_lineage1 <- data.frame(
    pseudotime = pseudotime1,
    expression = gene_expression_1,
    cell_id = colnames(expression_data),
    lineage = lineage1
  )
  # Create a data frame for Lineage 2
  df_lineage2 <- data.frame(
    pseudotime = pseudotime2,
    expression = gene_expression_2,
    cell_id = colnames(expression_data),
    lineage = lineage2
  )
  geom_vline(xintercept=25, linetype="dashed", color = "red", linewidth=1) +
  # Combine the data frames for both lineages
  df <- rbind(df_lineage1, df_lineage2)
  # Plot expression for the gene in both lineages with only smooth lines (no points)
  ggplot(df, aes(x = pseudotime, y = expression, colour = lineage)) +
    geom_smooth(method = "loess", aes(colour = lineage), linetype = "solid", size = 1.2) +  # Smoothing without points
    labs(title = paste("Expression of", gene_name, "in Lineage 1 and Lineage 2"),
         x = Pseudotime, y = Expression) +
    theme_minimal() +
    scale_colour_manual(values = c("purple", "red"))  # Color for the lineages
}
plot_gene_expression_diff(PF3D7-11330000, )
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
sc_annotated_data <- subset(sc_annotated_data, subset = nCount_RNA > 1000 & nCount_RNA < 10000 & nFeature_RNA > 300 & nFeature_RNA < 3000 & percent.mt < 2)
sc_annotated_data <- subset(sc_annotated_data, my_ann != "Unknown")
sc_annotated_data <- NormalizeData(sc_annotated_data)
sc_annotated_data <- FindVariableFeatures(sc_annotated_data, nfeatures = 1000)
sc_annotated_data <- ScaleData(sc_annotated_data)
sc_annotated_data <- RunPCA(sc_annotated_data, npcs = 50)
ElbowPlot(sc_annotated_data)
sc_annotated_data <- FindNeighbors(sc_annotated_data, dims = 1:20)
sc_annotated_data <- FindClusters(sc_annotated_data, resolution = 0.2) 
sc_annotated_data <- RunUMAP(
  sc_annotated_data,
  dims = 1:19,               
  n.components = 3,          
  min.dist = 0.3,            
  n.neighbors = 50,          
  spread = 1,                
  local.connectivity = 1,    
  a = 70,                    
  b = 1                      
)
my_stage_colors <- c(
  "Early Trophozoite" = "#FFEDA1",
  "Late Trophozoite" = "#FFAE2C",  
  "Early Schizont" = "#b69e00",
  "Late Schizont" = "#7aad00",
  "Early Ring" = "#f8746b",
  "Late Ring" = "#df8c00",
  "Stalk" = "#00c08b",
  "Branching" = "#00b3f0",
  "early male" = "#fa87ed",  
  "Early and Late Female" = "#c67bff",
  "Late Male" = "#ff64b0"
)

RK_umap_plot <- DimPlot(sc_annotated_data, reduction = "umap", group.by = "my_ann") +
  scale_color_manual(values = my_stage_colors) +
  ggplot2::ggtitle("UMAP with Automatic Clusters")
RK_umap_plot
ggplot2::ggsave("RK_whole_umap_plot_new.png", plot = RK_umap_plot, width = 8, height = 6, dpi = 300)
print(lineages)

late_stalk_cells <- subset(sc_annotated_data, subset = my_ann == "Stalk")

femgenes <- c("PF3D7-1457100", "PF3D7-1133000", 
              "PF3D7-1435200", "PF3D7-1246400","PF3D7-0827200", "PF3D7-0816800", "PF3D7-0323800")
malgenes <- c("PF3D7-1341500", "PF3D7-0930200", "PF3D7-0615500", "PF3D7-0602000")
FeaturePlot(sc_annotated_data, features = "PF3D7-1357600") 

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
Idents(sc_annotated_data) <- sc_annotated_data$fm
late_stalk_cells <- NormalizeData(late_stalk_cells)

defm <- FindMarkers(sc_annotated_data,ident.1 = mcellsz, ident.2 = fcellsz, only.pos = FALSE, test.use = "MAST")
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
sc_annotated_data <- FindClusters(sc_annotated_data, resolution = 0.25) 
sc_annotated_data <- RunUMAP(
  sc_annotated_data,
  dims = 1:4,              
  n.components = 2,       
  min.dist = 0.3,          
  n.neighbors = 50,        
  spread = 0.5,              
  local.connectivity = 2,  
  a = 50,                  
  b = 2                       
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

FeaturePlot(sc_annotated_data, features = "PF3D7-1202300") 

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

