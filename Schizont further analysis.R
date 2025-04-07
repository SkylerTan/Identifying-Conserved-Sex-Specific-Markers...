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
ElbowPlot(Dogga_data_seurat)
Dogga_data_seurat <- FindNeighbors(Dogga_data_seurat, dims = 1:20)
Dogga_data_seurat <- FindClusters(Dogga_data_seurat, resolution = 0.2)
Dogga_data_seurat <- FindNeighbors(Dogga_data_seurat, dims = 1:20)
Dogga_data_seurat <- FindClusters(Dogga_data_seurat, resolution = 0.25)
DG_schizonts <- subset(Dogga_data_seurat,
                       subset = stageLR %in% c("early schizont", "late schizont"))
DG_rings <- subset(Dogga_data_seurat,
                   subset = stageLR %in% c("early ring", "late ring"))
table(Dogga_data_seurat@meta.data$stageLR)

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

DG_committed_schizonts <- subset(
  Dogga_data_seurat,
  subset =
    (stageLR %in% c("early schizont", "late schizont") &
       (GetAssayData(Dogga_data_seurat, layer = "counts")["PF3D7-1222600", ] > mean_expr_1222600 |
          GetAssayData(Dogga_data_seurat, layer = "counts")["PF3D7-0935400", ] > mean_expr_0935400 |
          GetAssayData(Dogga_data_seurat, layer = "counts")["PF3D7-1328800", ] > mean_expr_1328800)))
table(DG_committed_schizonts@meta.data$stageLR)
Idents(DG_committed_schizonts) <- DG_committed_schizonts$stageLR

DG_committed_rings <- subset(
  Dogga_data_seurat,
  subset = 
    (stageLR %in% c("early ring", "late ring") &
       (GetAssayData(Dogga_data_seurat, layer = "counts")["PF3D7-1222600", ] > mean_expr_1222600_rings |
          GetAssayData(Dogga_data_seurat, layer = "counts")["PF3D7-0406500", ] > mean_expr_0406500 |
          GetAssayData(Dogga_data_seurat, layer = "counts")["PF3D7-0936500", ] > mean_expr_0936500)))
table(DG_committed_rings@meta.data$stageLR)
Idents(DG_committed_rings) <- DG_committed_rings$stageLR

#labels
#schizonts
DG_schizonts$cell_type <- "asexual"
cells_committed <- colnames(DG_committed_schizonts)
DG_schizonts$cell_type[cells_committed] <- "committed"
table(DG_schizonts@meta.data$cell_type)
Idents(DG_schizonts) <- DG_schizonts$cell_type
#rings
DG_rings$cell_type <- "asexual"
cells_committed <- colnames(DG_committed_rings)
DG_rings$cell_type[cells_committed] <- "committed"
table(DG_rings@meta.data$cell_type)
Idents(DG_rings) <- DG_rings$cell_type

schizont_deg_results <- FindMarkers(DG_schizonts, ident.1 = "committed", ident.2 = "asexual", 
                           test.use = "wilcox", logfc.threshold = 0.25, min.pct = 0.1)
ring_deg_results <- FindMarkers(DG_rings, ident.1 = "committed", ident.2 = "asexual", 
                                    test.use = "wilcox", logfc.threshold = 0.25, min.pct = 0.1)
#Schizontcano plot
schizont_deg_results$log_pval <- -log10(schizont_deg_results$p_val_adj)
logFC_threshold <- 0.5 
pval_threshold <- 0.05

schizont_deg_results$significance <- "Not Significant"
schizont_deg_results$significance[schizont_deg_results$p_val_adj < pval_threshold & schizont_deg_results$avg_log2FC > logFC_threshold] <- "Upregulated"
schizont_deg_results$significance[schizont_deg_results$p_val_adj < pval_threshold & schizont_deg_results$avg_log2FC < -logFC_threshold] <- "Downregulated"

# Plot
ggplot(schizont_deg_results, aes(x = avg_log2FC, y = log_pval, color = significance)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("Downregulated" = "blue", "Not Significant" = "gray", "Upregulated" = "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot of DEGs",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted p-value") +
  geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black")
s_upregulated_genes <- schizont_deg_results[schizont_deg_results$p_val_adj < pval_threshold & schizont_deg_results$avg_log2FC > logFC_threshold, ]
head(s_upregulated_genes)

s_upregulated_gene_names <- rownames(s_upregulated_genes)
print(s_upregulated_gene_names)
#Ringcano plot
ring_deg_results$log_pval <- -log10(ring_deg_results$p_val_adj)

ring_deg_results$significance <- "Not Significant"
ring_deg_results$significance[ring_deg_results$p_val_adj < pval_threshold & ring_deg_results$avg_log2FC > logFC_threshold] <- "Upregulated"
ring_deg_results$significance[ring_deg_results$p_val_adj < pval_threshold & ring_deg_results$avg_log2FC < -logFC_threshold] <- "Downregulated"

# Plot
ggplot(ring_deg_results, aes(x = avg_log2FC, y = log_pval, color = significance)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("Downregulated" = "blue", "Not Significant" = "gray", "Upregulated" = "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot of DEGs",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted p-value") +
  geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black")

r_upregulated_genes <- ring_deg_results[ring_deg_results$p_val_adj < pval_threshold & ring_deg_results$avg_log2FC > logFC_threshold, ]

head(r_upregulated_genes)

r_upregulated_gene_names <- rownames(r_upregulated_genes)
print(r_upregulated_gene_names)
DG_early_committed_schizonts <- subset(DG_committed_schizonts, subset= (stageHR %in% "early schizont"))
DG_early_committed_schizonts <- FindVariableFeatures(DG_early_committed_schizonts, nfeatures = 1000)
DG_early_committed_schizonts <- ScaleData(DG_early_committed_schizonts)
DG_early_committed_schizonts <- RunPCA(DG_early_committed_schizonts, npcs = 50)
ElbowPlot(DG_early_committed_schizonts)
DG_early_committed_schizonts <- FindNeighbors(DG_early_committed_schizonts, dims = 1:20)
DG_early_committed_schizonts <- FindClusters(DG_early_committed_schizonts, resolution = 0.2)
DG_early_committed_schizonts <- FindNeighbors(DG_early_committed_schizonts, dims = 1:20)
DG_early_committed_schizonts <- FindClusters(DG_early_committed_schizonts, resolution = 0.1)
DG_early_committed_schizonts <- RunUMAP(
  DG_early_committed_schizonts,
  dims = 1:20,               
  n.components = 3,          
  min.dist = 0.3,            
  n.neighbors = 50,          
  spread = 1,                
  local.connectivity = 1,    
  a = 70,                    
  b = 1                      
)
schizont_umap_plot <- DimPlot(DG_early_committed_schizonts, reduction = "umap", group.by = "seurat_clusters") +
  ggplot2::ggtitle("UMAP with Automatic Clusters")
schizont_ann_plot <- DimPlot(DG_early_committed_schizonts, reduction = "umap", group.by = "stageHR") +
  ggplot2::ggtitle("UMAP with Automatic Clusters")
schizont_umap_plot
table(DG_early_committed_schizonts@meta.data$seurat_clusters)
schizont_ann_plot

deg_0_vs_2 <- FindMarkers(DG_early_committed_schizonts, ident.1 = 0, ident.2 = 2)
deg_1_vs_2 <- FindMarkers(DG_early_committed_schizonts, ident.1 = 1, ident.2 = 2)

deg_0_vs_2$comparison <- "Cluster 0 vs 2"
deg_1_vs_2$comparison <- "Cluster 1 vs 2"

combined_deg <- rbind(deg_0_vs_2, deg_1_vs_2)
combined_deg$gene <- rownames(combined_deg)  # Add gene names as a column
combined_deg$Significance <- "Not Significant"
combined_deg$Significance[combined_deg$p_val_adj < 1e-15 & combined_deg$avg_log2FC > 1.5] <- "Upregulated"
combined_deg$Significance[combined_deg$p_val_adj < 1e-15 & combined_deg$avg_log2FC < -1.5] <- "Downregulated"

ggplot(combined_deg, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significance)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~ comparison) +  # Separate plots for each comparison
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
  theme_minimal() +
  labs(title = "Volcano Plot of DEGs", x = "Log2 Fold Change", y = "-Log10 Adjusted p-value") +
  theme(legend.title = element_blank())

sig_deg_0_vs_2 <- subset(deg_0_vs_2, p_val_adj < 1e-15)
sig_deg_2_vs_1 <- subset(deg_1_vs_2, p_val_adj < 1e-15)

common_genes <- intersect(rownames(sig_deg_0_vs_1), rownames(sig_deg_2_vs_1))
