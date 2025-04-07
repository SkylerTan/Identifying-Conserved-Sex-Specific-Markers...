install.packages("Seurat")
BiocManager::install("SingleCellExperiment")

library(Seurat)
library(ggplot2)
library(ggrepel)
library(SingleCellExperiment)
library(Matrix)
library(Cairo)


sc_annotated_data <- readRDS("/storage/work/skt5723/RKsc/PfE5_SingleR.rds")
Dogga_data <- readRDS("/storage/group/mul27/default/scRNAseq_share/v3lab.ref.sce.RDS")
assay_data <- as(counts(Dogga_data), "dgCMatrix")
assay_data <- as(logcounts(Dogga_data), "CsparseMatrix")
meta_data <- as.data.frame(colData(Dogga_data))
Dogga_data_seurat <- CreateSeuratObject(
  counts = assay_data,  # Raw counts or normalized counts
  meta.data = meta_data
)

committed_cells <- subset(sc_annotated_data, subset = SingleR_labels == "committed")
sc_data <- RunUMAP(committed_cells, dims = 1:10)
DimPlot(sc_data, reduction = "umap", group.by = "my_ann") +
  ggplot2::ggtitle("UMAP with my_ann Labels")
umap_plot <- DimPlot(committed_cells, reduction = "umap", group.by = "my_ann") +
  ggplot2::ggtitle("UMAP with my_ann Labels")
ggplot2::ggsave("umap_plot1.png", plot = umap_plot, width = 8, height = 6, dpi = 300)

committed_cells <- NormalizeData(committed_cells)
committed_cells <- FindVariableFeatures(committed_cells)
committed_cells <- ScaleData(committed_cells)
committed_cells <- RunPCA(committed_cells, npcs = 30)
ElbowPlot(committed_cells)
committed_cells <- FindNeighbors(committed_cells, dims = 1:10)  # Use the first 10 PCs
committed_cells <- FindClusters(committed_cells, resolution = 0.5)  # Adjust resolution for finer or coarser clusters
committed_cells <- RunUMAP(committed_cells, dims = 1:10)  # Use the same dimensions as for clustering
umap_plot <- DimPlot(committed_cells, reduction = "umap", group.by = "seurat_clusters") +
  ggplot2::ggtitle("UMAP with Automatic Clusters")

ggplot2::ggsave("umap_plot2.png", plot = umap_plot, width = 8, height = 6, dpi = 300)

Dogga_committed_cells <- subset(Dogga_data_seurat, subset = stageHR == "committed")
Dogga_committed_cells <- NormalizeData(Dogga_committed_cells)
Dogga_committed_cells <- FindVariableFeatures(Dogga_committed_cells)
Dogga_committed_cells <- ScaleData(Dogga_committed_cells)
Dogga_committed_cells <- RunPCA(Dogga_committed_cells, npcs = 30)
CairoPNG("elbow_plot.png")
ElbowPlot(Dogga_committed_cells)
dev.off() 
Dogga_committed_cells <- FindNeighbors(Dogga_committed_cells, dims = 1:15)
Dogga_committed_cells <- FindClusters(Dogga_committed_cells, resolution = 0.25)  # Adjust resolution for finer or coarser clusters
Dogga_committed_cells <- RunUMAP(Dogga_committed_cells, dims = 1:10)
umap_plot2 <- DimPlot(Dogga_committed_cells, reduction = "umap", group.by = "seurat_clusters") +
  ggplot2::ggtitle("UMAP with Automatic Clusters")

ggplot2::ggsave("umap_plot2.png", plot = umap_plot2, width = 8, height = 6, dpi = 300)

umap_plot3 <- DimPlot(Dogga_committed_cells, reduction = "umap", group.by = "stageLR") +
  ggplot2::ggtitle("UMAP with Automatic Clusters")

ggplot2::ggsave("umap_plot3.png", plot = umap_plot3, width = 8, height = 6, dpi = 300)

#early Stalk
Dogga_early_stalk <- subset(Dogga_data_seurat, subset = stageHR == "early stalk")
Dogga_early_stalk <- NormalizeData(Dogga_early_stalk)
Dogga_early_stalk <- FindVariableFeatures(Dogga_early_stalk)
Dogga_early_stalk <- ScaleData(Dogga_early_stalk)
Dogga_early_stalk <- RunPCA(Dogga_early_stalk, npcs = 30)
CairoPNG("elbow_plot.png")
ElbowPlot(Dogga_early_stalk)
dev.off() 
Dogga_early_stalk <- FindNeighbors(Dogga_early_stalk, dims = 1:15)
Dogga_early_stalk <- FindClusters(Dogga_early_stalk, resolution = 0.25)  # Adjust resolution for finer or coarser clusters
Dogga_early_stalk <- RunUMAP(Dogga_early_stalk, dims = 1:10)
CairoPNG("feature_plot.png")
FeaturePlot(Dogga_early_stalk, features = "PF3D7-0617900")
dev.off() 
umap_early_stalk <- DimPlot(Dogga_early_stalk, reduction = "umap", group.by = "seurat_clusters") +
  ggplot2::ggtitle("UMAP with Automatic Clusters")

ggplot2::ggsave("umap_early_stalk.png", plot = umap_early_stalk, width = 8, height = 6, dpi = 300)

deg_results <- FindMarkers(Dogga_early_stalk, ident.1 = "0", ident.2 = "2", 
                           logfc.threshold = 0.25, min.pct = 0.1, test.use = "wilcox")
head(deg_results)

top_up <- head(deg_results[order(-deg_results$avg_log2FC), ], 10)
print("Top Overexpressed Genes:")
print(top_up)

top_down <- head(deg_results[order(deg_results$avg_log2FC), ], 10)
print("Top Underexpressed Genes:")
print(top_down)

Dogga_early_gametocytes <- subset(Dogga_early_stalk, idents = c("0", "1"))
Dogga_early_gametocytes <- NormalizeData(Dogga_early_gametocytes)
Dogga_early_gametocytes <- FindVariableFeatures(Dogga_early_gametocytes)
Dogga_early_gametocytes <- ScaleData(Dogga_early_gametocytes)
Dogga_early_gametocytes <- RunPCA(Dogga_early_gametocytes, npcs = 30)
CairoPNG("elbow_plot.png")
ElbowPlot(Dogga_early_gametocytes)
dev.off() 
Dogga_early_gametocytes <- FindNeighbors(Dogga_early_gametocytes, dims = 1:15)
Dogga_early_gametocytes <- FindClusters(Dogga_early_gametocytes, resolution = 0.25) 
Dogga_early_gametocytes <- RunUMAP(Dogga_early_gametocytes, dims = 1:10)
CairoPNG("feature_plot2.png")
FeaturePlot(Dogga_early_gametocytes, features = "PF3D7-0729800")
dev.off() 
umap_early_gametocytes <- DimPlot(Dogga_early_gametocytes, reduction = "umap", group.by = "seurat_clusters") +
  ggplot2::ggtitle("UMAP with Automatic Clusters")
ggplot2::ggsave("umap_early_gametocytes.png", plot = umap_early_stalk, width = 8, height = 6, dpi = 300)
Idents(Dogga_early_gametocytes) <- "seurat_clusters"

deg_results <- FindMarkers(Dogga_early_gametocytes, ident.1 = 0, ident.2 = 1, 
                           logfc.threshold = 0.25, min.pct = 0.1, test.use = "wilcox")

head(deg_results)
deg_results$gene <- rownames(deg_results)
deg_results$significance <- "Not Significant"
deg_results$significance[deg_results$avg_log2FC > 0.25 & deg_results$p_val_adj < 0.05] <- "Upregulated"
deg_results$significance[deg_results$avg_log2FC < -0.25 & deg_results$p_val_adj < 0.05] <- "Downregulated"
volcano_plot <- ggplot(deg_results, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
  geom_point(alpha = 0.7, size = 2) +  
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differentially Expressed Genes",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-Value") +
  theme(legend.title = element_blank())

top_genes <- deg_results[deg_results$p_val_adj < 0.01 & abs(deg_results$avg_log2FC) > 1, ]
volcano_plot + geom_text_repel(data = top_genes, aes(label = gene), size = 3, max.overlaps = 15)
ggsave("volcano_plot.png", plot = volcano_plot, width = 8, height = 6, dpi = 300)

genes_of_interest <- c("PF3D7-0729800", "PF3D7-1007100", "PF3D7-1475700", "PF3D7-1466800", "PF3D7-1146800",
                       "PF3D7-0110000", "PF3D7-1403200", "PF3D7-0905300", "PF3D7-1403200", "PF3D7-0513000", 
                       "PF3D7-0421900", "PF3D7-1429300", "PF3D7-0517500", "PF3D7-0904200", "PF3D7-1228100", 
                       "PF3D7-1122900", "PF3D7-1362600", "PF3D7-1021600")
available_genes <- genes_of_interest[genes_of_interest %in% rownames(Dogga_early_gametocytes)]
missing_genes <- setdiff(genes_of_interest, rownames(Dogga_early_gametocytes))
print(available_genes)
print(missing_genes)
scaled_data <- GetAssayData(Dogga_early_gametocytes, layer = "scale.data")
dim(scaled_data)
all_genes <- rownames(Dogga_early_gametocytes)
print(available_genes[available_genes %in% all_genes])
# Two genes are missing
genes_in_scaled_data <- available_genes[available_genes %in% rownames(scaled_data)]
print(genes_in_scaled_data)
if (length(genes_in_scaled_data) > 0) {
  scaled_data[genes_in_scaled_data, ] <- scaled_data[genes_in_scaled_data, ] * 888
} else {
  warning("Uh oh")
}
Dogga_early_gametocytes <- SetAssayData(Dogga_early_gametocytes, layer = "scale.data", new.data = scaled_data)
if (length(genes_in_scaled_data) > 0) {
  gene_expression <- scaled_data[genes_in_scaled_data, ]
  # View a summary of gene expression to check how many cells express each gene
  summary(gene_expression)
} else {
  warning("Uh oh")
}
Dogga_early_gametocytes <- RunUMAP(Dogga_early_gametocytes, reduction = "pca", dims = 1:10)  # Adjust dims as needed
Dogga_early_gametocytes <- FindClusters(Dogga_early_gametocytes, resolution = 0.25)
umap_early_gametocytes <- DimPlot(Dogga_early_gametocytes, reduction = "umap", group.by = "seurat_clusters") +
  ggtitle("UMAP")
ggsave("umap_early_gametocytes.png", plot = umap_early_gametocytes, width = 8, height = 6, dpi = 300)
CairoPNG("feature_plot3.png")
FeaturePlot(Dogga_early_gametocytes, features = "PF3D7-0110000")
dev.off() 
Idents(Dogga_early_gametocytes) <- "seurat_clusters"
deg_results <- FindMarkers(Dogga_early_gametocytes, ident.1 = 0, ident.2 = 1, 
                           logfc.threshold = 0.25, min.pct = 0.1, test.use = "wilcox")
head(deg_results)
male_genes <- c("PF3D7-1122900", "PF3D7-1228100", "PF3D7-0905300", "PF3D7-0729800", 
                "PF3D7-1007100", "PF3D7-1475700")
female_genes <- c("PF3D7-1466800", "PF3D7-1146800", "PF3D7-0110000", "PF3D7-1403200",
                  "PF3D7-1403200", "PF3D7-0513000", "PF3D7-0421900", "PF3D7-1429300",
                  "PF3D7-0517500", "PF3D7-0904200", "PF3D7-1362600", "PF3D7-1021600")
# Cells that express late stage sex genes
Dogga_expression_data <- GetAssayData(Dogga_early_gametocytes, layer = "data")[genes_of_interest, ]
Dogga_cells_expressing_genes <- colnames(Dogga_expression_data)[colSums(Dogga_expression_data > 0) > 0]
Dogga_sex_genes_subset <- subset(Dogga_early_gametocytes, cells = Dogga_cells_expressing_genes)
Dogga_sex_genes_subset <- NormalizeData(Dogga_sex_genes_subset)
Dogga_sex_genes_subset <- FindVariableFeatures(Dogga_sex_genes_subset)
Dogga_sex_genes_subset <- ScaleData(Dogga_sex_genes_subset)
Dogga_sex_genes_subset <- RunPCA(Dogga_sex_genes_subset, npcs = 30)
CairoPNG("elbow_plot.png")
ElbowPlot(Dogga_sex_genes_subset)
dev.off() 
Dogga_sex_genes_subset <- FindNeighbors(Dogga_sex_genes_subset, dims = 1:15)
Dogga_sex_genes_subset <- FindClusters(Dogga_sex_genes_subset, resolution = 0.25)  # Adjust resolution for finer or coarser clusters
Dogga_sex_genes_subset <- RunUMAP(Dogga_sex_genes_subset, dims = 1:10)
CairoPNG("feature_plot2.png")
FeaturePlot(Dogga_sex_genes_subset, features = "PF3D7-0729800")
dev.off() 
umap_subset_gametocytes <- DimPlot(Dogga_sex_genes_subset, reduction = "umap", group.by = "seurat_clusters") +
  ggplot2::ggtitle("UMAP with Automatic Clusters")
ggplot2::ggsave("umap_subset_gametocytes.png", plot = umap_subset_gametocytes, width = 8, height = 6, dpi = 300)
CairoPNG("feature_plot_female.png")
FeaturePlot(Dogga_sex_genes_subset, features = "PF3D7-0110000")
dev.off() 
CairoPNG("feature_plot_male.png")
FeaturePlot(Dogga_sex_genes_subset, features = "PF3D7-1007100")
dev.off() 

male_expression <- colMeans(GetAssayData(Dogga_sex_genes_subset, layer = "data")[male_genes, ])
female_expression <- colMeans(GetAssayData(Dogga_sex_genes_subset, layer = "data")[female_genes, ])

Dogga_sex_genes_subset$male_expression <- male_expression
Dogga_sex_genes_subset$female_expression <- female_expression

Dogga_sex_genes_subset$sex_gene_cluster <- ifelse(
  male_expression > female_expression, "Male", 
  ifelse(female_expression > male_expression, "Female", "Unclassified")
)
DimPlot(Dogga_sex_genes_subset, reduction = "umap", group.by = "sex_gene_cluster") +
  ggtitle("UMAP Clustered by Male/Female Gene Expression")

FeaturePlot(Dogga_sex_genes_subset, features = male_genes, reduction = "umap", ncol = 3)

FeaturePlot(Dogga_sex_genes_subset, features = female_genes, reduction = "umap", ncol = 4)

umap_plot <- DimPlot(Dogga_sex_genes_subset, reduction = "umap", group.by = "sex_gene_cluster") +
  ggtitle("UMAP Clustered by Male/Female Gene Expression")
ggsave("umap_sex_gene_cluster.png", plot = umap_plot, width = 8, height = 6, dpi = 300)


Idents(Dogga_sex_genes_subset) <- "seurat_clusters"
deg_results <- FindMarkers(Dogga_sex_genes_subset, ident.1 = 0, ident.2 = 1, 
                           logfc.threshold = 0.25, min.pct = 0.1, test.use = "wilcox")
head(deg_results)
deg_results$gene <- rownames(deg_results)
deg_results$significance <- "Not Significant"
deg_results$significance[deg_results$avg_log2FC > 0.25 & deg_results$p_val_adj < 0.05] <- "Upregulated"
deg_results$significance[deg_results$avg_log2FC < -0.25 & deg_results$p_val_adj < 0.05] <- "Downregulated"
volcano_plot <- ggplot(deg_results, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
  geom_point(alpha = 0.7, size = 2) +  
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
  theme_minimal() +
  labs(title = "Volcano Plot of DE Subset Genes",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-Value") +
  theme(legend.title = element_blank())
ggsave("volcano_plot2.png", plot = volcano_plot, width = 8, height = 6, dpi = 300)
write.csv(deg_results, file = "deg_results.csv", row.names = TRUE)
Idents(Dogga_sex_genes_subset) <- "sex_gene_cluster"
deg_results <- FindMarkers(Dogga_sex_genes_subset, ident.1 = "Male", ident.2 = "Female", 
                           logfc.threshold = 0.25, min.pct = 0.1, test.use = "wilcox")
head(deg_results)
write.csv(deg_results, file = "deg_results_male_vs_female.csv", row.names = TRUE)
deg_results$neg_log10_pval <- -log10(deg_results$p_val)

volcano_plot <- ggplot(deg_results, aes(x = avg_log2FC, y = neg_log10_pval)) +
  geom_point(aes(color = ifelse(avg_log2FC > 0.25 & p_val_adj < 0.05, "Upregulated in Male",
                                ifelse(avg_log2FC < -0.25 & p_val_adj < 0.05, "Upregulated in Female", "Not significant"))), 
             alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Upregulated in Male" = "red", "Upregulated in Female" = "blue", "Not significant" = "gray")) +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change (avg_log2FC)",
    y = "-Log10 P-value",
    title = "Volcano Plot: Male vs Female DEGs",
    color = "Expression"
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  # Significance threshold
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "black")  # Fold-change thresholds

ggsave("volcano_plot_male_vs_female.png", plot = volcano_plot, width = 8, height = 6, dpi = 300)
