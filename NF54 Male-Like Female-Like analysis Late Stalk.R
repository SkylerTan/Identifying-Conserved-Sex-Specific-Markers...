BiocManager::install("MAST")
# Load necessary libraries
install.packages("Seurat")
install.packages("RANN", type="source")
BiocManager::install("DESeq2")
BiocManager::install("tradeSeq")
library(Seurat)
library(dplyr)
library(ggplot2)
library(MAST)
library(tradeSeq)

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
elbow_plot <- ElbowPlot(Dogga_data_seurat)
ggsave("elbow_plot.png", plot = elbow_plot, width = 6, height = 4, dpi = 300)
Dogga_data_seurat <- FindNeighbors(Dogga_data_seurat, dims = 1:20)
Dogga_data_seurat <- FindClusters(Dogga_data_seurat, resolution = 0.2)  
Dogga_data_seurat <- RunUMAP(
  Dogga_data_seurat,
  dims = 1:4,              
  n.components = 2,        
  min.dist = 0.3,          
  n.neighbors = 50,        
  spread = 1,              
  local.connectivity = 1,   
  a = 70,                  
  b = 1                   
)
dg_umap_plot <- DimPlot(Dogga_data_seurat, reduction = "umap", group.by = "stageHR") +
  ggplot2::ggtitle("UMAP with Automatic Clusters")
dg_umap_plot
ggplot2::ggsave("dg_whole_umap_plot.png", plot = dg_umap_plot, width = 8, height = 6, dpi = 300)

female_genes <- c("PF3D7-0110000", "PF3D7-1146800", "PF3D7-1416400", "PF3D7-1403200", "PF3D7-1466800")  
male_genes <- c("PF3D7-0422300", "PF3D7-0905300", "PF3D7-1122900", "PF3D7-1202300", "PF3D7-1227500")   

late_stalk_cells <- subset(Dogga_data_seurat, subset = stageHR == "late stalk")
late_stalk_cells <- NormalizeData(late_stalk_cells)
late_stalk_cells <- ScaleData(late_stalk_cells)

female_genes <- intersect(female_genes, rownames(late_stalk_cells))
male_genes <- intersect(male_genes, rownames(late_stalk_cells))

#Dogga script
Idents(Dogga_data_seurat) <- "stageHR"
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

combined_df <- data.frame(
  mcellsz = c(mcellsz, rep(NA, max(0, length(fcellsz) - length(mcellsz)))),
  fcellsz = c(fcellsz, rep(NA, max(0, length(mcellsz) - length(fcellsz))))
)
write.csv(combined_df, "mcellsz_fcellsz.csv", row.names = FALSE)


intersect(mcellsx, fcellsy) %>% length()
mcellsz <- setdiff(mcellsx,fcellsy)
fcellsz <- setdiff(fcellsy,mcellsx)
late_stalk_cells@meta.data$cellid <- rownames(late_stalk_cells@meta.data)

late_stalk_cells@meta.data$fm <- "other"
late_stalk_cells@meta.data[which(late_stalk_cells@meta.data$stageHR == "late stalk"),]$fm <- "other late stak"
late_stalk_cells@meta.data[which(late_stalk_cells@meta.data$cellid %in% mcellsz),]$fm <- "male-like"
late_stalk_cells@meta.data[which(late_stalk_cells@meta.data$cellid %in% fcellsz),]$fm <- "female-like"
View(late_stalk_cells@meta.data)
table(late_stalk_cells$fm)

Dogga_data_seurat@meta.data$cellid <- rownames(Dogga_data_seurat@meta.data)
Dogga_data_seurat@meta.data$fm <- "other"
Dogga_data_seurat@meta.data[which(Dogga_data_seurat@meta.data$stageHR == "late stalk"),]$fm <- "other late stak"
Dogga_data_seurat@meta.data[which(Dogga_data_seurat@meta.data$cellid %in% mcellsz),]$fm <- "male-like"
Dogga_data_seurat@meta.data[which(Dogga_data_seurat@meta.data$cellid %in% fcellsz),]$fm <- "female-like"
table(Dogga_data_seurat$fm)

saveRDS(late_stalk_cells, file = "classified_late_stalk_cells.rds")
saveRDS(Dogga_data_seurat, file = "Classified_Dogga.rds")
write.csv(late_stalk_cells@meta.data, "late_stalk_cells.csv", row.names = FALSE)

late_stalk_cells <- readRDS("/storage/work/skt5723/Single Cell Gametocyte stuff/classified_late_stalk_cells.rds")
Dogga_data_seurat <- readRDS("/storage/work/skt5723/Single Cell Gametocyte stuff/Classified_Dogga.rds")
Idents(late_stalk_cells) <- late_stalk_cells$Cell_Type
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
write.csv(defm, "/storage/work/skt5723/Single Cell Gametocyte stuff/DG_defm_results.csv", row.names = TRUE)
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
ggplot2::ggsave("fig_s9_c1.png", plot = fig_s9_C1, width = 8, height = 6, dpi = 300)
rm(p)

data <- cbind(Dogga_data_seurat@meta.data[,c("orig.ident","fm")], Dogga_data_seurat@reductions$umap@cell.embeddings)
data$label <- "select"
data[which(data$fm != "female-like"),]$label <- "other"
data$order <- ifelse(data$label=="other", 1, 2)
p <- ggplot(data) + geom_point(aes(x = umap_1, y = umap_2, color = label),size=3, shape=16) + 
  geom_point(data = subset(data, label == "select"), aes(x = umap_1, y = umap_2, color = label),size=3, shape=16) + 
  theme_classic() +
  ggtitle("Female Like Late Stalk") +
  theme(legend.title = element_text(colour="Black", size=12, face="bold"),
        legend.text = element_text(colour="Black", size=16))  +   scale_color_manual(values = c("other"="grey","select"="#a3767f")) + theme(legend.position="none") 

ggplot2::ggsave("fig_s9_c2.png", plot = p, width = 8, height = 6, dpi = 300)
rm(p)

data <- cbind(Dogga_data_seurat@meta.data[,c("orig.ident","fm")], Dogga_data_seurat@reductions$umap@cell.embeddings)
data$label <- "select"
data[which(data$fm != "male-like"),]$label <- "other"
data$order <- ifelse(data$label=="other", 1, 2)
p <- ggplot(data) + geom_point(aes(x = umap_1, y = umap_2, color = label),size=3, shape=16) + 
  geom_point(data = subset(data, label == "select"), aes(x = umap_1, y = umap_2, color = label),size=3, shape=16) + 
  theme_classic() +
  ggtitle("Male Like Late Stalk") +
  theme(legend.title = element_text(colour="Black", size=12, face="bold"),
        legend.text = element_text(colour="Black", size=16))  +   scale_color_manual(values = c("other"="grey","select"="#a3767f")) + theme(legend.position="none")

ggplot2::ggsave("fig_s9_c3.png", plot = p, width = 8, height = 6, dpi = 300)
rm(p)
