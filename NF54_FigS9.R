library(Seurat)
library(dplyr)
library(ggplot2)
library(MAST)
library(tradeSeq)

setwd("/storage/work/skt5723/Single Cell Gametocyte stuff")
sc_annotated_data <- readRDS("/storage/work/skt5723/RKsc/PfE5_SingleR.rds")

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
ggplot2::ggsave("RK_whole_umap_plot.png", plot = dg_umap_plot, width = 8, height = 6, dpi = 300)

#genes
female_genes <- c("PF3D7-0320200", "PF3D7-0904000", "PF3D7-1416400", "PF3D7-1357600")  
male_genes <- c("PF3D7-1202300", "PF3D7-0521600", "PF3D7-0530500", "PF3D7-1468200")   

late_stalk_cells <- subset(sc_annotated_data, subset = SingleR_labels == "late stalk")
late_stalk_cells <- NormalizeData(late_stalk_cells)
late_stalk_cells <- ScaleData(late_stalk_cells)

female_genes <- intersect(female_genes, rownames(late_stalk_cells))
male_genes <- intersect(male_genes, rownames(late_stalk_cells))

#Dogga script
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

late_stalk_cells@meta.data$fm <- "other"
late_stalk_cells@meta.data[which(late_stalk_cells@meta.data$SingleR_labels == "late stalk"),]$fm <- "other late stak"
late_stalk_cells@meta.data[which(late_stalk_cells@meta.data$cellid %in% mcellsz),]$fm <- "male-like"
late_stalk_cells@meta.data[which(late_stalk_cells@meta.data$cellid %in% fcellsz),]$fm <- "female-like"
View(late_stalk_cells@meta.data)
table(late_stalk_cells$fm)

sc_annotated_data@meta.data$cellid <- rownames(sc_annotated_data@meta.data)
sc_annotated_data@meta.data$fm <- "other"
sc_annotated_data@meta.data[which(sc_annotated_data@meta.data$SingleR_labels == "late stalk"),]$fm <- "other late stak"
sc_annotated_data@meta.data[which(sc_annotated_data@meta.data$cellid %in% mcellsz),]$fm <- "male-like"
sc_annotated_data@meta.data[which(sc_annotated_data@meta.data$cellid %in% fcellsz),]$fm <- "female-like"
table(sc_annotated_data$fm)

saveRDS(late_stalk_cells, file = "classified_late_stalk_cells.rds")
saveRDS(sc_annotated_data, file = "Classified_Dogga.rds")
write.csv(late_stalk_cells@meta.data, "late_stalk_cells.csv", row.names = FALSE)

late_stalk_cells <- readRDS("/storage/work/skt5723/Single Cell Gametocyte stuff/classified_late_stalk_cells.rds")
sc_annotated_data <- readRDS("/storage/work/skt5723/Single Cell Gametocyte stuff/Classified_Dogga.rds")
Idents(late_stalk_cells) <- late_stalk_cells$fm
late_stalk_cells <- NormalizeData(late_stalk_cells)

defm <- FindMarkers(late_stalk_cells,ident.1 = mcellsz, ident.2 = fcellsz, only.pos = FALSE, test.use = "MAST")
defm$geneid2 <- rownames(defm)

defm <- defm %>% mutate(genename2 =
                          case_when(genename == 'N/A' ~ geneid2,
                                    TRUE ~ genename))

defm$updown <- "other" 
defm[which(defm$avg_log2FC > 0.4 & defm$p_val_adj < 0.1),]$updown <- "up"
defm[which(defm$avg_log2FC < -0.4 & defm$p_val_adj < 0.1),]$updown <- "down"
table(defm$updown)

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

ggplot2::ggsave("rk_fig_s9_c3.png", plot = p, width = 8, height = 6, dpi = 300)
rm(p)
