library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(Matrix.utils)
# Load the syogren dataset
syogren_data <- Read10X(data.dir = "C:/Users/LSCS/Desktop/bbb")
# Initialize the Seurat object with the raw (non-normalized data).
#syogren <- CreateSeuratObject(counts = syogren_data, project = "syogren", min.cells = 3, min.features = 200)
syogren <- CreateSeuratObject(counts = syogren_data, project = "syogren")
cell_batch.tsv <- read.delim("C:/Users/LSCS/Desktop/bbb/cell_batch.tsv.gz")
cell_batch <- read.delim("C:/Users/LSCS/Desktop/bbb/cell_batch.tsv")
syogren <- AddMetaData(object = syogren, metadata = cell_batch.tsv$batch, col.name = "batch")
syogren <- AddMetaData(object = syogren, metadata = cell_batch$control, col.name = "control")
syogren

syogren[["percent.mt"]] <- PercentageFeatureSet(syogren, pattern = "^MT-")
VlnPlot(syogren, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3, pt.size=0)
syogren@active.ident <- as.factor(syogren$orig.ident)
syogren <- subset(syogren, subset = nFeature_RNA > 300 & nFeature_RNA < 3500 & percent.mt < 25)
syogren

plot1 <- FeatureScatter(syogren, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(syogren, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

syogren.list <- SplitObject(syogren, split.by = "control")

# normalize and identify variable features for each dataset independently
syogren.list <- lapply(X = syogren.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = syogren.list)
syogren.anchors <- FindIntegrationAnchors(object.list = syogren.list, anchor.features = features)

syogren.combined <- IntegrateData(anchorset = syogren.anchors)

DefaultAssay(syogren.combined) <- "integrated"



syogren.combined <- ScaleData(syogren.combined, verbose = FALSE)
syogren.combined <- RunPCA(syogren.combined, npcs = 100, verbose = FALSE)
syogren.combined <- RunUMAP(syogren.combined, reduction = "pca", dims = 1:50)
syogren.combined <- FindNeighbors(syogren.combined, reduction = "pca", dims = 1:50)
syogren.combined <- FindClusters(syogren.combined, resolution = 0.2)

syogren.combined <- JackStraw(syogren.combined, num.replicate = 100,dims = 50)
syogren.combined <- ScoreJackStraw(syogren.combined, dims = 1:50)

JackStrawPlot(syogren.combined, dims = 1:50)

# Visualization
p1 <- DimPlot(syogren.combined, reduction = "umap", group.by = "control")
p2 <- DimPlot(syogren.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

for (i in 30){
  
  syogren.combined <- FindNeighbors(syogren.combined, reduction = "pca", dims = 1:i)
  syogren.combined <- FindClusters(syogren.combined, resolution = 0.2)
  syogren.combined <- RunUMAP(syogren.combined, reduction = "pca", dims = 1:i)
  
  s1 <- DimPlot(syogren.combined, reduction = "umap", label = TRUE, split.by = "control")
  s2 <- FeaturePlot(syogren.combined, features = c("CD3D","CD247","CD79A","CD14","CD68","CD1C"))
  #png(paste0(i,"_new_plot4.png"),width=6000,height=4000,res=500)
  #print(s1+s2)
  #Sys.sleep(20)
  #s1 + s2
  #dev.off()
}


VlnPlot(syogren.combined, features = c("CD3D","CD247","CD79A","CD14","CD68","CD1C"), group.by = "control", pt.size = 0)
syogren.markers <- FindAllMarkers(syogren.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
syogren.markers %>%
  group_by(cluster) %>%
  top_n(n = 3, wt = avg_log2FC) -> top10

DoHeatmap(syogren.combined, features = top10$gene) + NoLegend()


syogren.combined <- RenameIdents(syogren.combined, `0` = "NKs", `1` = "T Cell", `2` = "T Cell",
                                `3` = "T Cell", `4` = "Monocyte", `5` = "T Cell", `6` = "B Cell",
                                `7` = "T Cell", `8`= "T Cell", `9` = "Macro", `10` = 'platelets', `11` = "T Cell",
                                `12` = "DC", `13` = "T Cell", `14` = "unknown"
                                  )

DimPlot(syogren.combined, label = TRUE)

syogren.combined@meta.data$dim_30 <- syogren.combined@active.ident


syogren.combined@active.ident

DimPlot(syogren.combined, group.by = "dim_30", label = FALSE)
write.table(syogren.combined$dim_30,"annotation.tsv",quote = F,sep = "\t")

saveRDS(syogren.combined, file = "syogren.combined.rds")

Tcell <- subset(syogren.combined, dim_30 == "T Cell")

Tcell  <- ScaleData(Tcell , verbose = FALSE)
Tcell  <- RunPCA(Tcell , npcs = 50, verbose = FALSE)
Tcell  <- RunUMAP(Tcell , reduction = "pca", dims = 1:40)
Tcell  <- FindNeighbors(Tcell , reduction = "pca", dims = 1:40)
Tcell  <- FindClusters(Tcell , resolution = 0.2)

p1 <- DimPlot(Tcell, reduction = "umap", group.by = "control")
p2 <- DimPlot(Tcell, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

s1 <- DimPlot(Tcell, reduction = "umap", label = TRUE, split.by = "control")
s2 <- FeaturePlot(Tcell, features = c("CCR7","CXCR4","IL17RA","GZMB","CD8A","FOXP3"))

s1 + s2

Tcell.markers <- FindAllMarkers(Tcell, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Tcell.markers %>%
  group_by(cluster) %>%
  top_n(n = 3, wt = avg_log2FC) -> top10_tcell

DoHeatmap(Tcell, features = top10_tcell$gene) + NoLegend()

Tcell <- RenameIdents(Tcell, `0` = "CD4_Naive", `1` = "CD4_effector_memory", `2` = "CD8_CTLs_1", `3` = "CD8_Naive",
                      `4` = "CD8_CTLs_2", `5` = "CD4_CTLs_1", `6` = "CD4_CTLs_2", `7` = "Treg",
                      `8`= "Hemoglobin", `9` = "TRAV13-2_CD8_Tcell", `10` = 'unknown', `11` = "TRBV7-9_T_CD8_cell",
                      `12`="unknown",`13` = "TRBV10-2_CD8_Tcell", `14`="Bcell", `15`="unknown")

DimPlot(Tcell, label = TRUE)

for (i in 40){
  
  Tcell <- FindNeighbors(Tcell, reduction = "pca", dims = 1:i)
  Tcell <- FindClusters(Tcell, resolution = 0.2)
  Tcell <- RunUMAP(Tcell, reduction = "pca", dims = 1:i)
  
  s1 <- DimPlot(Tcell, reduction = "umap", label = TRUE, split.by = "control")
  s2 <- DimPlot(Tcell, reduction = "umap", label = TRUE, group.by = "dim_40")
  png(paste0(i,"_Tcell_plot2.png"),width=6000,height=4000,res=500)
  print(s1+s2)
  #Sys.sleep(20)
  s1 + s2
  dev.off()
}


Tcell@meta.data$dim_40 <- Tcell@active.ident

DimPlot(Tcell, group.by = "dim_40", label = T)
write.table(Tcell$dim_40,"Tcell_annotation_40.tsv",quote = F,sep = "\t")

saveRDS(Tcell, file = "Tcell_40.rds")
DefaultAssay(Tcell) <- "RNA"
Tcell <- FindVariableFeatures(Tcell, selection.method = "vst", nfeatures = 2000)
write.table(Tcell@assays$integrated@var.features, quote = F, row.names = F, "gene_list.txt")

pSS_Tcell <- subset(Tcell, control == "pSS")
HC_Tcell <- subset(Tcell, control == "HC")
DefaultAssay(pSS_Tcell) <- "RNA"
DefaultAssay(HC_Tcell) <- "RNA"
pSS_Tcell <- FindVariableFeatures(pSS_Tcell, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(pSS_Tcell), 10)
HC_Tcell <- FindVariableFeatures(HC_Tcell, selection.method = "vst", nfeatures = 2000)
top10_HC <- head(VariableFeatures(HC_Tcell), 10)

pSS_cluster <- as.data.frame(pSS_Tcell@meta.data$dim_40)
pSS_cluster$`pSS_Tcell@meta.data$dim_40` <- gsub("CD4_Naive","0",pSS_cluster$`pSS_Tcell@meta.data$dim_40`)
pSS_cluster$`pSS_Tcell@meta.data$dim_40` <- gsub("CD4_effector_memory","1",pSS_cluster$`pSS_Tcell@meta.data$dim_40`)
pSS_cluster$`pSS_Tcell@meta.data$dim_40` <- gsub("CD8_CTLs_1","2",pSS_cluster$`pSS_Tcell@meta.data$dim_40`)
pSS_cluster$`pSS_Tcell@meta.data$dim_40` <- gsub("CD8_Naive","3",pSS_cluster$`pSS_Tcell@meta.data$dim_40`)
pSS_cluster$`pSS_Tcell@meta.data$dim_40` <- gsub("CD8_CTLs_2","4",pSS_cluster$`pSS_Tcell@meta.data$dim_40`)
pSS_cluster$`pSS_Tcell@meta.data$dim_40` <- gsub("CD4_CTLs_1","5",pSS_cluster$`pSS_Tcell@meta.data$dim_40`)
pSS_cluster$`pSS_Tcell@meta.data$dim_40` <- gsub("CD4_CTLs_2","6",pSS_cluster$`pSS_Tcell@meta.data$dim_40`)
pSS_cluster$`pSS_Tcell@meta.data$dim_40` <- gsub("Treg","7",pSS_cluster$`pSS_Tcell@meta.data$dim_40`)
pSS_cluster$`pSS_Tcell@meta.data$dim_40` <- gsub("Hemoglobin","8",pSS_cluster$`pSS_Tcell@meta.data$dim_40`)
pSS_cluster$`pSS_Tcell@meta.data$dim_40` <- gsub("TRAV13-2_CD8_Tcell","9",pSS_cluster$`pSS_Tcell@meta.data$dim_40`)
pSS_cluster$`pSS_Tcell@meta.data$dim_40` <- gsub("unknown","10",pSS_cluster$`pSS_Tcell@meta.data$dim_40`)
pSS_cluster$`pSS_Tcell@meta.data$dim_40` <- gsub("TRBV7-9_T_CD8_cell","11",pSS_cluster$`pSS_Tcell@meta.data$dim_40`)
pSS_cluster$`pSS_Tcell@meta.data$dim_40` <- gsub("unknown","12",pSS_cluster$`pSS_Tcell@meta.data$dim_40`)
pSS_cluster$`pSS_Tcell@meta.data$dim_40` <- gsub("TRBV10-2_CD8_Tcell","13",pSS_cluster$`pSS_Tcell@meta.data$dim_40`)
pSS_cluster$`pSS_Tcell@meta.data$dim_40` <- gsub("Bcell","14",pSS_cluster$`pSS_Tcell@meta.data$dim_40`)
pSS_cluster$`pSS_Tcell@meta.data$dim_40` <- gsub("unknown","15",pSS_cluster$`pSS_Tcell@meta.data$dim_40`)

pSS_UMAP <- DimPlot(pSS_Tcell, reduction="umap", label = TRUE)$data
pSS_UMAP$ident <- gsub("CD4_Naive","0",pSS_UMAP$ident)
pSS_UMAP$ident <- gsub("CD4_effector_memory","1",pSS_UMAP$ident)
pSS_UMAP$ident <- gsub("CD8_CTLs_1","2",pSS_UMAP$ident)
pSS_UMAP$ident <- gsub("CD8_Naive","3",pSS_UMAP$ident)
pSS_UMAP$ident <- gsub("CD8_CTLs_2","4",pSS_UMAP$ident)
pSS_UMAP$ident <- gsub("CD4_CTLs_1","5",pSS_UMAP$ident)
pSS_UMAP$ident <- gsub("CD4_CTLs_2","6",pSS_UMAP$ident)
pSS_UMAP$ident <- gsub("Treg","7",pSS_UMAP$ident)
pSS_UMAP$ident <- gsub("Hemoglobin","8",pSS_UMAP$ident)
pSS_UMAP$ident <- gsub("TRAV13-2_CD8_Tcell","9",pSS_UMAP$ident)
pSS_UMAP$ident <- gsub("unknown","10",pSS_UMAP$ident)
pSS_UMAP$ident <- gsub("TRBV7-9_T_CD8_cell","11",pSS_UMAP$ident)
pSS_UMAP$ident <- gsub("unknown","12",pSS_UMAP$ident)
pSS_UMAP$ident <- gsub("TRBV10-2_CD8_Tcell","13",pSS_UMAP$ident)
pSS_UMAP$ident <- gsub("Bcell","14",pSS_UMAP$ident)
pSS_UMAP$ident <- gsub("unknown","15",pSS_UMAP$ident)


write.table(colnames(pSS_Tcell@assays$integrated),"cell_id.txt")
write.table(pSS_Tcell@assays$integrated@var.features, quote = F, row.names = F, "gene_list.txt")
write.table(pSS_cluster$`pSS_Tcell@meta.data$dim_40`, "clustering2.txt")
write.table(pSS_UMAP,"umap.txt")

HC_cluster <- as.data.frame(HC_Tcell@meta.data$dim_40)
HC_cluster$`HC_Tcell@meta.data$dim_40` <- gsub("CD4_Naive","0",HC_cluster$`HC_Tcell@meta.data$dim_40`)
HC_cluster$`HC_Tcell@meta.data$dim_40` <- gsub("CD4_effector_memory","1",HC_cluster$`HC_Tcell@meta.data$dim_40`)
HC_cluster$`HC_Tcell@meta.data$dim_40` <- gsub("CD8_CTLs_1","2",HC_cluster$`HC_Tcell@meta.data$dim_40`)
HC_cluster$`HC_Tcell@meta.data$dim_40` <- gsub("CD8_Naive","3",HC_cluster$`HC_Tcell@meta.data$dim_40`)
HC_cluster$`HC_Tcell@meta.data$dim_40` <- gsub("CD8_CTLs_2","4",HC_cluster$`HC_Tcell@meta.data$dim_40`)
HC_cluster$`HC_Tcell@meta.data$dim_40` <- gsub("CD4_CTLs_1","5",HC_cluster$`HC_Tcell@meta.data$dim_40`)
HC_cluster$`HC_Tcell@meta.data$dim_40` <- gsub("CD4_CTLs_2","6",HC_cluster$`HC_Tcell@meta.data$dim_40`)
HC_cluster$`HC_Tcell@meta.data$dim_40` <- gsub("Treg","7",HC_cluster$`HC_Tcell@meta.data$dim_40`)
HC_cluster$`HC_Tcell@meta.data$dim_40` <- gsub("Hemoglobin","8",HC_cluster$`HC_Tcell@meta.data$dim_40`)
HC_cluster$`HC_Tcell@meta.data$dim_40` <- gsub("TRAV13-2_CD8_Tcell","9",HC_cluster$`HC_Tcell@meta.data$dim_40`)
HC_cluster$`HC_Tcell@meta.data$dim_40` <- gsub("unknown","10",HC_cluster$`HC_Tcell@meta.data$dim_40`)
HC_cluster$`HC_Tcell@meta.data$dim_40` <- gsub("TRBV7-9_T_CD8_cell","11",HC_cluster$`HC_Tcell@meta.data$dim_40`)
HC_cluster$`HC_Tcell@meta.data$dim_40` <- gsub("unknown","12",HC_cluster$`HC_Tcell@meta.data$dim_40`)
HC_cluster$`HC_Tcell@meta.data$dim_40` <- gsub("TRBV10-2_CD8_Tcell","13",HC_cluster$`HC_Tcell@meta.data$dim_40`)
HC_cluster$`HC_Tcell@meta.data$dim_40` <- gsub("Bcell","14",HC_cluster$`HC_Tcell@meta.data$dim_40`)
HC_cluster$`HC_Tcell@meta.data$dim_40` <- gsub("unknown","15",HC_cluster$`HC_Tcell@meta.data$dim_40`)

HC_UMAP <- DimPlot(HC_Tcell, reduction="umap", label = TRUE)$data
HC_UMAP$ident <- gsub("CD4_Naive","0",HC_UMAP$ident)
HC_UMAP$ident <- gsub("CD4_effector_memory","1",HC_UMAP$ident)
HC_UMAP$ident <- gsub("CD8_CTLs_1","2",HC_UMAP$ident)
HC_UMAP$ident <- gsub("CD8_Naive","3",HC_UMAP$ident)
HC_UMAP$ident <- gsub("CD8_CTLs_2","4",HC_UMAP$ident)
HC_UMAP$ident <- gsub("CD4_CTLs_1","5",HC_UMAP$ident)
HC_UMAP$ident <- gsub("CD4_CTLs_2","6",HC_UMAP$ident)
HC_UMAP$ident <- gsub("Treg","7",HC_UMAP$ident)
HC_UMAP$ident <- gsub("Hemoglobin","8",HC_UMAP$ident)
HC_UMAP$ident <- gsub("TRAV13-2_CD8_Tcell","9",HC_UMAP$ident)
HC_UMAP$ident <- gsub("unknown","10",HC_UMAP$ident)
HC_UMAP$ident <- gsub("TRBV7-9_T_CD8_cell","11",HC_UMAP$ident)
HC_UMAP$ident <- gsub("unknown","12",HC_UMAP$ident)
HC_UMAP$ident <- gsub("TRBV10-2_CD8_Tcell","13",HC_UMAP$ident)
HC_UMAP$ident <- gsub("Bcell","14",HC_UMAP$ident)
HC_UMAP$ident <- gsub("unknown","15",HC_UMAP$ident)

write.table(colnames(HC_Tcell@assays$integrated),"cell_id.txt")
write.table(HC_Tcell@assays$integrated@var.features, quote = F, row.names = F, "gene_list.txt")
write.table(HC_cluster$`HC_Tcell@meta.data$dim_40`, "clustering2.txt")
write.table(HC_UMAP,"umap.txt")


pSS_Tcell
HC_Tcell

LabelPoints(plot = VariableFeaturePlot(pSS_Tcell), points = top10, repel = TRUE)
LabelPoints(plot = VariableFeaturePlot(HC_Tcell), points = top10_HC, repel = TRUE)


VlnPlot(Tcell, features = filter(markers,markers$cluster=="CD4_CTLs_2")$gene[1:10])

Tcell@active.ident <- as.factor(Tcell$control)
Tcell.markers_control <- FindAllMarkers(Tcell, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

write.table(filter(Tcell.markers_control,Tcell.markers_control$cluster=="pSS")$gene,"gene_list.txt",quote = F,row.names = F,col.names = F)
Tcell@active.ident <- as.factor(Tcell$dim_22)
DotPlot(Tcell, features = Tcell.markers_control$gene[1:10],cols = c("red","blue"),  dot.scale = 8) +
  RotatedAxis()

HC_Tcell_markers <- FindAllMarkers(HC_Tcell, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pSS_Tcell_markers <- FindAllMarkers(pSS_Tcell, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(filter(HC_Tcell_markers,HC_Tcell_markers$cluster=="CD4_CTLs_2")$gene, quote = F, row.names = F, "gene_list.txt")
write.table(filter(pSS_Tcell_markers,pSS_Tcell_markers$cluster=="CD4_CTLs_2")$gene, quote = F, row.names = F, "gene_list.txt")

FeaturePlot(Tcell, features = c("JUN"))

write.table(colnames(Tcell@assays$integrated),"cell_id.txt")
write.table(Tcell@assays$integrated@var.features, quote = F, row.names = F, "gene_list.txt")
write.table(as.data.frame(Tcell@meta.data$dim_40), "clustering2.txt")
write.table(DimPlot(Tcell, reduction="umap", label = TRUE)$data,"umap.txt")

