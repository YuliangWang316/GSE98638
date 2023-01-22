HCC_TPM<-read.table("D:/Jmjd1c_Treg_Tumor/GSE98638/GSE98638_HCC.TCell.S5063.TPM.txt",sep = "\t",row.names = 1,header = TRUE)
HCC_count<-read.table("D:/Jmjd1c_Treg_Tumor/GSE98638/GSE98638_HCC.TCell.S5063.count.txt",sep = "\t",row.names = 1,header = TRUE)
HCC_normalize_centered<-read.table("D:/Jmjd1c_Treg_Tumor/GSE98638/GSE98638_HCC.TCell.S4070.norm.centered.txt",sep = "\t",row.names = 1,header = TRUE)

library(dplyr)
HCC_TPM_PTR<-select(HCC_TPM,starts_with("PTR"))
HCC_TPM_TTR<-select(HCC_TPM,starts_with("TTR"))
HCC_count_PTR<-select(HCC_count,starts_with("PTR"))
HCC_count_TTR<-select(HCC_count,starts_with("TTR"))
HCC_normalize_centered_PTR<-select(HCC_normalize_centered,starts_with("PTR"))
HCC_normalize_centered_TTR<-select(HCC_normalize_centered,starts_with("TTR"))
HCC_TPM_PTR_TTR<-select(HCC_TPM,starts_with("PTR"),starts_with("TTR"))
HCC_count_PTR_TTR<-select(HCC_count,starts_with("PTR"),starts_with("TTR"))
HCC_normalize_centered_PTR_TTR<-select(HCC_normalize_centered,starts_with("PTR"),starts_with("TTR"))

library(cowplot)
library(Seurat)
HCC_count_PTR_seurat <- CreateSeuratObject(counts = HCC_count_PTR, project = "IMMUNE_CTRL", min.cells = 5)
HCC_count_PTR_seurat$stim <- "CTRL"
HCC_count_PTR_seurat <- subset(HCC_count_PTR_seurat, subset = nFeature_RNA > 500)
HCC_count_PTR_seurat <- NormalizeData(HCC_count_PTR_seurat, verbose = FALSE)
HCC_count_PTR_seurat <- FindVariableFeatures(HCC_count_PTR_seurat, selection.method = "vst", nfeatures = 2000)

HCC_count_TTR_seurat <- CreateSeuratObject(counts = HCC_count_TTR, project = "IMMUNE_STIM", min.cells = 5)
HCC_count_TTR_seurat$stim <- "STIM"
HCC_count_TTR_seurat <- subset(HCC_count_TTR_seurat, subset = nFeature_RNA > 500)
HCC_count_TTR_seurat <- NormalizeData(HCC_count_TTR_seurat, verbose = FALSE)
HCC_count_TTR_seurat <- FindVariableFeatures(HCC_count_TTR_seurat, selection.method = "vst", nfeatures = 2000)

immune.anchors <- FindIntegrationAnchors(object.list = list(HCC_count_PTR_seurat, HCC_count_TTR_seurat), dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)

DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)

immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

DimPlot(immune.combined, reduction = "umap", split.by = "stim")
Cluster0.markers <- FindConservedMarkers(immune.combined, ident.1 = 0, grouping.var = "stim", verbose = FALSE)
Cluster1.markers <- FindConservedMarkers(immune.combined, ident.1 = 1, grouping.var = "stim", verbose = FALSE)
Cluster2.markers <- FindConservedMarkers(immune.combined, ident.1 = 2, grouping.var = "stim", verbose = FALSE)
write.table(Cluster0.markers,file = "c:/Users/Administrator/Desktop/GSE98638/Cluster0.conservemarkers.txt",sep = "\t",row.names = TRUE,col.names = TRUE)
write.table(Cluster1.markers,file = "c:/Users/Administrator/Desktop/GSE98638/Cluster1.conservemarkers.txt",sep = "\t",row.names = TRUE,col.names = TRUE)
write.table(Cluster2.markers,file = "c:/Users/Administrator/Desktop/GSE98638/Cluster2.conservemarkers.txt",sep = "\t",row.names = TRUE,col.names = TRUE)
immune.combined$celltype.stim <- paste(Idents(immune.combined), immune.combined$stim, sep = "_")
immune.combined$celltype <- Idents(immune.combined)
Idents(immune.combined) <- "celltype.stim"
Cluster0.response <- FindMarkers(immune.combined, ident.1 = "0_STIM", ident.2 = "0_CTRL", verbose = FALSE)
Cluster1.response <- FindMarkers(immune.combined, ident.1 = "1_STIM", ident.2 = "1_CTRL", verbose = FALSE)
Cluster2.response <- FindMarkers(immune.combined, ident.1 = "2_STIM", ident.2 = "2_CTRL", verbose = FALSE)
write.table(Cluster0.response,file = "c:/Users/Administrator/Desktop/GSE98638/Cluster0.response.txt",sep = "\t",row.names = TRUE,col.names = TRUE)
write.table(Cluster1.response,file = "c:/Users/Administrator/Desktop/GSE98638/Cluster1.response.txt",sep = "\t",row.names = TRUE,col.names = TRUE)
write.table(Cluster2.response,file = "c:/Users/Administrator/Desktop/GSE98638/Cluster2.response.txt",sep = "\t",row.names = TRUE,col.names = TRUE)

FeaturePlot(immune.combined, features = c("221037"), split.by = "stim", max.cutoff = 3, 
            label = TRUE)
VlnPlot(immune.combined, features = c("221037"), split.by = "stim", 
        pt.size = 0, combine = FALSE)

HCC_count_PTR_TTR_Seurat<-CreateSeuratObject(counts = HCC_count_PTR_TTR, project = "Treg", min.cells = 3, min.features = 200)
HCC_count_PTR_TTR_Seurat <- subset(HCC_count_PTR_TTR_Seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 6000)
HCC_count_PTR_TTR_Seurat <- NormalizeData(HCC_count_PTR_TTR_Seurat)
HCC_count_PTR_TTR_Seurat <- FindVariableFeatures(HCC_count_PTR_TTR_Seurat, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(HCC_count_PTR_TTR_Seurat)
HCC_count_PTR_TTR_Seurat <- ScaleData(HCC_count_PTR_TTR_Seurat, features = all.genes)
HCC_count_PTR_TTR_Seurat <- RunPCA(HCC_count_PTR_TTR_Seurat, features = VariableFeatures(object = HCC_count_PTR_TTR_Seurat))
ElbowPlot(HCC_count_PTR_TTR_Seurat)
HCC_count_PTR_TTR_Seurat <- FindNeighbors(HCC_count_PTR_TTR_Seurat, dims = 1:10)
HCC_count_PTR_TTR_Seurat <- FindClusters(HCC_count_PTR_TTR_Seurat, resolution = 0.5)

PTR<-colnames(HCC_count_PTR)
TTR<-colnames(HCC_count_TTR)

Idents(HCC_count_PTR_TTR_Seurat, cells = PTR)<-'PTR'
Idents(HCC_count_PTR_TTR_Seurat, cells = TTR)<-'TTR'

HCC_count_PTR_TTR_Seurat <- RunUMAP(HCC_count_PTR_TTR_Seurat, dims = 1:5)
HCC_count_PTR_TTR_Seurat <- RunTSNE(HCC_count_PTR_TTR_Seurat, dims = 1:10)

DimPlot(HCC_count_PTR_TTR_Seurat, reduction = "umap")
DimPlot(HCC_count_PTR_TTR_Seurat, reduction = "tsne")

HCC.markers <- FindAllMarkers(HCC_count_PTR_TTR_Seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(HCC.markers,file = "c:/Users/Administrator/Desktop/GSE98638/HCC.markers.txt",sep = "\t",row.names = TRUE,col.names = TRUE)

FeaturePlot(HCC_count_PTR_TTR_Seurat,features = "221037",label = TRUE)
VlnPlot(HCC_count_PTR_TTR_Seurat,features = "221037",pt.size = 0)

write.table(HCC_count_PTR_TTR,file = "c:/Users/Administrator/Desktop/GSE98638/HCC_count_PTR_TTR.txt",sep = "\t",row.names = TRUE,col.names = TRUE)
write.table(HCC_TPM_PTR_TTR,file = "c:/Users/Administrator/Desktop/GSE98638/HCC_TPM_PTR_TTR.txt",sep = "\t",row.names = TRUE,col.names = TRUE)
write.table(HCC_normalize_centered_PTR_TTR,file = "c:/Users/Administrator/Desktop/GSE98638/HCC_normlize_centered_PTR_TTR.txt",sep = "\t",row.names = TRUE,col.names = TRUE)

FeaturePlot(HCC_count_PTR_TTR_Seurat,features = "221037",cols = c("yellow","purple"),min.cutoff = "q40",max.cutoff = "q90",pt.size = 1.0)

HCC_normalize_centered_PTR_TTR_seurat<-CreateSeuratObject(counts = HCC_normalize_centered_PTR_TTR, project = "Treg", min.cells = 3, min.features = 200)
HCC_normalize_centered_PTR_TTR_seurat@assays[["RNA"]]@scale.data<-as.matrix(HCC_normalize_centered_PTR_TTR)
PTR<-colnames(HCC_normalize_centered_PTR)
TTR<-colnames(HCC_normalize_centered_TTR)

Idents(HCC_normalize_centered_PTR_TTR_seurat, cells = PTR)<-'PTR'
Idents(HCC_normalize_centered_PTR_TTR_seurat, cells = TTR)<-'TTR'
HCC_normalize_centered_PTR_TTR_seurat <- RunPCA(HCC_normalize_centered_PTR_TTR_seurat, features = rownames(HCC_normalize_centered_PTR_TTR))
HCC_normalize_centered_PTR_TTR_seurat <- RunUMAP(HCC_normalize_centered_PTR_TTR_seurat, dims = 1:10)
HCC_normalize_centered_PTR_TTR_seurat <- RunTSNE(HCC_normalize_centered_PTR_TTR_seurat, dims = 1:10)

DimPlot(HCC_normalize_centered_PTR_TTR_seurat, reduction = "umap")
DimPlot(HCC_normalize_centered_PTR_TTR_seurat, reduction = "tsne")

FeaturePlot(HCC_normalize_centered_PTR_TTR_seurat,features = "221037")
