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

IL2_stat5<-read.table("g:/GSE98638/IL2_STAT5_entrezID.txt",sep = "\t",header = TRUE)
IL6_stat3<-read.table("g:/GSE98638/IL6_STAT3_entrezID.txt",sep = "\t",header = TRUE)

IL2_stat5<-IL2_stat5[2:200,]
IL6_stat3<-IL6_stat3[2:88,]

HCC_count_PTR_TTR_Seurat <- CellCycleScoring(HCC_count_PTR_TTR_Seurat, s.features = IL2_stat5, g2m.features = IL6_stat3, set.ident = FALSE)
IL2_stat5_list<-list(IL2_stat5)
IL6_stat3_list<-list(IL6_stat3)
HCC_count_PTR_TTR_Seurat <- AddModuleScore(HCC_count_PTR_TTR_Seurat,features = IL2_stat5_list,name = "IL2_stat5")
HCC_count_PTR_TTR_Seurat <- AddModuleScore(HCC_count_PTR_TTR_Seurat,features = IL6_stat3_list,name = "IL6_stat3")

VlnPlot(HCC_count_PTR_TTR_Seurat,features = c("IL2_stat51","IL6_stat31","221037","8829","5133","3458"),pt.size = 0)
FeaturePlot(HCC_count_PTR_TTR_Seurat,features = c("IL2_stat51","IL6_stat31","221037","8829","5133","3458"))

HCC_count_PTR_TTR_markers<-FindMarkers(HCC_count_PTR_TTR_Seurat,ident.1 = "TTR",ident.2 = "PTR",logfc.threshold = 0,min.pct = 0,test.use = "DESeq2")
write.table(HCC_count_PTR_TTR_markers,file = "c:/Users/xjmik/Desktop/HCCDESeq2markers.txt",sep = "\t")
avg.HCC_count_PTR_TTR<-AverageExpression(HCC_count_PTR_TTR_Seurat)
write.table(avg.HCC_count_PTR_TTR,file = "G:/GSE98638/avg.HCC_count_PTR_TTR.txt",sep = "\t")
GSEA_TTR_PTR<-as.data.frame(HCC_count_PTR_TTR_Seurat@assays[["RNA"]]@scale.data)
write.table(GSEA_TTR_PTR,file = "G:/GSE98638/GSEA_TTR_PTR.txt",sep = "\t")


VlnPlot(HCC_count_PTR_TTR_Seurat,features = c("221037","5133"),pt.size = 0)
