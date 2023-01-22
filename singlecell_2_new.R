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

library(stringr)
library(Seurat)
library(clusterProfiler)
library(patchwork)
library(dplyr)
library(org.Hs.eg.db)
library(tidyverse)

epimarker<-read.table("D:/GSE139325_T_S_epigenetic_enzyme_for_volcano_plot_new.txt",sep = "\t",header = TRUE,row.names = 1)
for (i in c("wilcox","bimod","roc","t","negbinom","poisson","LR","MAST","DESeq2")) {
  Tregmarker<-FindMarkers(HCC_count_PTR_TTR_Seurat,ident.1 = "TTR",ident.2 = "PTR",logfc.threshold = 0,min.pct = 0,test.use = i)
  genename<-rownames(Tregmarker)
  g<-bitr(geneID = genename,fromType = "ENTREZID",toType = "SYMBOL",OrgDb = org.Hs.eg.db,drop = TRUE)
  g<-g[!duplicated(g$SYMBOL),]
  Tregmarker<-Tregmarker[g$ENTREZID,]
  rownames(Tregmarker)<-g$SYMBOL
  Newname<-intersect(toupper(rownames(epimarker)),rownames(Tregmarker))
  Tregmarker_new<-Tregmarker[Newname,]
  write.table(Tregmarker_new,file = paste0("c:/Users/xjmik/Desktop/HCCmarker",i,".txt"),sep = "\t")
  remove(Tregmarker,Tregmarker_new,Newname,g,genename)
}
