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
metadata<-read.table("D:/Jmjd1c_Treg_Tumor/GSE98638/GSE98638_family.txt",sep = "\t",header = TRUE)
PTR_metadata<-filter(metadata,sampleType == "PTR")
TTR_metadata<-filter(metadata,sampleType == "TTR")
PTR_TTR_metadata<-rbind(PTR_metadata,TTR_metadata)
rownames(PTR_TTR_metadata)<-PTR_TTR_metadata[,1]

HCC_count_PTR_TTR_Seurat<-CreateSeuratObject(counts = HCC_count_PTR_TTR, project = "Treg", min.cells = 3, min.features = 200,meta.data = PTR_TTR_metadata)
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

B_TTR_seurat<-subset(HCC_count_PTR_TTR_Seurat,idents = "TTR")

Idents(B_TTR_seurat)<-B_TTR_seurat@meta.data$Patient

for (i in unique(B_TTR_seurat@meta.data[["Patient"]])) {
  assign(i,subset(B_TTR_seurat,idents = i))
  assign("a",slot(get(i),"assays"))
  assign("b",as.data.frame(t(as.data.frame(a$RNA@scale.data))))
  assign("c",select(b,one_of("221037","3458")))
  colnames(c)<-c("JMJD1C","IFNG")
  assign(paste("scaldata_",i,"_TTR",sep = ""),data.frame(sum(c$JMJD1C),sum(c$IFNG)))
}
scaldata<-rbind(scaldata_P0322_TTR,scaldata_P0407_TTR,scaldata_P0508_TTR,scaldata_P1116_TTR)
colnames(scaldata)<-c("JMJD1C","IFNG")
library(ggplot2)
library(ggpubr)
ggplot(data = scaldata,aes(x=JMJD1C,y=IFNG)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=TRUE,size=1.5,color="red")+stat_cor(data = scaldata,method = "spearman") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 

