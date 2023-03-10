HCC_TPM<-read.table("G:/GSE98638/GSE98638_HCC.TCell.S5063.TPM.txt",sep = "\t",row.names = 1,header = TRUE)
HCC_count<-read.table("G:/GSE98638/GSE98638_HCC.TCell.S5063.count.txt",sep = "\t",row.names = 1,header = TRUE)
HCC_normalize_centered<-read.table("G:/GSE98638/GSE98638_HCC.TCell.S4070.norm.centered.txt",sep = "\t",row.names = 1,header = TRUE)

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
HCC_count_PTR_seurat  <- CreateSeuratObject(counts =HCC_count_PTR, project = "HCC_count_PTR", min.cells = 3, min.features = 200)
VlnPlot(HCC_count_PTR_seurat, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)
HCC_count_PTR_seurat <- subset(HCC_count_PTR_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 )
HCC_count_PTR_seurat <- NormalizeData(HCC_count_PTR_seurat, verbose = FALSE)
HCC_count_PTR_seurat <- FindVariableFeatures(HCC_count_PTR_seurat, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(HCC_count_PTR_seurat)
HCC_count_PTR_seurat <- ScaleData(HCC_count_PTR_seurat, features = all.genes)
HCC_count_PTR_seurat <- RunPCA(HCC_count_PTR_seurat, features = VariableFeatures(object = HCC_count_PTR_seurat))
ElbowPlot(HCC_count_PTR_seurat)
HCC_count_PTR_seurat <- FindNeighbors(HCC_count_PTR_seurat, dims = 1:10)
HCC_count_PTR_seurat <- FindClusters(HCC_count_PTR_seurat, resolution = 0.5)
HCC_count_PTR_seurat <- RunUMAP(HCC_count_PTR_seurat, dims = 1:10)
DimPlot(HCC_count_PTR_seurat, reduction = "umap")

HCC_count_PTR_nomrlize_scaledata<-HCC_count_PTR_seurat@assays[["RNA"]]@scale.data
HCC_count_PTR_nomrlize_scaledata<-as.data.frame(HCC_count_PTR_nomrlize_scaledata)

HCC_count_PTR_nomrlize_scaledata_t<-t(HCC_count_PTR_nomrlize_scaledata)
HCC_count_PTR_nomrlize_scaledata_t<-as.data.frame(HCC_count_PTR_nomrlize_scaledata_t)

IL2_stat5<-read.table("g:/GSE98638/IL2_STAT5_entrezID.txt",sep = "\t",header = TRUE)
IL6_stat3<-read.table("g:/GSE98638/IL6_STAT3_entrezID.txt",sep = "\t",header = TRUE)

IL2_stat5<-IL2_stat5[2:200,]
IL6_stat3<-IL6_stat3[2:88,]

HCC_count_PTR_seurat <- CellCycleScoring(HCC_count_PTR_seurat, s.features = IL2_stat5, g2m.features = IL6_stat3, set.ident = FALSE)
IL2_stat5_list<-list(IL2_stat5)
IL6_stat3_list<-list(IL6_stat3)
HCC_count_PTR_seurat <- AddModuleScore(HCC_count_PTR_seurat,features = IL2_stat5_list,name = "IL2_stat5")
HCC_count_PTR_seurat <- AddModuleScore(HCC_count_PTR_seurat,features = IL6_stat3_list,name = "IL6_stat3")
JMJD1C<-select(HCC_count_PTR_nomrlize_scaledata_t,starts_with("221037"))
NRP1<-select(HCC_count_PTR_nomrlize_scaledata_t,starts_with("8829"))
Pdcd1<-select(HCC_count_PTR_nomrlize_scaledata_t,starts_with("5133") & ends_with("5133"))
Ifng<-select(HCC_count_PTR_nomrlize_scaledata_t,starts_with("3458"))
IL2_stat5_score<-as.data.frame(HCC_count_PTR_seurat$IL2_stat51)
IL6_stat3_score<-as.data.frame(HCC_count_PTR_seurat$IL6_stat31)
data<-cbind(IL2_stat5_score,IL6_stat3_score,JMJD1C,Pdcd1,Ifng)
colnames(data)<-c("IL2_stat5","IL6_stat3","JMJD1C","Pdcd1","Ifng")
#library(corrplot)
#corrplot(matrix, method = "color")  

library(ggplot2)
library(ggpubr)
ggplot(data = data,aes(x=IL2_stat5,y=JMJD1C)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=FALSE,size=1.5,color="red")+stat_cor(data = data,method = "pearson") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 
data_new<-filter(data,IL2_stat5>0 & JMJD1C>0)
ggplot(data = data_new,aes(x=IL2_stat5,y=JMJD1C)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=FALSE,size=1.5,color="red")+stat_cor(data = data_new,method = "pearson") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 
ggplot(data = data,aes(x=IL6_stat3,y=JMJD1C)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=FALSE,size=1.5,color="red")+stat_cor(data = data,method = "pearson") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 
data_new2<-filter(data,IL6_stat3>0 & JMJD1C>0)
ggplot(data = data_new2,aes(x=IL6_stat3,y=JMJD1C)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=FALSE,size=1.5,color="red")+stat_cor(data = data_new2,method = "pearson") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 
ggplot(data = data,aes(x=Pdcd1,y=JMJD1C)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=FALSE,size=1.5,color="red")+stat_cor(data = data,method = "pearson") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 
data_new3<-filter(data,Pdcd1>0 & JMJD1C>0)
ggplot(data = data_new3,aes(x=Pdcd1,y=JMJD1C)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=FALSE,size=1.5,color="red")+stat_cor(data = data_new3,method = "pearson") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 
ggplot(data = data,aes(x=Ifng,y=JMJD1C)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=FALSE,size=1.5,color="red")+stat_cor(data = data,method = "pearson") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 
data_new4<-filter(data,Ifng>0 & JMJD1C>0)
ggplot(data = data_new4,aes(x=Ifng,y=JMJD1C)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=FALSE,size=1.5,color="red")+stat_cor(data = data_new4,method = "pearson") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 

HCC_count_TTR_seurat  <- CreateSeuratObject(counts =HCC_count_TTR, project = "HCC_count_TTR", min.cells = 3, min.features = 200)
VlnPlot(HCC_count_TTR_seurat, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)
HCC_count_TTR_seurat <- subset(HCC_count_TTR_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 )
HCC_count_TTR_seurat <- NormalizeData(HCC_count_TTR_seurat, verbose = FALSE)
HCC_count_TTR_seurat <- FindVariableFeatures(HCC_count_TTR_seurat, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(HCC_count_TTR_seurat)
HCC_count_TTR_seurat <- ScaleData(HCC_count_TTR_seurat, features = all.genes)
HCC_count_TTR_seurat <- RunPCA(HCC_count_TTR_seurat, features = VariableFeatures(object = HCC_count_TTR_seurat))
ElbowPlot(HCC_count_TTR_seurat)
HCC_count_TTR_seurat <- FindNeighbors(HCC_count_TTR_seurat, dims = 1:10)
HCC_count_TTR_seurat <- FindClusters(HCC_count_TTR_seurat, resolution = 0.5)
HCC_count_TTR_seurat <- RunUMAP(HCC_count_TTR_seurat, dims = 1:10)
DimPlot(HCC_count_TTR_seurat, reduction = "umap")

HCC_count_TTR_nomrlize_scaledata<-HCC_count_TTR_seurat@assays[["RNA"]]@scale.data
HCC_count_TTR_nomrlize_scaledata<-as.data.frame(HCC_count_TTR_nomrlize_scaledata)

HCC_count_TTR_nomrlize_scaledata_t<-t(HCC_count_TTR_nomrlize_scaledata)
HCC_count_TTR_nomrlize_scaledata_t<-as.data.frame(HCC_count_TTR_nomrlize_scaledata_t)

IL2_stat5<-read.table("g:/GSE98638/IL2_STAT5_entrezID.txt",sep = "\t",header = TRUE)
IL6_stat3<-read.table("g:/GSE98638/IL6_STAT3_entrezID.txt",sep = "\t",header = TRUE)

IL2_stat5<-IL2_stat5[2:200,]
IL6_stat3<-IL6_stat3[2:88,]

HCC_count_TTR_seurat <- CellCycleScoring(HCC_count_TTR_seurat, s.features = IL2_stat5, g2m.features = IL6_stat3, set.ident = FALSE)
IL2_stat5_list<-list(IL2_stat5)
IL6_stat3_list<-list(IL6_stat3)
HCC_count_TTR_seurat <- AddModuleScore(HCC_count_TTR_seurat,features = IL2_stat5_list,name = "IL2_stat5")
HCC_count_TTR_seurat <- AddModuleScore(HCC_count_TTR_seurat,features = IL6_stat3_list,name = "IL6_stat3")
JMJD1C<-select(HCC_count_TTR_nomrlize_scaledata_t,starts_with("221037"))
NRP1<-select(HCC_count_TTR_nomrlize_scaledata_t,starts_with("8829"))
Pdcd1<-select(HCC_count_TTR_nomrlize_scaledata_t,starts_with("5133") & ends_with("5133"))
Ifng<-select(HCC_count_TTR_nomrlize_scaledata_t,starts_with("3458") & ends_with("3458"))
IL2_stat5_score<-as.data.frame(HCC_count_TTR_seurat$IL2_stat51)
IL6_stat3_score<-as.data.frame(HCC_count_TTR_seurat$IL6_stat31)
data<-cbind(IL2_stat5_score,IL6_stat3_score,JMJD1C,Pdcd1,Ifng)
colnames(data)<-c("IL2_stat5","IL6_stat3","JMJD1C","Pdcd1","Ifng")
#library(corrplot)
#corrplot(matrix, method = "color")  

library(ggplot2)
library(ggpubr)
ggplot(data = data,aes(x=IL2_stat5,y=JMJD1C)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=FALSE,size=1.5,color="red")+stat_cor(data = data,method = "pearson") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 
data_new<-filter(data,IL2_stat5>0 & JMJD1C>0)
ggplot(data = data_new,aes(x=IL2_stat5,y=JMJD1C)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=FALSE,size=1.5,color="red")+stat_cor(data = data_new,method = "pearson") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 
ggplot(data = data,aes(x=IL6_stat3,y=JMJD1C)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=FALSE,size=1.5,color="red")+stat_cor(data = data,method = "pearson") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 
data_new2<-filter(data,IL6_stat3>0 & JMJD1C>0)
ggplot(data = data_new2,aes(x=IL6_stat3,y=JMJD1C)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=FALSE,size=1.5,color="red")+stat_cor(data = data_new2,method = "pearson") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 
ggplot(data = data,aes(x=Pdcd1,y=JMJD1C)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=FALSE,size=1.5,color="red")+stat_cor(data = data,method = "pearson") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 
data_new3<-filter(data,Pdcd1>0 & JMJD1C>0)
ggplot(data = data_new3,aes(x=Pdcd1,y=JMJD1C)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=FALSE,size=1.5,color="red")+stat_cor(data = data_new3,method = "pearson") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 
ggplot(data = data,aes(x=Ifng,y=JMJD1C)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=FALSE,size=1.5,color="red")+stat_cor(data = data,method = "pearson") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 
data_new4<-filter(data,Ifng>0 & JMJD1C>0)
ggplot(data = data_new4,aes(x=Ifng,y=JMJD1C)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=FALSE,size=1.5,color="red")+stat_cor(data = data_new4,method = "pearson") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 
