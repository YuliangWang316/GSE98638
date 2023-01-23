library(Seurat)
library(cowplot)
library(dplyr)
library(patchwork)

pbmc<-readRDS("c:/Users/xjmik/Desktop/GSE98638/pbmc.rds")
pbmc<-subset(pbmc,idents = "TR")
Idents(pbmc)<-pbmc@meta.data$Bulksample
pbmc<-subset(pbmc,idents = "T")
pbmc_JMJD1C<-FetchData(pbmc,vars = "JMJD1C")
pbmc_JMJD1C$barcode<-rownames(pbmc_JMJD1C)
pbmc_JMJD1C_zero<-rownames(pbmc_JMJD1C[which(pbmc_JMJD1C$JMJD1C == 0),])
pbmc_JMJD1C_nozero<-pbmc_JMJD1C[which(pbmc_JMJD1C$JMJD1C > 0),]
pbmc_JMJD1C_hi<-rownames(pbmc_JMJD1C_nozero[which(pbmc_JMJD1C_nozero$JMJD1C > median(pbmc_JMJD1C_nozero$JMJD1C)),])
pbmc_JMJD1C_lo<-rownames(pbmc_JMJD1C_nozero[which(pbmc_JMJD1C_nozero$JMJD1C <= median(pbmc_JMJD1C_nozero$JMJD1C)),])
remove(pbmc_JMJD1C,pbmc_JMJD1C_nozero)
pbmc_hilo<-subset(pbmc,cells = c(pbmc_JMJD1C_hi,pbmc_JMJD1C_lo))
pbmc_hilo@meta.data$sample_new<-""
for (i in 1:length(rownames(pbmc_hilo@meta.data))) {
  for (j in 1:length(pbmc_JMJD1C_hi)) {
    if(rownames(pbmc_hilo@meta.data)[i] == pbmc_JMJD1C_hi[j]){
      pbmc_hilo@meta.data$sample_new[i] <- "hi"      
    }
  }
  remove(j)
  for (j in 1:length(pbmc_JMJD1C_lo)) {
    if(rownames(pbmc_hilo@meta.data)[i] == pbmc_JMJD1C_lo[j]){
      pbmc_hilo@meta.data$sample_new[i] <- "lo"      
    }
  }
}
remove(i,j,pbmc_JMJD1C_hi,pbmc_JMJD1C_lo)
pbmc_zero<-subset(pbmc,cells = pbmc_JMJD1C_zero)
remove(pbmc_JMJD1C_zero)
pbmc_hilo<-RunLDA(pbmc_hilo,labels = pbmc_hilo$sample_new,reduction.name = "lda_hilo")
pbmc_hilo_new<-pbmc_hilo
a<-pbmc_hilo_new@reductions$lda_hilo@cell.embeddings
a<-as.data.frame(a)
a$ldahilo_2<-a$ldahilo_1
a<-as.matrix(a)
b<-pbmc_hilo_new@reductions$lda_hilo@feature.loadings
b<-as.data.frame(b)
b$LDA_2<-b$LDA_1
b<-as.matrix(b)
c<-pbmc_hilo_new@reductions$lda_hilo@feature.loadings.projected
c<-as.data.frame(c)
c$ldahilo_2<-c$ldahilo_1
c<-as.matrix(c)
pbmc_hilo_new@reductions$lda_hilo@cell.embeddings<-a
pbmc_hilo_new@reductions$lda_hilo@feature.loadings<-b
pbmc_hilo_new@reductions$lda_hilo@feature.loadings.projected<-c
remove(a,b,c)
Idents()
immune.anchors <- FindIntegrationAnchors(object.list = list(pbmc_hilo, pbmc_zero),reference = pbmc_hilo)

