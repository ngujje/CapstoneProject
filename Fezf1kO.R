
install.packages("Matrix")
canceinstall.packages('remotes')
install.packages("dendextend")

remotes::install_github(repo = 'satijalab/seurat', ref = 'develop')
library(Seurat)
library(dendextend)

install.packages("remotes") 
remotes::install_github("YuLab-SMU/ggtree")

setwd("/Users/namrathagujje/Downloads/NamrathaProject/Fezf1KO-scRNA-seq")
library(Seurat)
library(ggplot2)
library(cowplot)
library (aplot)
library(ggdendro)
library(glmGamPoi)
library(harmony)
library(stringr)
library(xgboost)
library(reshape)
library(tidyr)
library(Matrix)
library(remotes)
library(ggtree)

#Combine the two replicates to generate one large dgcMatrix
L71_Fezf1KO_S6=Read10X(data.dir = "/Users/namrathagujje/Downloads/NamrathaProject/Fezf1KO-scRNA-seq/L71_Fezf1KO_S6/filtered_feature_bc_matrix")

colnames(L71_Fezf1KO_S6) <- paste0("Fezf1KO_S6", colnames(L71_Fezf1KO_S6))
save(L71_Fezf1KO_S6, file="/Users/namrathagujje/Downloads/NamrathaProject/Fezf1KO-scRNA-seq/Fezf1KO_S6.R")


L70_WT_S5=Read10X(data.dir = "/Users/namrathagujje/Downloads/NamrathaProject/Fezf1KO-scRNA-seq/L70_WT_S5/outs/filtered_feature_bc_matrix")



#Might combine the two Rdata
Count.mat_Fezf1WTandKO= cbind(L71_Fezf1KO_S6, L70_WT_S5)
save(Count.mat_Fezf1WTandKO, file="/Users/namrathagujje/Downloads/NamrathaProject/Fezf1KO-scRNA-seq/Fezf1KO6andWT5.R")

#Set the output and figure directory
projectName = "Fezf1WT"
min.feature = 200
out_path = paste0("/Users/namrathagujje/Downloads/NamrathaProject/Fezf1KO-scRNA-seq", projectName, "_", min.feature, "_minfeature_", Sys.Date())
dir.create(out_path)
fig.dir = paste0(out_path,"/Figs_preprocess/")
dir.create(fig.dir)

#Create the Seurate object
KO<-CreateSeuratObject(counts=L71_Fezf1KO_S6, project ="Fezf1KO_S6", min.cells =10, min.feature =200)
head(KO@meta.data)

#Add percent.mt, percent.rp, and age to the metadata
KO$age <- c("8wks")
KO$Cre <- c("Fezf1")
KO$Genotype <- c("KO")
KO$percent.mt <- PercentageFeatureSet(object=KO, pattern ="^mt-")
KO$percent.rp <- PercentageFeatureSet(object=KO, pattern ="^Rp[sl]")

head(KO@meta.data)

#Check the distribution of nFeature_RNA, nCount_RNA, percent.mt, percent.rp to set up the threshold for filtering data.
pdf(paste0(fig.dir,"SampleMetricsBeforefilter.pdf"),w=12,h=12)
VlnPlot(KO, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rp"), ncol = 2)
dev.off()

KO <- subset(KO, subset = nFeature_RNA < 7500 & nCount_RNA < 30000 & percent.mt < 20 & percent.rp < 10)

pdf(paste0(fig.dir,"SampleMetricsAfterfilter.pdf"),w=12,h=12)
VlnPlot(KO, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rp"), ncol = 2)
dev.off()

#Normalization data with LogNormalize and a scale.factor 10000
KO <- NormalizeData(KO, normalization.method = "LogNormalize", scale.factor = 10000)
#Check normalized data using GetAssayData(WT,slot = "data") or head(WT@assays$RNA@data)

#Find highly variable features (Genes)
KO <- FindVariableFeatures(KO,selection.method = "vst")
#Check var.features using head(WT@assays$RNA@var.features)

#Scale data for PCA
KO <- ScaleData(KO, vars.to.regress = c("nCount_RNA","nFeature_RNA", "percent.mt","percent.rp"))
#Check Scaled data using GetAssayData(WT,slot = "scale.data") or head(WT@assays$RNA@scale.data)
KO <- RunPCA(KO, npcs = 100)

#Use elbow plot to check the PCs
pdf(paste0(fig.dir,"ElbowPlot100PCs.pdf"),w=12,h=12)
ElbowPlot(object = KO,ndims = 100)
dev.off()

#Cluster cells
KO <- FindNeighbors(KO, dim = 1:75)
KO <- FindClusters(KO, resolution = 1)
KO <- RunUMAP(KO, dims = 1:75)

pdf(paste0(fig.dir,"UmapKO.pdf"),w=12,h=12)
DimPlot(KO, reduction = "umap", label = T)
dev.off()

#Generate dotplot
genelist = c("Olig2","Bcas1","Apoe","Glul","Rlbp1","Gfap","Crip1","Id3","Actd2","C1qa","Kcnj8","Pecam1","Rcvrn","Ptprc","Cd53","Laptm5","Rho","Gngt1","Arr3","Gngt2","Opn1mw","Opn1sw","Vsx2","Otx2","Sox2","Isl1","Grm6","Scgn","Tacr3","Syt2","Neto1","Irx6","Prkar2b","Grik1","Kcng4","Cabp5","Vsx1","Prkca","Lhx1","Onecut1","Onecut2","Rbpms","Slc17a6","Thy1","Pou4f1","Opn4","Spp1", "Pax6","Tfap2a", "Tfap2b","Tfap2c", "Gad1", "Gad2","Slc32a1","Slc6a9","Slc17a8","Gjd2")
pdf(paste0(fig.dir,"DotPlotAllMarkers_KO.pdf"),w=12,h=12)
DotPlot(KO, features=genelist) + RotatedAxis() + coord_flip()
dev.off()

#Save the Seurat object
saveRDS(KO, file = "/Users/namrathagujje/Downloads/NamrathaProject/Fezf1KO-scRNA-seq/Fezf1KO6andWT5.R")

KO <- readRDS(file = "/Users/namrathagujje/Downloads/NamrathaProject/Fezf1KO-scRNA-seq/Fezf1KO6andWT5.R")

#Cleanup WT object
KO$percent.rbpms <- PercentageFeatureSet(object=KO, pattern ="Rbpms")
KO$percent.pou4f1 <- PercentageFeatureSet(object=KO, pattern ="Pou4f1")
KO$percent.slc17a6 <- PercentageFeatureSet(object=KO, pattern ="Slc17a6")

#remove cluster 1, 9, 13, 14, 21, 24, 35, 42, 43, 45, 46 from WT object, from 10,434 cells to 8848 cells
fig.dir2 <- c("/Users/namrathagujje/Downloads/NamrathaProject/Fezf1KO-scRNA-seq/Fezf1KO6andWT5.R")

table(KO@active.ident)
RGCcluster = c (0, 2:8, 10:12, 15:20, 22:23, 25:34, 36:41, 44)
KOclean1 <- subset(KO, idents = RGCcluster)
table(KOclean1@active.ident)

KOclean1 <- NormalizeData(KOclean1, normalization.method = "LogNormalize", scale.factor = 10000)
KOclean1 <- FindVariableFeatures(KOclean1,selection.method = "vst")
KOclean1 <- ScaleData(KOclean1, vars.to.regress = c("nCount_RNA","nFeature_RNA", "percent.mt","percent.rp"))
KOclean1 <- RunPCA(KOclean1, npcs = 100)

dev.off()
pdf(paste0(fig.dir2,"ElbowPlot100PCs_KOclean1.pdf"),w=12,h=12)
ElbowPlot(object = KOclean1,ndims = 100)
dev.off()

#Cluster cells
KOclean1 <- FindNeighbors(KOclean1, dim = 1:75)
KOclean1 <- FindClusters(KOclean1, resolution = 1)
KOclean1 <- RunUMAP(KOclean1, dims = 1:75)

pdf(paste0(fig.dir2,"UmapKOclean1.pdf"),w=12,h=12)
DimPlot(KOclean1, reduction = "umap", label = T)
dev.off()

#Generate dotplot
genelist = c("Olig2","Bcas1","Apoe","Glul","Rlbp1","Gfap","Crip1","Id3","Actd2","C1qa","Kcnj8","Pecam1","Rcvrn","Ptprc","Cd53","Laptm5","Rho","Gngt1","Arr3","Gngt2","Opn1mw","Opn1sw","Vsx2","Otx2","Sox2","Isl1","Grm6","Scgn","Tacr3","Syt2","Neto1","Irx6","Prkar2b","Grik1","Kcng4","Cabp5","Vsx1","Prkca","Lhx1","Onecut1","Onecut2","Rbpms","Slc17a6","Thy1","Pou4f1","Opn4","Spp1", "Pax6","Tfap2a", "Tfap2b","Tfap2c", "Gad1", "Gad2","Slc32a1","Slc6a9","Slc17a8","Gjd2")
pdf(paste0(fig.dir2,"DotPlotAllMarkers_KOclean1.pdf"),w=12,h=12)
DotPlot(KOclean1, features=genelist) + RotatedAxis() + coord_flip()
dev.off()

table(KOclean1@active.ident)
#remove cluster #8 to 8767 cells
RGCcluster2 = c (0:7, 9:38)
KOclean2 <- subset(KOclean1, idents = RGCcluster2)
dim(KOclean2)

KOclean2 <- NormalizeData(KOclean2, normalization.method = "LogNormalize", scale.factor = 10000)
KOclean2 <- FindVariableFeatures(KOclean2,selection.method = "vst")
KOclean2 <- ScaleData(KOclean2, vars.to.regress = c("nCount_RNA","nFeature_RNA", "percent.mt","percent.rp"))
KOclean2 <- RunPCA(KOclean2, npcs = 100)

dev.off()
pdf(paste0(fig.dir2,"ElbowPlot100PCs_KOclean2.pdf"),w=12,h=12)
ElbowPlot(object = KOclean2,ndims = 100)
dev.off()

#Cluster cells
KOclean2 <- FindNeighbors(KOclean2, dim = 1:75)
KOclean2 <- FindClusters(KOclean2, resolution = 1)
KOclean2 <- RunUMAP(KOclean2, dims = 1:75)

pdf(paste0(fig.dir2,"UmapKOclean2.pdf"),w=12,h=12)
DimPlot(KOclean2, reduction = "umap", label = T)
dev.off()

#Generate dotplot
genelist = c("Olig2","Bcas1","Apoe","Glul","Rlbp1","Gfap","Crip1","Id3","Actd2","C1qa","Kcnj8","Pecam1","Rcvrn","Ptprc","Cd53","Laptm5","Rho","Gngt1","Arr3","Gngt2","Opn1mw","Opn1sw","Vsx2","Otx2","Sox2","Isl1","Grm6","Scgn","Tacr3","Syt2","Neto1","Irx6","Prkar2b","Grik1","Kcng4","Cabp5","Vsx1","Prkca","Lhx1","Onecut1","Onecut2","Rbpms","Slc17a6","Thy1","Pou4f1","Opn4","Spp1", "Pax6","Tfap2a", "Tfap2b","Tfap2c", "Gad1", "Gad2","Slc32a1","Slc6a9","Slc17a8","Gjd2")
pdf(paste0(fig.dir2,"DotPlotAllMarkers_KOclean2.pdf"),w=12,h=12)
DotPlot(KOclean2, features=genelist) + RotatedAxis() + coord_flip()
dev.off()

saveRDS(KOclean1, file = "/Users/namrathagujje/Downloads/NamrathaProject/Fezf1KO-scRNA-seq/Fezf1KO_S6.RUmapKOclean1.pdf")
saveRDS(KOclean2, file = "/Users/namrathagujje/Downloads/NamrathaProject/Fezf1KO-scRNA-seq/Fezf1KO_S6.RUmapKOclean2.pdf")

table(KOclean2@active.ident)
#remove cluster #18 to 8767 cells
RGCcluster4 = c (0:17, 19:37)
KOclean4 <- subset(KOclean2, idents = RGCcluster4)
dim(KOclean4)

class(KOclean4)

KOclean4 <- NormalizeData(KOclean4, normalization.method = "LogNormalize", scale.factor = 10000)
KOclean4 <- FindVariableFeatures(KOclean4, selection.method = "vst")
KOclean4 <- ScaleData(KOclean4, vars.to.regress = c("nCount_RNA","nFeature_RNA", "percent.mt","percent.rp"))
KOclean4 <- RunPCA(KOclean4, npcs = 100)

dev.off()
pdf(paste0(fig.dir2,"ElbowPlot100PCs_KOclean4.pdf"),w=12,h=12)
ElbowPlot(object = KOclean4,ndims = 100)
dev.off()

#Cluster cells
KOclean4 <- FindNeighbors(KOclean4, dim = 1:75)
KOclean4 <- FindClusters(KOclean4, resolution = 1.25)
KOclean4 <- RunUMAP(KOclean4, dims = 1:75)

pdf(paste0(fig.dir2,"UmapKOclean4.pdf"),w=12,h=12)
DimPlot(KOclean4, reduction = "umap", label = T)
dev.off()

#Generate dotplot
genelist = c("Olig2","Bcas1","Apoe","Glul","Rlbp1","Gfap","Crip1","Id3","Actd2","C1qa","Kcnj8","Pecam1","Rcvrn","Ptprc","Cd53","Laptm5","Rho","Gngt1","Arr3","Gngt2","Opn1mw","Opn1sw","Vsx2","Otx2","Sox2","Isl1","Grm6","Scgn","Tacr3","Syt2","Neto1","Irx6","Prkar2b","Grik1","Kcng4","Cabp5","Vsx1","Prkca","Lhx1","Onecut1","Onecut2","Rbpms","Slc17a6","Thy1","Pou4f1","Opn4","Spp1", "Pax6","Tfap2a", "Tfap2b","Tfap2c", "Gad1", "Gad2","Slc32a1","Slc6a9","Slc17a8","Gjd2")
pdf(paste0(fig.dir2,"DotPlotAllMarkers_KOclean4.pdf"),w=12,h=12)
DotPlot(KOclean4, features=genelist) + RotatedAxis() + coord_flip()
dev.off()

table(KOclean4@active.ident)
#subset cells with percent.mt < 15 & percent.rp < 5 to 8379 cells
KOclean3 <- subset(KOclean4, subset = percent.mt < 15 & percent.rp < 5)
table(KOclean3@active.ident)

KOclean3 <- NormalizeData(KOclean3, normalization.method = "LogNormalize", scale.factor = 10000)
KOclean3 <- FindVariableFeatures(KOclean3,selection.method = "vst")
KOclean3 <- ScaleData(KOclean3, vars.to.regress = c("nCount_RNA","nFeature_RNA", "percent.mt","percent.rp"))
KOclean3 <- RunPCA(KOclean3, npcs = 100)

dev.off()
pdf(paste0(fig.dir2,"ElbowPlot100PCs_KOclean3.pdf"),w=12,h=12)
ElbowPlot(object = KOclean3,ndims = 100)
dev.off()

#Cluster cells
KOclean3 <- FindNeighbors(KOclean3, dim = 1:75)
KOclean3 <- FindClusters(KOclean3, resolution = 1)
KOclean3 <- RunUMAP(KOclean3, dims = 1:75)

pdf(paste0(fig.dir2,"UmapKOclean3.pdf"),w=12,h=12)
DimPlot(KOclean3, reduction = "umap", label = T)
dev.off()

#Generate dotplot
genelist = c("Olig2","Bcas1","Apoe","Glul","Rlbp1","Gfap","Crip1","Id3","Actd2","C1qa","Kcnj8","Pecam1","Rcvrn","Ptprc","Cd53","Laptm5","Rho","Gngt1","Arr3","Gngt2","Opn1mw","Opn1sw","Vsx2","Otx2","Sox2","Isl1","Grm6","Scgn","Tacr3","Syt2","Neto1","Irx6","Prkar2b","Grik1","Kcng4","Cabp5","Vsx1","Prkca","Lhx1","Onecut1","Onecut2","Rbpms","Slc17a6","Thy1","Pou4f1","Opn4","Spp1", "Pax6","Tfap2a", "Tfap2b","Tfap2c", "Gad1", "Gad2","Slc32a1","Slc6a9","Slc17a8","Gjd2")
pdf(paste0(fig.dir2,"DotPlotAllMarkers_KOclean3.pdf"),w=12,h=12)
DotPlot(KOclean3, features=genelist) + RotatedAxis() + coord_flip()
dev.off()

saveRDS(KOclean1, file = "/Users/namrathagujje/Downloads/NamrathaProject/Fezf1KO-scRNA-seq/Fezf1KO_S6.RUmapKOclean1.pdf")
saveRDS(KOclean2, file = "/Users/namrathagujje/Downloads/NamrathaProject/Fezf1KO-scRNA-seq/Fezf1KO_S6.RUmapKOclean2.pdf")


saveRDS(KOclean3, file = "/Users/namrathagujje/Downloads/NamrathaProject/Fezf1KO-scRNA-seq/Fezf1KO_S6.RUmapKOclean3.pdf")
KOclean3$Orig.cluster <- KOclean3@meta.data$seurat_clusters
KO$Orig.cluster <- KO@meta.data$seurat_clusters
saveRDS(KO, file = "/Users/namrathagujje/Downloads/NamrathaProject/Fezf1KO-scRNA-seq/Fezf1KO_S6.RUmapKOclean2.pdf")

#Fezf1KOMarkers=FindAllMarkers(KOclean3, test.use = "MAST", max.cells.per.ident = 2000, min.pct = 0.1, logfc.threshold = 0.5, min.cells.feature = 10, return.thresh = 0.001, only.pos = T)
write.csv(Fezf1KOMarkers, file="./Mouse/Fezf1koMapping/SeuratObject/Fezf1KO/Fezf1WTClean3AllMarkers.csv")


#Creating the Dendrogram


library(Seurat)

# Assuming 'your_data' is your gene expression matrix
KOclean3 <- subset(KOclean4, subset = percent.mt < 15 & percent.rp < 5)
KOclean3 <- NormalizeData(KOclean3)
KOclean3 <- ScaleData(KOclean3)
KOclean3 <- RunPCA(KOclean3)
KOclean3 <- FindNeighbors(KOclean3)
KOclean3 <- FindClusters(KOclean3)


DotPlot(KOclean3, features = c("Olig2","Bcas1","Apoe","Glul","Rlbp1","Gfap","Crip1","Id3","Actd2","C1qa","Kcnj8","Pecam1","Rcvrn","Ptprc","Cd53","Laptm5","Rho","Gngt1","Arr3","Gngt2","Opn1mw","Opn1sw","Vsx2","Otx2","Sox2","Isl1","Grm6","Scgn","Tacr3","Syt2","Neto1","Irx6","Prkar2b","Grik1","Kcng4","Cabp5","Vsx1","Prkca","Lhx1","Onecut1","Onecut2","Rbpms","Slc17a6","Thy1","Pou4f1","Opn4","Spp1", "Pax6","Tfap2a", "Tfap2b","Tfap2c", "Gad1", "Gad2","Slc32a1","Slc6a9","Slc17a8","Gjd2"))
dev.off()

# Visualize dot plot with annotation bar
DotPlot(KOclean3, features=genelist) + RotatedAxis() + coord_flip()
dev.off()

#Analyze the WT data
mat_Fezf1WTL70 = Read10X(data.dir = "/Users/namrathagujje/Downloads/NamrathaProject/Fezf1KO-scRNA-seq/L70_WT_S5/outs/filtered_feature_bc_matrix")
colnames(mat_Fezf1WTL70) = paste0("mat_Fezf1WTL70_", colnames(mat_Fezf1WTL70))
save(mat_Fezf1WTL70, file = "/Users/namrathagujje/Downloads/NamrathaProject/Fezf1KO-scRNA-seq/L70WT.R")

load("/Users/namrathagujje/Downloads/NamrathaProject/Fezf1KO-scRNA-seq/L70WT.R")

#Set the output and figure directory
projectName = "Fezf1WT"
min.feature = 200
out_path = paste0("/Users/namrathagujje/Downloads/NamrathaProject/Fezf1KO-scRNA-seq/L70WT.R", projectName, "_", min.feature, "_minfeature_", Sys.Date())
dir.create(out_path)
fig.dir = paste0(out_path,"/Figs_preprocess/")
dir.create(fig.dir)

#Create the Seurate object
WT<-CreateSeuratObject(counts=mat_Fezf1WTL70, project ="Fezf1WT70", min.cells =10, min.feature =200)
head(WT@meta.data)

#Add percent.mt, percent.rp, and age to the metadata
WT$age <- c("8wks")
WT$Cre <- c("Fezf1")
WT$Genotype <- c("WT")
WT$percent.mt <- PercentageFeatureSet(object=WT, pattern ="^mt-")
WT$percent.rp <- PercentageFeatureSet(object=WT, pattern ="^Rp[sl]")

head(WT@meta.data)

#Check the distribution of nFeature_RNA, nCount_RNA, percent.mt, percent.rp to set up the threshold for filtering data.
pdf(paste0(fig.dir,"SampleMetricsBeforefilter.pdf"),w=12,h=12)
VlnPlot(WT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"), ncol = 2)
dev.off()

WT <- subset(WT, subset = nFeature_RNA < 10000 & nCount_RNA < 45000 & percent.mt < 40 & percent.rp < 10)

pdf(paste0(fig.dir,"SampleMetricsAfterfilter.pdf"),w=12,h=12)
VlnPlot(WT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"), ncol = 2)
dev.off()

#Normalization data with LogNormalize and a scale.factor 10000
WT <- NormalizeData(WT, normalization.method = "LogNormalize", scale.factor = 10000)
#Check normalized data using GetAssayData(WT,slot = "data") or head(WT@assays$RNA@data)

#Find highly variable features (Genes)
WT <- FindVariableFeatures(WT,selection.method = "vst")
#Check var.features using head(WT@assays$RNA@var.features)

#Scale data for PCA
WT <- ScaleData(WT, vars.to.regress = c("nCount_RNA","nFeature_RNA", "percent.mt","percent.rp"))
#Check Scaled data using GetAssayData(WT,slot = "scale.data") or head(WT@assays$RNA@scale.data)
WT <- RunPCA(WT, npcs = 100)

#Use elbow plot to check the PCs
pdf(paste0(fig.dir,"ElbowPlot100PCs.pdf"),w=12,h=12)
ElbowPlot(object = WT,ndims = 100)
dev.off()

#Cluster cells
WT <- FindNeighbors(WT, dim = 1:75)
WT <- FindClusters(WT, resolution = 1.25)
WT <- RunUMAP(WT, dims = 1:75)

pdf(paste0(fig.dir,"UmapWT.pdf"),w=12,h=12)
DimPlot(WT, reduction = "umap", label = T)
dev.off()

#Generate dotplot
genelist = c("Olig2","Bcas1","Apoe","Glul","Rlbp1","Gfap","Crip1","Id3","Actd2","C1qa","Kcnj8","Pecam1","Rcvrn","Ptprc","Cd53","Laptm5","Rho","Gngt1","Arr3","Gngt2","Opn1mw","Opn1sw","Vsx2","Otx2","Sox2","Isl1","Grm6","Scgn","Tacr3","Syt2","Neto1","Irx6","Prkar2b","Grik1","Kcng4","Cabp5","Vsx1","Prkca","Lhx1","Onecut1","Onecut2","Rbpms","Slc17a6","Thy1","Pou4f1","Opn4","Spp1", "Pax6","Tfap2a", "Tfap2b","Tfap2c", "Gad1", "Gad2","Slc32a1","Slc6a9","Slc17a8","Gjd2")
pdf(paste0(fig.dir,"DotPlotAllMarkers_WT.pdf"),w=12,h=12)
DotPlot(WT, features=genelist) + RotatedAxis() + coord_flip()
dev.off()

#Save the Seurate object
saveRDS(WT, file = "/Users/namrathagujje/Downloads/NamrathaProject/Fezf1KO-scRNA-seq/L70WT.R")

#refine KO object
#Remove Clusters: #2, 5, 12, 19, 33, 41, 42, 43, 45, 47, from 7031 cells to 6025 cells
fig.dir <- c("/Users/namrathagujje/Downloads/NamrathaProject/Fezf1KO-scRNA-seq/L70WT.R")

RGCWTcluster = c (0:3, 5:6, 8:10, 12:17, 19:20, 22:33, 35:43, 46, 48, 50)
WTclean1 <- subset(WT, idents = RGCWTcluster)

WTclean1 <- NormalizeData(WTclean1, normalization.method = "LogNormalize", scale.factor = 10000)
WTclean1 <- FindVariableFeatures(WTclean1,selection.method = "vst")
WTclean1 <- ScaleData(WTclean1, vars.to.regress = c("nCount_RNA","nFeature_RNA", "percent.mt","percent.rp"))
WTclean1 <- RunPCA(WTclean1, npcs = 100)

dev.off()
pdf(paste0(fig.dir,"ElbowPlot100PCs_WTclean1.pdf"),w=12,h=12)
ElbowPlot(object = WTclean1,ndims = 100)
dev.off()

#Cluster cells
WTclean1 <- FindNeighbors(WTclean1, dim = 1:75)
WTclean1 <- FindClusters(WTclean1, resolution = 1.25)
WTclean1 <- RunUMAP(WTclean1, dims = 1:75)

pdf(paste0(fig.dir,"UmapWTclean1.pdf"),w=12,h=12)
DimPlot(WTclean1, reduction = "umap", label = T)
dev.off()

genelist = c("Olig2","Bcas1","Apoe","Glul","Rlbp1","Gfap","Crip1","Id3","Actd2","C1qa","Kcnj8","Pecam1","Rcvrn","Ptprc","Cd53","Laptm5","Rho","Gngt1","Arr3","Gngt2","Opn1mw","Opn1sw","Vsx2","Otx2","Sox2","Isl1","Grm6","Scgn","Tacr3","Syt2","Neto1","Irx6","Prkar2b","Grik1","Kcng4","Cabp5","Vsx1","Prkca","Lhx1","Onecut1","Onecut2","Rbpms","Slc17a6","Thy1","Pou4f1","Opn4","Spp1", "Pax6","Tfap2a", "Tfap2b","Tfap2c", "Gad1", "Gad2","Slc32a1","Slc6a9","Slc17a8","Gjd2")
pdf(paste0(fig.dir,"DotPlotAllMarkers_WTclean1.pdf"),w=12,h=12)
DotPlot(WTclean1, features=genelist) + RotatedAxis() + coord_flip()
dev.off()

dim(WTclean1)

#Further remove #5, 18 reduce to 5779
#RGCWTcluster2 = c (0:5, 6:17, 19:42)

#WTclean2 <- subset(WTclean1, idents = RGCWTcluster2)
#dim(WTclean2)

#WTclean2 <- NormalizeData(WTclean2, normalization.method = "LogNormalize", scale.factor = 10000)
#WTclean2 <- FindVariableFeatures(WTclean2,selection.method = "vst")
#WTclean2 <- ScaleData(WTclean2, vars.to.regress = c("nCount_RNA","nFeature_RNA", "percent.mt","percent.rp"))
#WTclean2 <- RunPCA(WTclean2, npcs = 100)

#dev.off()
#pdf(paste0(fig.dir,"ElbowPlot100PCs_WTclean2.pdf"),w=12,h=12)
#ElbowPlot(object = WTclean2,ndims = 100)
#dev.off()

#Cluster cells
#WTclean2 <- FindNeighbors(WTclean2, dim = 1:75)
#WTclean2 <- FindClusters(WTclean2, resolution = 1)
#WTclean2 <- RunUMAP(WTclean2, dims = 1:75)

#pdf(paste0(fig.dir,"UmapWTclean2.pdf"),w=12,h=12)
#DimPlot(WTclean2, reduction = "umap", label = T)
#dev.off()

#genelist = c("Olig2","Bcas1","Apoe","Glul","Rlbp1","Gfap","Crip1","Id3","Actd2","C1qa","Kcnj8","Pecam1","Rcvrn","Ptprc","Cd53","Laptm5","Rho","Gngt1","Arr3","Gngt2","Opn1mw","Opn1sw","Vsx2","Otx2","Sox2","Isl1","Grm6","Scgn","Tacr3","Syt2","Neto1","Irx6","Prkar2b","Grik1","Kcng4","Cabp5","Vsx1","Prkca","Lhx1","Onecut1","Onecut2","Rbpms","Slc17a6","Thy1","Pou4f1","Opn4","Spp1", "Pax6","Tfap2a", "Tfap2b","Tfap2c", "Gad1", "Gad2","Slc32a1","Slc6a9","Slc17a8","Gjd2")
#pdf(paste0(fig.dir,"DotPlotAllMarkers_WTcleanup2.pdf"),w=12,h=12)
#DotPlot(WTclean2, features=genelist) + RotatedAxis() + coord_flip()
#dev.off()

#Further remove #38 reduce to 5764
#RGCWTcluster3 = c (0:37)

#WTclean3 <- subset(WTclean2, idents = RGCWTcluster3)
#dim(WTclean3)

#WTclean3 <- NormalizeData(WTclean3, normalization.method = "LogNormalize", scale.factor = 10000)
#WTclean3 <- FindVariableFeatures(WTclean3,selection.method = "vst")
#WTclean3 <- ScaleData(WTclean3, vars.to.regress = c("nCount_RNA","nFeature_RNA", "percent.mt","percent.rp"))
#WTclean3 <- RunPCA(WTclean3, npcs = 100)

#dev.off()
#pdf(paste0(fig.dir,"ElbowPlot100PCs_WTclean3.pdf"),w=12,h=12)
#ElbowPlot(object = WTclean3,ndims = 100)
#dev.off()

#Cluster cells
#WTclean3 <- FindNeighbors(WTclean3, dim = 1:75)
#WTclean3 <- FindClusters(WTclean3, resolution = 1)
#WTclean3 <- RunUMAP(WTclean3, dims = 1:75)

#pdf(paste0(fig.dir,"UmapWTclean3.pdf"),w=12,h=12)
#DimPlot(WTclean3, reduction = "umap", label = T)
#dev.off()

#genelist = c("Olig2","Bcas1","Apoe","Glul","Rlbp1","Gfap","Crip1","Id3","Actd2","C1qa","Kcnj8","Pecam1","Rcvrn","Ptprc","Cd53","Laptm5","Rho","Gngt1","Arr3","Gngt2","Opn1mw","Opn1sw","Vsx2","Otx2","Sox2","Isl1","Grm6","Scgn","Tacr3","Syt2","Neto1","Irx6","Prkar2b","Grik1","Kcng4","Cabp5","Vsx1","Prkca","Lhx1","Onecut1","Onecut2","Rbpms","Slc17a6","Thy1","Pou4f1","Opn4","Spp1", "Pax6","Tfap2a", "Tfap2b","Tfap2c", "Gad1", "Gad2","Slc32a1","Slc6a9","Slc17a8","Gjd2")
#pdf(paste0(fig.dir,"DotPlotAllMarkers_WTcleanup3.pdf"),w=12,h=12)
##DotPlot(WTclean3, features=genelist) + RotatedAxis() + coord_flip()
#dev.off()

#saveRDS(WTclean3, file = "./Mouse/Fezf1koMapping/SeuratObject/Fezf1WTclean3.rds")
saveRDS(WTclean1, file = "/Users/namrathagujje/Downloads/NamrathaProject/Fezf1KO-scRNA-seq/L70WT.R")
#saveRDS(WTclean2, file = "./Mouse/Fezf1koMapping/SeuratObject/Fezf1WTclean2.rds")

WTclean1$percent.rbpms <- PercentageFeatureSet(object=WTclean1, pattern ="Rbpms")
WTclean1$percent.pou4f1 <- PercentageFeatureSet(object=WTclean1, pattern ="Pou4f1")
WTclean1$percent.slc17a6 <- PercentageFeatureSet(object=WTclean1, pattern ="Slc17a6")

install.packages("MAST")

#Find marker genes from individual cluster
#Fezf1WTMarkers=FindAllMarkers(WTclean4, test.use = "MAST", max.cells.per.ident = 2000, min.pct = 0.1, logfc.threshold = 0.25, min.cells.feature = 10, return.thresh = 0.05, only.pos = T)
#write.csv(Fezf1WTMarkers, file="./Mouse/Fezf1koMapping/SeuratObject/Fezf1KO/Fezf1KOClean3AllMarkers.csv")

#subset cells with percent.mt < 18 & percent.rp < 18 to 5457 cells
WTclean4 <- subset(WTclean1, subset = percent.mt < 18 & percent.rp < 5)

WTclean4 <- NormalizeData(WTclean4, normalization.method = "LogNormalize", scale.factor = 10000)
WTclean4 <- FindVariableFeatures(WTclean4,selection.method = "vst")
WTclean4 <- ScaleData(WTclean4, vars.to.regress = c("nCount_RNA","nFeature_RNA", "percent.mt","percent.rp"))
WTclean4 <- RunPCA(WTclean4, npcs = 100)
WTclean4 <- FindNeighbors(WTclean4, dim = 1:75)
WTclean4 <- FindClusters(WTclean4, resolution = 1.75)
WTclean4 <- RunUMAP(WTclean4, dims = 1:75)

dev.off()
pdf(paste0(fig.dir,"ElbowPlot100PCs_WTclean4.pdf"),w=12,h=12)
ElbowPlot(object = WTclean4, ndims = 100)
dev.off()

pdf(paste0(fig.dir,"UmapWTclean4.pdf"),w=12,h=12)
DimPlot(WTclean4, reduction = "umap", label = T)
dev.off()

genelist = c("Olig2","Bcas1","Apoe","Glul","Rlbp1","Gfap","Crip1","Id3","Actd2","C1qa","Kcnj8","Pecam1","Rcvrn","Ptprc","Cd53","Laptm5","Rho","Gngt1","Arr3","Gngt2","Opn1mw","Opn1sw","Vsx2","Otx2","Sox2","Isl1","Grm6","Scgn","Tacr3","Syt2","Neto1","Irx6","Prkar2b","Grik1","Kcng4","Cabp5","Vsx1","Prkca","Lhx1","Onecut1","Onecut2","Rbpms","Slc17a6","Thy1","Pou4f1","Opn4","Spp1", "Pax6","Tfap2a", "Tfap2b","Tfap2c", "Gad1", "Gad2","Slc32a1","Slc6a9","Slc17a8","Gjd2")
pdf(paste0(fig.dir,"DotPlotAllMarkers_KOcleanup4.pdf"),w=12,h=12)
DotPlot(WTclean4, features=genelist) + RotatedAxis() + coord_flip()
dev.off()

saveRDS(WTclean4, file = "/Users/namrathagujje/Downloads/NamrathaProject/Fezf1KO-scRNA-seq/L70WT.R")
#Fezf1WTMarkers2=FindAllMarkers(WTclean4, test.use = "MAST", max.cells.per.ident = 2000, min.pct = 0.1, logfc.threshold = 0.5, min.cells.feature = 10, return.thresh = 0.001, only.pos = T)
#write.csv(Fezf1WTMarkers2, file="/Users/namrathagujje/Downloads/NamrathaProject/Fezf1KO-scRNA-seq/L70WT.R")

WTclean4$Orig.cluster <- WTclean4@meta.data$seurat_clusters
saveRDS(WTclean4, file = "/Users/namrathagujje/Downloads/NamrathaProject/Fezf1KO-scRNA-seq/L70WT.R")

#Add suffix to a column
WT$Orig.cluster <- paste0("Fezf1WT_",WT$Orig.cluster)
saveRDS(WT, file = "/Users/namrathagujje/Downloads/NamrathaProject/Fezf1KO-scRNA-seq/L70WT.R")

install.packages("MAST")

#Find marker genes from individual cluster
#Fezf1WTMarkers=FindAllMarkers(WTclean4, test.use = "MAST", max.cells.per.ident = 2000, min.pct = 0.1, logfc.threshold = 0.25, min.cells.feature = 10, return.thresh = 0.05, only.pos = T)
#write.csv(Fezf1WTMarkers, file="./Mouse/Fezf1koMapping/SeuratObject/Fezf1KO/Fezf1KOClean3AllMarkers.csv")


#Perform transcriptomic mapping
setwd("/Volumes/LABS01/Peng/PengYR/DataAnalysis/Scripts/Yirong/Transcriptional_Mapping_Demo/")
source("./SingleCell_pipeline_Lin_Zhang_Nov_2021_2.R")
Prepare()
Mouse8wksRGC <- readRDS("./Mouse8wksSubsetRGC.rds")
Mouse63wksRGC <- readRDS("./Mouse63wksRGC.rds")
DimPlot(Mouse8wksRGC, reduction ="umap", label = T, pt.size = 0.8)
DimPlot(Mouse63wksRGC, reduction ="umap", label = T, pt.size = 0.8)
TranCorrespon_original(train = Mouse8wksRGC,test = Mouse63wksRGC, file1 = "Mouse8wksRGC",file2 = "MouseRGC63wks",xlab.use = "MouseRGC8wks",ylab.use = "MouseRGC63wks", hp = 800,wp = 800, hs=8, ws=10)
TranCorrespon_original(train = Mouse63wksRGC,test = Mouse8wksRGC, file1 = "MouseRGC63wks",file2 = "MouseRGC8wks",xlab.use = "MouseRGC63wks",ylab.use = "MouseRGC8wks", hp = 800,wp = 800, hs=8, ws=10)

#Run XgBoost Transcriptomic Mapping
source("./Scripts/Yirong/SingleCellSupervisedMapping-master/xgboost_script.R")
train <- Seurat::AddMetaData(KOclean3,Idents(KOclean3),col.name = "celltype")
train_data <- as.matrix(Seurat::GetAssayData(KOclean3, layer = "data"))
train_id <- train$celltype
train_vargenes <- VariableFeatures(train)

test <- Seurat::AddMetaData(WTclean4,Idents(WTclean4),col.name = "celltype")
test_data  <-  as.matrix(Seurat::GetAssayData(WTclean4, slot = "data"))
test_id <- test$celltype
test_vargenes <- VariableFeatures(test)

var.genes_train = intersect(train_vargenes, test_vargenes)
train_data = train_data[var.genes_train,]; test_data = test_data[var.genes_train,]

XGBoost_train(train_data, var.genes = var.genes_train, do.scale = TRUE)

XGBoost_train = function(train_object, var.genes = NULL, do.scale=FALSE, max.cells.per.ident = 1000, train.frac = 0.9, weights = TRUE, plot.validation=FALSE){
  
  # Subsetting and Scaling
  train_Data = as.matrix(train_object@assays$RNA[var.genes,])
  scale.mean = rowMeans(train_Data)
  scale.var = apply(train_Data, 1, sd)
  train_Data = t(scale(t(train_Data), center=scale.mean, scale=scale.var))
  
  training.set = c(); validation.set=c()
  training.label = c(); validation.label=c();
  print(paste0("Using either ", train.frac*100, " percent cells or ", max.cells.per.ident, " cells per cluster for training, whichever is smaller"))
  
  for (cc in c(1:length(levels(Idents(train_object))))){
    i = levels(Idents(train_object))[cc]
    cells.in.clust = WhichCells(train_object,idents=i);
    
    n = min(max.cells.per.ident, round(length(cells.in.clust)*train.frac))
    train.temp = cells.in.clust[sample(length(cells.in.clust))][1:n]
    validation.temp = setdiff(cells.in.clust, train.temp)
    training.set = c(training.set,train.temp); validation.set=c(validation.set,validation.temp)
    training.label = c(training.label, rep(cc-1,length(train.temp))); validation.label = c(validation.label, rep(cc-1, length(validation.temp)));
  }
  print(table(training.label))
  train_matrix <- xgb.DMatrix(data = t(train_Data[,training.set]), label=training.label)
  validation_matrix <- xgb.DMatrix(data = t(train_Data[,validation.set]), label=validation.label)
  
  numberOfClasses <- length(unique(training.label))
  print(numberOfClasses)
  xgb_params <- list("objective" = "multi:softprob",
                     "eval_metric" = "mlogloss",
                     "num_class" = numberOfClasses,
                     "eta" = 0.2,"max_depth"=6, subsample = 0.6)
  nround    <- 200 # number of XGBoost rounds
  if (weights==TRUE){
    freqs = table(training.label)
    weights1 = c()
    minf = min(freqs)
    maxf = max(freqs)
    for (i in 1:length(freqs)){
      f = freqs[i]
      weights1 = c(weights1, rep(0.6 + 0.4*(f-minf)/(maxf-minf), f))
    }
    bst_model <- xgb.train(params = xgb_params,
                           data = train_matrix,
                           nrounds = nround,weight = weights1)
  } else {
    bst_model <- xgb.train(params = xgb_params,
                           data = train_matrix,
                           nrounds = nround)
  }
  
  
  # Predict hold-out validation set
  validation_pred <- predict(bst_model, newdata = validation_matrix)
  validation_prediction <- matrix(validation_pred, nrow = numberOfClasses,
                                  ncol=length(validation_pred)/numberOfClasses)
  valid_predlabels=apply(validation_prediction,2,which.max)-1
  
  # Confusion matrix
  if (plot.validation){
    A = table(validation.label, valid_predlabels)
    colnames(A) = levels(Idents(train_object)); rownames(A) = colnames(A)
    plotConfusionMatrix(A, order="Row", xlab.use = "True", ylab.use = "Predicted")
  }
  
  
  to.return = list()
  to.return$bst_model = bst_model
  to.return$scale_mean = scale.mean
  to.return$scale_var = scale.var
  
  return(to.return)
  
}

class(cluster_data)
bst_model = XGBoost_train(train_data, train_labels = train_id, var.genes = var.genes_train, do.scale = TRUE)



numberOfClasses = length(levels(train_id))
test_data_scale = xgb.DMatrix(scale(t(test_data), scale=bst_model$scale_var, center = bst_model$scale_mean ))
test_pred <- predict(bst_model$bst_model, newdata = test_data_scale)
test_prediction <- matrix(test_pred, nrow = numberOfClasses,
                          ncol=length(test_pred)/numberOfClasses)

# Find best class for each cell
test_pred_margins = apply(test_prediction,2,max)
test_predlabels=apply(test_prediction,2,which.max)
names(test_predlabels) = colnames(test_data)
test_pred_names = levels(train_id)[test_predlabels]
names(test_pred_names) = names(test_predlabels)

# Plot confusion matrix
C = table(test_id, test_pred_names)
plotConfusionMatrix(C, order="Row", xlab.use = "Training clusters", ylab.use = "Test clusters", title.use = "Performance Confusion Matrix")

#Merge WT and KO cluster
Fezf1Six <- merge(WTclean4, y = KOclean3, project = "Fezf1Six")
Fezf1Six <- NormalizeData(Fezf1Six, normalization.method = "LogNormalize", scale.factor = 10000)

Fezf1Six <- FindVariableFeatures(Fezf1Six,selection.method = "vst")
Fezf1Six <- ScaleData(Fezf1Six, vars.to.regress = c("nCount_RNA","nFeature_RNA", "percent.mt","percent.rp"))
Fezf1Six <- RunPCA(Fezf1Six, npcs = 100)

fig.dir <- c("/Users/namrathagujje/Downloads/NamrathaProject/Fezf1KO-scRNA-seq/L70WT.R")
dev.off()
pdf(paste0(fig.dir,"ElbowPlot100PCs_Fezf1All.pdf"),w=12,h=12)
ElbowPlot(object = Fezf1Six,ndims = 100)
dev.off()

#Cluster cells
Fezf1Six <- FindNeighbors(Fezf1Six, dim = 1:75)
Fezf1Six <- FindClusters(Fezf1Six, resolution = 1)
Fezf1Six <- RunUMAP(Fezf1Six, dims = 1:75)

pdf(paste0(fig.dir,"Fezf1Six.pdf"),w=12,h=12)
DimPlot(Fezf1Six, reduction = "umap", label = T)
dev.off()

pdf(paste0(fig.dir,"Fezf1SixwithGenotype.pdf"),w=12,h=12)
DimPlot(Fezf1Six, reduction = "umap", group.by="Genotype",label = T)
dev.off()

genelist = c("Olig2","Bcas1","Apoe","Glul","Rlbp1","Gfap","Crip1","Id3","Actd2","C1qa","Kcnj8","Pecam1","Rcvrn","Ptprc","Cd53","Laptm5","Rho","Gngt1","Arr3","Gngt2","Opn1mw","Opn1sw","Vsx2","Otx2","Sox2","Isl1","Grm6","Scgn","Tacr3","Syt2","Neto1","Irx6","Prkar2b","Grik1","Kcng4","Cabp5","Vsx1","Prkca","Lhx1","Onecut1","Onecut2","Rbpms","Slc17a6","Thy1","Pou4f1","Opn4","Spp1", "Pax6","Tfap2a", "Tfap2b","Tfap2c", "Gad1", "Gad2","Slc32a1","Slc6a9","Slc17a8","Gjd2")
pdf(paste0(fig.dir,"DotPlotAllMarkers_Fezf1All.pdf"),w=12,h=12)
DotPlot(Fezf1Six, features=genelist) + RotatedAxis() + coord_flip()
dev.off()

saveRDS(Fezf1Six, file = "/Users/namrathagujje/Downloads/NamrathaProject/Fezf1KO-scRNA-seq/L70WT.R")

Fezf1Six$Merge.cluster <- Fezf1Six@meta.data$seurat_clusters
Fezf1Six$Merge.cluster <- paste0("SixMerge_", Fezf1Six$Merge.cluster)

#load published mouse RGC types
load("./Lamprey/Lamprey/AnalysisFromJunqiang/rgc.published.45celltypes.rda")
TranRGC <- seu
head(TranRGC@meta.data)
TranRGC<-SetIdent(TranRGC, value = TranRGC@meta.data$cell_type)

train <- Seurat::AddMetaData(TranRGC,Idents(TranRGC),col.name = "celltype")
train_data <- as.matrix(Seurat::GetAssayData(TranRGC, slot = "data"))
train_id <- train$celltype
train_vargenes <- VariableFeatures(train)

test <- Seurat::AddMetaData(Fezf1Six,Idents(Fezf1Six),col.name = "celltype")
test_data  <-  as.matrix(Seurat::GetAssayData(Fezf1Six, slot = "data"))
test_id <- test$celltype
test_vargenes <- VariableFeatures(test)

var.genes_train = intersect(train_vargenes, test_vargenes)
train_data = train_data[var.genes_train,]; test_data = test_data[var.genes_train,]

bst_model = XGBoost_train(train_data, train_labels = train_id, var.genes = var.genes_train, do.scale = TRUE)

numberOfClasses = length(levels(train_id))
test_data_scale = xgb.DMatrix(scale(t(test_data), scale=bst_model$scale_var, center = bst_model$scale_mean ))
test_pred <- predict(bst_model$bst_model, newdata = test_data_scale)
test_prediction <- matrix(test_pred, nrow = numberOfClasses,
                          ncol=length(test_pred)/numberOfClasses)

# Find best class for each cell
test_pred_margins = apply(test_prediction,2,max)
test_predlabels=apply(test_prediction,2,which.max)
names(test_predlabels) = colnames(test_data)
test_pred_names = levels(train_id)[test_predlabels]
names(test_pred_names) = names(test_predlabels)

# Plot confusion matrix
C = table(test_id, test_pred_names)
plotConfusionMatrix(C, order="Row", xlab.use = "Training clusters", ylab.use = "Test clusters", x.lab.rot = T,title.use = "Performance Confusion Matrix")

#Plot cell by dendrogram
getLeftDecendants=function(tree,node) {
  daughters<-tree$edge[which(tree$edge[,1]==node),2]
  if (daughters[1] <= (tree$Nnode+1)) return(daughters[1])
  daughter.use=getDescendants(tree,daughters[1])
  daughter.use=daughter.use[daughter.use<=(tree$Nnode+1)]
  return(daughter.use)
}

getRightDecendants=function(tree,node) {
  daughters<-tree$edge[which(tree$edge[,1]==node),2]
  if (daughters[2] <= (tree$Nnode+1)) return(daughters[2])
  daughter.use=getDescendants(tree,daughters[2])
  daughter.use=daughter.use[daughter.use<=(tree$Nnode+1)]
  return(daughter.use)
}

getDescendants<-function(tree,node,curr=NULL){
  if(is.null(curr)) curr<-vector()
  daughters<-tree$edge[which(tree$edge[,1]==node),2]
  curr<-c(curr,daughters)
  w<-which(daughters>=length(tree$tip))
  if(length(w)>0) for(i in 1:length(w)) 
    curr<-getDescendants(tree,daughters[w[i]],curr)
  return(curr)
}

Fezf1Six <- BuildClusterTree(object = Fezf1Six)
Fezf1Six@tools$BuildClusterTree
tree = Fezf1Six@tools$BuildClusterTree

ident.order=tree$tip.label
nodes.1=ident.order[getLeftDecendants(tree,tree$Nnode+2)]
nodes.2=ident.order[getRightDecendants(tree,tree$Nnode+2)]

ident.order = c(nodes.1,nodes.2)
myLevels <-ident.order
Idents(Fezf1Six) <- myLevels

table(Idents(Fezf1Six))


RGCgenes <- c("Rbpms", "Pou4f1", "Pou4f2", "Pou4f3", "Spp1","Calb1","Opn4","Eomes","Tbx20","Nmb", "Tbr1", "Jam2", "Calb2","Foxp2", "Foxp1","Irx4","Satb2", "Satb1","Cartpt","Mmp17","Col25a1","Gpr88","Fam19a4","Penk","Neurod2","Irx3", "Tusc5")

RGCgenes2 <- c("Rbpms", "Pou4f1", "Pou4f2", "Pou4f3","Calb1","Opn4","Eomes","Tbx20","Tbr1", "Jam2", "Foxp2", "Foxp1","Bnc2","Satb2", "Satb1","Cartpt","Mmp17","Col25a1","Gpr88","Penk","Neurod2","Irx3")

dotplot <- DotPlot(Fezf1Six, features = RGCgenes2, dot.scale = 8, dot.min = 0.1) + RotatedAxis() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab('') + xlab ('')

tree = Fezf1Six@tools$BuildClusterTree
data = as.hclust(tree)
data.tree=as.dendrogram(data)

dev.off()

plot(data.tree) 

plot(dotplot) 

ddata <- dendro_data(data.tree, type = "rectangle")
pdendro <- ggplot(segment(ddata)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  theme_dendro() + coord_flip() + scale_y_reverse()

ggtree_plot_yset <- (pdendro + dotplot)

plot(ggtree_plot_yset)

ggtree_plot_yset <- pdendro + (dotplot)

plot_grid(ggtree_plot_yset, NULL, dotplot, nrow = 1, rel_widths = c(0.5,-0.09, 2), align = 'h')

plot_grid(ggtree_plot_yset, dotplot, nrow = 1, rel_widths = c(0.5,-0.09, 2), align = 'h')

Fezf1Six <-SetIdent(Fezf1Six, value = Fezf1Six@meta.data$seurat_clusters)
Number <- table(Fezf1Six@active.ident, Fezf1Six@meta.data$Genotype)
totals_per_cluster <- rowSums(Number)
fractions <- data.frame(matrix(ncol = ncol(Number), nrow = nrow(Number)))
colnames(fractions) <- colnames(Number)
rownames(fractions) <- rownames(Number)

for (i in 1:nrow(Number)) {
  for (j in 1:ncol(Number)) {
    fractions[i, j] <- Number[i, j] / totals_per_cluster[i]
  }
}

print(fractions)

colnames(fractions) <- fractions[,1]
rownames(fractions) <- fractions[1,]



ggplot(long_data, aes(x = rownames(long_data), y = Fraction, fill = Genotype)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "Cluster", y = "Fraction", fill = "Genotype") +
  ggtitle("Stacked Plot of Fractions per Cluster")

ggplot(long_data, aes(x = rownames(long_data), y = Fraction, fill = Genotype)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal() +
  labs(x = "Cluster", y = "Fraction", fill = "Genotype") +
  ggtitle("Bar Plot of Fractions per Cluster by Genotype")

ggplot(long_data, aes(x = rownames(long_data), y = Fraction, fill = Genotype)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(x = "Cluster", y = "Fraction", fill = "Genotype") +
  ggtitle("Stacked Bar Plot of Fractions per Cluster")

#Subset cluster DSGCs

DSGC <- subset(Fezf1Six, idents = c(23, 18, 11, 22, 9, 14, 29, 37))

DSGC <- NormalizeData(DSGC, normalization.method = "LogNormalize", scale.factor = 10000)
DSGC <- FindVariableFeatures(DSGC,selection.method = "vst")
DSGC <- ScaleData(DSGC, vars.to.regress = c("nCount_RNA","nFeature_RNA", "percent.mt","percent.rp"))
DSGC <- RunPCA(DSGC, npcs = 50)
ElbowPlot(DSGC, ndims = 50)
DSGC <- FindNeighbors(DSGC, dim = 1:20)
DSGC <- FindClusters(DSGC, resolution = 1)
DSGC <- RunUMAP(DSGC, dims = 1:20)
DimPlot(DSGC, group.by = "Genotype")
DSGC<-SetIdent(DSGC, value = DSGC@meta.data$Genotype)

ONDSGCmarker=FindAllMarkers(ONDSGC, test.use = "MAST", max.cells.per.ident = 2000, min.pct = 0.1, logfc.threshold = 0.25, min.cells.feature = 10, return.thresh = 0.05, only.pos = T)
write.csv(ONDSGCmarker, file="./Mouse/Fezf1koMapping/SeuratObject/ONDSGCMarkers.csv")















library(ggplot2)
library(ggdendro)

hcWT <- hclust(dist(WT), "ave")
ggdendrogram(hcWT, rotate = FALSE, size = 2)

# Example data
#data_matrix <- matrix(rnorm(100), ncol = 10)  # Generate random data matrix

rownames(WT) <- paste0("Obs", 1:nrow(WT))  # Assign row names

# Calculate distance matrix
distance_matrix <- dist(WT)

# Perform hierarchical clustering
hclust_result <- hclust(distance_matrix, method = "complete")

# Plot the dendrogram
plot(hclust_result, main = "Dendrogram")

install.packages("dendextend")
library(dendextend)

class(WT)

AggregateExpression(WT)


dendrogram <- BuildClusterTree(WT)

ist_matrix <- as.dist(1 - seurat_obj@assays$RNA@scale.data)  # Compute distance matrix
hclust_tree <- hclust(dist_matrix, method = "ward.D2")  # Perform hierarchical clustering

# Plot the dendrogram
plot(hclust_tree, hang = -1, main = "Dendrogram")

plot_dendrogram(dendrogram)

dev.off()

dim(WTclean4) 

data("WTclean4")
WTclean4 <- BuildClusterTree(object = WTclean4)
windows(width = 18915, height = 6075)
PlotClusterTree(object = WTclean4)

dev.off()

# Install and load Seurat package
install.packages("Seurat")
library(Seurat)
library(dplyr)
library(patchwork)

if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("HGC")

if(!require(devtools))
  install.packages("devtools")
devtools::install_github("XuegongLab/HGC")

library(HGC)
dim(WTclean4)
dim(WT)

table(WTclean4)

WTclean4.PCs <- WTclean4[["PCs"]]






dendrogram <- BuildClusterTree(seurat_obj)

# Plot the dendrogram
plot(dendrogram)



# Assuming you have a Seurat object named seurat_obj
# Perform clustering
WTclean4 <- FindNeighbors(WTclean4, dims = 1:75)
WTclean4 <- FindClusteringTree(WTclean4, graph.type = "SNN")

# Build cluster tree
cluster_tree <- BuildClusterTree(WTclean4)

par(mar = c(2, 2, 2, 2))  # Set margins to 2 lines on each side
plot(cluster_tree)

windows(width = 10, height = 8)  # Open a new Windows device with specified dimensions
plot(cluster_tree)

# Plot dendrogram
plot(cluster_tree@meta.data$seurat_clusters)


WTclean4 <- subset(WTclean1, subset = percent.mt < 18 & percent.rp < 5)

WTclean4 <- NormalizeData(WTclean4, normalization.method = "LogNormalize", scale.factor = 10000)
WTclean4 <- FindVariableFeatures(WTclean4,selection.method = "vst")
WTclean4 <- ScaleData(WTclean4, vars.to.regress = c("nCount_RNA","nFeature_RNA", "percent.mt","percent.rp"))
WTclean4 <- RunPCA(WTclean4, npcs = 100)
WTclean4 <- FindNeighbors(WTclean4, dim = 1:75)
WTclean4 <- FindClusters(WTclean4, resolution = 1.75)
WTclean4 <- RunUMAP(WTclean4, dims = 1:75)




