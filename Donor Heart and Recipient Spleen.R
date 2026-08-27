# load the library
library(Seurat)
library(tidyverse)
# Read barcodes, features(genes) and matrix files
Heart_raw<-Read10X(data.dir = "Heart_raw_feature_bc_matrix")
Spleen_raw<-Read10X(data.dir = "Spleen_raw_feature_bc_matrix")
gc()

#Create Seurat Objects
Heart_raw<-CreateSeuratObject(counts= Heart_raw, project = "Heart_raw", min.cells = 5, min.features = 200)
Spleen_raw<-CreateSeuratObject(counts= Spleen_raw, project = "Spleen_raw", min.cells = 5, min.features = 200)

# Store mitochondrial percentage in object meta data
Heart_raw <-PercentageFeatureSet(Heart_raw, pattern = "^mt-", col.name = "percent.mt")
view(Heart_raw@meta.data)
Spleen_raw <-PercentageFeatureSet(Spleen_raw, pattern = "^mt-", col.name = "percent.mt")
view(Spleen_raw@meta.data)

# visualation
VlnPlot(Heart_raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Spleen_raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Cut-offs subset data
Heart_raw<- subset(Heart_raw, subset=nFeature_RNA>200 & nFeature_RNA< 6000 & nCount_RNA<32000 & percent.mt < 10)
Spleen_raw<- subset(Spleen_raw, subset=nFeature_RNA>500 & nFeature_RNA< 6000 & nCount_RNA<25000 & percent.mt < 5)

# visualation after cut-off
VlnPlot(Heart_raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Spleen_raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

######################################################################################################
#QC and clustering for Heart_raw
Heart_raw <- NormalizeData(Heart_raw, normalization.method = "LogNormalize", scale.factor = 10000)
Heart_raw <- FindVariableFeatures(Heart_raw, selection.method = "vst", nfeatures = 2000)
Heart_raw<- ScaleData(Heart_raw)
all.genes <- rownames(Heart_raw)
Heart_raw <- ScaleData(Heart_raw, features = all.genes)
Heart_raw <-RunPCA(Heart_raw, verbose=FALSE)
ElbowPlot(Heart_raw, 50)
Heart_raw<- RunUMAP(Heart_raw, reduction = "pca", dims = 1:15)
Heart_raw <- FindNeighbors(Heart_raw, dims = 1:15)
Heart_raw <- FindClusters(Heart_raw, resolution = 0.1)

DimPlot(Heart_raw, reduction = "umap", label = T)

# find all markers for annotaion
Heart_raw <- JoinLayers(Heart_raw)
Heart_raw_markers <- FindAllMarkers(Heart_raw, only.pos = TRUE,
                                        min.pct = 0.25,
                                        logfc.threshold = 0.25)
head(Heart_raw_markers)
write.csv(Heart_raw_markers, "Heart_raw_markers.csv")

# removing cluster 7 (endothelial cells)
Heart_raw_New <- subset(x = Heart_raw, idents= c("0", "1", "2", "3", "4", "5", "6"))

# annotation
Heart_raw_New <-RenameIdents(Heart_raw_New, '0'="Mac 1", '1'="Mac 2", '2'= "Dendritic Cells", 
                             '3'="NK Cells", '4'="Proliferating", '5'="T Cells", '6'="B Cells")
Idents(Heart_raw_New)
DimPlot(Heart_raw_New, reduction = "umap", label = T)
DimPlot(Heart_raw_New, reduction = "umap", label = F)
##################################################################################
#QC and clustering for Spleen_raw
Spleen_raw <- NormalizeData(Spleen_raw, normalization.method = "LogNormalize", scale.factor = 10000)
Spleen_raw <- FindVariableFeatures(Spleen_raw, selection.method = "vst", nfeatures = 2000)
Spleen_raw<- ScaleData(Spleen_raw)
all.genes <- rownames(Spleen_raw)
Spleen_raw <- ScaleData(Spleen_raw, features = all.genes)
Spleen_raw <-RunPCA(Spleen_raw, verbose=FALSE)
ElbowPlot(Spleen_raw, 50)
Spleen_raw<- RunUMAP(Spleen_raw, reduction = "pca", dims = 1:6)
Spleen_raw <- FindNeighbors(Spleen_raw, dims = 1:6)
Spleen_raw <- FindClusters(Spleen_raw, resolution = 0.2)
DimPlot(Spleen_raw, reduction = "umap", label = T)

# find all markers for annotaion
Spleen_raw <- JoinLayers(Spleen_raw)
Spleen_raw_markers <- FindAllMarkers(Spleen_raw, only.pos = TRUE,
                                    min.pct = 0.25,
                                    logfc.threshold = 0.25)
head(Spleen_raw_markers)
write.csv(Spleen_raw_markers, "Spleen_raw_markers.csv")

# annotation
Spleen_raw <-RenameIdents(Spleen_raw, '0'="NK cells", '1'="Non-Classical Mon", '2'= "Dendritic Cells", '3'="Classical Mon", '4'="B Cells", '5'="MHC-II Myeloid", '6'="Regulatory Myeloid", '7'="Proliferating" )
Idents(Spleen_raw)

DimPlot(Spleen_raw, reduction = "umap", label = T)
DimPlot(Spleen_raw, reduction = "umap", label = F)
