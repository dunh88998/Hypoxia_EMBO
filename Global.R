# load the library
library(Seurat)
library(tidyverse)
library(harmony)
library(RColorBrewer)
library(dplyr)
library(ggplot2)

################################################################################
# Read barcodes, features(genes) and matrix files
Denver_raw<-Read10X(data.dir="Denver_raw")
Sea_level_raw<-Read10X(data.dir = "Sea_level_raw")
Post_HTx_raw<-Read10X(data.dir = "raw_feature_bc_matrix")
gc()

#Create Seurat Objects
Denver_raw<-CreateSeuratObject(counts = Denver_raw, project = "Denver_raw", min.cells = 5, min.features = 200)
Sea_level_raw<-CreateSeuratObject(counts= Sea_level_raw, project = "Sea_level_raw", min.cells = 5, min.features = 200)
Post_HTx_raw<-CreateSeuratObject(counts= Post_HTx_raw, project = "Post_HTx_raw", min.cells = 5, min.features = 200)

# Create a list of Seurat objects
pbmc_list <- list(Denver_raw, Sea_level_raw, Post_HTx_raw)

# Add metadata to each object in the list
for (i in seq_along(pbmc_list)) {
  pbmc_list[[i]] <- AddMetaData(pbmc_list[[i]], paste0("sample", i), col.name = "sample")
}

# Merge all objects in the list
pbmc_merged <- merge(x = pbmc_list[[1]], y = pbmc_list[2:length(pbmc_list)], add.cell.ids = paste0("sample", 1:length(pbmc_list)))

# Check meta.data
view(pbmc_merged@meta.data)

# Store mitochondrial percentage in object meta data
Merge_Data <-PercentageFeatureSet(pbmc_merged, pattern = "^mt-", col.name = "percent.mt")
view(Merge_Data@meta.data)

# check Merged_data
VlnPlot(Merge_Data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# filter out low-quality single cells
Merge_Data<- subset(Merge_Data, subset=nFeature_RNA>500 & nFeature_RNA< 6500 & nCount_RNA<32000 & percent.mt < 5)

# check Merged_data after cutoff
VlnPlot(Merge_Data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
################################################################################
#normalization and find variable features
Merge_Data <- NormalizeData(Merge_Data, normalization.method = "LogNormalize", scale.factor = 10000)
Merge_Data <- FindVariableFeatures(Merge_Data, selection.method = "vst", nfeatures = 2000)

# Top features/feature plot
plot1 <- VariableFeaturePlot(Merge_Data)
plot1 + ggtitle("Variable Feature Plot")

# Identify the 10, 20 and 30 most highly variable genes
top10 <- head(VariableFeatures(Merge_Data), 10)
top20 <- head(VariableFeatures(Merge_Data), 20)
top30 <- head(VariableFeatures(Merge_Data), 30)

plot10 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot10 + ggtitle("The 10 most highly variable genes")

plot20 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
plot20 + ggtitle("The 20 most highly variable genes")

plot30 <- LabelPoints(plot = plot1, points = top30,  repel = TRUE)
plot30 + ggtitle("The 30 most highly variable genes")
################################################################################
# scale data
Merge_Data<- ScaleData(Merge_Data) # 2000 identified variable features
all.genes <- rownames(Merge_Data)
Merge_Data <- ScaleData(Merge_Data, features = all.genes)

# run PCA
Merge_Data <-RunPCA(Merge_Data, verbose=FALSE)
ElbowPlot(Merge_Data, 50)

# run Harmony
Merge_Data <-RunHarmony(Merge_Data, group.by.vars="orig.ident")

# clustering
Merge_Data <- RunUMAP(Merge_Data, reduction = "harmony", dims = 1:16)
Merge_Data <- FindNeighbors(Merge_Data, dims = 1:16)
Merge_Data <- FindClusters(Merge_Data, resolution = 0.054)

# visualizing UMAP
P <- DimPlot(Merge_Data, group.by="seurat_clusters", label=T)
P
################################################################################
#Find markers for cell type annotation
Merge_Data <- JoinLayers(Merge_Data)
Merge_Data_all_markers <- FindAllMarkers(Merge_Data, only.pos = TRUE,
                                         min.pct = 0.25,
                                         logfc.threshold = 0.25)
head(Merge_Data_all_markers)
write.csv(Merge_Data_all_markers, "Merge_Data_all_markers.csv")

#Annotations
Merge_Data <-RenameIdents(Merge_Data, '0'="NK cells", '1'="Myeloid 1", 
                          '2'= "B cells", '3'="Myeloid 2", '4'="Proliferating", 
                          '5'="Myeloid 3", '6'="Granulocytes", '7'="Myeloid 4" )

Idents(Merge_Data)
view(Merge_Data@meta.data)
#################################################################################
# global cluster composition
pt <- table(Idents(Merge_Data), Merge_Data$orig.ident)
pt <- as.data.frame(pt)
pt$Freq<- as.numeric(pt$Freq)

# Calculate total pt and add percentage column
pt_with_percentage <- pt%>%
  group_by(Var2) %>%
  mutate(Percentage = round(Freq /sum(Freq)*100, 1)) %>%
  ungroup()
colnames(pt_with_percentage)
print(pt_with_percentage)

#rename(Cell_Types = Var1, Groups = Var2)
colnames(pt_with_percentage)[1:2] <- c("Cell_Types", "Groups")
pt_with_percentage <- as.data.frame(pt_with_percentage)

ggplot(pt_with_percentage, aes(x = Groups, y = Percentage, fill = Cell_Types))+
  geom_bar(stat = "identity", width = 0.5)+
  geom_text(aes(label = paste("")), position = position_stack(vjust = 0.5))+
  theme(legend.title = element_blank())

###############################################################################
# Heatmap showing the top differentially expressed genes across annotated global clusters. 
# Calculate markers for all clusters
#Identify DEGs for every cluster compared to all other clusters.
Merge_Data<- JoinLayers(Merge_Data)
all_markers <- FindAllMarkers(Merge_Data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#Filter Top 5 Genes per Cluster
#Use dplyr to group the results by cluster and 
#select the top 5 genes based on the average log-fold change.
top5_genes <- all_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

# Generate the Heatmap
# use scaled data for visualization.
# Ensure data is scaled for the selected features

Merge_Data <- ScaleData(Merge_Data, features = top5_genes$gene)

# Plot heatmap
DoHeatmap(Merge_Data, features = top5_genes$gene) + NoLegend()



