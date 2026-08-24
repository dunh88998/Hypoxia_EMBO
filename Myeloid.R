library(ggplot2)
library(patchwork)
install.packages("viridis")
library(viridis)
library(tidyverse)
library(RColorBrewer)


# global clusters
Merge_Data <-RenameIdents(Merge_Data, '0'="NK cells", '1'="Myeloid 1", '2'= "B cells", 
                          '3'="Myeloid 2", '4'="Proliferating", '5'="Myeloid 3", 
                          '6'="Granulocytes", '7'="Myeloid 4" )
Idents(Merge_Data)
view(Merge_Data@meta.data)
Merge_Data$Celltypes<-Idents(Merge_Data)

# subset of Myeloid 1 from global clusters
Myeloid_1 <- subset(x = Merge_Data, idents= c("Myeloid 1"))

# clustering
Myeloid_1 <- NormalizeData(Myeloid_1, normalization.method = "LogNormalize", scale.factor = 10000)
Myeloid_1 <- FindVariableFeatures(Myeloid_1, selection.method = "vst", nfeatures = 2000)
Myeloid_1<- ScaleData(Myeloid_1)
all.genes <- rownames(Myeloid_1)
Myeloid_1 <- ScaleData(Myeloid_1, features = all.genes)
Myeloid_1 <-RunPCA(Myeloid_1, verbose=FALSE)
ElbowPlot(Myeloid_1)
Myeloid_1 <-RunHarmony(Myeloid_1, group.by.vars="orig.ident")
Myeloid_1 <- RunUMAP(Myeloid_1, reduction = "harmony", dims = 1:16)
Myeloid_1 <- FindNeighbors(Myeloid_1, dims = 1:16)
Myeloid_1 <- FindClusters(Myeloid_1, resolution = 0.4)

# visualizing UMAP
P<- DimPlot(Myeloid_1, group.by="seurat_clusters", label=TRUE)
P

# find markers for cell type annotation
Myeloid_1 <- JoinLayers(Myeloid_1)
Myeloid_1_markers <- FindAllMarkers(Myeloid_1, only.pos = TRUE,
                                            min.pct = 0.25,
                                            logfc.threshold = 0.25)
head(Myeloid_1_markers)
write.csv(Myeloid_1_markers, "Myeloid_1_markers.csv")

# remove cluster 4 and 5 
Myeloid_1_New <- subset(x = Myeloid_1, idents = c("0", "1", "2", "3", "6"))

# rename identities
Myeloid_1_New <-RenameIdents(Myeloid_1_New, '0'="CM1", '1'="NCM1", '2'= "CM2", '3'="NM2", '6'="IFN Mo")
Myeloid_1_New$celltypes1 <-Idents(Myeloid_1_New)
# visualizing UMAP
P<- DimPlot(Myeloid_1_New, group.by="celltypes1", label=TRUE)
P
##############################################################################################
# antigen-presentation related signatures
my_gene_list_1 <- list(
  module1_genes = c("H2-Ea", "H2-Ob", "H2-DMb1", "H2-Eb1", "H2-Aa", "H2-Ab1", 
                    "H2-DMb2", "Cd74", "Cd86"))

# join layers
Myeloid_1_New <- JoinLayers(Myeloid_1_New)

# Add module scores to Myeloid_1_New

Myeloid_1_New <- AddModuleScore(
  object = Myeloid_1_New,
  features = my_gene_list_1,
  name = "ModuleName" # Prefix for the metadata column names
)

# View the new metadata columns
head(Myeloid_1_New@meta.data)

# Calculate the Z-score for the module score column
Myeloid_1_New$Immune_Module_Score_Zscore<- scale(Myeloid_1_New$ModuleName1)

#Visualization
# Visualize the z-scores using a FeaturePlot (UMAP/tSNE)
plot_list <- FeaturePlot(object = Myeloid_1_New, 
                         features = "Immune_Module_Score_Zscore", 
                         keep.scale = "feature",
                         combine = FALSE)

modified_plot_list <- lapply(plot_list, function(x) {
  x + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
    labs(title = NULL)
})

# Combine the modified plots using patchwork syntax
combined_plot <- wrap_plots(modified_plot_list) + plot_annotation (title = 'Module Score of Antigen Presentation Signatures')
# Print the final plot
print(combined_plot)

VlnPlot(Myeloid_1_New, group.by="celltypes1", features = "Immune_Module_Score_Zscore") 

###############################################################################
# tolerance-associated signatures
my_gene_list_2 <- list(
  module1_genes = c("Cd209a", "Cd274", "Csf1", "Ido1", "Hmox1", 
                    "Tgfb1", "Tgfb2", "Tgfbi", "Il10", "Tgfbr1", "Tgfbr2", "Tgfbr3","Il4"))

Myeloid_1_New <- AddModuleScore(
  object = Myeloid_1_New,
  features = my_gene_list_2,
  name = "ModuleName" # Prefix for the metadata column names
)


# View the new metadata columns
head(Myeloid_1_New@meta.data)

# Calculate the Z-score for the module score column
Myeloid_1_New$Immune_Module_Score_Zscore<- scale(Myeloid_1_New$ModuleName1)

#Visualization
# Visualize the z-scores using a FeaturePlot (UMAP/tSNE)
plot_list <- FeaturePlot(object = Myeloid_1_New, 
                         features = "Immune_Module_Score_Zscore", 
                         keep.scale = "feature",
                         combine = FALSE)

modified_plot_list <- lapply(plot_list, function(x) {
  x + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
    labs(title = NULL)
})

combined_plot <- wrap_plots(modified_plot_list) + plot_annotation (title = 'Module Score of Tolerance Signatures', 
                                            theme = theme(plot.title = element_text(hjust = 0.5)))
# Print the final plot
print(combined_plot)

VlnPlot(Myeloid_1_New, group.by="celltypes1", features = "Immune_Module_Score_Zscore") 

##################################################################################
# Heat Shock signatures
my_gene_list_3 <- list(
  module1_genes = c("Hspa5", "Hsp90aa1", "Hsp90ab1", "Hspa8", "Ptges3")
)

Myeloid_1_New <- AddModuleScore(
  object = Myeloid_1_New,
  features = my_gene_list_3,
  name = "ModuleName" # Prefix for the metadata column names
)


# View the new metadata columns
head(Myeloid_1_New@meta.data)

# Calculate the Z-score for the module score column
Myeloid_1_New$Immune_Module_Score_Zscore<- scale(Myeloid_1_New$ModuleName1)

#Visualization of heat shock protein-associated module score
VlnPlot(Myeloid_1_New , split.by="orig.ident", group.by="celltypes1", features = "Immune_Module_Score_Zscore") 

#####################################################################################
# vlnplot for relative genes
VlnPlot(Myeloid_1_New, split.by="orig.ident", group.by="celltypes1", features = c("Ccr2"))
VlnPlot(Myeloid_1_New, split.by="orig.ident", group.by="celltypes1", features = c("Ccl24"))
VlnPlot(Myeloid_1_New, split.by="orig.ident", group.by="celltypes1", features = c("Hspa8"))
VlnPlot(Myeloid_1_New, split.by="orig.ident", group.by="celltypes1", features = c("Hsp90ab1"))
VlnPlot(Myeloid_1_New, split.by="orig.ident", group.by="celltypes1", features = c("Hsp90aa1"))
VlnPlot(Myeloid_1_New, split.by="orig.ident", group.by="celltypes1", features = c("Hspa5"))
