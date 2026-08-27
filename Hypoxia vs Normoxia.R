# load the library
library(Seurat)
library(ggplot2)
library(patchwork)
library(viridis)
library(tidyverse)
library(RColorBrewer)
library(dplyr)
library(clusterProfiler)
library(harmony)

# Read barcodes, features(genes) and matrix files
Denver_raw<-Read10X(data.dir="//data.ucdenver.pvt/dept/SOM/CARD/CARD/KopeckyShared/Hao/Single Cell/new single cell project/Denver_raw")
Sea_level_raw<-Read10X(data.dir = "//data.ucdenver.pvt/dept/SOM/CARD/CARD/KopeckyShared/Hao/Single Cell/new single cell project/Sea_level_raw")
gc()

#Create Seurat Objects
Denver_raw<-CreateSeuratObject(counts = Denver_raw, project = "Denver_raw", min.cells = 5, min.features = 200)
Sea_level_raw<-CreateSeuratObject(counts= Sea_level_raw, project = "Sea_level_raw", min.cells = 5, min.features = 200)

# Merge two of Seurat Objects
Merge_Two<-merge(Denver_raw, y=Sea_level_raw, add.cell.ids = c(x="Denver_raw", y="Sea_level_raw"), project="Merge_Data")
view(Merge_Two@meta.data)

# Store mitochondrial percentage in object meta data
Merge_Two <-PercentageFeatureSet(Merge_Two, pattern = "^mt-", col.name = "percent.mt")
view(Merge_Two@meta.data)

# check Merged data
VlnPlot(Merge_Two, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Cut-offs
Merge_Two<- subset(Merge_Two, subset=nFeature_RNA>500 & nFeature_RNA< 6500 & nCount_RNA<32000 & percent.mt < 5)

# check Merged data after cut-off
VlnPlot(Merge_Data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# QC and clustering
Merge_Two <- NormalizeData(Merge_Two, normalization.method = "LogNormalize", scale.factor = 10000)
Merge_Two <- FindVariableFeatures(Merge_Two, selection.method = "vst", nfeatures = 2000)
Merge_Two<- ScaleData(Merge_Two)
all.genes <- rownames(Merge_Two)
Merge_Two <- ScaleData(Merge_Two, features = all.genes)
Merge_Two <-RunPCA(Merge_Two, verbose=FALSE)
ElbowPlot(Merge_Two, 50)
Merge_Two <-RunHarmony(Merge_Two, group.by.vars="orig.ident")
Merge_Two <- RunUMAP(Merge_Two, reduction = "harmony", dims = 1:16, min.dist = 0.2)
Merge_Two <- FindNeighbors(Merge_Two, dims = 1:16)
Merge_Two <- FindClusters(Merge_Two, resolution = 0.06)

# visualizing UMAP
P<- DimPlot(Merge_Two, group.by="seurat_clusters", label=T)
P

#Find markers for cell type annotation
Merge_Two <- JoinLayers(Merge_Two)
Merge_Two_markers <- FindAllMarkers(Merge_Two, only.pos = TRUE,
                                         min.pct = 0.25,
                                         logfc.threshold = 0.25)

write.csv(Merge_Two_markers, "Merge_Two_markers.csv")

# rename identities
Merge_Two <-RenameIdents(Merge_Two, '0'="NK cells", '1'="Myeloid 1", '2'= "B cells", 
                         '3'="Myeloid 2", '4'="Myeloid 3", '5'="Proliferating", 
                         '6'="Granulocytes", '7'="Myeloid 4" )

# add a new column
Merge_Two$Celltypes <-Idents(Merge_Two)

################################################################################
# global cluster composition
pt <- table(Idents(Merge_Two), Merge_Two$orig.ident)
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

#Visualization
ggplot(pt_with_percentage, aes(x = Groups, y = Percentage, fill = Cell_Types))+
  geom_bar(stat = "identity", width = 0.5)+
  geom_text(aes(label = paste("")), position = position_stack(vjust = 0.5))+
  theme(legend.title = element_blank())

################################################################################
# antigen-presentation related signatures
my_gene_list_1 <- list(
  module1_genes = c("H2-Ea", "H2-Ob", "H2-DMb1", "H2-Eb1", "H2-Aa", "H2-Ab1", 
                    "H2-DMb2", "Cd74", "Cd86"))

# join layers
Merge_Two <- JoinLayers(Merge_Two)

# Add module scores
Merge_Two <- AddModuleScore(
  object = Merge_Two,
  features = my_gene_list_1,
  name = "ModuleName" # Prefix for the metadata column names
)

# View the new metadata columns
head(Merge_Two@meta.data)

# Calculate the Z-score for the module score column
Merge_Two$Immune_Module_Score_Zscore<- scale(Merge_Two$ModuleName1)

#Visualization
# Visualize the z-scores using a FeaturePlot (UMAP/tSNE)
plot_list <- FeaturePlot(object = Merge_Two, 
                         features = "Immune_Module_Score_Zscore", 
                         split.by = "orig.ident",
                         ncol = 2, keep.scale = "feature",
                         combine = FALSE)

modified_plot_list <- lapply(plot_list, function(x) {
  x + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
    labs(title = NULL)
})

# Combine the modified plots using patchwork syntax
combined_plot <- wrap_plots(modified_plot_list) + 
  plot_annotation (title = 'Module Score of Antigen Presentation Signatures')
# Print the final plot
print(combined_plot)

VlnPlot(Merge_Two, split.by="orig.ident", group.by="Celltypes", 
        features = "Immune_Module_Score_Zscore") 

################################################################################
# tolerance-associated signatures
my_gene_list_2 <- list(
  module1_genes = c("Cd209a", "Cd209b", "Cd274", "Ido1", "Hmox1", 
                    "Tgfb1", "Tgfb2", "Tgfbr1", "Tgfbr2")
)

Merge_Two <- AddModuleScore(
  object = Merge_Two,
  features = my_gene_list_2,
  name = "ModuleName" # Prefix for the metadata column names
)


# View the new metadata columns
head(Merge_Two@meta.data)

# Calculate the Z-score for the module score column
Merge_Two$Immune_Module_Score_Zscore<- scale(Merge_Two$ModuleName1)

#Visualization
# Visualize the z-scores using a FeaturePlot (UMAP/tSNE)
plot_list <- FeaturePlot(object = Merge_Two, 
                         features = "Immune_Module_Score_Zscore", 
                         split.by = "orig.ident",
                         ncol = 2, keep.scale = "feature",
                         combine = FALSE)

modified_plot_list <- lapply(plot_list, function(x) {
  x + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
    labs(title = NULL)
})

combined_plot <- wrap_plots(modified_plot_list) + 
  plot_annotation (title = 'Module Score of Tolerance Signatures', 
  theme = theme(plot.title = element_text(hjust = 0.5)))

# Print the final plot
print(combined_plot)

VlnPlot(Merge_Two, split.by="orig.ident", group.by="Celltypes", 
        features = "Immune_Module_Score_Zscore") 
#################################################################################
# Gene Ontology (GO) enrichment analysis
# add condition into metadata
Merge_Two$Condition <- Merge_Two$orig.ident 
view(Merge_Two@meta.data)

# Set the identity to the 'condition' variable
Idents(Merge_Two) <- Merge_Two$Condition

Merge_Two <- JoinLayers(Merge_Two)
degs <- FindMarkers(Merge_Two, 
                    ident.1 = "Denver_raw", 
                    ident.2 = "Sea_level_raw",
                    group.by = "Condition",
                    min.pct = 0.25) # Optional: filter out genes expressed in < 25% of cells

write.csv(degs, "degs.csv")

significant_degs <- degs %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25)

# Identify up-regulated genes in Denver_raw compared to Sea_level_raw
upregulated_in_group_Denver <- subset(significant_degs, avg_log2FC > 0)

# Identify up-regulated genes in Sea_level_raw compared to Denver_raw
upregulated_in_group_Sea_Level <- subset(significant_degs, avg_log2FC < 0)

# select top 50 DEGs from each groups
top50_upregulated_in_group_Denver <- upregulated_in_group_Denver %>%
  arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
  head(n = 50)
top50_upregulated_in_group_Denver_list <- rownames(top50_upregulated_in_group_Denver)


top50_upregulated_in_group_Sea_Level <- upregulated_in_group_Sea_Level %>%
  arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
  head(n = 50)
top50_upregulated_in_group_Sea_Level_list <- rownames(top50_upregulated_in_group_Sea_Level)

# Convert gene symbols to Entrez IDs
top50_denver_entrez_ids <- bitr(top50_upregulated_in_group_Denver_list,
                                fromType = "SYMBOL",
                                toType = "ENTREZID",
                                OrgDb = "org.Mm.eg.db")

top50_sea_level_entrez_ids <- bitr(top50_upregulated_in_group_Sea_Level_list,
                                   fromType = "SYMBOL",
                                   toType = "ENTREZID",
                                   OrgDb = "org.Mm.eg.db")


#Run Gene Ontology (GO) enrichment analysis. 
denver_go_enrichment <- enrichGO(gene = top50_denver_entrez_ids$ENTREZID,
                                    OrgDb = "org.Mm.eg.db",
                                    ont = "BP", # 'BP' for Biological Process
                                    pvalueCutoff = 0.05,
                                    readable = TRUE)

dotplot(denver_go_enrichment, showCategory = 10, 
        title = "Top 10 BP Enriched Pathways in Denver")


sea_level_go_enrichment <- enrichGO(gene = top50_sea_level_entrez_ids$ENTREZID,
                                       OrgDb = "org.Mm.eg.db",
                                       ont = "BP",            
                                       pvalueCutoff = 0.05,
                                       readable = TRUE)

dotplot(sea_level_go_enrichment, showCategory = 10, 
        title = "Top 10 BP Enriched Pathways in Sea Level")
#################################################################################
# set Merge_Two identity to Celltypes
Idents(Merge_Two)<-Merge_Two$Celltypes

# select Myeloid cells
Myeloid_cells<- subset(x = Merge_Two,idents= c("Myeloid 1", "Myeloid 2", "Myeloid 3", "Myeloid 4"))

# Check the extraction result 
DimPlot(Myeloid_cells, reduction = "umap", label = T)

# clustering Myeloid Cells
Myeloid_cells <- NormalizeData(Myeloid_cells, normalization.method = "LogNormalize", scale.factor = 10000)
Myeloid_cells <- FindVariableFeatures(Myeloid_cells, selection.method = "vst", nfeatures = 2000)
Myeloid_cells<- ScaleData(Myeloid_cells) # 2000 identified variable features
all.genes <- rownames(Myeloid_cells)
Myeloid_cells <- ScaleData(Myeloid_cells, features = all.genes)
Myeloid_cells <-RunPCA(Myeloid_cells, verbose=FALSE)
ElbowPlot(Myeloid_cells, 50)
Myeloid_cells <-RunHarmony(Myeloid_cells, group.by.vars="orig.ident")
Myeloid_cells <- RunUMAP(Myeloid_cells, reduction = "harmony", dims = 1:16)
Myeloid_cells <- FindNeighbors(Myeloid_cells, dims = 1:16)
Myeloid_cells <- FindClusters(Myeloid_cells, resolution = 0.1)

# visualize UMAP
P<- DimPlot(Myeloid_cells, group.by="seurat_clusters", label=TRUE)
P
################################################################################
#Find markers for cell type annotation
Myeloid_cells <- JoinLayers(Myeloid_cells)
Myeloid_cells_markers <- FindAllMarkers(Myeloid_cells, only.pos = TRUE,
                                    min.pct = 0.25,
                                    logfc.threshold = 0.25)

write.csv(Myeloid_cells_markers, "Myeloid_cells_markers.csv")

# removing cluster 3 and 4
New_Myeloid <- subset(x = Myeloid_cells,idents= c("0", "1", "2", "5", "6"))
Idents(New_Myeloid) <- New_Myeloid$seurat_clusters
Idents(New_Myeloid)

# rename clusters
New_Myeloid <-RenameIdents(New_Myeloid, '0'="Classical Monocyte", 
                           '1'="NCM", '2'= "Mo-DCs", '5'="cDC1", '6'="cDc2")

# add new meta data
New_Myeloid$celltypes<-Idents(New_Myeloid)
view(New_Myeloid@meta.data)
#################################################################################
# Myeloid cells composition
pt <- table(Idents(New_Myeloid), New_Myeloid$orig.ident)
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

#Visualization
ggplot(pt_with_percentage, aes(x = Groups, y = Percentage, fill = Cell_Types))+
  geom_bar(stat = "identity", width = 0.5)+
  geom_text(aes(label = paste("")), position = position_stack(vjust = 0.5))+
  theme(legend.title = element_blank())
################################################################################
# Gene Ontology (GO) enrichment analysis
# Set the identity to the 'condition' variable
Idents(New_Myeloid) <- New_Myeloid$Condition
Idents(New_Myeloid)
# join layers
New_Myeloid <- JoinLayers(New_Myeloid)
New_Myeloid_degs <- FindMarkers(New_Myeloid, 
                    ident.1 = "Denver_raw", 
                    ident.2 = "Sea_level_raw",
                    group.by = "Condition",
                    min.pct = 0.25) 

write.csv(New_Myeloid_degs, "New_Myeloid_degs.csv")

significant_degs <- New_Myeloid_degs %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25)

# Identify up-regulated genes in Denver_raw compared to Sea_level_raw
upregulated_in_group_Denver <- subset(significant_degs, avg_log2FC > 0)

# Identify up-regulated genes in Sea_level_raw compared to Denver_raw
upregulated_in_group_Sea_Level <- subset(significant_degs, avg_log2FC < 0)

# select top 50 DEGs from each groups
top50_upregulated_in_group_Denver <- upregulated_in_group_Denver %>%
  arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
  head(n = 50)
top50_upregulated_in_group_Denver_list <- rownames(top50_upregulated_in_group_Denver)


top50_upregulated_in_group_Sea_Level <- upregulated_in_group_Sea_Level %>%
  arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
  head(n = 50)
top50_upregulated_in_group_Sea_Level_list <- rownames(top50_upregulated_in_group_Sea_Level)

# Convert gene symbols to Entrez IDs
top50_denver_entrez_ids <- bitr(top50_upregulated_in_group_Denver_list,
                                fromType = "SYMBOL",
                                toType = "ENTREZID",
                                OrgDb = "org.Mm.eg.db")

top50_sea_level_entrez_ids <- bitr(top50_upregulated_in_group_Sea_Level_list,
                                   fromType = "SYMBOL",
                                   toType = "ENTREZID",
                                   OrgDb = "org.Mm.eg.db")


#Run Gene Ontology (GO) enrichment analysis. 
denver_go_enrichment <- enrichGO(gene = top50_denver_entrez_ids$ENTREZID,
                                 OrgDb = "org.Mm.eg.db",
                                 ont = "BP",             # 'BP' for Biological Process
                                 pvalueCutoff = 0.05,
                                 readable = TRUE)

dotplot(denver_go_enrichment, showCategory = 10, title = "Top 10 BP Enriched Pathways in Denver")


sea_level_go_enrichment <- enrichGO(gene = top50_sea_level_entrez_ids$ENTREZID,
                                    OrgDb = "org.Mm.eg.db",
                                    ont = "BP",            
                                    pvalueCutoff = 0.05,
                                    readable = TRUE)

dotplot(sea_level_go_enrichment, showCategory = 10, 
        title = "Top 10 BP Enriched Pathways in Sea Level")

################################################################################
# antigen-presentation related signatures in Myeloid Cells
my_gene_list_1 <- list(
  module1_genes = c("H2-Ea", "H2-Ob", "H2-DMb1", "H2-Eb1", "H2-Aa", "H2-Ab1", 
                    "H2-DMb2", "Cd74", "Cd86"))

# join layers
New_Myeloid <- JoinLayers(New_Myeloid)

# Add module scores
New_Myeloid <- AddModuleScore(
  object = New_Myeloid,
  features = my_gene_list_1,
  name = "ModuleName" # Prefix for the metadata column names
)
view(New_Myeloid@meta.data)


# View the new metadata columns
head(New_Myeloid@meta.data)

# Calculate the Z-score for the module score column
New_Myeloid$Immune_Module_Score_Zscore<- scale(New_Myeloid$ModuleName1)

#Visualization
# Visualize the z-scores using a FeaturePlot (UMAP/tSNE)
plot_list <- FeaturePlot(object = New_Myeloid, 
                         features = "Immune_Module_Score_Zscore", 
                         split.by = "orig.ident",
                         ncol = 2, keep.scale = "feature",
                         #min.cutoff = "q15",
                         #max.cutoff = "q70",
                         combine = FALSE)

modified_plot_list <- lapply(plot_list, function(x) {
  x + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
    labs(title = NULL)
})

# Combine the modified plots using patchwork syntax
combined_plot <- wrap_plots(modified_plot_list) + 
  plot_annotation (title = 'Module Score of Antigen Presentation Signatures')
# Print the final plot
print(combined_plot)

VlnPlot(New_Myeloid, split.by="orig.ident", group.by="celltypes", 
        features = "Immune_Module_Score_Zscore") 

################################################################################
# tolerance-associated signatures in Myeloid cells
my_gene_list_2 <- list(
  module1_genes = c("Cd209a", "Cd209b", "Cd274", "Ido1", "Hmox1", 
                    "Tgfb1", "Tgfb2", "Tgfbr1", "Tgfbr2")
)

New_Myeloid <- AddModuleScore(
  object = New_Myeloid,
  features = my_gene_list_2,
  name = "ModuleName" # Prefix for the metadata column names
)


# View the new metadata columns
head(New_Myeloid@meta.data)

# Calculate the Z-score for the module score column
New_Myeloid$Immune_Module_Score_Zscore<- scale(New_Myeloid$ModuleName1)

#Visualization
# Visualize the z-scores using a FeaturePlot (UMAP/tSNE)
plot_list <- FeaturePlot(object = New_Myeloid, 
                         features = "Immune_Module_Score_Zscore", 
                         split.by = "orig.ident",
                         ncol = 2, keep.scale = "feature",
                         combine = FALSE)

modified_plot_list <- lapply(plot_list, function(x) {
  x + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
    labs(title = NULL)
})

combined_plot <- wrap_plots(modified_plot_list) + 
  plot_annotation (title = 'Module Score of Tolerance Signatures', 
  theme = theme(plot.title = element_text(hjust = 0.5)))
# Print the final plot
print(combined_plot)

VlnPlot(New_Myeloid, split.by="orig.ident", group.by="celltypes", 
        features = "Immune_Module_Score_Zscore") 

################################################################################
# Hif1a-associated signatures in Myeloid cells
my_gene_list_3 <- list(
  module1_genes = c("Hif1a", "Aldoa", "Aldoc", "Cxcr4", "Nrp1", "Eno1", "Gpi1", 
                    "Hk1", "Ldha", "Nfe2l2", "Pdk1", "Pfkp")
)
New_Myeloid <- AddModuleScore(
  object = New_Myeloid,
  features = my_gene_list_3,
  name = "ModuleName" # Prefix for the metadata column names
)


# View the new metadata columns
head(New_Myeloid@meta.data)

# Calculate the Z-score for the module score column
New_Myeloid$Immune_Module_Score_Zscore<- scale(New_Myeloid$ModuleName1)

VlnPlot(New_Myeloid, split.by="orig.ident", group.by="celltypes", 
        features = "Immune_Module_Score_Zscore") 
################################################################################
# Hif2-associated signatures in Myeloid cells
my_gene_list_4 <- list(
  module1_genes = c("Vegfc", "Abl2", "Flt1", "Plin2", "Epo", "Epas1")
)
New_Myeloid <- AddModuleScore(
  object = New_Myeloid,
  features = my_gene_list_4,
  name = "ModuleName" # Prefix for the metadata column names
)


# View the new metadata columns
head(New_Myeloid@meta.data)

# Calculate the Z-score for the module score column
New_Myeloid$Immune_Module_Score_Zscore<- scale(New_Myeloid$ModuleName1)

VlnPlot(New_Myeloid, split.by="orig.ident", group.by="celltypes", 
        features = "Immune_Module_Score_Zscore") 


