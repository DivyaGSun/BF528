# Analyst 
# Load Packages
library(tidyverse)
library(Seurat)
set.seed(20)

# Load Data
RDS <- readRDS("/projectnb2/bf528/users/lava_lamp/project_4/pbmc_program.rds")

# Cluster biomarkers
cells_markers <- FindAllMarkers(RDS, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top_markers <- cells_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

# Marker genes used in original analysis
marker_genes <- read_csv("marker_genes.csv")
allgenes <- unlist(strsplit(marker_genes$Genes, ","))

# Visualize Clusters, by Cluster and by Biomarker
png("umap.png")
DimPlot(RDS, reduction = "umap")
dev.off()
png("features.png", width = 1000, height = 1000)
FeaturePlot(RDS, features = allgenes)
dev.off()

# Heatmap
png("heatmap_clusters.png", height = 800, width = 1000)
DoHeatmap(RDS, features = top_markers$gene) + NoLegend()
dev.off()

# Violin Plots Using Features Provided in Paper
pdf("violins.pdf")
VlnPlot(RDS, features = c("GCG")) # Alpha
VlnPlot(RDS, features = c("INS")) # Beta
VlnPlot(RDS, features = c("SST")) # Delta
VlnPlot(RDS, features = c("KRT19")) # Ductal 
VlnPlot(RDS, features = c("PPY")) # Gamma
VlnPlot(RDS, features = c("CPA1")) #Acinar
VlnPlot(RDS, features = c("PDGFRB")) # Stellate
VlnPlot(RDS, features = c("VWF", "PECAM1","CD34")) # Vascular
VlnPlot(RDS, features = c("SDS", "CD163")) # Macrophage
dev.off()

# Assign Cell Type to Clusters
new_cluster_ids <- c("Delta", "Alpha", "Alpha", "Acinar", "Ductal", "Stellate", "Beta", "Gamma", "Beta", "Acinar", "Vascular", "Macrophage")
names(new_cluster_ids) <- levels(RDS)
genes <- RenameIdents(RDS, new_cluster_ids)

# Clustered UMAP
png("umap_celltype.png")
DimPlot(genes, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

# Clustered Heatmap
png("heatmap_celltype.png", height = 800, width = 1000)
DoHeatmap(genes, features = top_markers$gene) + NoLegend()
dev.off()

# Save Seurat Object
saveRDS(genes, file = "analyst_output.rds")

# Save Marker Genes
write_csv(cells_markers, "Marker_Genes_Found.csv")
