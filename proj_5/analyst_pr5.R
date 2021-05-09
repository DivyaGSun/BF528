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

# Clustering
# From heatmap, we know clusters (1,2), (3,9), and (6,8) are likely to overlap
panglao <- read_tsv("https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz") %>%  
  filter(str_detect(species,"Hs")) %>% 
  select(c("official gene symbol", "cell type", "organ"))


# Delta (SST) = 0
# Alpha (GCG) = 1, 2
# Acinar (CPA1) - 3, 9
# Ductal (KRT19) - 4
# Stellate (PDGFRB) - 5
# Beta(INS) = 6, 8
# Gamma (PPY) = 7
# Vascular (VWF, PECAM1, CD34) - 10
# Macrophage (CD163, CD58, IgG) - 11
# Epsilon (GHRL) - None Found
# Cytotoxic T (CD3, CD8, TRAC) - None Found
# Mast (TPSAB1, KIT, CPA3) - None Found

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