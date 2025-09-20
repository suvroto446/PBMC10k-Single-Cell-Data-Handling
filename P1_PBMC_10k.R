
install.packages("hdf5r")
library(hdf5r)
library(dplyr)
library(Seurat)
library(patchwork)

data_h5 <- "C:/Users/SUVROTO/Documents/R scripts/Dataset/PBMC_10k_raw_feature_bc_matrix.h5"
pbmc10k.data <- Read10X_h5(filename = data_h5)

names(pbmc10k.data)

#separate the gene expression data
rna_counts <- pbmc10k.data[["Gene Expression"]]

# Initialize the Seurat object with the raw (non-normalized data).
pbmc10k <- CreateSeuratObject(counts = rna_counts, project = "pbmc10k", min.cells = 3, min.features = 200)

# Lets examine a few genes in the first thirty cells
rna_counts[c("CD3D", "CD19", "MS4A1"), 1:100]

# To add a new column with percentage of mitochondrial genes
pbmc10k[["percent.mt"]] <- PercentageFeatureSet(pbmc10k, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc10k, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, the number at the
#top shows pearson correlation
plot1 <- FeatureScatter(pbmc10k, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc10k, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter out low-quality cells (e.g., high mitochondrial %)
pbmc10k_sub <- subset(pbmc10k, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Global-scaling normalization method “LogNormalize”
pbmc10k_sub <- NormalizeData(pbmc10k_sub, normalization.method = "LogNormalize", scale.factor = 10000)
#The same outcome can be gained by pbmc10k_sub <- NormalizeData(pbmc10k_sub)

#Finding features(genes) that exhibit high cell-to-cell variation in the data set
pbmc10k_sub <- FindVariableFeatures(pbmc10k_sub, selection.method = "vst", nfeatures = 2000)
#2000 features will be returned

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc10k_sub), 10)

# plot variable features with and without labels
plot3 <- VariableFeaturePlot(pbmc10k_sub)
plot4 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) #repel is true so that the labels do not overlap
plot3 + plot4

#Scaling the data
all.genes <- rownames(pbmc10k_sub)
pbmc10k_sub <- ScaleData(pbmc10k_sub, features = all.genes)

##Performing linear dimension reduction
#PCA
pbmc10k_sub <- RunPCA(pbmc10k_sub, features = VariableFeatures(object = pbmc10k_sub))

# Examine and visualize PCA results a few different ways
print(pbmc10k_sub[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc10k_sub, dims = 1:2, reduction = "pca")
DimPlot(pbmc10k_sub, reduction = "pca") + NoLegend()
DimHeatmap(pbmc10k_sub, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc10k_sub, dims = 1:9, cells = 500, balanced = TRUE)

#Determine the ‘dimensionality’ of the dataset
ElbowPlot(pbmc10k_sub) #To visualize which PC's has the most variance

#Clustering the cells
pbmc10k_sub <- FindNeighbors(pbmc10k_sub, dims = 1:12)
pbmc10k_sub <- FindClusters(pbmc10k_sub, resolution = 0.5)

Idents(pbmc10k_sub) #To check the distribution of cells in 12 clusters (0-11)
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc10k_sub), 5)

#Run non-linear dimensional reduction
pbmc10k_sub <- RunUMAP(pbmc10k_sub, dims = 1:12)

#Visualizing the graph
DimPlot(pbmc10k_sub, reduction = "umap", label = TRUE)

#Finding differentially expressed features (cluster biomarkers)
# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc10k_sub, ident.1 = 2)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc10k_sub, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(pbmc10k_sub, only.pos = TRUE)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

VlnPlot(pbmc10k_sub, features = c("CD79A", "CD19"))
FeaturePlot(pbmc10k_sub, features = c("MS4A1", "CD19", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))

###Identifying the clusters
VlnPlot(pbmc10k_sub, features = c("CD3E","CD4","CD8A","NKG7","GNLY","LYZ","CD14","MS4A1","CD79A","CD19","PPBP","S100A4"))

#Cell markers used here are canonical here
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "Naive CD8 T","B cell","CD14+ Mono",
                     "NK","CD14+ Mono","Cytotoxic CD8 T","Cytotoxic CD8 T", "B cell", "Platelets")
names(new.cluster.ids) <- levels(pbmc10k_sub)
pbmc10k_sub <- RenameIdents(pbmc10k_sub, new.cluster.ids)
DimPlot(pbmc10k_sub, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

#0 = CD4 T cell (Naive)
#1 = CD14⁺ classical monocytes
#2 = CD4 T cell (Memory)
#3 = Naive CD8
#4 = B cell
#5 = CD14+ monocytes
#6 = NK
#7 = CD14+ monocytes
#8 = Cytotoxic CD8
#9 = Cytotoxic CD8
#10 = B cell
#11 = Platelets




