## --------------------------------------------------------------------------------------------------------------------------------------------------
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)

data <- Read10X(data.dir = "C:/Users/willd/OneDrive/Desktop/scRNA_seq/experiment/outs/count/filtered_feature_bc_matrix/")

mesc <- CreateSeuratObject(counts = data, project = "mESC_diff", min.cells = 3, min.features = 200)



## --------------------------------------------------------------------------------------------------------------------------------------------------
mesc[["percent.mt"]] <- PercentageFeatureSet(mesc, pattern = "^mt-")

#Show QC metrics for the first 5 cells
head(mesc@meta.data, 5)


## --------------------------------------------------------------------------------------------------------------------------------------------------
# Visualize QC metrics as a violin plot
VlnPlot(mesc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


## ----fig.height=8, fig.width=12--------------------------------------------------------------------------------------------------------------------
# We can observe the relationship of the features and read counts
plot1 <- FeatureScatter(mesc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2 <- FeatureScatter(mesc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot1 + plot2


## --------------------------------------------------------------------------------------------------------------------------------------------------
# keep cells with features between 200-3500.
mesc <- subset(mesc, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 5)


## --------------------------------------------------------------------------------------------------------------------------------------------------
mesc <- SCTransform(mesc, vars.to.regress = "percent.mt")


## --------------------------------------------------------------------------------------------------------------------------------------------------
mesc <- RunPCA(mesc, features = VariableFeatures(object = mesc))


## --------------------------------------------------------------------------------------------------------------------------------------------------
mesc <- FindNeighbors(mesc, dims = 1:30)
mesc <- FindClusters(mesc, resolution = 0.5)


## ----echo=FALSE------------------------------------------------------------------------------------------------------------------------------------
# load reticulate to interface with python and make sure you download umap-learn into the python environment. If you haven't downloaded it run the below conda_install line. You need this to run the UMAP function.

library(reticulate)
use_virtualenv("bioinfo")
#conda_install(envname='bioinfo', packages='umap-learn')


## ----message=FALSE---------------------------------------------------------------------------------------------------------------------------------
mesc <- RunUMAP(mesc, dims = 1:30, umap.method = 'umap-learn', metric = 'correlation' )


## ----echo=TRUE-------------------------------------------------------------------------------------------------------------------------------------
DimPlot(mesc, reduction = "umap", label = TRUE)


## ----echo=FALSE, message=FALSE---------------------------------------------------------------------------------------------------------------------
# find all markers of cluster 1, this compares the cells in cluster 1 versus all other cells.
cluster1.markers <- FindMarkers(mesc, ident.1 = 1, min.pct = 0.25, test.use = "MAST")
cluster0.markers <- FindMarkers(mesc, ident.1 = 0, min.pct = 0.25, test.use = "MAST")
cluster2.markers <- FindMarkers(mesc, ident.1 = 2, min.pct = 0.25, test.use = "MAST")
cluster3.markers <- FindMarkers(mesc, ident.1 = 3, min.pct = 0.25, test.use = "MAST")
cluster4.markers <- FindMarkers(mesc, ident.1 = 4, min.pct = 0.25, test.use = "MAST")
cluster5.markers <- FindMarkers(mesc, ident.1 = 5, min.pct = 0.25, test.use = "MAST")
cluster6.markers <- FindMarkers(mesc, ident.1 = 6, min.pct = 0.25, test.use = "MAST")
cluster7.markers <- FindMarkers(mesc, ident.1 = 7, min.pct = 0.25, test.use = "MAST")

head(arrange(cluster4.markers, desc(avg_logFC)), n=15)


## ----echo=FALSE------------------------------------------------------------------------------------------------------------------------------------
# compare between cluster 4 and 9
cluster4_9.markers <- FindMarkers(mesc, ident.1 = 4, ident.2 = 9, min.pct = 0.25)
head(arrange(cluster4_9.markers, desc(avg_logFC)), n=20)


## ----fig.height=8, fig.width=8---------------------------------------------------------------------------------------------------------------------
VlnPlot(mesc, features = c("Nanog", "Klf4", "Tbx3", "Sox1", "Dnmt3b", "Otx2"), 
    pt.size = 0.2, ncol = 3)


## ----fig.height=8, fig.width=8---------------------------------------------------------------------------------------------------------------------
# Visualize canonical marker genes on the sctransform embedding.
FeaturePlot(mesc, features = c("Nanog", "Klf4", "Tbx3", "Dnmt3b"), pt.size = 0.2, 
    ncol =2)


## ----fig.height=20, fig.width=8, message=FALSE, warning=FALSE--------------------------------------------------------------------------------------
# find markers for every cluster compared to all remaining cells, report only the positive ones
mesc.markers <- FindAllMarkers(mesc, test.use = "MAST")
top20 <- mesc.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
DoHeatmap(mesc, features = top20$gene) #+ NoLegend()


## --------------------------------------------------------------------------------------------------------------------------------------------------
new.cluster.ids <- c("0","1", "2", "Early Epiblast-like", "mESC", "5", "Early Epiblast-like", "mESC", "8", "Epiblast-like", "10", "11", "12")
names(new.cluster.ids) <- levels(mesc)
mesc <- RenameIdents(mesc, new.cluster.ids)
DimPlot(mesc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

