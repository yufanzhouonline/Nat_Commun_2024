
library(dplyr)
library(Seurat)
library(patchwork)

readdata <- function(readsdir, projectname, mincells, minfeatures, percentmt){
    print(paste("Reading", readsdir))
    pbmc.data <- Read10X(data.dir = readsdir)
    pbmc <- CreateSeuratObject(counts = pbmc.data, project = projectname, min.cells = mincells, min.features = minfeatures)
    pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
    pbmc <- subset(pbmc, subset = percent.mt < percentmt)
    return(pbmc)
}


###MCF7R1
mcf7r1path <- "/data/yufan/data20210603_scrnaseq_mcf7/Run210602_JA21193_UT_Austin_GSAF_SC/CellRanger_Counts/MCF7R1_Count_outs/filtered_feature_bc_matrix"
###MCF7R2
mcf7r2path <- "/data/yufan/data20210603_scrnaseq_mcf7/Run210602_JA21193_UT_Austin_GSAF_SC/CellRanger_Counts/MCF7R2_Count_outs/filtered_feature_bc_matrix"
###TRR1
trr1path <- "/data/yufan/data20210603_scrnaseq_mcf7/Run210602_JA21193_UT_Austin_GSAF_SC/CellRanger_Counts/TRR1_Count_outs/filtered_feature_bc_matrix"
###TRR2
trr2path <- "/data/yufan/data20210603_scrnaseq_mcf7/Run210602_JA21193_UT_Austin_GSAF_SC/CellRanger_Counts/TRR2_Count_outs/filtered_feature_bc_matrix"
###M1R1
m1r1path <- "/data/yufan/data20210709_scrnaseq_mcf7/Run210708_NB_SC/CellRanger_Counts/M1R1_Count_outs/filtered_feature_bc_matrix"
###M1R2
m1r2path <- "/data/yufan/data20210709_scrnaseq_mcf7/Run210708_NB_SC/CellRanger_Counts/M1R2_Count_outs/filtered_feature_bc_matrix"


mcf7r1 <- readdata(mcf7r1path, 'MCF7_R1', 3, 200, 30)
mcf7r2 <- readdata(mcf7r2path, 'MCF7_R2', 3, 200, 30)
trr1 <- readdata(trr1path, 'MCF7TR_R1', 3, 200, 30)
trr2 <- readdata(trr2path, 'MCF7TR_R2', 3, 200, 30)
m1r1 <- readdata(m1r1path, 'MCF7M1_R1', 3, 200, 30)
m1r2 <- readdata(m1r2path, 'MCF7M1_R2', 3, 200, 30)


mcf7all <- merge(mcf7r1, y = c(mcf7r2, m1r1, m1r2, trr1, trr2), add.cell.ids = c("MCF7_R1", "MCF7_R2", "MCF7M1_R1", "MCF7M1_R2", "MCF7TR_R1", "MCF7TR_R2"), project = "MCF7 Single cells")
mysample <- mcf7all



hist(colSums(mysample$RNA@data),
     breaks = 100,
     main = "Total expression before normalization",
     xlab = "Sum of expression")

mysample <- NormalizeData(mysample, normalization.method = "LogNormalize", scale.factor = 10000)


hist(colSums(mysample$RNA@data),
     breaks = 100,
     main = "Total expression after normalization",
     xlab = "Sum of expression")

mysample <- FindVariableFeatures(mysample, selection.method = "vst", nfeatures = 2000)



top10 <- head(VariableFeatures(mysample), 10) #head(mysample$RNA@var.features,10)
# "PPBP"   "LYZ"    "S100A9" "IGLL5"  "GNLY"   "FTL"    "PF4"    "FTH1"   "GNG11"  "S100A8"

###varible genes with vst.variance
head(mysample@assays$RNA@meta.features)
vargene <- mysample@assays$RNA@meta.features
vargene[top10,]
vargene <- vargene[order(-vargene$vst.variance.standardized),]
vargene <- vargene[which(vargene$vst.variable==TRUE),]
vargene$gene <- rownames(vargene)
setwd("/data/yufan/schic2/scrnaseq")
getwd()
#write.table(vargene[,c("gene", "vst.variance.standardized")],file="vargene.rnk",row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
#write.table(vargene[,c("gene", "vst.variance.standardized")],file="vargenenew.rnk",row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)



plot1 <- VariableFeaturePlot(mysample)

plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2


all.genes <- rownames(mysample)
mysample <- ScaleData(mysample, features = all.genes)


mysample <- ScaleData(mysample, vars.to.regress = "percent.mt")


mysample <- RunPCA(mysample, features = VariableFeatures(object = mysample))



print(mysample[["pca"]], dims = 1:5, nfeatures = 5)



VizDimLoadings(mysample, dims = 1:2, reduction = "pca")



DimPlot(mysample, reduction = "pca",split.by = 'ident')


DimHeatmap(mysample, dims = 1, cells = 500, balanced = TRUE)


DimHeatmap(mysample, dims = 1:15, cells = 500, balanced = TRUE)

mysample <- JackStraw(mysample, num.replicate = 100)
mysample <- ScoreJackStraw(mysample, dims = 1:20)
JackStrawPlot(mysample, dims = 1:15)
ElbowPlot(mysample)


mysample <- FindNeighbors(mysample, dims = 1:10)
#mysample <- FindClusters(mysample, resolution = 0.5)
#mysample <- FindClusters(mysample, resolution = 2)
#mysample <- FindClusters(mysample, resolution = 1)
mysample <- FindClusters(mysample, resolution = 0.75)
#The FindClusters function implements the procedure, 
#and contains a resolution parameter that sets the ‘granularity’ of the downstream clustering, 
#with increased values leading to a greater number of clusters. 
#We find that setting this parameter between 0.6-1.2 typically returns good results
#for single cell datasets of around 3K cells. Optimal resolution often increases for larger datasets.
#The clusters are saved in the object@ident slot.


head(Idents(mysample), 5)

# Save old identity classes (the cluster labels) for reference.
mysample[["old.ident"]] <- Idents(object = mysample)

# Rename classes.
mysample <- RenameIdents(object = mysample, `0` = "D1", `1` = "D2", `2` = "D3", `3` = "D4", `4` = "D5", `5` = "D6", `6` = "D7", `7` = "D8", `8` = "D9", `9` = "D10", `10` = "D11", `11` = "D12", `12` = "D13")


# install UMAP： reticulate::py_install(packages ='umap-learn')
mysample <- RunUMAP(mysample, dims = 1:10)


DimPlot(mysample, reduction = "umap")


DimPlot(mysample, reduction = "umap",label = TRUE)
LabelClusters(DimPlot(mysample, reduction = "umap"),id = 'ident')

###TSNE
mysample <- RunTSNE(mysample, dims = 1:10)
DimPlot(mysample, reduction = "tsne")
DimPlot(mysample, reduction = "tsne",label = TRUE)
LabelClusters(DimPlot(mysample, reduction = "tsne"),id = 'ident')

###Get the tsne coordinates
tsne1 <- mysample@reductions$tsne@cell.embeddings[,1]
tsne2 <- mysample@reductions$tsne@cell.embeddings[,2]

###Pull number of cells in cluster from seurat object
cluster <- mysample@meta.data
table(cluster$seurat_clusters)
#write.table(table(cluster$seurat_clusters),file="clusternumber.txt",row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

###get the tsne coordinates to cluster
cluster$tsne1 = tsne1
cluster$tsne2 = tsne2

clusterlist <- min(as.numeric(levels(cluster$seurat_clusters))):max(as.numeric(levels(cluster$seurat_clusters)))
for (i in clusterlist)
{
    cat(i)
    cat("\t")
    cat(mean(cluster[which(cluster$seurat_clusters==i),]$tsne1))
    cat("\t")
    cat(mean(cluster[which(cluster$seurat_clusters==i),]$tsne2))
    cat("\n")
}

euclidean <- function(a, b) sqrt(sum((a - b)^2))
apply(celltype, 2, function(x) sum(x))
celltypep <- apply(celltype, 2, function(x) x/sum(x)*100)


######PCA
mysample <- RunPCA(mysample, dims = 1:10)
DimPlot(mysample, reduction = "pca")
DimPlot(mysample, reduction = "pca",label = TRUE)
LabelClusters(DimPlot(mysample, reduction = "pca"),id = 'ident')


#saveRDS(mysample, file = "../output/mysample_tutorial.rds")

# find all markers of cluster 1
cluster1.markers <- FindMarkers(mysample, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(mysample, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
mysample.markers <- FindAllMarkers(mysample, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
mysample.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
mysample.markers$gene <- rownames(mysample.markers)


write.table(mysample.markers, file="clusterDEGsnewp05.txt",row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)


nrow(mysample.markers)

degp05 <- read.csv("clusterDEGsnewp05.txt", header=TRUE, sep="\t")
degp10 <- read.csv("clusterDEGsnewp10.txt", header=TRUE, sep="\t")
degp25 <- read.csv("clusterDEGsnewp25.txt", header=TRUE, sep="\t")
degp50 <- read.csv("clusterDEGsnewp50.txt", header=TRUE, sep="\t")
degp75 <- read.csv("clusterDEGsnewp75.txt", header=TRUE, sep="\t")
degp100 <- read.csv("clusterDEGsnew.txt", header=TRUE, sep="\t")


nrow(degp75[which(degp75[which(degp75$cluster == 5),]$gene %in% degp100[which(degp100$cluster == 4),]$gene),])

overlapcount <- data.frame("SN" = 1:2, "Age" = c(21,15), "Name" = c("John","Dora"))

count <- c()
for (i in 0:12)
{
    count = c(count, nrow(degp100[which(degp100$cluster == i),]))
}


# find markers for every cluster compared to all remaining cells, report only the positive ones
mysample.markers <- FindAllMarkers(mysample, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0)
mysample.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
mysample.markers$gene <- rownames(mysample.markers)
write.table(mysample.markers, file="/data/yufan/schic2/scrnaseq/nodegonlypos.txt",row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

mysample.markers <- FindAllMarkers(mysample, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0)
mysample.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
mysample.markers$gene <- rownames(mysample.markers)
write.table(mysample.markers, file="/data/yufan/schic2/scrnaseq/nodegposneg.txt",row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

### 75% cells
overlap75 <- data.frame("D1"=c(), "D2"=c(), "D3"=c(), "D4"=c(), "D5"=c(), "D6"=c(), "D7"=c(), "D8"=c(), "D9"=c(), "D10"=c(), "D11"=c(), "D12"=c(), "D13"=c())
#overlap75 <- data.frame()

for (i in 0:12)
{
    overlaprow <- c()
    for (j in 0:12)
    {
        overlapcount = nrow(degp75[which(degp75[which(degp75$cluster == i),]$gene %in% degp100[which(degp100$cluster == j),]$gene),])
        overlaprow <- c(overlaprow, overlapcount)
    }
    overlap75 <- rbind(overlap75, overlaprow / count * 100)
}

colnames(overlap75) <- c("D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10", "D11", "D12", "D13")
rownames(overlap75) <- colnames(overlap75)

library("pheatmap")
par(mar = c(5, 10, 4, 2) + 0.1)

pheatmap(overlap75, cluster_rows=FALSE, cluster_cols=FALSE, color = colorRampPalette(c("#0088FF", "white", "#FF5500"))(100), fontsize=20, fontsize_row=20, fontsize_col=20)

### 50% cells
overlap75 <- data.frame("D1"=c(), "D2"=c(), "D3"=c(), "D4"=c(), "D5"=c(), "D6"=c(), "D7"=c(), "D8"=c(), "D9"=c(), "D10"=c(), "D11"=c(), "D12"=c(), "D13"=c())
#overlap75 <- data.frame()

for (i in 0:12)
{
    overlaprow <- c()
    for (j in 0:12)
    {
        overlapcount = nrow(degp50[which(degp50[which(degp50$cluster == i),]$gene %in% degp100[which(degp100$cluster == j),]$gene),])
        overlaprow <- c(overlaprow, overlapcount)
    }
    overlap75 <- rbind(overlap75, overlaprow / count * 100)
}

colnames(overlap75) <- c("D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10", "D11", "D12", "D13")
rownames(overlap75) <- colnames(overlap75)

library("pheatmap")
par(mar = c(5, 10, 4, 2) + 0.1)

pheatmap(overlap75, cluster_rows=FALSE, cluster_cols=FALSE, color = colorRampPalette(c("#0088FF", "white", "#FF5500"))(100), fontsize=20, fontsize_row=20, fontsize_col=20)

### 25% cells
overlap75 <- data.frame("D1"=c(), "D2"=c(), "D3"=c(), "D4"=c(), "D5"=c(), "D6"=c(), "D7"=c(), "D8"=c(), "D9"=c(), "D10"=c(), "D11"=c(), "D12"=c(), "D13"=c())
#overlap75 <- data.frame()

for (i in 0:12)
{
    overlaprow <- c()
    for (j in 0:12)
    {
        overlapcount = nrow(degp25[which(degp25[which(degp25$cluster == i),]$gene %in% degp100[which(degp100$cluster == j),]$gene),])
        overlaprow <- c(overlaprow, overlapcount)
    }
    overlap75 <- rbind(overlap75, overlaprow / count * 100)
}

colnames(overlap75) <- c("D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10", "D11", "D12", "D13")
rownames(overlap75) <- colnames(overlap75)

library("pheatmap")
par(mar = c(5, 10, 4, 2) + 0.1)

pheatmap(overlap75, cluster_rows=FALSE, cluster_cols=FALSE, color = colorRampPalette(c("#0088FF", "white", "#FF5500"))(100), fontsize=20, fontsize_row=20, fontsize_col=20)

### 10% cells
overlap75 <- data.frame("D1"=c(), "D2"=c(), "D3"=c(), "D4"=c(), "D5"=c(), "D6"=c(), "D7"=c(), "D8"=c(), "D9"=c(), "D10"=c(), "D11"=c(), "D12"=c(), "D13"=c())
#overlap75 <- data.frame()

for (i in 0:12)
{
    overlaprow <- c()
    for (j in 0:12)
    {
        overlapcount = nrow(degp10[which(degp10[which(degp10$cluster == i),]$gene %in% degp100[which(degp100$cluster == j),]$gene),])
        overlaprow <- c(overlaprow, overlapcount)
    }
    overlap75 <- rbind(overlap75, overlaprow / count * 100)
}

colnames(overlap75) <- c("D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10", "D11", "D12", "D13")
rownames(overlap75) <- colnames(overlap75)

library("pheatmap")
par(mar = c(5, 10, 4, 2) + 0.1)

pheatmap(overlap75, cluster_rows=FALSE, cluster_cols=FALSE, color = colorRampPalette(c("#0088FF", "white", "#FF5500"))(100), fontsize=20, fontsize_row=20, fontsize_col=20)


cluster1.markers <- FindMarkers(mysample, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster1.markers, n = 5)

