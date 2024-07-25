###single cell rnaseq data for WTC11 cells and MCF7 cells

library(Seurat)
library(dplyr)
library(Matrix)


readdata <- function(readsdir, projectname, mincells, minfeatures, percentmt){
    print(paste("Reading", readsdir))
    #mysample.data <- Read10X(data.dir = readsdir)
    mysample.data <- read.table(file=readsdir,sep="\t")
    mysample <- CreateSeuratObject(counts = mysample.data, project = projectname, min.cells = mincells, min.features = minfeatures)
    mysample[["percent.mt"]] <- PercentageFeatureSet(mysample, pattern = "^MT-")
    mysample <- subset(mysample, subset = percent.mt < percentmt)
    return(mysample)
}

readfolder <- function(readsdir, projectname, mincells, minfeatures, percentmt){
    print(paste("Reading", readsdir))
    mysample.data <- Read10X(data.dir = readsdir)
    #mysample.data <- read.table(file=readsdir,sep="\t")
    mysample <- CreateSeuratObject(counts = mysample.data, project = projectname, min.cells = mincells, min.features = minfeatures)
    mysample[["percent.mt"]] <- PercentageFeatureSet(mysample, pattern = "^MT-")
    mysample <- subset(mysample, subset = percent.mt < percentmt)
    return(mysample)
}

#read scRNA-seq data
###wtcd0r1
wtcd0r1path <- "/data/yufan/schic2/wtccells/scrnaseq/fastq/Day0Rep1_ReadMappedCount.txt"
###wtcd0r2
wtcd0r2path <- "/data/yufan/schic2/wtccells/scrnaseq/fastq/Day0Rep2_ReadMappedCount.txt"
###wtcd2r1
wtcd2r1path <- "/data/yufan/schic2/wtccells/scrnaseq/fastq/Day2Rep1_ReadMappedCount.txt"
###wtcd2r2
wtcd2r2path <- "/data/yufan/schic2/wtccells/scrnaseq/fastq/Day2Rep2_ReadMappedCount.txt"
###wtcd5r1
wtcd5r1path <- "/data/yufan/schic2/wtccells/scrnaseq/fastq/Day5Rep1_ReadMappedCount.txt"
###wtcd5r2
wtcd5r2path <- "/data/yufan/schic2/wtccells/scrnaseq/fastq/Day5Rep2_ReadMappedCount.txt"
###wtcd15r1
wtcd15r1path <- "/data/yufan/schic2/wtccells/scrnaseq/fastq/Day15Rep1_ReadMappedCount.txt"
###wtcd15r2
wtcd15r2path <- "/data/yufan/schic2/wtccells/scrnaseq/fastq/Day15Rep2_ReadMappedCount.txt"
###wtcd30r1
wtcd30r1path <- "/data/yufan/schic2/wtccells/scrnaseq/fastq/Day30Rep1_ReadMappedCount.txt"
###wtcd30r2
wtcd30r2path <- "/data/yufan/schic2/wtccells/scrnaseq/fastq/Day30Rep2_ReadMappedCount.txt"


wtcd0r1 <- readdata(wtcd0r1path, 'WTC_R1', 3, 200, 30)
wtcd0r2 <- readdata(wtcd0r2path, 'WTC_R2', 3, 200, 30)
#wtcd2r1 <- readdata(wtcd2r1path, 'WTC_Day2_R1', 3, 200, 30)
#wtcd2r2 <- readdata(wtcd2r2path, 'WTC_Day2_R2', 3, 200, 30)
#wtcd5r1 <- readdata(wtcd5r1path, 'WTC_Day5_R1', 3, 200, 30)
#wtcd5r2 <- readdata(wtcd5r2path, 'WTC_Day5_R2', 3, 200, 30)
#wtcd15r1 <- readdata(wtcd15r1path, 'WTC_Day15_R1', 3, 200, 30)
#wtcd15r2 <- readdata(wtcd15r2path, 'WTC_Day15_R2', 3, 200, 30)
#wtcd30r1 <- readdata(wtcd30r1path, 'WTC_Day30_R1', 3, 200, 30)
#wtcd30r2 <- readdata(wtcd30r2path, 'WTC_Day30_R2', 3, 200, 30)


###MCF7R1
mcf7r1path <- "/data/yufan/data20210603_scrnaseq_mcf7/Run210602_JA21193_UT_Austin_GSAF_SC/CellRanger_Counts/MCF7R1_Count_outs/filtered_feature_bc_matrix"
###MCF7R2
mcf7r2path <- "/data/yufan/data20210603_scrnaseq_mcf7/Run210602_JA21193_UT_Austin_GSAF_SC/CellRanger_Counts/MCF7R2_Count_outs/filtered_feature_bc_matrix"

mcf7r1 <- readfolder(mcf7r1path, 'MCF7_R1', 3, 200, 30)
mcf7r2 <- readfolder(mcf7r2path, 'MCF7_R2', 3, 200, 30)


#mcf7all <- merge(wtcd0r1, y = c(wtcd0r2, wtcd2r1, wtcd2r2, wtcd5r1, wtcd5r2, wtcd15r1, wtcd15r2, wtcd30r1, wtcd30r2, mcf7r1, mcf7r2), add.cell.ids = c("WTC_Day0_R1", "WTC_Day0_R2", "WTC_Day2_R1", "WTC_Day2_R2", "WTC_Day5_R1", "WTC_Day5_R2", "WTC_Day15_R1", "WTC_Day15_R2", "WTC_Day30_R1", "WTC_Day30_R2", "MCF7_R1", "MCF7_R2"), project = "WTC/MCF7 Single cells")
mcf7all <- merge(wtcd0r1, y = c(wtcd0r2, mcf7r1, mcf7r2), add.cell.ids = c("WTC_R1", "WTC_R2", "MCF7_R1", "MCF7_R2"), project = "WTC/MCF7 Single cells")

mysample <- mcf7all


#total expression per cell before normalization
hist(colSums(mysample$RNA@data),
     breaks = 100,
     main = "Total expression before normalization",
     xlab = "Sum of expression")

mysample <- NormalizeData(mysample, normalization.method = "LogNormalize", scale.factor = 10000)

#total expression per cell after normalization
hist(colSums(mysample$RNA@data),
     breaks = 100,
     main = "Total expression after normalization",
     xlab = "Sum of expression")

mysample <- FindVariableFeatures(mysample, selection.method = "vst", nfeatures = 2000)


# top10 highly variable genes
top10 <- head(VariableFeatures(mysample), 10) #head(mysample$RNA@var.features,10)
# "PPBP"   "LYZ"    "S100A9" "IGLL5"  "GNLY"   "FTL"    "PF4"    "FTH1"   "GNG11"  "S100A8"

###varible genes with vst.variance
head(mysample@assays$RNA@meta.features)
vargene <- mysample@assays$RNA@meta.features
vargene[top10,]
vargene <- vargene[order(-vargene$vst.variance.standardized),]
vargene <- vargene[which(vargene$vst.variable==TRUE),]
vargene$gene <- rownames(vargene)
#setwd("/data/yufan/schic2/scrnaseq")
#getwd()
#write.table(vargene[,c("gene", "vst.variance.standardized")],file="vargene.rnk",row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
write.table(vargene[,c("gene", "vst.variance.standardized")],file="/data/yufan/schic2/wtccells/varGeneNew.rnk",row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)


# plot genes and observe the distribution
plot1 <- VariableFeaturePlot(mysample)
# plot variable genes and label top 10
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2


### data scale to average at 0 and standard deviation at 1
### results store at mysample[["RNA"]]@scale.data

all.genes <- rownames(mysample)
mysample <- ScaleData(mysample, features = all.genes)

#move factors affecting standard deviation
mysample <- ScaleData(mysample, vars.to.regress = "percent.mt")

#PCA analysis

mysample <- RunPCA(mysample, features = VariableFeatures(object = mysample))

#print the major components 

print(mysample[["pca"]], dims = 1:5, nfeatures = 5)

#visualize the major components 

VizDimLoadings(mysample, dims = 1:2, reduction = "pca")

#plot two major components

DimPlot(mysample, reduction = "pca",split.by = 'ident')

#DimHeatmap draws heatmap of major components

DimHeatmap(mysample, dims = 1, cells = 500, balanced = TRUE)

#DimHeatmap draws heatmap of multiple major components
DimHeatmap(mysample, dims = 1:15, cells = 500, balanced = TRUE)

#data dimension
mysample <- JackStraw(mysample, num.replicate = 100)
mysample <- ScoreJackStraw(mysample, dims = 1:20)
JackStrawPlot(mysample, dims = 1:15)
ElbowPlot(mysample)

#cell cluster
mysample <- FindNeighbors(mysample, dims = 1:10)
mysample <- FindClusters(mysample, resolution = 0.5)
#mysample <- FindClusters(mysample, resolution = 2)
#mysample <- FindClusters(mysample, resolution = 1)
#mysample <- FindClusters(mysample, resolution = 0.75)
#The FindClusters function implements the procedure, 
#and contains a resolution parameter that sets the ‘granularity’ of the downstream clustering, 
#with increased values leading to a greater number of clusters. 
#We find that setting this parameter between 0.6-1.2 typically returns good results
#for single cell datasets of around 3K cells. Optimal resolution often increases for larger datasets.
#The clusters are saved in the object@ident slot.

#which cluster are cells in
head(Idents(mysample), 5)

# Save old identity classes (the cluster labels) for reference.
mysample[["old.ident"]] <- Idents(object = mysample)

# Rename classes.
#mysample <- RenameIdents(object = mysample, `0` = "D1", `1` = "D2", `2` = "D3", `3` = "D4", `4` = "D5", `5` = "D6", `6` = "D7", `7` = "D8", `8` = "D9", `9` = "D10", `10` = "D11", `11` = "D12", `12` = "D13")
#mysample <- RenameIdents(object = mysample, `0` = "D1", `1` = "D2", `2` = "D3", `3` = "D4", `4` = "D5", `5` = "D6", `6` = "D7", `7` = "D8", `8` = "D9", `9` = "D10", `10` = "D11", `11` = "D12", `12` = "D13", `13` = "D14")
mysample <- RenameIdents(object = mysample, `0` = "DD1", `1` = "DD2", `2` = "DD3", `3` = "DD4", `4` = "DD5", `5` = "DD6", `6` = "DD7", `7` = "DD8", `8` = "DD9", `9` = "DD10")

#dimensionality reduction
# install UMAP： reticulate::py_install(packages ='umap-learn')
mysample <- RunUMAP(mysample, dims = 1:10)

#plot
DimPlot(mysample, reduction = "umap")

#add label of cell cluster
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
write.table(table(cluster$seurat_clusters),file="/data/yufan/schic2/wtccells/clusterNumber.txt",row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

table(cluster$orig.ident)

celltype <- table(cluster$orig.ident, cluster$seurat_clusters)
print(celltype)
#write.table(celltype,file="clustercelltype.txt",row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
write.table(celltype,file="/data/yufan/schic2/wtccells/clusterCellType.txt",row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
