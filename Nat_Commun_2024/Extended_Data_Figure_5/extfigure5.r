###SCCNV

chrlist = paste0("chr", 1:22)
sccnv <- read.csv("result.sccnv.txt", header=TRUE, sep="\t")

sccnvnoxy <- sccnv[which((sccnv$chr != "chrX") & (sccnv$chr != "chrY")),]

sccnvnoxy[is.na(sccnvnoxy)] <- 2

sccnvnoxy[,8:112]

library("umap")

umap.defaults
custom.settings = umap.defaults
custom.settings$n_components = 2

preumap <- t(sccnvnoxy[,8:112])

sccnv.umap <- umap(preumap, config = custom.settings)

plot(sccnv.umap$layout, col = rep(1:3, c(33, 38, 34)), pch = 19, xlab="UMAP1", ylab="UMAP2")
legend("bottomright", legend = c("MCF7", "MCF7TR", "MCF7M1"), col = 1:3, pch = 19, bty = "n")

#################################################################################
#################################################################################

chrlist = paste0("chr", 1:22)
sccnv <- read.csv("result.sccnv.txt", header=TRUE, sep="\t")

sccnvnoxy <- sccnv[which((sccnv$chr != "chrX") & (sccnv$chr != "chrY")),]

sccnvnoxy[is.na(sccnvnoxy)] <- 2
sccnvnoxy <- sccnvnoxy[complete.cases(sccnvnoxy),]

### hit cnv in any cells
hitcnv <- data.frame()
for (i in 1:nrow(sccnvnoxy))
{
    for (j in 8:112)
    {
        if (sccnvnoxy[i, j] !=2 )
        {
            print(sccnvnoxy[i, 1:3])
            hitcnv <- rbind(hitcnv, sccnvnoxy[i, 1:3])
            break
        }
    }
}

rownames(hitcnv) <- 1:nrow(hitcnv)

### hit cnv in 10% of cells
hitcnv <- data.frame()
for (i in 1:nrow(sccnvnoxy))
{
    rowsum <- 0
    for (j in 8:112)
    {
        if (sccnvnoxy[i, j] !=2 )
        {
            rowsum <- rowsum + 1
        }
    }
    if (rowsum > 11)
    {
        print(i)
        hitcnv <- rbind(hitcnv, sccnvnoxy[i, 1:3])
    }
}

rownames(hitcnv) <- 1:nrow(hitcnv)

###
contact <- read.csv("data.txt", header=TRUE, sep="\t")

cell <- read.csv("labelname.txt", header=FALSE, sep="\t")

colnames(cell) <- c("celltype")
cell$cellid <- 0:292

### all cell
hitindex <- c()
for (i in 1:nrow(contact))
{
    overlap <- hitcnv[which(((contact$chrom1[i] == hitcnv$chr) & (contact$pos1[i] >= hitcnv$pos1) & (contact$pos1[i] <= hitcnv$pos2)) | ((contact$chrom2[i] == hitcnv$chr) & (contact$pos2[i] >= hitcnv$pos1) & (contact$pos2[i] <= hitcnv$pos2))),]
    if (nrow(overlap) > 0)
    {
        print(i)
        #print(overlap)
        hitindex <- c(hitindex, i)
    }
}

### one cell
contacta001 <- contact[which(contact$cell_name=="A001"),]

hitindex <- c()
for (i in 1:nrow(contacta001))
{
    overlap <- hitcnv[which(((contacta001$chrom1[i] == hitcnv$chr) & (contacta001$pos1[i] >= hitcnv$pos1) & (contacta001$pos1[i] <= hitcnv$pos2)) | ((contacta001$chrom2[i] == hitcnv$chr) & (contacta001$pos2[i] >= hitcnv$pos1) & (contacta001$pos2[i] <= hitcnv$pos2))),]
    if (nrow(overlap) > 0)
    {
        print(i)
        #print(overlap)
        hitindex <- c(hitindex, i)
    }
}

### change overlap from hitcnv
hitdf <- data.frame()
for (i in 1:nrow(hitcnv))
{
    overlap <- contact[which(((contact$chrom1 == hitcnv$chr[1]) & (contact$pos1 >= hitcnv$pos1[1]) & (contact$pos1 <= hitcnv$pos2[1])) | ((contact$chrom2 == hitcnv$chr[1]) & (contact$pos2 >= hitcnv$pos1[1]) & (contact$pos2 <= hitcnv$pos2[1]))),]
    print(paste(i, "/", nrow(hitcnv)))
    hitdf <- rbind(hitdf, overlap)
}

### MCF7
mcf7id <- cell[which(cell$celltype=="MCF7"),]$cellid

mcf7ratio <- c()
for (i in mcf7id)
{
    print(i)
    ratio <- nrow(hitdf[which(hitdf$cell_id == mcf7id[i]),]) / nrow(contact[which(contact$cell_id == mcf7id[i]),])
    mcf7ratio <- c(mcf7ratio, ratio)
}

ratiolist1 <- mcf7ratio[complete.cases(mcf7ratio)]

### MCF7M1
mcf7id <- cell[which(cell$celltype=="MCF7M1"),]$cellid

mcf7ratio <- c()
for (i in mcf7id)
{
    print(i)
    ratio <- nrow(hitdf[which(hitdf$cell_id == mcf7id[i]),]) / nrow(contact[which(contact$cell_id == mcf7id[i]),])
    mcf7ratio <- c(mcf7ratio, ratio)
}

ratiolist2 <- mcf7ratio[complete.cases(mcf7ratio)]

### MCF7TR
mcf7id <- cell[which(cell$celltype=="MCF7TR"),]$cellid

mcf7ratio <- c()
for (i in mcf7id)
{
    print(i)
    ratio <- nrow(hitdf[which(hitdf$cell_id == mcf7id[i]),]) / nrow(contact[which(contact$cell_id == mcf7id[i]),])
    mcf7ratio <- c(mcf7ratio, ratio)
}

ratiolist3 <- mcf7ratio[complete.cases(mcf7ratio)]

wilcox.test(ratiolist1, ratiolist3)
#p-value = 0.1792
par(mar = c(5, 10, 4, 2) + 0.1)
boxplot(ratiolist1*100, ratiolist3*100, medcol="red", xaxt="n", ylab="Contacts in CNVs (%)", ylim=c(0, 20), cex.names=2, cex.axis=2, cex.lab=2)
