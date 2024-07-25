###Survival analysis

##################
###Loi et al
library(survival)
library(pheatmap)

cohortinfo <- read.csv("/data/yufan/schic/cohort/LoiGeneProfilingBCaTam/cohortinfo1.txt", header=TRUE, sep="\t")
cohortmat <- read.csv("/data/yufan/schic/cohort/LoiGeneProfilingBCaTam/cohortmatrix1.txt", header=TRUE, sep="\t")
recurrent = cohortinfo[which((cohortinfo$tamoxifen==1) & (cohortinfo$radioTherapy==0) & (cohortinfo$chemoTherapy==0) & (cohortinfo$ER.==1)),]
nrow(recurrent)
#[1] 349
rownames(recurrent) <- 1:nrow(recurrent)

survobj <- with(recurrent, Surv(rfs.time, rfs)) 
fit0 <- survfit(survobj~1, data=recurrent) 
summary(fit0) 
plot(fit0, xlab="Survival Time in Days",  
     ylab="% Surviving", yscale=100,   
     main="Survival Distribution (Overall)") 

recurrent$old = recurrent$age %/% 60
fit1 <- survfit(survobj~old,data=recurrent) 
plot(fit1, xlab="Survival time in days",
  ylab="Recurrence-free survival (%)", yscale=100,col=c("blue", "red"),
  main="Survival distributions by age")
  legend("topright", title="Age",c("<60", ">=60"),
  fill=c("blue", "red")) 

survdiff(survobj~old, data=recurrent)

chromatingene <- c('BRWD1', 'ATXN7', 'KDM5B', 'KMT2E', 'ENY2', 'HMG20B', 'MBIP', 'PRMT6', 'CCND1', 'ELP2', 'KMT5A', 'SMARCB1', 'TADA3', 'JADE1', 'MORF4L1')
rnapol2 <- c('PABPN1', 'PPM1D', 'BNIP3L', 'NELFA', 'MED1', 'YEATS4', 'SRSF1', 'RPRD1A', 'CEBPB', 'UBE2I', 'THOC7', 'DYRK2', 'COX7A2L', 'BTG2', 'TXNRD1', 'EAF1', 'CNOT6', 'RPS27A', 'TIGAR', 'ZNF223', 'ZNF221')

#genename <- chromatingene[9]
genename <- rnapol2[17]
esr1 <- cohortmat[which(cohortmat$symbol==genename),]
#esr1 <- cohortmat[which(cohortmat$symbol==rnapol2[21]),]
esr1$symbol <- NULL
geneexp <- apply(esr1, 2, function(x) mean(x))
geneexp <- as.data.frame(geneexp)
colnames(geneexp) <- c("exp")
geneexp$samplename <- rownames(geneexp)

recurrent["esr1"] <- 0
for (i in rownames(recurrent))
{
    print(paste(i, "/", nrow(recurrent)))
    recurrent[i, "esr1"] <- geneexp[which(geneexp$samplename == as.character(recurrent[i, "samplename"])),][1, "exp"]
}

top25 <- quantile(recurrent$esr1, 0.75)

recurrent$esr1b <- ifelse(recurrent$esr1 > top25, 1, 0)

fit1 <- survfit(survobj~esr1b,data=recurrent) 
survdiff(survobj~esr1b, data=recurrent) 

plot(fit1, xlab="Time (Days)",
  ylab="Recurrence-free survival (%)", yscale=100,col=c("blue", "red"),
  main="Survival distributions by gene expression")

legend("bottomleft", inset=0.02, title=paste(genename, "gene expression"),c("Low", "High"),
  fill=c("blue", "red")) # 绘图

#############################
####heatmap
#cohortgroup1 <- cohortmat[which(cohortmat$symbol %in% chromatingene),]
cohortgroup1 <- cohortmat[which(cohortmat$symbol %in% rnapol2),]
rfmap1 <- aggregate(cohortgroup1[, as.character(recurrent[which(recurrent$rfs==0),]$samplename)], by=list(cohortgroup1$symbol), FUN=mean)
relapsemap1 <- aggregate(cohortgroup1[, as.character(recurrent[which(recurrent$rfs==1),]$samplename)], by=list(cohortgroup1$symbol), FUN=mean)
rownames(rfmap1) <- rfmap1[,"Group.1"]
rfmap1[,"Group.1"] <- NULL
rownames(relapsemap1) <- relapsemap1[,"Group.1"]
relapsemap1[,"Group.1"] <- NULL
#library(pheatmap)
#library(RColorBrewer)
breaksList = seq(-3, 3, by = 0.01)
#pheatmap(rfmap1, border_color=NA, color = colorRampPalette(c("green", "black", "red"))(300), cluster_rows=FALSE, fontsize_row=10)
#pheatmap(relapsemap1, border_color=NA, color = colorRampPalette(c("green", "black", "red"))(300), cluster_rows=FALSE, fontsize_row=10)
pheatmap(rfmap1,
    border_color=NA,
    color = colorRampPalette(c("#00FFFF", "#FFFFFF", "#FF00FF"))(length(breaksList)),
    ###color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
    cluster_rows=FALSE,
    breaks = breaksList,
    fontsize_row=10)
pheatmap(relapsemap1,
    border_color=NA,
    color = colorRampPalette(c("#00FFFF", "#FFFFFF", "#FF00FF"))(length(breaksList)),
    cluster_rows=FALSE,
    breaks = breaksList,
    fontsize_row=10)


####################################################
###Sotiriou et al.
###

cohortinfo <- read.csv("/data/yufan/schic/cohort/Sotiriou_Christos/cohortinfo3.txt", header=TRUE, sep="\t")
cohortmat <- read.csv("/data/yufan/schic/cohort/Sotiriou_Christos/cohortmatrix3.txt", header=TRUE, sep="\t")
recurrent = cohortinfo[which((cohortinfo$tamoxifen==1) & (cohortinfo$radioTherapy==0) & (cohortinfo$chemoTherapy==0) & (cohortinfo$ER.==1)),]
nrow(recurrent)
#[1] 64
rownames(recurrent) <- 1:nrow(recurrent)

cohortmat[,"symbol"] <- ""
for (i in 1:nrow(cohortmat))
{
    print(paste(i, "/", nrow(cohortmat)))
    cohortmat[i, "symbol"] <- strsplit(as.character(cohortmat[,"Gene.Symbol"][i]), "///")[[1]][1]
}

chromatingene <- c('BRWD1', 'ATXN7', 'KDM5B', 'KMT2E', 'ENY2', 'HMG20B', 'MBIP', 'PRMT6', 'CCND1', 'ELP2', 'KMT5A', 'SMARCB1', 'TADA3', 'JADE1', 'MORF4L1')
rnapol2 <- c('PABPN1', 'PPM1D', 'BNIP3L', 'NELFA', 'MED1', 'YEATS4', 'SRSF1', 'RPRD1A', 'CEBPB', 'UBE2I', 'THOC7', 'DYRK2', 'COX7A2L', 'BTG2', 'TXNRD1', 'EAF1', 'CNOT6', 'RPS27A', 'TIGAR', 'ZNF223', 'ZNF221')

####heatmap
#cohortgroup1 <- cohortmat[which(cohortmat[,"symbol"] %in% chromatingene),]
cohortgroup1 <- cohortmat[which(cohortmat[,"symbol"] %in% rnapol2),]
rfmap1 <- aggregate(cohortgroup1[, as.character(recurrent[which(recurrent$rfs==0),]$geo_accession)], by=list(cohortgroup1$symbol), FUN=mean)
relapsemap1 <- aggregate(cohortgroup1[, as.character(recurrent[which(recurrent$rfs==1),]$geo_accession)], by=list(cohortgroup1$symbol), FUN=mean)
rownames(rfmap1) <- rfmap1[,"Group.1"]
rfmap1[,"Group.1"] <- NULL
rownames(relapsemap1) <- relapsemap1[,"Group.1"]
relapsemap1[,"Group.1"] <- NULL
#library(pheatmap)
#library(RColorBrewer)
breaksList = seq(4, 10, by = 0.01)
#pheatmap(rfmap1, border_color=NA, color = colorRampPalette(c("green", "black", "red"))(300), cluster_rows=FALSE, fontsize_row=10)
#pheatmap(relapsemap1, border_color=NA, color = colorRampPalette(c("green", "black", "red"))(300), cluster_rows=FALSE, fontsize_row=10)
pheatmap(rfmap1,
    border_color=NA,
    color = colorRampPalette(c("green", "black", "red"))(length(breaksList)),
    ###color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
    cluster_rows=FALSE,
    breaks = breaksList,
    fontsize_row=10)
pheatmap(relapsemap1,
    border_color=NA,
    color = colorRampPalette(c("green", "black", "red"))(length(breaksList)),
    cluster_rows=FALSE,
    breaks = breaksList,
    fontsize_row=10)
