withintscore <- read.csv("notindeg_intscore_c4_dataframe.txt", header=TRUE, sep="\t")
withintscore <- withintscore[order(-withintscore$intscore),]
rownames(withintscore) <- withintscore$gene
#newscore <- withintscore[1:50,]
newscore <- withintscore
newscore <- newscore[-c(1:7)]
library("pheatmap")
pheatmap(newscore, cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=FALSE, show_colnames=FALSE, color = colorRampPalette(c("white", "blue", "black"))(100))

### figure for enrichment pathway

library(ggplot2)

### david/go signaling pathway enrichment
davidgo <- read.csv("hifreq153genesdavidgo.txt", header=TRUE, sep="\t")
#davidgo <- davidgo[1:10,]

nes <- davidgo$Fold.Enrichment
fdr <- davidgo$PValue
genenumber <- davidgo$Count
pathway <- davidgo$Term

a3adata <- data.frame("PATHWAY" = pathway, "NES" = nes, "FDR" = fdr, "Gene_Number" = genenumber)

ggplot(a3adata,aes(x=FDR,y=PATHWAY,size=Gene_Number,color=NES)) + 
  geom_point() + theme_bw() +
  scale_colour_gradient(low="red",high="green") + 
  theme(axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 20)) +
  theme(axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 20)) +
  labs(x="p value",y="Pathway",title="Top20 enriched pathways",
       colour="Enrichment fold",size="Gene number") +
  theme(plot.title = element_text(hjust = 0.5))
