library(SummarizedExperiment)
library(DESeq2)
library(parallel)
library(ComplexHeatmap)
library(biomaRt)
library(xml2)
library(XML)
library(pheatmap)
library(ggplot2)


####generate pheatmap  Day5 #####



setwd("C:/Users/mariu/Documents/newChRIP_analysis")

dataOI <- read.table('annotated_RUS_Peaks_10012022_corr3.csv', header = TRUE, sep = ";")
genes <- dataOI$gene
loci <- dataOI$locus

se_d5 <- readRDS("seD5.rds")

se <- se_d5[,colData(se_d5)$RNAi %in% c("ctrl","SH1")]
#data <- data[data$undesirable == 0,]
seD5 <- se[apply(assays(se_d5)$tpm, 1, function(row) all(row !=0 )), ]


se_d7 <- readRDS("seD7.rds")

se <- se_d7[,colData(se_d7)$RNAi %in% c("ctrl","SH1")]
#data <- data[data$undesirable == 0,]
seD7 <- se[apply(assays(se)$tpm, 1, function(row) all(row !=0 )), ]



MAT_D5 <- seD5[match(genes, rowData(seD5)$symbol),]
MAT_D5 <- MAT_D5[na.omit(rownames(MAT_D5)),]

MAT_D7 <- seD7[match(genes, rowData(seD7)$symbol),]
MAT_D7 <- MAT_D7[na.omit(rownames(MAT_D7)),]


MAT_D5 <- MAT_D5[match(rownames(MAT_D7), rownames(MAT_D5)),]
MAT_D5 <- MAT_D5[na.omit(rownames(MAT_D5)),]

MAT_D7 <- MAT_D7[match(rownames(MAT_D5), rownames(MAT_D7)),]
MAT_D7 <- MAT_D7[na.omit(rownames(MAT_D7)),]

mat_d5 <- assays(MAT_D5)$tpm
#mat_d5 <- log2(mat_d5) +0.1
rownames(mat_d5) <- rowData(MAT_D5)$symbol

mat_d7 <- assays(MAT_D7)$tpm
#mat_d7 <- log2(mat_d7) +0.1
rownames(mat_d7) <- rowData(MAT_D7)$symbol


###ctrl normalization
ctrl_D5 <- apply(mat_d5[,grep("ctrl", colnames(mat_d5))],1,mean)
mat_d5_norm <- mat_d5 / ctrl_D5
ctrl_D7 <- apply(mat_d7[,grep("ctrl", colnames(mat_d7))],1,mean)
mat_d7_norm <- mat_d7 / ctrl_D7

mat <- cbind(mat_d5_norm, mat_d7_norm)

###clustering

#mat <- cbind(mat_d5_norm, mat_d7_norm)
#distances <- dist(mat, method ="euclidean", diag = TRUE, upper = TRUE, p = 1)
#clustering <- hclust(distances, method ="complete")
#setEPS()
#postscript('dendogramm_cluster.eps')
#, width = length(rownames(mat))*0.015, height=1)

#plot(clustering)
#dev.off()
#plot(clustering)

#ordering <- clustering$order

col_breaks = c(seq(-2,-0.3,length=100),  # for red
               seq(-0.2,0.2,length=50),              # for white
               seq(0.3,2,length=100))              # for blue

ann <- data.frame(rnai=se_d5$RNAi)
ann
rownames(ann) <- colnames(se_d5)
ann2 <- data.frame(rnai=se_d7$RNAi)
ann2
rownames(ann2) <- colnames(se_d7)
ann3 <- rbind(ann,ann2)





ann_colors = list(
  rnai = c(ctrl = "blue", SH1 = "red")
)
rownames(ann2) <- colnames(se_d7)
ann_colors = list(
  rnai = c(ctrl = "blue", SH1 = "red"),
  locus = c(promoter = "green", Intron = "cyan" ,  Intergenic = "darkred")
)



hm_1<-pheatmap(mat, scale = "row", cluster_rows=T, cluster_cols=F) 
ordering2 <- hm_1$tree_row$order
dev.off()
#ordering
### reorder matrixes

mat_d5_ordered <- mat_d5[ordering2,] 
mat_d7_ordered <- mat_d7[ordering2,] 

lociOI <- dataOI$locus[match(rownames(mat_d5_ordered), dataOI$gene)]

lociOI <- data.frame(locus = lociOI)
rownames(lociOI) <- rownames(mat_d5_ordered)
annRow <- lociOI
my_palette <- colorRampPalette(c("darkgreen", "yellow", "darkorange"))(n = 249)

# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(-2,-0.6,length=100),  # for red
               seq(-0.5,0.5,length=50),              # for white
               seq(0.6,2,length=100))              # for blue

setEPS()
postscript('pheatmap_ChIRP_genes_Day5_hclustered_withLocus_10012022_3.eps', width = 5, height = length(rownames(mat))*0.15)


pheatmap(log2(mat_d5_ordered+0.1), scale = "row", annotation_col = ann, annotation_row = annRow, annotation_colors = ann_colors, color = my_palette, breaks = col_breaks , cluster_rows=F, cluster_cols=F)


dev.off()



setEPS()
postscript('pheatmap_ChIRP_genes_Day7_hclustered_withLocus_10012022_3.eps', width = 5, height = length(rownames(mat))*0.15)


pheatmap(log2(mat_d7_ordered+0.1), scale = "row", annotation_col = ann2, annotation_row = annRow, annotation_colors = ann_colors, color = my_palette, breaks = col_breaks , cluster_rows=F, cluster_cols=F)


dev.off()

setEPS()
postscript('pheatmap_ChIRP_log2genes_merged_timepoints_hclustered_withLocus_10012022.eps', width = 5, height = length(rownames(mat))*0.15)

pheatmap(mat, scale = "row", annotation_col = ann2, annotation_row = annRow, annotation_colors = ann_colors2, color = my_palette, breaks = col_breaks , cluster_rows=T, cluster_cols=T)
dev.off()


###make correlation plots for the side
my_palette2 <- colorRampPalette(c("darkgreen", "white", "purple"))(n = 13)
col_breaks2 = c(seq(-1,-0.2,length=6),  # for red
               seq(-0.1,0.1,length=2),              # for white
               seq(0.2,1,length=6))  


genes <- rownames(mat)[ordering2]
genes2 <- c("Gm20754",genes)


#d5
se <-readRDS("seD5.rds")
se <- se[,colData(se)$RNAi %in% c("ctrl","SH1")]
matD5 <-assays(se)$tpm
rownames(matD5) = rowData(se)$symbol 
x2 <- matD5[rownames(matD5) %in% genes2,]


x2_corr <- cor(t(x2), method = 'pearson') 
x2_ref <- x2_corr[colnames(x2_corr)=="Gm20754",]
x2_ref <- data.frame(x2_ref)
x2_ref <- data.frame(gene = colnames(x2_corr), correlation = x2_ref)


x2_ref_red <- x2_ref[x2_ref$gene != "Gm20754", ]

x2_ref_red_reord <- x2_ref_red[match(genes,x2_ref_red$gene), ]
x2_ref_red_reord <- data.frame(x2_ref_red_reord, refgene="Gm20754")
x2_ref_red_reord$gene <- factor(x2_ref_red_reord$gene , levels = rev(x2_ref_red_reord$gene))
colnames(x2_ref_red_reord) <- c("gene","correlation","refgene")

cor5 <- x2_ref_red_reord$correlation


setEPS()
postscript('corrMatrixD5_10012022.eps', width = 2.5, height = length(genes)*0.15)

ggplot(x2_ref_red_reord, aes(y=gene,x=refgene, fill= correlation)) + 
  geom_tile()  +   scale_fill_gradient2(low = "purple",
                                        mid = "white",
                                        high = "darkgreen",
                                        midpoint = 0,
                                        space = "Lab",
                                        na.value = "grey50",
                                        guide = "colourbar",
                                        aesthetics = "fill",,limits = c(-1,1)
  ) +  theme_classic()


dev.off()

###d7


se <-readRDS("seD7.rds")
se <- se[,colData(se)$RNAi %in% c("ctrl","SH1")]
matD7 <-assays(se)$tpm
rownames(matD7) = rowData(se)$symbol 
x2 <- matD7[rownames(matD7) %in% genes2,]


x2_corr <- cor(t(x2), method = 'pearson') 
x2_ref <- x2_corr[colnames(x2_corr)=="Gm20754",]
x2_ref <- data.frame(x2_ref)
x2_ref <- data.frame(gene = colnames(x2_corr), correlation = x2_ref)


x2_ref_red <- x2_ref[x2_ref$gene != "Gm20754", ]

x2_ref_red_reord <- x2_ref_red[match(genes,x2_ref_red$gene), ]
x2_ref_red_reord <- data.frame(x2_ref_red_reord, refgene="Gm20754")
x2_ref_red_reord$gene <- factor(x2_ref_red_reord$gene , levels = rev(x2_ref_red_reord$gene))
colnames(x2_ref_red_reord) <- c("gene","correlation","refgene")


cor7 <- x2_ref_red_reord$correlation


setEPS()
postscript('corrMatrixD7_10012022.eps', width = 2.5, height = length(genes)*0.15)

ggplot(x2_ref_red_reord, aes(y=gene,x=refgene, fill= correlation)) + 
  geom_tile()  +   scale_fill_gradient2(low = "purple",
    mid = "white",
    high = "darkgreen",
    midpoint = 0,
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "fill",limits = c(-1,1)
  ) +  theme_classic()


dev.off()



se_in <- readRDS("seD5.rds")

se <- se_in[,colData(se_in)$RNAi %in% c("ctrl","SH1")]
#data <- data[data$undesirable == 0,]
#se <- se[apply(assays(se)$tpm, 1, function(row) all(row !=0 )), ]
colData(se)
dds <- DESeqDataSet(se, design = ~ condition)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds, contrast=c("condition","kd","ctrl"))
summary(res)
res$symbol <- rowData(se)$symbol

dataOI_5 <- res[res$symbol %in%genes,]
dataOI_5 <- data.frame(dataOI_5)
dataOI_5 <-dataOI_5[match(genes, dataOI_5$symbol),]
SE <- assays(se)$tpm
SE <- data.frame(SE)
SE$symbol <- rowData(se)$symbol
dataOI_5b <- SE[SE$symbol %in%genes,]
dataOI_5b <-dataOI_5b[match(genes, dataOI_5b$symbol),]


dataOI_5 <- cbind(dataOI_5,dataOI_5b, cor5 )
write.table(dataOI_5, 'targetGenes_expressionD5_10012022_expanded2.csv', sep = ";", row.names= F, quote= F)




se_in <- readRDS("seD7.rds")

se <- se_in[,colData(se_in)$RNAi %in% c("ctrl","SH1")]
#data <- data[data$undesirable == 0,]
#se <- se[apply(assays(se)$tpm, 1, function(row) all(row !=0 )), ]
colData(se)
dds <- DESeqDataSet(se, design = ~ condition)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds, contrast=c("condition","kd","ctrl"))
summary(res)
res$symbol <- rowData(se)$symbol

dataOI_7 <- res[res$symbol %in%genes,]
dataOI_7 <- data.frame(dataOI_7)
dataOI_7 <-dataOI_7[match(genes, dataOI_7$symbol),]


SE <- assays(se)$tpm
SE <- data.frame(SE)
SE$symbol <- rowData(se)$symbol
dataOI_7b <- SE[SE$symbol %in%genes,]
dataOI_7b <-dataOI_7b[match(genes, dataOI_7b$symbol),]


dataOI_7 <- cbind(dataOI_7,dataOI_7b, cor7 )
write.table(dataOI_7, 'targetGenes_expressionD7_10012022_expanded.csv', sep = ";", row.names= F, quote= F)