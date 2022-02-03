library(SummarizedExperiment)
library(DESeq2)
library(parallel)
library(ComplexHeatmap)
library(biomaRt)
library(xml2)
library(XML)
library(pheatmap)
library(ggplot2)


####generate pheatmap#####
GO_up <- c('apoptotic_signaling_pathway')
color_up <- c('purple')
GO_down <- c("cell_cycle","microtubule.based_process","neurogenesis")
color_down <- c('green','navy','red')
go_OI <- c(GO_down, GO_up)
go_OI <- unique(go_OI)
colors_GO_pos <- c(color_down, color_up)
colors_GO_pos <- unique(colors_GO_pos)
print(go_OI)
print(colors_GO_pos)



setwd("C:/Users/mariu/Documents/newRNASeq_analysis")

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
print(res[1:55,])


go_up <- read.table("analysis_d5_go_up.csv", sep= ";", header = TRUE)
go_down <- read.table("analysis_d5_go_do.csv", sep = ";", header = TRUE)
go_down <- go_down[go_down$label %in% GO_down,]

go_up <- go_up[ go_up$label %in% GO_up,]
#head(go_up)
#head(go_down)

merged_GOs <- rbind(go_up, go_down)


GOs <- data.frame(merged_GOs$label, merged_GOs$mapped_id,1)
colnames(GOs) <- c('label','symbol','numbre')


head(GOs)

selected_GOs <- GOs[GOs$label %in% go_OI,]

head(selected_GOs)
'generate 2 go tables: one with yes or no for plot, second for filtering'
GO_table <- tapply(selected_GOs$numbre,list(selected_GOs$symbol ,selected_GOs$label)  ,sum)
GO_table2 <- tapply(selected_GOs$numbre,list(selected_GOs$symbol ,selected_GOs$label)  ,sum, na.rm=TRUE)
GO_table[!is.na(GO_table)] <- 'yes'
GO_table[is.na(GO_table)] <- 'no'



GO_table2[is.na(GO_table2)] <- 0
GO_table <- data.frame(GO_table)
GO_table$symbol <- rownames(GO_table)
head(GO_table)

###avoid using genes with 3 go labels
row_sum_table <- rowSums(GO_table2) 
GO_table2 <- GO_table2[row_sum_table <=2,]
GO_table <- GO_table[row_sum_table <=2,]
head(GO_table)
genes <- rownames(GO_table)





res_GO <- res[row.names(res) %in% GENES,]
mat <- assays(se[GENES])$tpm
mat <- log2(mat+0.1)
rownames(mat) <- GO_table$symbol

GO_table_final <- data.frame(GO_table ,res_GO, mat,GO_table2)
head(GO_table_final)

write.table(GO_table_final, 'GO_oi_table_D5.csv', sep = ';', row.names=FALSE, quote = FALSE)
GO_table_final <- read.table('GO_oi_table_D5.csv', sep = ';', header = TRUE)



#dynamic filling of ann colors
  
ann_colors = list(
    rnai = c(ctrl = "blue", SH1 = "red"),
    condition = c(ctrl = "blue", kd = "red")

)
  
for(i in 1:length(go_OI)){
  list_row <- list(condition = c(yes = colors_GO_pos[i], no ='white'))

  names(list_row) <- as.character(go_OI[i])

  
  ann_colors <- append(ann_colors,list_row) 
}

my_palette <- colorRampPalette(c("darkgreen", "yellow", "darkorange"))(n = 249)

# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(-2,-0.6,length=100),  # for red
               seq(-0.5,0.5,length=50),              # for white
               seq(0.6,2,length=100))              # for blue
ann <- data.frame(condition=se_in$condition, rnai=se_in$RNAi)

rownames(ann) <- colnames(se_in)
setEPS()
### 2 strategies to fill the annRow , both worked
annRow <-GO_table[,1:length(go_OI)]


for(i in 1:length(length(go_OI))){
  df_or <- data.frame(annRow[,i])
  names(df_or) <- colnames(annRow)[i]
  print(head(df_or))
  annRow2 <- data.frame( df_or,annRow)
}
annRow2 <-annRow2[,1:length(go_OI)]
rownames(annRow2) <- rownames(GO_table)
annRow <- data.frame(annRow)
colnames(annRow)
head(annRow)
head(annRow2)

#selecting data for plot
genesOI_table <- read.table('selected_GenesD5.csv',sep=";", header = TRUE)
genesOI <- c("")
GO_table_final <- data.frame(GO_table ,res_GO, mat,GO_table2)
select_vector <- GO_table_final$padj
select_vector <- select_vector[!GO_table_final$symbol %in% genesOI]

GO_table_final <- GO_table_final[order(select_vector),]

rows_oi <- c(genesOI,GO_table_final$symbol[1:100])

rows_oi <- rows_oi[1:100]
head(rows_oi)

library(textshape)

library("scales")

mat_ex <-GO_table2[rownames(GO_table2) %in% rows_oi,]
head(mat_ex)
length(mat_ex[,1])
mat_ex_sort <- cluster_matrix(mat_ex, dim = 'row')
dist_mat <- dist(mat_ex, 'euclidean')
clust_mat <- hclust(dist_mat, method = "complete")
sort_vec <- rownames(mat_ex_sort )
head(mat_ex_sort)
annRow2 <- annRow[rownames(GO_table2) %in% rows_oi,]
annRow2 <- annRow2[order(match(rownames(annRow2), sort_vec)),1:length(go_OI)]
head(annRow2)
mat2 <-mat[rownames(GO_table2) %in% rows_oi,]
mat2 <- mat2[order(match(rownames(mat2), sort_vec)),]
#annRow <- data.frame(annRow)

setEPS()
postscript('pheatmap_GO_RNASeq-D5Vb2_05012022.eps', width = 7, height = length(rows_oi)*0.2)
par(lwd=3,cex=3)
pheatmap(mat2, annotation_col = ann, annotation_row = annRow2, scale = "row", annotation_colors = ann_colors, color = my_palette, breaks = col_breaks, ,cluster_rows= F, cluster_cols=T)

#
#pheatmap(mat, annotation_col = ann,annotation_row = annRow, scale = "row", annotation_colors = ann_colors, color = my_palette, cluster_rows = F)
dev.off()
postscript('pheatmap_GO_RNASeq-D5V1b_05012022.eps', width = 7, height = length(rows_oi)*0.2)

pheatmap(mat2, annotation_col = ann, annotation_row = annRow2, scale = "row", annotation_colors = ann_colors, color = my_palette, breaks = col_breaks ,cluster_rows= F, cluster_cols=F)
dev.off()
  
  
  
  
  
  
  
  
  
  
  
  ####
  
  