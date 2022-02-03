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
GO_up <- c('apoptotic_signaling_pathway',"neurogenesis")
color_up <- c('purple','red')
GO_down <- c("cell_cycle","microtubule.based_process")
color_down <- c('green','navy')
go_OI <- c(GO_up, GO_down)
go_OI <- unique(go_OI)
colors_GO_pos <- c(color_up, color_down)
colors_GO_pos <- unique(colors_GO_pos)
print(go_OI)
print(colors_GO_pos)



setwd("C:/Users/mariu/Documents/newRNASeq_analysis")

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




mat <- assays(se)$tpm
#mat <- log2(mat+0.1)
Expression_df <- data.frame(res, symbol= rowData(se)$symbol,mat)
print(Expression_df[1:5,])


go_up <- read.table("analysis_upregulated_FD2FDR005.csv", sep= ";", header = TRUE)
go_down <- read.table("analysis_downregulated_FD2FDR005.csv", sep = ";", header = TRUE)








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

GO_table2b <- GO_table2
GO_table2b <- data.frame(GO_table2b)
GO_table2b$symbol <- rownames(GO_table2b)
write.table(Expression_df, 'expression_D5_05012021.csv', sep = ';', row.names=FALSE, quote = FALSE)

Expression_df_selected <- Expression_df[Expression_df$symbol %in% genes,]
Expression_df_selected <- Expression_df_selected[order(match(Expression_df_selected$symbol, genes)),]
Expression_df_selected <- Expression_df_selected[!duplicated(Expression_df_selected$symbol), ]
GO_table_final <- cbind(GO_table, Expression_df_selected, GO_table2)


# data.frame(GO_table ,res_GO, mat,GO_table2)
head(GO_table_final)

write.table(GO_table_final, 'GO_oi_table_D7_05012021.csv', sep = ';', row.names=FALSE, quote = FALSE)
#GO_table_final <- read.table('GO_oi_table_D5_05012021.csv', sep = ';', header = TRUE)



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
col_breaks = c(seq(-2,-0.4,length=100),  # for red
               seq(-0.5,0.5,length=50),              # for white
               seq(0.6,2,length=100))              # for blue
ann <- data.frame(condition=se$condition, rnai=se$RNAi)

rownames(ann) <- colnames(se)
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
rows_oi <- 50

GO_table_final <-  GO_table_final[GO_table_final$symbol !='Fgf2',]
select_vector <- GO_table_final$padj
GO_table_sort <- GO_table_final[order(select_vector),]
head(GO_table_sort)
GO_table_red <- GO_table_sort[1:50,]

library(textshape)

library("scales")

mat_ex <-GO_table_red[,(ncol(GO_table_red)-length(go_OI)+1):ncol(GO_table_red) ]
rownames(mat_ex) <- GO_table_red$symbol
head(mat_ex)
length(mat_ex[,1])
mat_ex_sort <- cluster_matrix(mat_ex, dim = 'row')
dist_mat <- dist(mat_ex, 'euclidean')
clust_mat <- hclust(dist_mat, method = "complete")
sort_vec <- rownames(mat_ex_sort )
head(mat_ex_sort)

GO_table_selected <- GO_table_red[order(match(GO_table_red$symbol, sort_vec)),]

#annRow2 <- annRow[rownames(GO_table2) %in% rows_oi,]
#annRow2 <- annRow2[order(match(rownames(annRow2), sort_vec)),1:length(go_OI)]
annRow2 <- GO_table_selected[,1:(length(go_OI))]
rownames(annRow2) <- GO_table_selected$symbol
mat2 <- GO_table_selected[,13:20]
rownames(mat2) <- GO_table_selected$symbol
head(annRow2)
#mat2 <-mat[rownames(GO_table2) %in% rows_oi,]
#mat2 <- mat2[order(match(rownames(mat2), sort_vec)),]
#annRow <- data.frame(annRow)
rownames(ann) <- colnames(mat2)
mat2 <- log2(mat2+0.1)
setEPS()
postscript('pheatmap_GO_RNASeq-D7Vb2_05012022.eps', width = 7, height = 50*0.2)
par(lwd=1,cex=1)
pheatmap(mat2, annotation_col = ann, annotation_row = annRow2, scale = "row", annotation_colors = ann_colors, color = my_palette, breaks = col_breaks ,cluster_rows= F, cluster_cols=T)

#
#pheatmap(mat, annotation_col = ann,annotation_row = annRow, scale = "row", annotation_colors = ann_colors, color = my_palette, cluster_rows = F)
dev.off()
postscript('pheatmap_GO_RNASeq-D7V1b_05012022.eps', width = 7, height = 50*0.2)

pheatmap(mat2, annotation_col = ann, annotation_row = annRow2, scale = "row", annotation_colors = ann_colors, color = my_palette, breaks = col_breaks ,cluster_rows= F, cluster_cols=F)
dev.off()
  
  
  
  
  
  
  
  
  
  
  
  ####
  
  