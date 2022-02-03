library(SummarizedExperiment)
library(DESeq2)
#library(parallel)
library(pheatmap)
library(biomaRt)

cuttOff_FD <- 2
cuttOff_FDR <- 0.05

setwd("C:/Users/mariu/Documents/newRNASeq_analysis")

se_in <- readRDS("seD7.rds")

se <- se_in[,colData(se_in)$RNAi %in% c("ctrl","SH1")]
#data <- data[data$undesirable == 0,]
se <- se[apply(assays(se)$tpm, 1, function(row) sum(row) >1), ]
colData(se)
dds <- DESeqDataSet(se, design = ~ condition)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds, contrast=c("condition","kd","ctrl"))
summary(res)


summary_df <-  data.frame(res,rowData(se)$symbol,assays(se)$tpm)
write.table(summary_df, paste(paste(cuttOff_FD,cuttOff_FDR, collapse= "_"),"Expression_genes_Day7_31012022.csv",collapse=""),sep =";" ,quote = FALSE, row.names = FALSE)


upregulated2 <- res[!is.na(res$padj) & res$padj < cuttOff_FDR & res$stat  >  cuttOff_FD,]
#head(upregulated2)
upregulated3 <- data.frame(upregulated2,rowData(se[rownames(upregulated2),])$symbol,assays(se[rownames(upregulated2),])$tpm)
write.table(upregulated3, paste(paste(cuttOff_FD,cuttOff_FDR, collapse= "_"),"upregulated_gm20754_genes_Day7_31012022.csv",collapse=""),sep =";" ,quote = FALSE, row.names = FALSE)
downregulated2 <- res[!is.na(res$padj) & res$padj < cuttOff_FDR & res$stat  <  (cuttOff_FD)*-1,]
#head(downregulated2)
downregulated3 <- data.frame(downregulated2,rowData(se[rownames(downregulated2),])$symbol,assays(se[rownames(downregulated2),])$tpm)
write.table(downregulated3, paste(paste(cuttOff_FD,cuttOff_FDR, collapse= "_"),"downregulated_gm20754_genes_Day7_31012022.csv",collapse=""),sep =";" ,quote = FALSE, row.names = FALSE)


se_in <- readRDS("seD5.rds")

se <- se_in[,colData(se_in)$RNAi %in% c("ctrl","SH1")]
#data <- data[data$undesirable == 0,]
se <- se[apply(assays(se)$tpm, 1, function(row) sum(row) >1), ]
colData(se)
dds <- DESeqDataSet(se, design = ~ condition)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds, contrast=c("condition","kd","ctrl"))
summary(res)


summary_df <-  data.frame(res,rowData(se)$symbol,assays(se)$tpm)
write.table(summary_df, paste(paste(cuttOff_FD,cuttOff_FDR, collapse= "_"),"Expression_genes_Day5_31012022.csv",collapse=""),sep =";" ,quote = FALSE, row.names = FALSE)


upregulated2 <- res[!is.na(res$padj) & res$padj < cuttOff_FDR & res$stat  >  cuttOff_FD,]
#head(upregulated2)
upregulated3 <- data.frame(upregulated2,rowData(se[rownames(upregulated2),])$symbol,assays(se[rownames(upregulated2),])$tpm)
write.table(upregulated3, paste(paste(cuttOff_FD,cuttOff_FDR, collapse= "_"),"upregulated_gm20754_genes_Day5_31012022.csv",collapse=""),sep =";" ,quote = FALSE, row.names = FALSE)
downregulated2 <- res[!is.na(res$padj) & res$padj < cuttOff_FDR & res$stat  <  (cuttOff_FD)*-1,]
#head(downregulated2)
downregulated3 <- data.frame(downregulated2,rowData(se[rownames(downregulated2),])$symbol,assays(se[rownames(downregulated2),])$tpm)
write.table(downregulated3, paste(paste(cuttOff_FD,cuttOff_FDR, collapse= "_"),"downregulated_gm20754_genes_Day5_31012022.csv",collapse=""),sep =";" ,quote = FALSE, row.names = FALSE)

