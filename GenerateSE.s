BiocManager::install("SummarizedExperiment")
BiocManager::install("biomaRt")
BiocManager::install("summarytools")
library(SummarizedExperiment)
library(biomaRt)
library(summarytools)

samples <- read.delim("samples.txt", header=F) 

  cmat <- c()
  tmat <- c()

  for (sample in samples$V2) {
    sample.name <-
      starcounts <- read.delim(paste(sample,".ReadsPerGene.out.tab",sep=""), header=F)
    cmat <- cbind(cmat, starcounts[,4])
    if (!sum(starcounts[,4])==0) {
      rsemtpm <- read.delim(paste(sample,"_rsem.genes.results",sep=""), header=T)
      tmat <- cbind(tmat, rsemtpm$TPM)
    } else {
      tmat <- cbind(tmat, rep(NA, nrow(tmat)))
    }

  }
  rownames(cmat) <- starcounts[,1]
  rownames(tmat) <- rsemtpm$gene_id

  flt <- !(grepl("^N_", rownames(cmat)))
  cmat <- cmat[flt,]

  colnames(cmat) <- samples$Sample.ID
  colnames(tmat) <- samples$Sample.ID
  tmat <- tmat[match(rownames(cmat), rownames(tmat)),]

  samplesM <- data.frame(samples, row.names = samples$V2)
  samplesM$RNAi <- sapply(as.character(samplesM$V1), function(x){strsplit(x,"-")[[1]][1]})
  samplesM$condition <- ifelse(samplesM$RNAi=="ctrl", "ctrl", "kd")
  samplesM$replicate <- sapply(as.character(samplesM$V2), function(x){strsplit(x,"_")[[1]][2]})
  
  ## get Annotation data
  ml <- listDatasets(useMart("ensembl"))
  print(ml[ml$dataset=="mmusculus_gene_ensembl",])
  ensembl <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
  bm <- getBM(filters = "ensembl_gene_id",
              attributes = c("ensembl_gene_id","entrezgene_id","external_gene_name","gene_biotype","description"),
              values = rownames(cmat),
              mart =  useEnsembl(biomart = "ensembl", 
                   dataset = "mmusculus_gene_ensembl", 
                   mirror = "useast"))  #useast, uswest, asia

  annot <- data.frame(id=rownames(cmat))
  annot$symbol <- bm$external_gene_name[match(rownames(cmat), bm$ensembl_gene_id)]
  annot$entrez <- bm$entrezgene[match(rownames(cmat), bm$ensembl_gene_id)]
  annot$biotype <- bm$gene_biotype[match(rownames(cmat), bm$ensembl_gene_id)]
  annot$description <- bm$description[match(rownames(cmat), bm$ensembl_gene_id)]


  se <- SummarizedExperiment(assays=list(counts=cmat, tpms=tmat), colData=samplesM, rowData=annot)
  dim(se)
  saveRDS(se, file="se.rds")