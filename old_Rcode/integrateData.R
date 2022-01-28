library(Matrix)
library(Seurat)
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/scQC.R')
library(ggplot2)
library(clustermole)
library(fgsea)

readData <- function(X){
  D <- readMM(paste0(X,'matrix.mtx.gz'))
  C <- readLines(paste0(X,'barcodes.tsv.gz'))
  C <- unlist(lapply(strsplit(C, '-'), function(X){X[1]}))
  G <- read.table(paste0(X, 'features.tsv.gz'), stringsAsFactors = FALSE)
  G <- G[,2]
  rownames(D) <- G
  colnames(D) <- C
  return(D)
}

WT1 <- readData('GSE140511/GSM4173504_WT_1_')
WT2 <- readData('GSE140511/GSM4173505_WT_2_')
WT3 <- readData('GSE140511/GSM4173506_WT_3_')
KO1 <- readData('GSE140511/GSM4173507_Trem2_KO_1_')
KO2 <- readData('GSE140511/GSM4173508_Trem2_KO_2_')
KO3 <- readData('GSE140511/GSM4173509_Trem2_KO_3_')

WT1 <- CreateSeuratObject(WT1, project = 'WT1')
WT2 <- CreateSeuratObject(WT2, project = 'WT2')
WT3 <- CreateSeuratObject(WT3, project = 'WT3')
KO1 <- CreateSeuratObject(KO1, project = 'KO1')
KO2 <- CreateSeuratObject(KO2, project = 'KO2')
KO3 <- CreateSeuratObject(KO3, project = 'KO3')

ALL <- merge(WT1, WT2)
ALL <- merge(ALL, WT3)
ALL <- merge(ALL, KO1)
ALL <- merge(ALL, KO2)
ALL <- merge(ALL, KO3)

rm(WT1,WT2,WT3,KO1,KO2,KO3)

ALL <- scQC(ALL, mtThreshold = 0.05)

ALL <- NormalizeData(ALL)
ALL <- FindVariableFeatures(ALL)
ALL <- ScaleData(ALL)
ALL <- RunPCA(ALL, verbose = FALSE)
ALL <- RunUMAP(ALL, dims = 1:50)
png('allDataInt.png', width = 1500, height = 1200, res = 300)
UMAPPlot(ALL) + theme_bw() + xlab('UMAP 1') + ylab('UMAP 2')
dev.off()

set.seed(1)
ALL <- FindNeighbors(ALL, reduction = 'umap', dims = 1:2)
set.seed(1)
ALL <- FindClusters(ALL, resolution = 0.02)

png('allClusters.png', width = 1300, height = 1200, res = 300)
UMAPPlot(ALL, label = TRUE) + theme_bw() + xlab('UMAP 1') + ylab('UMAP 2') + theme(legend.position = 'none')
dev.off()

ALL <- subset(ALL, idents = names(table(Idents(ALL))[table(Idents(ALL)) >= 500]))
UMAPPlot(ALL, label = TRUE) + theme_bw() + xlab('UMAP 1') + ylab('UMAP 2') + theme(legend.position = 'none')

DE <- FindAllMarkers(ALL)
write.csv(DE, 'ctDE.csv')
DE <- read.csv('ctDE.csv', row.names = 1)

ctList <- clustermole_markers(species = 'mm')
ctList <- ctList[ctList$db %in% c('PanglaoDB'),]
ctList <- ctList[ctList$organ %in% 'Brain',]
ctNames <- unique(ctList$celltype)
ctList <- lapply(ctNames, function(X){
  ctList$gene[ctList$celltype %in% X]
})
names(ctList) <- ctNames

Idents(ALL) <- ALL$RNA_snn_res.0.02
ctID <- pbapply::pbsapply(levels(Idents(ALL)), function(C){
  FC <- DE$avg_logFC[DE$cluster %in% C]
  names(FC) <- DE$gene[DE$cluster %in% C]
  E <- fgseaMultilevel(ctList, FC)
  E <- E[E$NES > 0 & E$padj < 0.05,]
  E <- E[order(E$pval,1/E$NES),]
  E$leadingEdge <- unlist(lapply(E$leadingEdge, function(X){paste0(X, collapse = ';')}))
  return(E[1,])
})
ctID <- t(ctID)
ctID <- as.data.frame(ctID)
ctID[,-8]

levels(Idents(ALL)) <- unlist(ctID[,1])

png('cellTypes.png',width = 1300*1.4, height = 1200*1.4, res = 300)
UMAPPlot(ALL, label = TRUE) + theme_bw() + xlab('UMAP 1') + ylab('UMAP 2') + theme(legend.position = 'none') + xlim(c(-15,15))
dev.off()

t(table(ALL$orig.ident, Idents(ALL)))

ALL$CellTypes <- Idents(ALL)
Idents(ALL) <- ALL$RNA_snn_res.0.02

Neurons <- subset(ALL, idents = '1')
Astrocytes <- subset(ALL, idents = '5')
Microglia <- subset(ALL, idents = '11')

save(Neurons, file = 'Neurons.RData')
save(Astrocytes, file = 'Astrocytes.RData')
save(Microglia, file = 'Microglia.RData')
