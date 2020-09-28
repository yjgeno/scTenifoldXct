library(fgsea)
GOBP <- gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2018')
BIOP <- gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=BioPlanet_2019')

load('na_Xct.RData')
DR <- scTenifoldNet::dRegulation(testOut$manifoldAlignment[,1:5])
X <- DR[grepl('^X_',DR$gene),]
Y <- DR[grepl('^Y_',DR$gene),]

ZX <- X$distance
BC <- MASS::boxcox(ZX~1)
BC <- abs(BC$x[which.max(BC$y)])
ZX <- as.numeric(scale(ZX^BC))
names(ZX) <- toupper(gsub('^X_','',X$gene))
EX <- fgseaMultilevel(GOBP, ZX)
EX <- EX[EX$NES > 0 & EX$pval < 0.05,]

ZY <- Y$distance
BC <- MASS::boxcox(ZY~1)
BC <- abs(BC$x[which.max(BC$y)])
ZY <- as.numeric(scale(ZY^BC))
names(ZY) <- toupper(gsub('^Y_','',Y$gene))
EY <- fgseaMultilevel(GOBP, ZY)
EY <- EY[EY$NES > 0 & EY$pval < 0.05,]

intersect(EX$pathway, EY$pathway)

plotEnrichment(GOBP$`long-term memory (GO:0007616)`, ZX)
plotEnrichment(GOBP$`long-term memory (GO:0007616)`, ZY)

source('https://raw.githubusercontent.com/dosorio/utilities/master/idConvert/hsa2mmu_SYMBOL.R')
N <- hsa2mmu_SYMBOL(unlist(EX$leadingEdge[EX$pathway == 'long-term memory (GO:0007616)']))
A <- hsa2mmu_SYMBOL(unlist(EY$leadingEdge[EY$pathway == 'long-term memory (GO:0007616)']))

ComplexHeatmap::Heatmap(as.matrix(testOut$xNet[paste0("X_",N),paste0("Y_",A)]), row_order = seq_along(N), column_order = seq_along(A), name = 'LTM WT') +
ComplexHeatmap::Heatmap(as.matrix(testOut$yNet[paste0("X_",N),paste0("Y_",A)]), row_order = seq_along(N), column_order = seq_along(A), name = 'LTM KO')

N <- hsa2mmu_SYMBOL(unlist(EX$leadingEdge[EX$pathway == 'regulation of neuronal synaptic plasticity (GO:0048168)']))
A <- hsa2mmu_SYMBOL(unlist(EY$leadingEdge[EY$pathway == 'regulation of neuronal synaptic plasticity (GO:0048168)']))

ComplexHeatmap::Heatmap(as.matrix(testOut$xNet[paste0("X_",N),paste0("Y_",A)]), row_order = seq_along(N), column_order = seq_along(A), name = 'LTM WT') +
ComplexHeatmap::Heatmap(as.matrix(testOut$yNet[paste0("X_",N),paste0("Y_",A)]), row_order = seq_along(N), column_order = seq_along(A), name = 'LTM KO')


load('Astrocytes.RData')
load('Neurons.RData')
library(Seurat)
library(Rmagic)
library(Matrix)
iNWT <- Neurons@assays$RNA@counts[,grepl('WT',Neurons$orig.ident)]
iNWT <- iNWT[rowMeans(iNWT) > 0.05,]
iNWT <- magic(as.matrix(t(iNWT)))
iNWT <- t(iNWT$result)

iNKO <- Neurons@assays$RNA@counts[,grepl('KO',Neurons$orig.ident)]
iNKO <- iNKO[rowMeans(iNKO) > 0.05,]
iNKO <- magic(as.matrix(t(iNKO)))
iNKO <- t(iNKO$result)


iAWT <- Astrocytes@assays$RNA@counts[,grepl('WT',Astrocytes$orig.ident)]
iAWT <- iAWT[rowMeans(iAWT) > 0.05,]
iAWT <- magic(as.matrix(t(iAWT)))
iAWT <- t(iAWT$result)

iAKO <- Astrocytes@assays$RNA@counts[,grepl('KO',Astrocytes$orig.ident)]
iAKO <- iAKO[rowMeans(iAKO) > 0.05,]
iAKO <- magic(as.matrix(t(iAKO)))
iAKO <- t(iAKO$result)

boxplot(list(iNWT['Camk4',],iNKO['Camk4',]))

boxplot(list(iNWT['Syngr1',],iNKO['Syngr1',]))
boxplot(list(iNWT['App',],iNKO['App',]))
boxplot(list(iAWT['Syp',], iAKO['Syp',]))
boxplot(list(iAWT['Apoe',],iAKO['Apoe',]))

Idents(Neurons) <- ifelse(grepl('WT',Neurons$orig.ident), 'WT', 'KO')
Idents(Astrocytes) <- ifelse(grepl('WT',Astrocytes$orig.ident), 'WT', 'KO')

DEN <- FindMarkers(Neurons, ident.1 = 'KO', ident.2 = 'WT')
DEA <- FindMarkers(Astrocytes, ident.1 = 'KO', ident.2 = 'WT')

X <- iNWT['Prnp',]
Y <- iAWT['Arc',]
nCells <- min(c(length(X), length(Y)))
set.seed(1)
dWT <- sapply(1:1000, function(B){
  sX <- sample(X, nCells, replace = TRUE)
  sY <- sample(Y, nCells, replace = TRUE)
  cor(sX,sY, method = 'sp')
})

X <- iNKO['Prnp',]
Y <- iAKO['Arc',]
nCells <- min(c(length(X), length(Y)))
set.seed(1)
dKO <- sapply(1:1000, function(B){
  sX <- sample(X, nCells, replace = TRUE)
  sY <- sample(Y, nCells, replace = TRUE)
  cor(sX,sY, method = 'sp')
})

plot(ecdf(dWT))
lines(ecdf(dKO), col = 'red')

load('Microglia.RData')
Idents(Microglia) <- ifelse(grepl('WT',Microglia$orig.ident), 'WT', 'KO')
DEM <- FindMarkers(Microglia, ident.1 = 'KO', ident.2 = 'WT', logfc.threshold = 0, test.use = 'MAST')
plot(DEM$avg_logFC, -log10(DEM$p_val))       
