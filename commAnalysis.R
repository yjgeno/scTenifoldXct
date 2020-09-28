library(fgsea)
library(Matrix)
library(scTenifoldNet)
library(igraph)
source('https://raw.githubusercontent.com/dosorio/utilities/master/idConvert/hsa2mmu_SYMBOL.R')
KEGG <- gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2019_Mouse')
BIOP <- gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=BioPlanet_2019')
GOBP <- gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2018')

load('na_Xct.RData')
DR <- scTenifoldNet::dRegulation(testOut$manifoldAlignment)
X <- DR[grepl('^X', DR$gene),]
Y <- DR[grepl('^Y', DR$gene),]

X$gene <- gsub('X_', '', X$gene)
Y$gene <- gsub('Y_', '', Y$gene)

X$FC <- (X$distance^2)/mean(X$distance^2)
Y$FC <- (Y$distance^2)/mean(Y$distance^2)

X$p.value <- pchisq(X$FC, df = 1, lower.tail = FALSE)
Y$p.value <- pchisq(Y$FC, df = 1, lower.tail = FALSE)

X$p.adj <- p.adjust(X$p.value)
Y$p.adj <- p.adjust(Y$p.value)
 
zX <- X$distance
names(zX) <- toupper(X$gene)
BC <- MASS::boxcox(zX~1)
zX <- zX ^ abs(BC$x[which.max(BC$y)])
plot(density(zX))

zY <- Y$distance
names(zY) <- toupper(Y$gene)
BC <- MASS::boxcox(zY~1)
zY <- zY ^ abs(BC$x[which.max(BC$y)])
plot(density(zX))


zX <- (zX - mean(zX))/sd(zX)
zY <- (zY - mean(zY))/sd(zY)

set.seed(1)
eX <- fgseaMultilevel(KEGG, zX)
set.seed(1)
eY <- fgseaMultilevel(KEGG, zY)

eX <- eX[eX$NES > 0 & eX$padj < 0.05,]
eX <- eX[order(eX$NES, decreasing = TRUE),]
eY <- eY[eY$NES > 0 & eY$padj < 0.05,]
eY <- eY[order(eY$NES, decreasing = TRUE),]


#writeLines(c(X$gene[X$p.value < 0.05], Y$gene[Y$p.value < 0.05]))

intersect(eX$pathway, eY$pathway)

load('Neurons.RData')
library(Seurat)
library(ggplot2)
Neurons$G <- ifelse(grepl('WT', Neurons$orig.ident), 'WT', 'KO')
Idents(Neurons) <- Neurons$G
# deNeurons <- FindMarkers(Neurons, ident.1 = 'WT', ident.2 = 'KO', logfc.threshold = 0)
# write.csv(deNeurons, 'deNeurons.csv')
deNeurons <- read.csv('deNeurons.csv', row.names = 1)
plot(deNeurons$avg_logFC, -log10(deNeurons$p_val))
FC <- deNeurons$avg_logFC
names(FC) <- toupper(rownames(deNeurons))
eNeurons <- fgseaMultilevel(KEGG, FC)
eNeurons <- eNeurons[eNeurons$padj < 0.05,]




load('Microglia.RData')
library(Seurat)
library(ggplot2)
Microglia$G <- ifelse(grepl('WT', Microglia$orig.ident), 'WT', 'KO')
Idents(Microglia) <- Microglia$G
# deMicroglia <- FindMarkers(Microglia, ident.1 = 'WT', ident.2 = 'KO', logfc.threshold = 0)
# write.csv(deMicroglia, 'deMicroglia.csv')
deMicroglia <- read.csv('deMicroglia.csv', row.names = 1)
plot(deMicroglia$avg_logFC, -log10(deMicroglia$p_val))
FC <- deMicroglia$avg_logFC
names(FC) <- toupper(rownames(deMicroglia))
eMicroglia <- fgseaMultilevel(KEGG, FC)
eMicroglia <- eMicroglia[eMicroglia$padj < 0.05,]

load('Astrocytes.RData')
library(Seurat)
library(ggplot2)
Astrocytes$G <- ifelse(grepl('WT', Astrocytes$orig.ident), 'WT', 'KO')
Idents(Astrocytes) <- Astrocytes$G
# deAstrocytes <- FindMarkers(Astrocytes, ident.1 = 'WT', ident.2 = 'KO', logfc.threshold = 0)
# write.csv(deAstrocytes, 'deAstrocytes.csv')
deAstrocytes <- read.csv('deAstrocytes.csv', row.names = 1)
plot(deAstrocytes$avg_logFC, -log10(deAstrocytes$p_val))
FC <- deAstrocytes$avg_logFC
names(FC) <- toupper(rownames(deAstrocytes))
eAstrocytes <- fgseaMultilevel(KEGG, FC)
eAstrocytes <- eAstrocytes[eAstrocytes$padj < 0.05,]

# deNeurons <- deNeurons[toupper(rownames(deNeurons)) %in% unlist(eX$leadingEdge),]
# DotPlot(Neurons, features = X$gene[1:20]) + coord_flip()
# 
# load('Microglia.RData')
# Microglia$G <- ifelse(grepl('WT', Microglia$orig.ident), 'WT', 'KO')
# Idents(Microglia) <- Microglia$G
# deMicroglia <- FindMarkers(Microglia, ident.1 = 'WT', ident.2 = 'KO', min.diff.pct = 0.05, logfc.threshold = 0.2)
# 
# 
# load('Astrocytes.RData')
# Astrocytes$G <- ifelse(grepl('WT', Astrocytes$orig.ident), 'WT', 'KO')
# Idents(Astrocytes) <- Astrocytes$G
# deAstrocytes <- FindMarkers(Astrocytes, ident.1 = 'WT', ident.2 = 'KO', min.diff.pct = 0.05, logfc.threshold = 0.15)
# deAstrocytes <- deAstrocytes[toupper(rownames(deAstrocytes)) %in% unlist(eY$leadingEdge),]

# load('Neurons.RData')
# N1 <- Neurons@assays$RNA@counts[,grepl('WT',Neurons$orig.ident)]
# N2 <- Neurons@assays$RNA@counts[,grepl('KO',Neurons$orig.ident)]
# N1 <- N1[rowMeans(N1!=0) > 0.05,]
# N2 <- N2[rowMeans(N2!=0) > 0.05,]
# library(Rmagic)
# N1 <- magic(as.matrix(t(N1)))
# N1 <- t(N1$result)
# N2 <- magic(as.matrix(t(N2)))
# N2 <- t(N2$result)

