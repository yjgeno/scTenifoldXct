library(RSpectra)
library(Matrix)
library(igraph)
library(scTenifoldNet)
load('Astrocytes.RData')
load('Neurons.RData')

GT <- read.csv('GT.txt', sep = '\t', stringsAsFactors = FALSE)
GT <- GT$Gene.name[GT$Gene.type %in% 'protein_coding']
onlyProtein <- function(X){
  X[rownames(X) %in% GT,]
}

A1 <- onlyProtein(Astrocytes@assays$RNA@counts[,grepl('WT',Astrocytes$orig.ident)])
A2 <- onlyProtein(Astrocytes@assays$RNA@counts[,grepl('KO',Astrocytes$orig.ident)])

N1 <- onlyProtein(Neurons@assays$RNA@counts[,grepl('WT',Neurons$orig.ident)])
N2 <- onlyProtein(Neurons@assays$RNA@counts[,grepl('KO',Neurons$orig.ident)])

A1 <- A1[rowMeans(A1 != 0) > 0.05,]
A2 <- A2[rowMeans(A2 != 0) > 0.05,]
N1 <- N1[rowMeans(N1 != 0) > 0.05,]
N2 <- N2[rowMeans(N2 != 0) > 0.05,]

X1 <- N1[1:1000,]
X2 <- N2[1:1000,]
Y1 <- A1[1:1000,]
Y2 <- A2[1:1000,]

scTenifoldXct <- function(X1,Y1, X2, Y2, nNet = 5){
  iNet <- function(X,Y){
    X <- X[rowMeans(X != 0) > 0.05,]
    Y <- Y[rowMeans(Y != 0) > 0.05,]
    
    maxCells <- max(c(ncol(X),ncol(Y)))
    
    tX <- X[,sample(seq_len(ncol(X)), maxCells, replace = TRUE)]
    tY <- Y[,sample(seq_len(ncol(Y)), maxCells, replace = TRUE)]
    
    tX <- (scale(t(tX)))
    tY <- (scale(t(tY)))
    
    SVD <- svds(tX, 5)$v
    
    iNet <- pbapply::pbapply(tY,2,function(y){
      score <- tX %*% SVD
      # Step 2: Regress the observed vector of outcomes on the selected principal components as covariates, using ordinary least squares regression to get a vector of estimated regression coefficients.
      score <-  Matrix::t(Matrix::t(score) / (apply(score, 2, function(X) {sqrt(sum(X ^ 2))}) ^ 2))
      # Step 3: Transform this vector back to the scale of the actual covariates, using the eigenvectors corresponding to the selected principal components to get the final PCR estimator for estimating the regression coefficients characterizing the original model.
      Beta <- colSums(y * score)
      Beta <- SVD %*% (Beta)
      return(Beta)
    })
    
    absN <- abs(iNet)
    iNet[absN < quantile(absN, 0.95)] <- 0
    iNet <- as(iNet, 'dgCMatrix')
    colnames(iNet) <- paste0('Y_', rownames(Y))
    rownames(iNet) <- paste0('X_', rownames(X))
    return(iNet)
  }
  i1 <- lapply(seq_len(nNet), function(X){iNet(X1,Y1)})
  i2 <- lapply(seq_len(nNet), function(X){iNet(X2,Y2)})
  
  save(i1,i2, file='tempI.RData')
  load('tempI.RData')
  
  gX <- unique(unlist(lapply(c(i1,i2), function(X){rownames(X)})))
  gY <- unique(unlist(lapply(c(i1,i2), function(X){colnames(X)})))
  
  i1 <- lapply(i1, function(X){
    O <- matrix(0, nrow = length(gX), ncol = length(gY))
    rownames(O) <- gX
    colnames(O) <- gY
    O[rownames(X),colnames(X)] <- as.matrix(X)
    O <- as(O, 'dgCMatrix')
    #O <- graph_from_incidence_matrix(O, weighted = TRUE)[]
    return(O)
  })
  
  i2 <- lapply(i2, function(X){
    O <- matrix(0, nrow = length(gX), ncol = length(gY))
    rownames(O) <- gX
    colnames(O) <- gY
    O[rownames(X),colnames(X)] <- as.matrix(X)
    O <- as(O, 'dgCMatrix')
    #O <- graph_from_incidence_matrix(O, weighted = TRUE)[]
    return(O)
  })
  
  t1 <- array(0, dim = c(length(gX), length(gY), nNet))
  t2 <- array(0, dim = c(length(gX), length(gY), nNet))
  
  for(i in seq_len(nNet)){
    t1[,,i] <- as.matrix(i1[[i]])
    t2[,,i] <- as.matrix(i2[[i]])
  }
  
  t1 <- rTensor::as.tensor(t1)
  t1 <- rTensor::cp(t1, 5, max_iter = 1e3)
  t1 <- t1$est@data
  
  t2 <- rTensor::as.tensor(t2)
  t2 <- rTensor::cp(t2, 5, max_iter = 1e3)
  t2 <- t2$est@data
  
  
  for(i in seq_len(nNet)[-1]){
    t1[,,1] <- t1[,,1] + t1[,,i]
    t2[,,1] <- t2[,,1] + t2[,,i]
  }
  
  t1 <- t1[,,1]
  t2 <- t2[,,1]
  
  t1 <- t1/nNet
  t2 <- t2/nNet
  
  t1 <- t1/max(abs(t1))
  t2 <- t2/max(abs(t2))
  
  rownames(t1) <- rownames(t2) <- gX
  colnames(t1) <- colnames(t2) <- gY
  
  i1 <- igraph::graph_from_incidence_matrix(t1, weighted = TRUE)
  i2 <- igraph::graph_from_incidence_matrix(t2, weighted = TRUE)
  
  i1 <- i1[]
  i2 <- i2[]
  
  i1 <- (i1 + t(i1))/2
  i2 <- (i2 + t(i2))/2
  
  MA <- manifoldAlignment(i1, i2)
  DR <- dRegulation(MA)
  
  O <- list(xNet = t1, yNet = t2, manifoldAlignment = MA, diffRegulation = DR)
  return(O)
}

testOut <- scTenifoldXct(Y1,X1,Y2,X2, nNet = 10)
DR <- dRegulation(testOut$manifoldAlignment[,1:30])
writeLines(gsub('X_|Y_','',DR$gene[DR$p.adj < 0.05]))

testOut <- scTenifoldXct(N1,A1,N2,A2)
save(testOut, file = 'na_Xct.RData')

# load('na_Xct.RData')
# library(scTenifoldNet)
# DR1 <- dRegulation(testOut$manifoldAlignment[,1:2])
#load('tigss_Test.RData')
#DR2 <- dRegulation(testOut$manifoldAlignment[,1:2])

# writeLines(gsub('X_|Y_','',testOut$diffRegulation$gene[testOut$diffRegulation$Z > 0]))
# 
# 
# o1 <- as(o1, 'dgCMatrix')
# o2 <- as(o2, 'dgCMatrix')
# 
# #save(o1,o2, file = 'O.RData')
# MA <- scTenifoldNet::manifoldAlignment(o1,o2)
# DR <- scTenifoldNet::dRegulation(MA[,1:5])
# 
# load('Otigss.RData')
# DR1 <- DR
# load('O.RData')
# DR2 <- DR
# 
D1 <- DR1$Z
names(D1) <- DR1$gene
#D2 <- DR2$Z
#names(D2) <- DR2$gene
# 
# gNames <- intersect(names(D1), names(D2))
# 
# cor(D1[gNames], D2[gNames], method = 'sp')
# plot(D1[gNames], D2[gNames])
# 
# boxplot(list(N1['Scg2',],N2['Scg2',]))
# 

library(fgsea)
BIOP <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=BioPlanet_2019')
GOBP <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2018')
load('O.RData')

drX <- DR[grepl('X_', DR$gene),]
zX <- drX$Z
names(zX) <- toupper(gsub('X_','',drX$gene))
set.seed(1)
eX <- fgseaMultilevel(GOBP, zX)
eX[eX$NES > 0 & eX$padj < 0.05,]

drY <- DR[grepl('Y_', DR$gene),]
zY <- drY$Z
names(zY) <- toupper(gsub('Y_','',drY$gene))
set.seed(1)
eY <- fgseaMultilevel(GOBP, zY)
eY[eY$NES > 0 & eY$padj < 0.05,]
# 
# 
# which.max(o1['X_Scd2',] - o2['X_Scd2',])

