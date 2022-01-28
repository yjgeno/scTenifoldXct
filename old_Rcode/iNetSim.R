# Simulating of a dataset following a negative binomial distribution with high sparcity (~67%)
nCells = 2000
nGenes = 100
set.seed(1)
X <- rnbinom(n = nGenes * nCells, size = 20, prob = 0.98)
X <- round(X)
X <- matrix(X, ncol = nCells)
rownames(X) <- c(paste0('ng', 1:90), paste0('mt-', 1:10))


iNet <- function(X, Y, 
                 nComp = 3,
                 scaleScores = TRUE,
                 symmetric = FALSE,
                 q = 0, verbose = TRUE) {
  
  pcCoefficients <- function(K) {
    # Selecting the gene to be regressed out
    y <- Y[,K]
    score <- Xi %*% coeff
    # Step 2: Regress the observed vector of outcomes on the selected principal components as covariates, using ordinary least squares regression to get a vector of estimated regression coefficients.
    score <-
      Matrix::t(Matrix::t(score) / (apply(score, 2, function(X) {
        sqrt(sum(X ^ 2))
      }) ^ 2))
    # Step 3: Transform this vector back to the scale of the actual covariates, using the eigenvectors corresponding to the selected principal components to get the final PCR estimator for estimating the regression coefficients characterizing the original model.
    Beta <- colSums(y * score)
    Beta <- coeff %*% (Beta)
    
    return(Beta)
  }
  
  X <- (scale(Matrix::t(X)))
  Y <- (scale(Matrix::t(Y)))
  
  # Identify the number of rows in the input matrix
  nX <- ncol(X)
  nY <- ncol(Y)
  
  nComp <- 5
  # Apply the principal component regression for each gene
  #if(verbose){
  Xi <- X
  # Step 1: Perform PCA on the observed covariates data matrix to obtain $n$ number of the principal components.
  coeff <- RSpectra::svds(Xi, nComp)$v
  B <- pbapply::pbsapply(seq_len(nY), pcCoefficients)  
  #} else {
  #  B <- sapply(seq_len(nY), pcCoefficients)  
  #}
  
  # Transposition of the Beta coefficient matrix
  B <- t(B)
  
  # Absolute values for scaling and filtering
  absB <- abs(B)
  
  # Scaling the output matrix
  #if (isTRUE(scaleScores)) {
  B <- (B / max(absB))
  #}
  
  # Filtering the output matrix
  q <- 0.95
  B[absB < quantile(absB, q)] <- 0
  
  rownames(B) <- paste0('Y_g',seq_len(ncol(Y)))
  colnames(B) <- paste0('X_g',seq_len(ncol(X)))
  
  return(B)
}

Y <- X[1:10,]
A <- iNet(X,Y)

set.seed(1)
Y[1,] <- runif(ncol(X))
B <- iNet(X,Y)

library(igraph)
library(Matrix)
ANet <- graph_from_incidence_matrix(A, weighted = TRUE)
BNet <- graph_from_incidence_matrix(B, weighted = TRUE)

ANet <- ANet[]
BNet <- BNet[]

ANet <- (ANet + t(ANet))/2
BNet <- (BNet + t(BNet))/2

oNet <- scTenifoldNet::manifoldAlignment(ANet,BNet, d = 3)
dRegulation <- function (manifoldOutput, eps = 1e-16){
  geneList <- rownames(manifoldOutput)
  geneList <- geneList[grepl("^X_", geneList)]
  geneList <- gsub("^X_", "", geneList)
  nGenes <- length(geneList)
  eGenes <- nrow(manifoldOutput)/2
  eGeneList <- rownames(manifoldOutput)
  eGeneList <- eGeneList[grepl("^Y_", eGeneList)]
  eGeneList <- gsub("^Y_", "", eGeneList)
  if (nGenes != eGenes) {
    stop("Number of identified and expected genes are not the same")
  }
  if (!all(eGeneList == geneList)) {
    stop("Genes are not ordered as expected. X_ genes should be followed by Y_ genes in the same order")
  }
  dMetric <- sapply(seq_len(nGenes), function(G) {
    X <- manifoldOutput[G, ]
    Y <- manifoldOutput[(G + nGenes), ]
    I <- rbind(X, Y)
    O <- dist(I)
    O <- as.numeric(O)
    return(O)
  })
    dMetric[dMetric < eps] <- 0
    lambdaValues <- seq(-2, 2, length.out = 1000)
    lambdaValues <- lambdaValues[lambdaValues != 0]
    if(any(dMetric != 0)){
      BC <- MASS::boxcox(dMetric[dMetric != 0] ~ 1, plot = FALSE, lambda = lambdaValues)
      BC <- BC$x[which.max(BC$y)]
      BC <- abs(BC)
    } else {
      BC <- 1
    }
    nD <- dMetric ^ BC
    Z <- scale(nD)
    E <- mean(dMetric^2)
    FC <- dMetric^2/E
    pValues <- pchisq(q = FC, df = 1, lower.tail = FALSE)
    pAdjusted <- p.adjust(pValues, method = "fdr")
    dOut <- data.frame(gene = geneList, distance = dMetric, Z = Z, 
                       FC = FC, p.value = pValues, p.adj = pAdjusted)
    dOut <- dOut[order(dOut$p.value), ]
    dOut <- as.data.frame.array(dOut)
    return(dOut)
}
oNet <- dRegulation(oNet)

ANet <- graph_from_incidence_matrix(A, weighted = TRUE)
BNet <- graph_from_incidence_matrix(B, weighted = TRUE)

ANet <- bipartite_projection(ANet)
BNet <- bipartite_projection(BNet)

AOut <- dRegulation(scTenifoldNet::manifoldAlignment(ANet[[1]][], BNet[[1]][], d = 3))
head(AOut)
BOut <- dRegulation(scTenifoldNet::manifoldAlignment(ANet[[2]][], BNet[[2]][], d = 3))
head(BOut)

library(UpSetR)
O <- list(iNet = oNet$gene[oNet$p.adj < 0.05], yNet = AOut$gene[AOut$p.value < 0.05], XNet = BOut$gene[BOut$p.adj < 0.05])
upset(fromList(O))

A <- union(ANet[[1]],BNet[[1]])
L <- layout.auto(A)
rownames(L) <- rownames(A[])
par(mfrow=c(1,2))
plot(ANet[[1]], layout=L[rownames(ANet[[1]][]),])     
plot(BNet[[1]], layout=L[rownames(BNet[[1]][]),])     

A <- union(ANet[[2]],BNet[[2]])
A <- A[rowSums(A[]) > 0,colSums(A[]) > 0]
A <- graph_from_incidence_matrix(A)
set.seed(1)
L <- layout.auto(A)
L <- scale(L)
L <- t(t(L)/apply(abs(L),2,max))
rownames(L) <- rownames(A[])

B1 <- ANet[[2]][]
B1 <- B1[rowSums(B1) > 0, colSums(B1) > 0]
B1 <- graph_from_incidence_matrix(B1, weighted = TRUE)
B2 <- BNet[[2]][]
B2 <- B2[rowSums(B2) > 0, colSums(B2) > 0]
B2 <- graph_from_incidence_matrix(B2, weighted = TRUE)

par(mfrow=c(1,2), mar=c(1,1,1,1))
xLim <- c(min(L[,1]), max(L[,1]))
yLim <- c(min(L[,2]), max(L[,2]))
gCol <- (hcl.colors(nrow(L), palette = 'RdGy'))
#gCol <- densCols(BOut$Z, colramp = function(X){hcl.colors(X, palette = 'Reds 2')})
#plot(seq_along(gCol), col = gCol, pch = 16)
names(gCol) <- BOut$gene[BOut$gene %in% rownames(L)]
plot(B1, layout=L[rownames(B1[]),], rescale = FALSE, xlim = xLim, ylim = yLim, vertex.color = gCol[names(V(B1))], vertex.label.family= 'Arial', vertex.label.color='black', vertex.label.cex= 1, vertex.frame.color = NA, vertex.size = 5)     
plot(B2, layout=L[rownames(B2[]),], rescale = FALSE, xlim = xLim, ylim = yLim, vertex.color = gCol[names(V(B2))], vertex.label.family= 'Arial', vertex.label.color='black', vertex.label.cex= 1, vertex.frame.color = NA, vertex.size = 5)     

