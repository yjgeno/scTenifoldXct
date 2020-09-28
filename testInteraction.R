testInteraction <- function(inputNetwork){
  geneDB <- read.csv('https://raw.githubusercontent.com/Teichlab/cellphonedb-data/master/data/gene_input.csv', stringsAsFactors = FALSE)
  geneID <- geneDB$hgnc_symbol
  names(geneID) <- geneDB$uniprot
  
  complexDB <- read.csv('https://raw.githubusercontent.com/Teichlab/cellphonedb-data/master/data/sources/complex_curated.csv', stringsAsFactors = FALSE)
  complexDB$uniprot_1 <- geneID[complexDB$uniprot_1]
  complexDB$uniprot_2 <- geneID[complexDB$uniprot_2]
  complexDB$uniprot_3 <- geneID[complexDB$uniprot_3]
  complexDB$uniprot_4 <- geneID[complexDB$uniprot_4]
  complexDB$geneID <- apply(complexDB[,2:5],1,function(X){paste0(X[!is.na(X)], collapse = ';')})
  complexDB$geneID <- paste(complexDB$geneID, complexDB$stoichiometry, sep = ';')
  complexDB <- complexDB[complexDB$geneID != ';',]
  complexDB$geneID <- unlist(lapply(strsplit(complexDB$geneID, ';'), function(X){paste0(unique(X[X!='']), collapse = ';')}))
  complexID <- complexDB$geneID
  names(complexID) <- complexDB$complex_name
  
  allID <- c(geneID, complexID)
  
  intDB <- read.csv('https://raw.githubusercontent.com/Teichlab/cellphonedb-data/master/data/sources/interaction_curated.csv', stringsAsFactors = FALSE)
  intDB$partner_a <- allID[intDB$partner_a]
  intDB$partner_b <- allID[intDB$partner_b]
  
  intGene <- apply(intDB, 1, function(X){
    A <- unlist(strsplit(X[2], split = ';'))
    B <- unlist(strsplit(X[3], split = ';'))
    expand.grid(A,B)
  })
  
  intGene <- do.call(rbind.data.frame, intGene)
  intGene <- unique(intGene)
  intGene <- as.data.frame.array(intGene)
  intGene <- intGene[intGene[,1] != intGene[,2],]
  
  load('nm_Xct.RData')
  inputNetwork <- testOut$xNet
  
  rownames(inputNetwork) <- gsub('X_', '', toupper(rownames(inputNetwork)))
  colnames(inputNetwork) <- gsub('Y_', '', toupper(colnames(inputNetwork)))
  
  intGene <- intGene[intGene$Var1 %in% rownames(inputNetwork),]
  intGene <- intGene[intGene$Var2 %in% colnames(inputNetwork),]
  
  png('weightComparison.png', width = 2000, height = 2000, res = 300)
  inputNetwork <- abs(inputNetwork)
  par(mar=c(3,3,1,1), mgp = c(1.5,0.5,0))
  plot(ecdf(inputNetwork[intGene$Var1, intGene$Var2]), col = 'red', main = 'Cell-Cell Communication | Gene - Gene Association Weight by PCR', xlab = 'Abs(Weight)', ylab = 'Empirical Cumulative Density Function')
  lines(ecdf(inputNetwork))
  ks.test(c(as.matrix(inputNetwork[intGene$Var1, intGene$Var2])),c(as.matrix(inputNetwork)), alternative = 'less')
  legend('bottomright', legend = c('CellPhoneDB Interactions','Full Network by PCR'), bty = 'n', col = c('red', 'black'), lty = c(1,1))
  dev.off()
}