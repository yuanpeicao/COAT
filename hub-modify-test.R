require(ROCR)
require(glmnet)
require(MCMCpack)
require(parallel)
require(MASS)
require(corpcor)

source("R/cclasso_modify.R");
source("R/REBACCA_modify.R")
source("R/coat.R");
source("simulation.R")

#---------------------------------------------------------------------------------------------------------
#  Correlation estimation
#  Input:
#     dataCell ------ 2 x nRep Cell, the first rows: composition X, the second rows: count W
#  Output:
#   mat.coat.s ------ p^2 x nRep matrix, the esimated correlation using COAT with soft thresholding
#    mat.ora.s ------ p^2  x nRep matrix, the esimated correlation using oracle procedure with soft thresholding
#  time.coat.s ------ nRep x 1 vector, execution time using COAt with soft thresholding
#   time.ora.s ------ nRep x 1 vector, execution time using oracle procedure with soft thresholding
#---------------------------------------------------------------------------------------------------------
summaryCorrEst <- function(dataCell){
  nRep <- ncol(dataCell)
  p <- ncol(dataCell[[1,1]])
  mat.coat.s <- matrix(0, p*p, nRep)
  mat.ora.s <- matrix(0, p*p, nRep)
  
  time.coat.s <- vector(length = nRep)
  time.ora.s <- vector(length = nRep)
  
  for (i in 1:nRep){
    X <- dataCell[[1,i]]
    res.coat.s <- coat(X, soft = 1)
    corr.coat.s <- res.coat.s$corr
    mat.coat.s[,i] <- as.vector(corr.coat.s)
    time.coat.s[i] <- res.coat.s$time
    
    W <- dataCell[[2,i]]
    logW <- log(W)
    startTime <- proc.time()
    res.ora.s <- adaptThresoldCov(logW, soft = 1)
    exeTimeClass <- proc.time() - startTime
    res.ora.s.time <- as.numeric(exeTimeClass[3])    
    corr.ora.s <- res.ora.s$corr
    mat.ora.s[,i] <- as.vector(corr.ora.s)
    time.ora.s[i] <- res.ora.s.time

  }
  return(list(mat.coat.s = mat.coat.s, mat.ora.s = mat.ora.s,  
              time.coat.s = time.coat.s, time.ora.s = time.ora.s))
}
#---------------------------------------------------------------------------------------------------------
# Compute d1, df, auc and plot ROC for correlation matrix
#  Input:
#     dataCell ------ 2 x nRep cell, the first rows: composition X, the second rows: count W
#   mat.coat.s ------ p^2 x nRep matrix, the esimated correlation using COAT with soft thresholding
#    mat.ora.s ------ p^2  x nRep matrix, the esimated correlation using oracle procedure with soft thresholding
#      mat.rho ------ p^2 x 1 vector, the true correlation
#  time.coat.s ------ nRep x 1 vector, execution time using COAt with soft thresholding
#   time.ora.s ------ nRep x 1 vector, execution time using oracle procedure with soft thresholding
#---------------------------------------------------------------------------------------------------------
dist.auc.roc <- function(dataCell, mat.coat.s, mat.ora.s,
                         mat.rho, time.coat.s, time.ora.s, 
                         main = "", eps = 1e-6) {
  R <- ncol(mat.coat.s)
  res <- matrix(0, 2, 7)
  rownames(res) <- c("mean", "sd")
  colnames(res) <- c("d1", "df", "d2", "auc", "TPR", "FPR", "time")
  #
  mat.rho <- matrix(rep(mat.rho, R), ncol = R)
  dmat.coat.s <- mat.coat.s - mat.rho
  dmat.ora.s <- mat.ora.s - mat.rho
  
  # COAT with soft thresholding
  pred.coat.s <- prediction(abs(mat.coat.s), abs(mat.rho) > eps + 0)
  d1.coat.s <- calL1NormLoss(mat.coat.s, mat.rho)
  df.coat.s <- sqrt(colSums(dmat.coat.s^2))
  d2.coat.s <- calSpectNormLoss(mat.coat.s, mat.rho)
  auc.coat.s <- unlist(performance(pred.coat.s, measure = "auc")@y.values)
  supp.coat.s <- calMatTprFpr(mat.coat.s, mat.rho)
  tpr.coat.s <- supp.coat.s$tpr
  fpr.coat.s <- supp.coat.s$fpr
  d1.df.auc.coat.s <- cbind(d1.coat.s, df.coat.s, d2.coat.s, auc.coat.s, tpr.coat.s, fpr.coat.s, time.coat.s)
  res.coat.s <- apply(d1.df.auc.coat.s, 2, function(x) {c(mean(x), sd(x))})
  print("COAT, soft thresholding")
  print(res.coat.s)
  
  # Oracle with soft thresholding
  pred.ora.s <- prediction(abs(mat.ora.s), abs(mat.rho) > eps + 0)
  d1.ora.s <- calL1NormLoss(mat.ora.s, mat.rho)
  df.ora.s <- sqrt(colSums(dmat.ora.s^2))
  d2.ora.s <- calSpectNormLoss(mat.ora.s, mat.rho)
  auc.ora.s <- unlist(performance(pred.ora.s, measure = "auc")@y.values)
  supp.ora.s <- calMatTprFpr(mat.ora.s, mat.rho)
  tpr.ora.s <- supp.ora.s$tpr
  fpr.ora.s <- supp.ora.s$fpr
  d1.df.auc.ora.s <- cbind(d1.ora.s, df.ora.s, d2.ora.s, auc.ora.s, tpr.ora.s, fpr.ora.s, time.ora.s)
  res.ora.s <- apply(d1.df.auc.ora.s, 2, function(x) {c(mean(x), sd(x))})
  print("Oracle, soft thresholding")
  print(res.ora.s)
  
  # ROC
  p <- ncol(dataCell[[1,1]])
  perf.coat <- calCoatROC(dataCell, matrix(mat.rho[,1], nrow = p), nPlotPoint = 76, nGrid = 75)
  perf.ora <- calOracleROC(dataCell, matrix(mat.rho[,1], nrow = p), nPlotPoint = 76, nGrid = 75)

  setEPS(width=8, height=6)
  postscript(file="plot/hub-modify/hub-Normal-p200-test.eps", family="Times")
  par(mar = c(4.5, 5, 3, 0.5))

  plot(perf.coat$fprGrid, perf.coat$tprGrid, ylab = "True positive rate", xlab ="False positive rate", cex.lab = 3, type = 'l',
       lty = 1, lwd = 3, xlim = c(0.0,1.0), ylim = c(0.0,1.0), mgp = c(3, 1, 0), cex.axis = 2)
  title(main = substitute(paste(italic('p'), " = 200")), cex.main = 3)
  par(new=TRUE)
  plot(perf.ora$fprGrid, perf.ora$tprGrid, type = 'l', lty = 2, lwd = 3, axes = FALSE, xlab = "", ylab = "", xlim = c(0.0,1.0), ylim = c(0.0,1.0))
  legend(x = "bottomright", bty = "n", cex = 2.5,
         legend = c("COAT","Oracle"), lty = c(1,2), lwd = c(3,3))
  abline(a = 0, b = 1, lty = 1, col = "grey")
  dev.off()
  
  return(NULL)
}

#---------------------------------------------------------------------------------------------------------
# Compute spectral norm loss for correlation matrix
#  Input:
#       mat.est  ------ p^2 x nRep matrix, the esimated correlation
#       mat.rho  ------ p^2 x nRep matrix, the true correlation
#  Output:
#     spectLoss  ------ spectrum norm loss
#---------------------------------------------------------------------------------------------------------
calSpectNormLoss <- function(mat.est, mat.rho){
  p <- sqrt(nrow(mat.est))
  nRep <- ncol(mat.est)
  corr.true <- matrix(mat.rho[,1], nrow = p)
  spectLoss <- vector(length = nRep)
  for (i in 1:nRep){
    corr.est <- matrix(mat.est[,i],nrow = p)
    res.svd <- fast.svd(corr.est - corr.true)
    spectLoss[i] <- max(res.svd$d)
  }
  return(spectLoss)
}

#-----------------------------------------------------------------------------------------------------------
#  Calculate the true positive rate and false positive rate
#  Input:
#       mat.est  ------ p^2 x nRep matrix, the esimated correlation
#       mat.rho  ------ p^2 x nRep matrix, the true correlation
#  Output:
#     tpr        ------- true positive rate
#     fpr        ------- false positive rate
#-----------------------------------------------------------------------------------------------------------
calMatTprFpr <- function(mat.est, mat.rho){
  p <- sqrt(nrow(mat.est))
  nRep <- ncol(mat.est)
  corr.true <- matrix(mat.rho[,1], nrow = p)
  tpr <- vector(length = nRep)
  fpr <- vector(length = nRep)
  for (i in 1:nRep){
    corr.est <- matrix(mat.est[,i],nrow = p)
    res.supp <- calTprFpr(sigmaTrue = corr.true, sigmaHat = corr.est)
    tpr[i] <- res.supp$tpr
    fpr[i] <- res.supp$fpr
  }
  return(list(tpr = tpr, fpr = fpr))
}

#---------------------------------------------------------------------------------------------------------
# Compute matrix L-1 norm loss for correlation matrix
#  Input:
#       mat.est  ------ p^2 x nRep matrix, the esimated correlation
#       mat.rho  ------ p^2 x nRep matrix, the true correlation
#  Output:
#     L1Loss  ------ spectrum norm loss
#---------------------------------------------------------------------------------------------------------
calL1NormLoss <- function(mat.est, mat.rho){
  p <- sqrt(nrow(mat.est))
  nRep <- ncol(mat.est)
  corr.true <- matrix(mat.rho[,1], nrow = p)
  L1Loss <- vector(length = nRep)
  for (i in 1:nRep){
    corr.est <- matrix(mat.est[,i],nrow = p)
    L1Loss[i] <- norm(corr.est - corr.true, "o")
  }
  return(L1Loss)
}

#-----------------------------------------------------------------------------------------------------------
#  Generate single data cell
#  Input:
#         n         -------- sample size
#         modelCov  -------- modelCov = 1: hub cov; modelCov = 2: block cov; modelCov = 3: sparse cov
#         p         -------- number of variables
#         nRep      -------- number of replications
#  Output:
#     dataCellN     -------- data cell (normal, n x p)
#     dataCellG     -------- data cell (gamma, n x p)
#     sigma         -------- true covariance matrix p x p
#-----------------------------------------------------------------------------------------------------------
generateSingleDataCell <- function(n = 500, modelCov, p = 500, nRep = 100){
  dataCellN <- matrix(list(), 2, nRep)
  dataCellG <- matrix(list(), 2, nRep)
  
  # Hegerogeneous mean
  a = 0
  b = 10
  mu <- a + (b-a) * runif(p,0,1) 
  
  if (modelCov == 1){
    # Create hub covariance
    sigma <- hub(p)
  }
  if (modelCov == 2){
    # Create block covariance
    sigma <- block(p)
  }
  if (modelCov == 3){
    # Create sparse covariance
    sigma <- sparse(p)
  }
  
  for (i in 1:nRep){
    tmpN <- generateCompNormal(n, mu, sigma)
    dataCellN[[1,i]] <- tmpN$X
    dataCellN[[2,i]] <- tmpN$W
    
    tmpG<- generateCompGamma(n, mu, sigma)
    dataCellG[[1,i]] <- tmpG$X
    dataCellG[[2,i]] <- tmpG$W
    
  }
  return(list(dataCellN = dataCellN, dataCellG = dataCellG, sigma = sigma))
}

#---------------------------------------------------------------------------------------------------------
set.seed(1000)

# Generate the date: the first rows are the composition X and second rows are the count W; number of the column/replication: 100 
dataHubCov <- generateDataCell(n = 200, modelCov = 1, nRep = 3)
dataCell <- dataHubCov$dataCellN2
mat.rho <- as.vector(dataHubCov$sigma2)
remove(dataHubCov)
res <- summaryCorrEst(dataCell)
t <- dist.auc.roc(dataCell = dataCell, mat.coat.s = res$mat.coat.s, mat.ora.s = res$mat.ora.s,
                  mat.rho = mat.rho, time.coat.s = res$time.coat.s, time.ora.s = res$time.ora.s)


#save.image("plot/hub-modify/hub-Normal-p100.RData")

