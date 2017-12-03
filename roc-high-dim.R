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
# Compute d1, df, auc and plot ROC for correlation matrix
#  Input:
#    dataCell1 ------ 2 x nRep cell, the first rows: composition X, the second rows: count W
#    dataCell2 ------ 2 x nRep cell, the first rows: composition X, the second rows: count W
#    dataCell3 ------ 2 x nRep cell, the first rows: composition X, the second rows: count W
#    mat.coat1 ------ p^2 x nRep matrix, the esimated correlation using COAT with soft thresholding
#     mat.ora1 ------ p^2  x nRep matrix, the esimated correlation using oracle procedure with soft thresholding
#    mat.coat2 ------ p^2 x nRep matrix, the esimated correlation using COAT with soft thresholding
#     mat.ora2 ------ p^2  x nRep matrix, the esimated correlation using oracle procedure with soft thresholding
#    mat.coat3 ------ p^2 x nRep matrix, the esimated correlation using COAT with soft thresholding
#     mat.ora3 ------ p^2  x nRep matrix, the esimated correlation using oracle procedure with soft thresholding
#     mat.rho1 ------ p^2 x 1 vector, the true correlation
#     mat.rho2 ------ p^2 x 1 vector, the true correlation
#     mat.rho3 ------ p^2 x 1 vector, the true correlation
#---------------------------------------------------------------------------------------------------------
dist.auc.roc.coat.ora.highD <- function(dataCell1, dataCell2, dataCell3, mat.coat1, mat.coat2, mat.coat3, mat.ora1, mat.ora2, mat.ora3,
                                  mat.rho1, mat.rho2, mat.rho3, main = "", eps = 1e-6) {
  R <- ncol(mat.coat1)
  
  # true covariance
  mat.rho1 <- matrix(rep(mat.rho1, R), ncol = R)
  mat.rho2 <- matrix(rep(mat.rho2, R), ncol = R)
  mat.rho3 <- matrix(rep(mat.rho3, R), ncol = R)

  
  # ROC
  p <- ncol(dataCell1[[1,1]])
  
  perf.coat1 <- calCoatROC(dataCell1, matrix(mat.rho1[,1], nrow = p), nPlotPoint = 76, nGrid = 75)
  perf.ora1 <- calOracleROC(dataCell1, matrix(mat.rho1[,1], nrow = p), nPlotPoint = 76, nGrid = 75)
  
  perf.coat2 <- calCoatROC(dataCell2, matrix(mat.rho2[,1], nrow = p), nPlotPoint = 76, nGrid = 75)
  perf.ora2 <- calOracleROC(dataCell2, matrix(mat.rho2[,1], nrow = p), nPlotPoint = 76, nGrid = 75)
  
  perf.coat3 <- calCoatROC(dataCell3, matrix(mat.rho3[,1], nrow = p), nPlotPoint = 76, nGrid = 75)
  perf.ora3 <- calOracleROC(dataCell3, matrix(mat.rho3[,1], nrow = p), nPlotPoint = 76, nGrid = 75)


  # setEPS(width=8, height=6)
  # postscript(file="plot/simulation/hub-Normal-p200.eps", family="Times")
  par(mar = c(4.5, 5, 3, 0.5))

  plot(perf.coat1$fprGrid, perf.coat1$tprGrid, ylab = "True positive rate", xlab ="False positive rate", cex.lab = 3, type = 'l',
       lty = 1, lwd = 3, xlim = c(0.0,1.0), ylim = c(0.0,1.0), mgp = c(3, 1, 0), cex.axis = 2)
  title(main = "Model 1", cex.main = 3)
  par(new=TRUE)
  plot(perf.ora1$fprGrid, perf.ora1$tprGrid, type = 'l', lty = 2, lwd = 3, axes = FALSE, xlab = "", ylab = "", xlim = c(0.0,1.0), ylim = c(0.0,1.0))
  par(new=TRUE)
  plot(perf.coat2$fprGrid, perf.coat2$tprGrid, type = 'l', lty = 3, lwd = 3, axes = FALSE, xlab = "", ylab = "", xlim = c(0.0,1.0), ylim = c(0.0,1.0))
  par(new=TRUE)
  plot(perf.ora2$fprGrid, perf.ora1$tprGrid, type = 'l', lty = 4, lwd = 3, axes = FALSE, xlab = "", ylab = "", xlim = c(0.0,1.0), ylim = c(0.0,1.0))
  par(new=TRUE)
  plot(perf.coat3$fprGrid, perf.coat3$tprGrid, type = 'l', lty = 5, lwd = 3, axes = FALSE, xlab = "", ylab = "", xlim = c(0.0,1.0), ylim = c(0.0,1.0))
  par(new=TRUE)
  plot(perf.ora3$fprGrid, perf.ora1$tprGrid, type = 'l', lty = 6, lwd = 3, axes = FALSE, xlab = "", ylab = "", xlim = c(0.0,1.0), ylim = c(0.0,1.0))
  par(new=TRUE)
  legend(x = "bottomright", bty = "n", cex = 2.5,
         legend = c("COAT, n = 500","Oracle, n = 500","COAT, n = 1000","Oracle, n = 1000","COAT, n = 2000","Oracle, n = 2000"), lty = c(1,2,3,4,5,6), lwd = c(3,3,3,3,3,3))
  abline(a = 0, b = 1, lty = 1, col = "grey")
  # dev.off()
  
  return(NULL)
}

#-----------------------------------------------------------------------------------------------------------------------------------------------

load("result/Hub-Normal-p500-n500.RData")
dataCell1 <- dataCell
mat.rho1 <- mat.rho
mat.coat1 <- res$mat.coat.s
mat.ora1 <- res$mat.ora.s

load("result/Hub-Normal-p500-n1000.RData")
dataCell2 <- dataCell
mat.rho2 <- mat.rho
mat.coat2 <- res$mat.coat.s
mat.ora2 <- res$mat.ora.s

load("result/Hub-Normal-p500-n2000.RData")
dataCell3 <- dataCell
mat.rho3 <- mat.rho
mat.coat3 <- res$mat.coat.s
mat.ora3 <- res$mat.ora.s

t <- dist.auc.roc.coat.ora.highD(dataCell1, dataCell2, dataCell3, mat.coat1, mat.coat2, mat.coat3, mat.ora1, mat.ora2, mat.ora3,
                                             mat.rho1, mat.rho2, mat.rho3)