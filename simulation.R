require(MASS)
#-----------------------------------------------------------------------------------------------------------
#  Generate correlation matrix using hub model
#  Input:
#     p        ------ number of the column (dimension of the variable)
#     nHub     ------ number of the hubs
#     prHub    ------ probability of the connection bewteen different hubs
#     prCommon ------ probability of the connection bewteen common points
#     edgeStr  ------ edge strength
#  Output:
#     corr     ------ correlation matrix
#-----------------------------------------------------------------------------------------------------------
hub <- function(p, nHub = 3, edgeStr = 0.3, prHub = 0.7, prCommon = 0.2){
  sigma <- matrix(0, p, p)
  for (i in 1 : (nHub - 1)){
    for (j in (i + 1) : p){
      sigma[i, j] <- edgeStr * (runif(1) <= prHub) 
      if (runif(1) <= 0.5){
        sigma[i, j] = sigma[i, j] * (-1)
      }
      sigma[j, i] <- sigma[i, j]
    }
  }
  for (i in (nHub + 1) : (p - 1)){
    for (j in (i + 1) : p){
      sigma[i, j] <- edgeStr * (runif(1) <= prCommon)
      if (runif(1) <= 0.5){
        sigma[i, j] = sigma[i, j] * (-1)
      }
      sigma[j, i] <- sigma[i, j]
    }
  }
  eigenSys <- eigen(sigma)
  sigma <- (abs(min(eigenSys$values)) + 0.01) * diag(p) + sigma
  corr <- diag(diag(sigma)^(-0.5))%*%sigma%*%diag(diag(sigma)^(-0.5))
  return(corr)
}

#-----------------------------------------------------------------------------------------------------------
#  Generate correlation matrix using block model
#  Input:
#     p            ------ number of the column (dimension of the variable)
#     nBlock       ------ number of the blocks
#     prInBk       ------ probability of the connection inside the same block
#     prCrsBk      ------ probability of the connection bewteen the different block
#     edgeStrInBk  ------ edge strength in the same block
#     edgeStrCrsBk ------ edge strength accross different blocks
#  Output:
#     corr         ------ correlation matrix
#-----------------------------------------------------------------------------------------------------------
block <- function(p, nBlock = 10, prInBlock = 0.5, prCrsBk = 0.2, edgeStrInBk = 0.3, edgeStrCrsBk = 0.3){
  sigma <- matrix(0, p, p)
  blockId <- c(sample(1:p, size = p, replace = FALSE) %% nBlock + 1)
  for (i in 1 : (p - 1)){
    for (j in (i + 1) : p){
      if (blockId[i] == blockId[j]){
        sigma[i, j] <- edgeStrInBk * (runif(1) <= prInBlock)
        if (runif(1) <= 0.5){
          sigma[i, j] = sigma[i, j] * (-1)
        }
        sigma[j, i] <- sigma[i, j]
      }else{
        sigma[i, j] <- edgeStrCrsBk * (runif(1) <= prCrsBk)
        if (runif(1) <= 0.5){
          sigma[i, j] = sigma[i, j] * (-1)
        }
        sigma[j, i] <- sigma[i, j]
      }
    }
  }
  eigenSys <- eigen(sigma)
  sigma <- (abs(min(eigenSys$values)) + 0.01) * diag(p) + sigma
  corr <- diag(diag(sigma)^(-0.5))%*%sigma%*%diag(diag(sigma)^(-0.5))  
  return(corr)
}

#-----------------------------------------------------------------------------------------------------------
#  Generate correlation matrix using sparse model
#  Input:
#     p       ------ number of the column (dimension of the variable)
#     p1      ------ number of the column in first block
#     strLb   ------ edge strength lower bound in the first block
#     strUb   ------ edge strength upper bound in the first block
#     prZero  ------ probability of being 0 in the first block
#  Output:
#     corr    ------ correlation matrix
#-----------------------------------------------------------------------------------------------------------
sparse <- function(p, p1 = floor(3*p^0.5), strLb = 0.5, strUb = 1, prZero = 0.7){
  sigma1 <- matrix(0, p1, p1)
  for (i in 1 : (p1 - 1)){
    for (j in (i + 1) : p1){
      sigma1[i, j] <- runif(1, min = strLb, max = strUb) * (runif(1) > prZero)
      if (runif(1) <= 0.5){
        sigma1[i, j] <- sigma1[i, j] * (-1)
      }
      sigma1[j, i] <- sigma1[i, j]
    }
  }
  eigenSys <- eigen(sigma1)
  sigma1 <- (abs(min(eigenSys$values)) + 0.01) * diag(p1) + sigma1
  sigma1 <- diag(diag(sigma1)^(-0.5))%*%sigma1%*%diag(diag(sigma1)^(-0.5)) 
  sigma2 <- diag(p - p1)
  corr <- rbind(cbind(sigma1,matrix(0, p1, p - p1)), cbind(matrix(0, p - p1, p1), sigma2))
  return(corr)
}

#-----------------------------------------------------------------------------------------------------------
#  Generate composition from multivariate normal distribution (log W is multivariate normal)
#  Input:
#     n      ------ sample size
#     mu     ------ mean of log W
#     corr   ------ covariance of log W
#  Output:
#     x      ------ composition matrix
#     w      ------ count matrix
#-----------------------------------------------------------------------------------------------------------
generateCompNormal <- function(n, mu, corr){
  logW <- mvrnorm(n = n, mu, corr)
  W <- exp(logW)
  X <- W / rowSums(W)
  return(list(X = X, W = W))
}

#-----------------------------------------------------------------------------------------------------------
#  Generate composition from multivariate Gamma distribution (log W is multivariate Gamma)
#  Input:
#     n      ------ sample size
#     mu     ------ mean of log W
#     corr   ------ covariance of log W
#     shape  ------ shape parameter
#     scale  ------ scale parameter
#  Output:
#     x      ------ composition matrix
#     w      ------ count matrix
#-----------------------------------------------------------------------------------------------------------
generateCompGamma <- function(n, mu, corr, shape = 10, scale = 1 ){
  s <- svd(corr)
  p <- ncol(corr)
  T <- s$u %*% diag(s$d)^0.5
  Z <- matrix(rgamma(n * p, shape = shape, scale = scale), nrow = p, ncol = n, byrow = TRUE)
  Z <- Z / (shape / scale^2)
  logW <- t(T %*% Z + as.matrix(mu,ncol = 1) %*% matrix(1,1,n))
  W <- exp(logW)
  X <- W / rowSums(W)
  return(list(X = X, W = W))
}

#-----------------------------------------------------------------------------------------------------------
#  Calculate the true positive rate and false positive rate
#  Input:
#     sigmaTrue  ------- true correlation
#     sigmaHat   ------- estimated correlation
#     eps        ------- if the strength is small than eps, then we think the strength is zero
#  Output:
#     tpr        ------- true positive rate
#     fpr        ------- false positive rate
#-----------------------------------------------------------------------------------------------------------
calTprFpr <- function(sigmaTrue, sigmaHat, eps = 1e-10){
  indTrueZero <- abs(sigmaTrue) < eps
  indTrueNonzero <- abs(sigmaTrue) >= eps
  indTestNonzero <- abs(sigmaHat) >= eps
  nTrueSparse <- sum(indTrueNonzero & indTestNonzero)
  tpr <- nTrueSparse/sum(indTrueNonzero)
  nFalseSparse <- sum(indTrueZero & indTestNonzero)
  fpr <- nFalseSparse/sum(indTrueZero)
  return(list(tpr = tpr, fpr = fpr))
}

#-----------------------------------------------------------------------------------------------------------
#  Calculate the true positive rate and false positive rate along ROC curve for COAT
#  Input:
#     nGrid       -------- number of grid
#     nPlotPoint  -------- number of the points plotted on ROC curve
#     sigmaTrue   -------- true covariance
#     dataCell    -------- cell, each element includes count W and composition X
#     soft        -------- soft = 1: soft thresholding; soft = 0: hard thrsholding
#  Output:
#     tprGrid     -------- true positive rate along ROC curve
#     fprGrid     -------- false positive rate along ROC curve
#-----------------------------------------------------------------------------------------------------------
calCoatROC <- function(dataCell, sigmaTrue, nPlotPoint = 21, nGrid = 20, soft = 1){
  nRep <- ncol(dataCell)
  plotLength <- 1/(nPlotPoint - 1)
  p <- ncol(dataCell[[1,1]])
  TprMat <- matrix(0, nRep, nPlotPoint)
  FprMat <- matrix(0, nRep, nPlotPoint)
  for (i in 1:nRep){
    W <- dataCell[[2,i]]
    logW <- log(W)
    clrX <- logW - rowSums(logW) %*%matrix(1,1,p) / p
    gridInfo <- adaptThresholdRange(clrX)
#    grid <- gridInfo$lower + (gridInfo$upper - gridInfo$lower)*rep(1:nGrid)/nGrid
    grid <- gridInfo$lower * (gridInfo$upper / gridInfo$lower)^(rep(1:nGrid)/nGrid)
    TPR <- matrix(0, nRep, 1)
    FPR <- matrix(0, nRep, 1)
    for (j in 1:nGrid){
      sigmaHat <- adaptThreshold(gridInfo$cov,gridInfo$theta,grid[j],soft)
      corrHat <- diag(diag(sigmaHat)^(-0.5))%*%sigmaHat%*%diag(diag(sigmaHat)^(-0.5)) 
      tprFpr <- calTprFpr(sigmaTrue, corrHat)
      TPR[j] <- tprFpr$tpr
      FPR[j] <- tprFpr$fpr
    }
    for (k in 0:(nPlotPoint-1)){
      TprMat[i, k+1] <- max(TPR[which(FPR < (plotLength/2 + k*plotLength))])
    }
  }
  tprGrid <- colMeans(TprMat)
  fprGrid <- rep(0:(nPlotPoint-1))/(nPlotPoint-1)
  return(list(tprGrid = tprGrid, fprGrid = fprGrid))
}

#-----------------------------------------------------------------------------------------------------------
#  Calculate the true positive rate and false positive rate along ROC curve for Oracle method
#  Input:
#     nGrid       -------- number of grid
#     nPlotPoint  -------- number of the points plotted on ROC curve
#     sigmaTrue   -------- true covariance
#     dataCell    -------- cell, each element includes count W and composition X
#     soft        -------- soft = 1: soft thresholding; soft = 0: hard thrsholding
#  Output:
#     tprGrid     -------- true positive rate along ROC curve
#     fprGrid     -------- false positive rate along ROC curve
#-----------------------------------------------------------------------------------------------------------
calOracleROC <- function(dataCell, sigmaTrue, nPlotPoint = 21, nGrid = 20, soft = 1){
  nRep <- ncol(dataCell)
  plotLength <- 1/(nPlotPoint - 1)
  p <- ncol(dataCell[[1,1]])
  TprMat <- matrix(0, nRep, nPlotPoint)
  FprMat <- matrix(0, nRep, nPlotPoint)
  for (i in 1:nRep){
    W <- dataCell[[2,i]]
    logW <- log(W)
    gridInfo <- adaptThresholdRange(logW)
#    grid <- gridInfo$lower + (gridInfo$upper - gridInfo$lower)*rep(1:nGrid)/nGrid
    grid <- gridInfo$lower * (gridInfo$upper / gridInfo$lower)^(rep(1:nGrid)/nGrid)
    TPR <- matrix(0, nRep, 1)
    FPR <- matrix(0, nRep, 1)
    for (j in 1:nGrid){
      sigmaHat <- adaptThreshold(gridInfo$cov,gridInfo$theta,grid[j],soft)
      corrHat <- diag(diag(sigmaHat)^(-0.5))%*%sigmaHat%*%diag(diag(sigmaHat)^(-0.5)) 
      tprFpr <- calTprFpr(sigmaTrue, corrHat)
      TPR[j] <- tprFpr$tpr
      FPR[j] <- tprFpr$fpr
    }
    for (k in 0:(nPlotPoint-1)){
      TprMat[i, k+1] <- max(TPR[which(FPR < (plotLength/2 + k*plotLength))])
    }
  }
  tprGrid <- colMeans(TprMat)
  fprGrid <- rep(0:(nPlotPoint-1))/(nPlotPoint-1)
  return(list(tprGrid = tprGrid, fprGrid = fprGrid))
}

#-----------------------------------------------------------------------------------------------------------
#  Generate data cell
#  Input:
#         n         -------- sample size
#         modelCov  -------- modelCov = 1: hub cov; modelCov = 2: block cov; modelCov = 3: sparse cov
#         p1        -------- number of variables in setting 1
#         p2        -------- number of variables in setting 2
#         p3        -------- number of variables in setting 3
#         nRep      -------- number of replications
#  Output:
#     dataCellN1    -------- data cell (normal, n x p1)
#     dataCellN2    -------- data cell (normal, n x p2)
#     dataCellN3    -------- data cell (normal, n x p3)
#     dataCellG1    -------- data cell (gamma, n x p1)
#     dataCellG2    -------- data cell (gamma, n x p2)
#     dataCellG3    -------- data cell (gamma, n x p3)
#     sigma1        -------- true covariance matrix p1 x p1
#     sigma2        -------- true covariance matrix p2 x p2
#     sigma3        -------- true covariance matrix p3 x p3
#-----------------------------------------------------------------------------------------------------------
generateDataCell <- function(n, modelCov, p1 = 50, p2 = 100, p3 = 200, nRep = 100){
  dataCellN1 <- matrix(list(), 2, nRep)
  dataCellN2 <- matrix(list(), 2, nRep)
  dataCellN3 <- matrix(list(), 2, nRep)
  dataCellG1 <- matrix(list(), 2, nRep)
  dataCellG2 <- matrix(list(), 2, nRep)
  dataCellG3 <- matrix(list(), 2, nRep)
  
  # Hegerogeneous mean
  a = 0
  b = 10
  mu1 <- a + (b-a) * runif(p1,0,1) 
  mu2 <- a + (b-a) * runif(p2,0,1) 
  mu3 <- a + (b-a) * runif(p3,0,1) 
  
  if (modelCov == 1){
    # Create hub covariance
    sigma1 <- hub(p1)
    sigma2 <- hub(p2)
    sigma3 <- hub(p3)
  }
  if (modelCov == 2){
    # Create block covariance
    sigma1 <- block(p1)
    sigma2 <- block(p2)
    sigma3 <- block(p3)
  }
  if (modelCov == 3){
    # Create sparse covariance
    sigma1 <- sparse(p1)
    sigma2 <- sparse(p2)
    sigma3 <- sparse(p3)
  }

  for (i in 1:nRep){
    tmpN1 <- generateCompNormal(n, mu1, sigma1)
    dataCellN1[[1,i]] <- tmpN1$X
    dataCellN1[[2,i]] <- tmpN1$W
    tmpN2 <- generateCompNormal(n, mu2, sigma2)
    dataCellN2[[1,i]] <- tmpN2$X
    dataCellN2[[2,i]] <- tmpN2$W
    tmpN3 <- generateCompNormal(n, mu3, sigma3)
    dataCellN3[[1,i]] <- tmpN3$X
    dataCellN3[[2,i]] <- tmpN3$W
    
    tmpG1 <- generateCompGamma(n, mu1, sigma1)
    dataCellG1[[1,i]] <- tmpG1$X
    dataCellG1[[2,i]] <- tmpG1$W
    tmpG2 <- generateCompGamma(n, mu2, sigma2)
    dataCellG2[[1,i]] <- tmpG2$X
    dataCellG2[[2,i]] <- tmpG2$W
    tmpG3 <- generateCompGamma(n, mu3, sigma3)
    dataCellG3[[1,i]] <- tmpG3$X
    dataCellG3[[2,i]] <- tmpG3$W
  }
  return(list(dataCellN1 = dataCellN1, dataCellN2 = dataCellN2, dataCellN3 = dataCellN3, 
              dataCellG1 = dataCellG1, dataCellG2 = dataCellG2, dataCellG3 = dataCellG3, 
              sigma1 = sigma1, sigma2 = sigma2, sigma3 = sigma3))
}
