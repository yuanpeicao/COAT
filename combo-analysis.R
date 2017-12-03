require(ROCR)
require(glmnet)
require(MCMCpack)
require(parallel)
require(MASS)
require(corpcor)
require(tictoc)

source("R/cclasso_modify.R");
source("R/coat.R");
source("simulation.R")
source("R/REBACCA_modify.R")

tic()

set.seed(300)

eps <- 1e-100

# Load the data
df <- read.table("data/BMI_and_Counts.txt", header = TRUE)

# split the data into obese group (BMI >= 25) and lean group (BMI < 25)
df_obese <- df[which(df$BMI >= 25),]
df_lean <- df[which(df$BMI < 25),]

w_obese <- as.matrix(df_obese[,-c(1,2)])
w_lean <- as.matrix(df_lean[,-c(1,2)])

# delete rare taxa
w_obese.binary <- w_obese > 0
w_lean.binary <- w_lean > 0
freq_obese <- colSums(w_obese.binary)
freq_lean <- colSums(w_lean.binary)
w_obese <- w_obese[, -which(freq_obese < 4 | freq_lean < 4)]
w_lean <- w_lean[, -which(freq_obese < 4 | freq_lean < 4)]

# replace zero counts by 0.5
w_obese[which(w_obese == 0)] <- 0.5
w_lean[which(w_lean == 0)] <- 0.5

# convert the count into composition
p <- ncol(w_obese)
x_obese <- w_obese/(rowSums(w_obese)%*%matrix(1,1,p))
x_lean <- w_lean/(rowSums(w_lean)%*%matrix(1,1,p))

# COAT
res.coat.obese <- coat(x_obese, nFoler = 10)
res.coat.lean <- coat(x_lean, nFoler = 10)
coat.ind_true_nonzero.obese <- abs(res.coat.obese$corr - diag(diag(res.coat.obese$corr))) >= eps
coat.ind_true_nonzero.lean <- abs(res.coat.lean$corr - diag(diag(res.coat.lean$corr))) >= eps

# # CCLasso
# res.ccl.obese <- cclasso(x_obese)
# res.ccl.lean <- cclasso(x_lean)
# ccl.ind_true_nonzero.obese <- abs(res.ccl.obese$cor_w - diag(p)) >= eps
# ccl.ind_true_nonzero.lean <- abs(res.ccl.lean$cor_w - diag(p)) >= eps
# 
# # REBECCA
# res.reb.obese <- rebacca_main(x_obese)
# res.reb.lean <- rebacca_main(x_lean)
# reb.ind_true_nonzero.obese <- abs(res.reb.obese$corr - diag(diag(res.reb.obese$corr))) >= eps
# reb.ind_true_nonzero.lean <- abs(res.reb.lean$corr - diag(diag(res.reb.lean$corr))) >= eps

#res.reb.obese <- rebacca(t(x_obese), nbootstrap = 1)
#res.reb.obese.tau <- stability_cutoff(res.reb.obese$Stability, res.reb.obese$q, B = 50, FWER = 0.05)
#res.reb.obese.adj <- sscore2adjmatrix(res.reb.obese$Stability, res.reb.obese.tau)
#res.reb.obese.est <- rebacca_adjm2corr(t(x_obese), res.reb.obese.adj)
#res.reb.lean <- rebacca(t(x_lean), nbootstrap = 1)
#res.reb.lean.tau <- stability_cutoff(res.reb.lean$Stability, res.reb.lean$q, B = 50, FWER = 0.05)
#res.reb.lean.adj <- sscore2adjmatrix(res.reb.lean$Stability, res.reb.lean.tau)
#res.reb.lean.est <- rebacca_adjm2corr(t(x_lean), res.reb.lean.adj)
#reb.ind_true_nonzero.obese <- abs(res.reb.obese.est$corr - diag(diag(res.reb.obese.est$corr))) >= eps
#reb.ind_true_nonzero.lean <- abs(res.reb.lean.est$corr - diag(diag(res.reb.lean.est$corr))) >= eps

# bootstraping stability
nRep <- 100
n_obese <- nrow(x_obese)
n_lean <- nrow(x_lean)
coat.reprod_obese <- vector(length = nRep)
coat.reprod_lean <- vector(length = nRep)
# ccl.reprod_obese <- vector(length = nRep)
# ccl.reprod_lean <- vector(length = nRep)
# reb.reprod_obese <- vector(length = nRep)
# reb.reprod_lean <- vector(length = nRep)

coat.corr.obese.boot.binary <- matrix(0, p, p)
coat.corr.lean.boot.binary <- matrix(0, p, p)
# ccl.corr.obese.boot.binary <- matrix(0, p, p)
# ccl.corr.lean.boot.binary <- matrix(0, p, p)
# reb.corr.obese.boot.binary <- matrix(0, p, p)
# reb.corr.lean.boot.binary <- matrix(0, p, p)

coat.time.obese <- 0
coat.time.lean <- 0
# ccl.time.obese <- 0
# ccl.time.lean <- 0
# reb.time.obese <- 0
# reb.time.lean <- 0

for (i in 1:nRep){
  # bootstrap resampling
  x_obese.boot <- x_obese[sample(n_obese, n_obese, replace = TRUE),]
  x_lean.boot <- x_lean[sample(n_lean, n_lean, replace = TRUE),]
  
  # # CCLasso  
  # res.ccl.obese.boot <-  cclasso(x_obese.boot)
  # ccl.ind_nonzero.obese <- abs(res.ccl.obese.boot$cor_w - diag(diag(res.ccl.obese.boot$cor_w))) >= eps
  # ccl.nsparse.obese <- sum(ccl.ind_true_nonzero.obese & ccl.ind_nonzero.obese)
  # ccl.reprod_obese[i] <- ccl.nsparse.obese / sum(ccl.ind_true_nonzero.obese)
  # 
  # res.ccl.lean.boot <-  cclasso(x_lean.boot)
  # ccl.ind_nonzero.lean <- abs(res.ccl.lean.boot$cor_w - diag(diag(res.ccl.lean.boot$cor_w))) >= eps
  # ccl.nsparse.lean <- sum(ccl.ind_true_nonzero.lean & ccl.ind_nonzero.lean)
  # ccl.reprod_lean[i] <- ccl.nsparse.lean / sum(ccl.ind_true_nonzero.lean)
  
  # while (is.na(ccl.reprod_obese[i]) | is.na(ccl.reprod_lean[i])){
  #   x_obese.boot <- x_obese[sample(n_obese, n_obese, replace = TRUE),]
  #   x_lean.boot <- x_lean[sample(n_lean, n_lean, replace = TRUE),]
  #   
  #   res.ccl.obese.boot <-  cclasso(x_obese.boot)
  #   ccl.ind_nonzero.obese <- abs(res.ccl.obese.boot$cor_w - diag(diag(res.ccl.obese.boot$cor_w))) >= eps
  #   ccl.nsparse.obese <- sum(ccl.ind_true_nonzero.obese & ccl.ind_nonzero.obese)
  #   ccl.reprod_obese[i] <- ccl.nsparse.obese / sum(ccl.ind_true_nonzero.obese)
  #   
  #   res.ccl.lean.boot <-  cclasso(x_lean.boot)
  #   ccl.ind_nonzero.lean <- abs(res.ccl.lean.boot$cor_w - diag(diag(res.ccl.lean.boot$cor_w))) >= eps
  #   ccl.nsparse.lean <- sum(ccl.ind_true_nonzero.lean & ccl.ind_nonzero.lean)
  #   ccl.reprod_lean[i] <- ccl.nsparse.lean / sum(ccl.ind_true_nonzero.lean)
  # }
  # ccl.corr.obese.boot.binary <- ccl.corr.obese.boot.binary + ccl.ind_nonzero.obese + diag(p)
  # ccl.corr.lean.boot.binary <- ccl.corr.lean.boot.binary + ccl.ind_nonzero.lean + diag(p)
  # ccl.time.obese <- ccl.time.obese + res.ccl.obese.boot$time
  # ccl.time.lean <- ccl.time.lean + res.ccl.lean.boot$time
  
  # COAT  
  res.coat.obese.boot <-  coat(x_obese.boot, nFoler = 10)
  coat.ind_nonzero.obese <- abs(res.coat.obese.boot$corr - diag(diag(res.coat.obese.boot$corr))) >= eps
  coat.nsparse.obese <- sum(coat.ind_true_nonzero.obese & coat.ind_nonzero.obese)
  coat.reprod_obese[i] <- coat.nsparse.obese / sum(coat.ind_true_nonzero.obese)
  
  res.coat.lean.boot <-  coat(x_lean.boot, nFoler = 10)
  coat.ind_nonzero.lean <- abs(res.coat.lean.boot$corr - diag(diag(res.coat.lean.boot$corr))) >= eps
  coat.nsparse.lean <- sum(coat.ind_true_nonzero.lean & coat.ind_nonzero.lean)
  coat.reprod_lean[i] <- coat.nsparse.lean / sum(coat.ind_true_nonzero.lean)
  
  coat.corr.obese.boot.binary <- coat.corr.obese.boot.binary + coat.ind_nonzero.obese + diag(p)
  coat.corr.lean.boot.binary <- coat.corr.lean.boot.binary + coat.ind_nonzero.lean + diag(p)
  coat.time.obese <- coat.time.obese + res.coat.obese.boot$time
  coat.time.lean <- coat.time.lean + res.coat.lean.boot$time
  
  # REBECCA
#  res.reb.obese.boot <-  rebacca(t(x_obese.boot), nbootstrap = 1)
#  res.reb.obese.boot.tau <- stability_cutoff(res.reb.obese.boot$Stability, res.reb.obese.boot$q, B = 50, FWER = 0.05)
#  res.reb.obese.boot.adj <- sscore2adjmatrix(res.reb.obese.boot$Stability, res.reb.obese.boot.tau)
#  res.reb.obese.boot.est <- rebacca_adjm2corr(t(x_obese), res.reb.obese.boot.adj)
#  reb.ind_nonzero.obese <- abs(res.reb.obese.boot.est$corr - diag(diag(res.reb.obese.boot.est$corr))) >= eps

  # res.reb.obese.boot <- rebacca_main(x_obese.boot)
  # reb.ind_nonzero.obese <- abs(res.reb.obese.boot$corr - diag(diag(res.reb.obese.boot$corr))) >= eps
  # reb.nsparse.obese <- sum(reb.ind_true_nonzero.obese & reb.ind_nonzero.obese)
  # reb.reprod_obese[i] <- reb.nsparse.obese / sum(reb.ind_true_nonzero.obese)

#  res.reb.lean.boot <-  rebacca(t(x_lean.boot), nbootstrap = 1)
#  res.reb.lean.boot.tau <- stability_cutoff(res.reb.lean.boot$Stability, res.reb.lean.boot$q, B = 50, FWER = 0.05)
#  res.reb.lean.boot.adj <- sscore2adjmatrix(res.reb.lean.boot$Stability, res.reb.lean.boot.tau)
#  res.reb.lean.boot.est <- rebacca_adjm2corr(t(x_lean), res.reb.lean.boot.adj)
#  reb.ind_nonzero.lean <- abs(res.reb.lean.boot.est$corr - diag(diag(res.reb.lean.boot.est$corr))) >= eps  

  # res.reb.lean.boot <- rebacca_main(x_lean.boot)
  # reb.ind_nonzero.lean <- abs(res.reb.lean.boot$corr - diag(diag(res.reb.lean.boot$corr))) >= eps
  # reb.nsparse.lean <- sum(reb.ind_true_nonzero.lean & reb.ind_nonzero.lean)
  # reb.reprod_lean[i] <- reb.nsparse.lean / sum(reb.ind_true_nonzero.lean)
  # 
  # reb.corr.obese.boot.binary <- reb.corr.obese.boot.binary + reb.ind_nonzero.obese + diag(p)
  # reb.corr.lean.boot.binary <- reb.corr.lean.boot.binary + reb.ind_nonzero.lean + diag(p)
  # reb.time.obese <- reb.time.obese + res.reb.obese.boot$time
  # reb.time.lean <- reb.time.lean + res.reb.lean.boot$time
  
}
stab.coat.obese <- mean(coat.reprod_obese)
stab.coat.lean <- mean(coat.reprod_lean)
# stab.ccl.obese <- mean(ccl.reprod_obese)
# stab.ccl.lean <- mean(ccl.reprod_lean)
# stab.reb.obese <- mean(reb.reprod_obese)
# stab.reb.lean <- mean(reb.reprod_lean)

coat.corr.final.obese <- res.coat.obese$corr
coat.corr.final.obese[which(coat.corr.obese.boot.binary < 85)] <- 0
coat.corr.final.lean <- res.coat.lean$corr
coat.corr.final.lean[which(coat.corr.lean.boot.binary < 85)] <- 0

# ccl.corr.final.obese <- res.ccl.obese$cor_w
# ccl.corr.final.obese[which(ccl.corr.obese.boot.binary < 0.8 * nRep)] <- 0
# ccl.corr.final.lean <- res.ccl.lean$cor_w
# ccl.corr.final.lean[which(ccl.corr.lean.boot.binary < 0.8 * nRep)] <- 0
# 
# reb.corr.final.obese <- res.reb.obese.boot$corr
# reb.corr.final.obese[which(reb.corr.obese.boot.binary < 0.8 * nRep)] <- 0
# reb.corr.final.lean <- res.reb.lean.boot$corr
# reb.corr.final.lean[which(reb.corr.lean.boot.binary < 0.8 * nRep)] <- 0

# record the number of edges and average execution time
# COAT
coat.edgeP.obese <- length(which(coat.corr.final.obese - diag(diag(coat.corr.final.obese)) >= eps))/2
coat.edgeP.lean <- length(which(coat.corr.final.lean - diag(diag(coat.corr.final.lean)) >= eps))/2
coat.edgeN.obese <- length(which(coat.corr.final.obese - diag(diag(coat.corr.final.obese)) <= -eps))/2
coat.edgeN.lean <- length(which(coat.corr.final.lean - diag(diag(coat.corr.final.lean)) <= -eps))/2
coat.avg.time.obese <- coat.time.obese / nRep
coat.avg.time.lean <- coat.time.lean / nRep

# # CCLasso
# ccl.edgeP.obese <- length(which(ccl.corr.final.obese - diag(diag(ccl.corr.final.obese)) >= eps))/2
# ccl.edgeP.lean <- length(which(ccl.corr.final.lean - diag(diag(ccl.corr.final.lean)) >= eps))/2
# ccl.edgeN.obese <- length(which(ccl.corr.final.obese - diag(diag(ccl.corr.final.obese)) <= -eps))/2
# ccl.edgeN.lean <- length(which(ccl.corr.final.lean - diag(diag(ccl.corr.final.lean)) <= -eps))/2
# ccl.avg.time.obese <- ccl.time.obese / nRep
# ccl.avg.time.lean <- ccl.time.lean / nRep
# 
# # REBECCA
# reb.edgeP.obese <- length(which(reb.corr.final.obese - diag(diag(reb.corr.final.obese)) >= eps))/2
# reb.edgeP.lean <- length(which(reb.corr.final.lean - diag(diag(reb.corr.final.lean)) >= eps))/2
# reb.edgeN.obese <- length(which(reb.corr.final.obese - diag(diag(reb.corr.final.obese)) <= -eps))/2
# reb.edgeN.lean <- length(which(reb.corr.final.lean - diag(diag(reb.corr.final.lean)) <= -eps))/2
# reb.avg.time.obese <- reb.time.obese / nRep
# reb.avg.time.lean <- reb.time.lean / nRep

toc()

# write.csv(coat.corr.final.lean, file = "coat-lean.csv")
# write.csv(coat.corr.final.obese, file = "coat-obese.csv")

lean.net.coat <- coat.corr.final.lean
obese.net.coat <- coat.corr.final.obese