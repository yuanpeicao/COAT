require('latex2exp')
require('ggplot2')
setwd("~/coat")
source('simulation.R')

maxRowSumQ <- function(cov, q){
  return(val = max(rowSums((abs(cov))^q)))
}

p1 <- 50
p2 <- 100
p3 <- 150
p4 <- 200
p5 <- 250
p6 <- 300
p7 <- 350
p8 <- 400
p9 <- 450
p10 <- 500

nRep <- 5
q <- 1

################################################
# hub covariance
################################################

hub_val1 <- 0
hub_val2 <- 0
hub_val3 <- 0
hub_val4 <- 0
hub_val5 <- 0
hub_val6 <- 0
hub_val7 <- 0
hub_val8 <- 0
hub_val9 <- 0
hub_val10 <- 0

for (i in 1:nRep){
  hubCov1 <- hub(p1)
  val1 <- maxRowSumQ(hubCov1, q)
  hub_val1 <- hub_val1 + val1
}
hub_val1 <- hub_val1/nRep

for (i in 1:nRep){
  hubCov2 <- hub(p2)
  val2 <- maxRowSumQ(hubCov2, q)
  hub_val2 <- hub_val2 + val2
}
hub_val2 <- hub_val2/nRep

for (i in 1:nRep){
  hubCov3 <- hub(p3)
  val3 <- maxRowSumQ(hubCov3, q)
  hub_val3 <- hub_val3 + val3
}
hub_val3 <- hub_val3/nRep

for (i in 1:nRep){
  hubCov4 <- hub(p4)
  val4 <- maxRowSumQ(hubCov4, q)
  hub_val4 <- hub_val4 + val4
}
hub_val4 <- hub_val4/nRep

for (i in 1:nRep){
  hubCov5 <- hub(p5)
  val5 <- maxRowSumQ(hubCov5, q)
  hub_val5 <- hub_val5 + val5
}
hub_val5 <- hub_val5/nRep

for (i in 1:nRep){
  hubCov6 <- hub(p6)
  val6 <- maxRowSumQ(hubCov6, q)
  hub_val6 <- hub_val6 + val6
}
hub_val6 <- hub_val6/nRep

for (i in 1:nRep){
  hubCov7 <- hub(p7)
  val7 <- maxRowSumQ(hubCov7, q)
  hub_val7 <- hub_val7 + val7
}
hub_val7 <- hub_val7/nRep

for (i in 1:nRep){
  hubCov8 <- hub(p8)
  val8 <- maxRowSumQ(hubCov8, q)
  hub_val8 <- hub_val8 + val8
}
hub_val8 <- hub_val8/nRep

for (i in 1:nRep){
  hubCov9 <- hub(p9)
  val9 <- maxRowSumQ(hubCov9, q)
  hub_val9 <- hub_val9 + val9
}
hub_val9 <- hub_val9/nRep

for (i in 1:nRep){
  hubCov10 <- hub(p10)
  val10 <- maxRowSumQ(hubCov10, q)
  hub_val10 <- hub_val10 + val10
}
hub_val10 <- hub_val10/nRep

Y_hub <- c(hub_val1/p1, hub_val2/p2, hub_val3/p3, hub_val4/p4, hub_val5/p5, hub_val6/p6, hub_val7/p7, hub_val8/p8, hub_val9/p9, hub_val10/p10)

################################################
# block covariance
################################################

blk_val1 <- 0
blk_val2 <- 0
blk_val3 <- 0
blk_val4 <- 0
blk_val5 <- 0
blk_val6 <- 0
blk_val7 <- 0
blk_val8 <- 0
blk_val9 <- 0
blk_val10 <- 0

for (i in 1:nRep){
  blockCov1 <- block(p1)
  val1 <- maxRowSumQ(blockCov1, q)
  blk_val1 <- blk_val1 + val1
}
blk_val1 <- blk_val1/nRep

for (i in 1:nRep){
  blockCov2 <- block(p2)
  val2 <- maxRowSumQ(blockCov2, q)
  blk_val2 <- blk_val2 + val2
}
blk_val2 <- blk_val2/nRep

for (i in 1:nRep){
  blockCov3 <- block(p3)
  val3 <- maxRowSumQ(blockCov3, q)
  blk_val3 <- blk_val3 + val3
}
blk_val3 <- blk_val3/nRep

for (i in 1:nRep){
  blockCov4 <- block(p4)
  val4 <- maxRowSumQ(blockCov4, q)
  blk_val4 <- blk_val4 + val4
}
blk_val4 <- blk_val4/nRep

for (i in 1:nRep){
  blockCov5 <- block(p5)
  val5 <- maxRowSumQ(blockCov5, q)
  blk_val5 <- blk_val5 + val5
}
blk_val5 <- blk_val5/nRep

for (i in 1:nRep){
  blockCov6 <- block(p6)
  val6 <- maxRowSumQ(blockCov6, q)
  blk_val6 <- blk_val6 + val6
}
blk_val6 <- blk_val6/nRep

for (i in 1:nRep){
  blockCov7 <- block(p7)
  val7 <- maxRowSumQ(blockCov7, q)
  blk_val7 <- blk_val7 + val7
}
blk_val7 <- blk_val7/nRep

for (i in 1:nRep){
  blockCov8 <- block(p8)
  val8 <- maxRowSumQ(blockCov8, q)
  blk_val8 <- blk_val8 + val8
}
blk_val8 <- blk_val8/nRep

for (i in 1:nRep){
  blockCov9 <- block(p9)
  val9 <- maxRowSumQ(blockCov9, q)
  blk_val9 <- blk_val9 + val9
}
blk_val9 <- blk_val9/nRep

for (i in 1:nRep){
  blockCov10 <- block(p10)
  val10 <- maxRowSumQ(blockCov10, q)
  blk_val10 <- blk_val10 + val10
}
blk_val10 <- blk_val10/nRep

Y_blk <- c(blk_val1/p1, blk_val2/p2, blk_val3/p3, blk_val4/p4, blk_val5/p5, blk_val6/p6, blk_val7/p7, blk_val8/p8, blk_val9/p9, blk_val10/p10)

################################################
# sparse covariance
################################################

sp_val1 <- 0
sp_val2 <- 0
sp_val3 <- 0
sp_val4 <- 0
sp_val5 <- 0
sp_val6 <- 0
sp_val7 <- 0
sp_val8 <- 0
sp_val9 <- 0
sp_val10 <- 0

for (i in 1:nRep){
  spCov1 <- sparse(p1)
  val1 <- maxRowSumQ(spCov1, q)
  sp_val1 <- sp_val1 + val1
}
sp_val1 <- sp_val1/nRep

for (i in 1:nRep){
  spCov2 <- sparse(p2)
  val2 <- maxRowSumQ(spCov2, q)
  sp_val2 <- sp_val2 + val2
}
sp_val2 <- sp_val2/nRep

for (i in 1:nRep){
  spCov3 <- sparse(p3)
  val3 <- maxRowSumQ(spCov3, q)
  sp_val3 <- sp_val3 + val3
}
sp_val3 <- sp_val3/nRep

for (i in 1:nRep){
  spCov4 <- sparse(p4)
  val4 <- maxRowSumQ(spCov4, q)
  sp_val4 <- sp_val4 + val4
}
sp_val4 <- sp_val4/nRep

for (i in 1:nRep){
  spCov5 <- sparse(p5)
  val5 <- maxRowSumQ(spCov5, q)
  sp_val5 <- sp_val5 + val5
}
sp_val5 <- sp_val5/nRep

for (i in 1:nRep){
  spCov6 <- sparse(p6)
  val6 <- maxRowSumQ(spCov6, q)
  sp_val6 <- sp_val6 + val6
}
sp_val6 <- sp_val6/nRep

for (i in 1:nRep){
  spCov7 <- sparse(p7)
  val7 <- maxRowSumQ(spCov7, q)
  sp_val7 <- sp_val7 + val7
}
sp_val7 <- sp_val7/nRep

for (i in 1:nRep){
  spCov8 <- sparse(p8)
  val8 <- maxRowSumQ(spCov8, q)
  sp_val8 <- sp_val8 + val8
}
sp_val8 <- sp_val8/nRep

for (i in 1:nRep){
  spCov9 <- sparse(p9)
  val9 <- maxRowSumQ(spCov9, q)
  sp_val9 <- sp_val9 + val9
}
sp_val9 <- sp_val9/nRep

for (i in 1:nRep){
  spCov10 <- sparse(p10)
  val10 <- maxRowSumQ(spCov10, q)
  sp_val10 <- sp_val10 + val10
}
sp_val10 <- sp_val10/nRep

Y_sp <- c(sp_val1/p1, sp_val2/p2, sp_val3/p3, sp_val4/p4, sp_val5/p5, sp_val6/p6, sp_val7/p7, sp_val8/p8, sp_val9/p9, sp_val10/p10)

## create data frame
plot_data <- data.frame(hub = Y_hub, blk = Y_blk, sp = Y_sp, p_dim = c(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10))

setEPS(width = 4.5, height = 3)
postscript(file="q1.eps")
ggplot(plot_data, aes(p_dim)) + ggtitle("q = 1") + 
  geom_line(aes(y = hub))  + geom_point(aes(y = hub, shape = "model 2")) +
  geom_line(aes(y = blk)) + geom_point(aes(y = blk, shape = "model 3")) +
  geom_line(aes(y = sp)) + geom_point(aes(y = sp, shape = "model 4")) +
  ylab(TeX('s_0(p)/p')) + xlab("p") + scale_x_continuous(breaks=seq(50, 500, 50)) +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 20, face="bold"), legend.title = element_blank(), 
        legend.text = element_text(size = 15), axis.title = element_text(size = 15, face = "italic"),
        axis.text = element_text(size = 8))
dev.off()
