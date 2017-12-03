require("qgraph")
require("RColorBrewer")
require(ggplot2)

lean.net.coat <- as.matrix(read.csv("plot/combo/coat-lean.csv"))
obese.net.coat <- as.matrix(read.csv("plot/combo/coat-obese.csv")) 

taxa <- matrix(unlist(strsplit(colnames(lean.net.coat), "\\.")), 6)
taxa <- matrix(unlist(strsplit(colnames(obese.net.coat), "\\.")), 6)
phyla <- taxa[2, ]
genera <- taxa[6, ]

lean.q.coat <- qgraph(lean.net.coat, fade = FALSE,layout="spring", groups=phyla, color=brewer.pal(nlevels(as.factor(phyla)), "Set3"), vsize=4, labels=genera,
                 label.cex=.7, label.scale=FALSE, legend.cex=.4, GLratio=6, mar=c(1,2,1,3), height=4.5, width=6, filetype="eps", filename="plot/coat-lean-net")
dev.off()
obese.q.coat <- qgraph(obese.net.coat, fade = FALSE, layout="spring", groups=phyla, color=brewer.pal(nlevels(as.factor(phyla)), "Set3"), vsize=4, labels=genera,
                  label.cex=.7, label.scale=FALSE, legend.cex=.4, GLratio=6, mar=c(1,2,1,3), height=4.5, width=6, filetype="eps", filename="plot/coat-obese-net")
dev.off()

p <- length(genera)

corr <- as.numeric(lean.net.coat[upper.tri(lean.net.coat)])
g1 <- matrix(genera, p, p, byrow=TRUE)[upper.tri(lean.net.coat)]
g2 <- matrix(genera, p, p)[upper.tri(lean.net.coat)]
abs.corr <-abs(corr) 
ord <- order(abs.corr, decreasing=TRUE)
lean.rank.coat <- data.frame(Correlation=corr[ord], Genus.1=g1[ord], Genus.2=g2[ord])
write.csv(lean.rank.coat, file="coat-lean-rank.csv", row.names=FALSE)

corr <- as.numeric(obese.net.coat[upper.tri(obese.net.coat)])
g1 <- matrix(genera, p, p, byrow=TRUE)[upper.tri(obese.net.coat)]
g2 <- matrix(genera, p, p)[upper.tri(obese.net.coat)]
abs.corr <-abs(corr) 
ord <- order(abs.corr, decreasing=TRUE)
obese.rank.coat <- data.frame(Correlation=corr[ord], Genus.1=g1[ord], Genus.2=g2[ord])
write.csv(obese.rank.coat, file="coat-obese-rank.csv", row.names=FALSE)

# heatmap

colnames(lean.net.coat) <- genera
colnames(obese.net.coat) <- genera

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

lean.net.coat2 <- reorder_cormat(lean.net.coat)
obese.net.coat2 <- reorder_cormat(obese.net.coat)

setEPS(width=9.6, height=7.2)
postscript(file="plot/combo/coat-lean-heatmap.eps", family="Times")
# par(mar = c(5, 1, 1, 0.5))
ggplot(data = melt(lean.net.coat2), aes(x = Var1, y = Var2, fill = value), mar=c(1,1,1,1)) +  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", name="") + theme_minimal() + 
  theme(axis.title.x = element_blank(), axis.text.x=  element_blank(),
        axis.ticks.x = element_blank(), axis.title.y = element_blank()) + coord_fixed()
dev.off()

setEPS(width=9.6, height=7.2)
postscript(file="plot/combo/coat-obese-heatmap.eps", family="Times")
# par(mar = c(5, 1, 1, 0.5))
ggplot(data = melt(obese.net.coat2), aes(x = Var1, y = Var2, fill = value), mar=c(1,1,1,1)) +  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", name="") + theme_minimal() + 
  theme(axis.title.x = element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), axis.title.y = element_blank()) + coord_fixed()
dev.off()

# network: power law decay

lean.net.coat.binary <- lean.net.coat -diag(p) > 0
lean.net.coat.connect <- rowSums(lean.net.coat.binary)
count.lean.net.coat <- unique(lean.net.coat.connect)
sum(lean.net.coat.connect == 13)
# x.lean.coat <- c(1,2,3,4,5,6,7,8,10,11,12,13,14,15)
# y.lean.coat <- c(2,4,3,4,7,2,2,2,2,1,4,2,4,1)
x.lean.coat <- c(1,2,3,4,5,6,7,8,9,11,12,13)
y.lean.coat <- c(4,5,6,5,2,2,1,2,3,3,1,5)

lean.coat.degree <- data.frame(degree = x.lean.coat, degreeP = y.lean.coat)
fit <- lm(log(degreeP)~log(degree), data = lean.coat.degree)
ggplot(lean.coat.degree, aes(x = log(degree), y = log(degreeP))) + geom_point() + stat_smooth(method = "lm", col = "red") +
  labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                     "Intercept =",signif(fit$coef[[1]],5 ),
                     " Slope =",signif(fit$coef[[2]], 5),
                     " P =",signif(summary(fit)$coef[2,4], 5)))


obese.net.coat.binary <- obese.net.coat - diag(p) > 0
obese.net.coat.connect <- rowSums(obese.net.coat.binary)
count.obese.net.coat <- unique(obese.net.coat.connect)
sum(obese.net.coat.connect == 1)
x.obese.coat <- c(1,2,3,4,5,6,7,8,9,11)
y.obese.coat <- c(4,7,9,5,3,1,1,1,2,1)
obese.coat.degree <- data.frame(degreeP = x.obese.coat, degree = y.obese.coat)
fit <- lm(log(degreeP)~log(degree), data = lean.coat.degree)
ggplot(obese.coat.degree, aes(x = log(degree), y = log(degreeP))) + geom_point() + stat_smooth(method = "lm", col = "red") +
  labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                     "Intercept =",signif(fit$coef[[1]],5 ),
                     " Slope =",signif(fit$coef[[2]], 5),
                     " P =",signif(summary(fit)$coef[2,4], 5)))