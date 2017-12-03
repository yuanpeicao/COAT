set.seed(0)
## Setup the parameters
n <- 200
p1 <- 50
p2 <- 100
p3 <- 200

## Setup the mean and covariance for normal distribution
Sigma1_true <- diag(p1)
Sigma2_true <- diag(p2)
Sigma3_true <- diag(p3)

## Setup the mean for Gamma distribution
mu1_g <- runif(p1,0,10)
mu2_g <- runif(p2,0,10)
mu3_g <- runif(p3,0,10)

## Generate 100 replications
corr1_n <- c()
corr2_n <- c()
corr3_n <- c()

corr1_g <- c()
corr2_g <- c()
corr3_g <- c()

corr1_clr_n <- c()
corr2_clr_n <- c()
corr3_clr_n <- c()

corr1_clr_g <- c()
corr2_clr_g <- c()
corr3_clr_g <- c()

corr1_logx_n <- c()
corr2_logx_n <- c()
corr3_logx_n <- c()

corr1_logx_g <- c()
corr2_logx_g <- c()
corr3_logx_g <- c()

corr1_logw_n <- c()
corr2_logw_n <- c()
corr3_logw_n <- c()

corr1_logw_g <- c()
corr2_logw_g <- c()
corr3_logw_g <- c()

N_replication <- 100
for (i in 1:N_replication) {
  
  # Generate the compositions from Normal distribution under 4 different parameters  
  
  log_w1_n <- MASS::mvrnorm(n, mu1_g, Sigma1_true)
  log_w2_n <- MASS::mvrnorm(n, mu2_g, Sigma2_true)
  log_w3_n <- MASS::mvrnorm(n, mu3_g, Sigma3_true)
  
  w1_n <- exp(log_w1_n)
  w2_n <- exp(log_w2_n)
  w3_n <- exp(log_w3_n)
  
  x1_n <- w1_n/rowSums(w1_n)
  x2_n <- w2_n/rowSums(w2_n)
  x3_n <- w3_n/rowSums(w3_n)
  
  # Generate the compositions from Gamma distribution under 4 different parameters 
  # log_w=U + mu; U is Gamma(10,1), mu is the mean
  log_w1_g <- matrix(c(rgamma(n*p1, shape=10, scale=1)),nrow=n)/sqrt(10) + matrix(c(rep(1,n)),ncol=1)%*%matrix(c(mu1_g),nrow=1)
  log_w2_g <- matrix(c(rgamma(n*p2, shape=10, scale=1)),nrow=n)/sqrt(10) + matrix(c(rep(1,n)),ncol=1)%*%matrix(c(mu2_g),nrow=1)
  log_w3_g <- matrix(c(rgamma(n*p3, shape=10, scale=1)),nrow=n)/sqrt(10) + matrix(c(rep(1,n)),ncol=1)%*%matrix(c(mu3_g),nrow=1)
  
  w1_g <- exp(log_w1_g)
  w2_g <- exp(log_w2_g)
  w3_g <- exp(log_w3_g)
  
  x1_g <- w1_g/rowSums(w1_g)
  x2_g <- w2_g/rowSums(w2_g)
  x3_g <- w3_g/rowSums(w3_g)
  
  ## Calculate the sample correlation
  sigma1_n <- cor(x1_n)
  sigma2_n <- cor(x2_n)
  sigma3_n <- cor(x3_n)
  
  sigma1_g <- cor(x1_g)
  sigma2_g <- cor(x2_g)
  sigma3_g <- cor(x3_g)
  
  corr1_n <- c(corr1_n,c(sigma1_n[lower.tri(sigma1_n)]))
  corr1_g <- c(corr1_g,c(sigma1_g[lower.tri(sigma1_g)]))
  
  corr2_n <- c(corr2_n,c(sigma2_n[lower.tri(sigma2_n)]))
  corr2_g <- c(corr2_g,c(sigma2_g[lower.tri(sigma2_g)]))
  
  corr3_n <- c(corr3_n,c(sigma3_n[lower.tri(sigma3_n)]))
  corr3_g <- c(corr3_g,c(sigma3_g[lower.tri(sigma3_g)]))
  
  ## Correlation on clrX
  clrX1_n <- log(x1_n)-matrix(rowMeans(log(x1_n)),ncol=1)%*%matrix(rep(1,p1),nrow=1)
  clrX2_n <- log(x2_n)-matrix(rowMeans(log(x2_n)),ncol=1)%*%matrix(rep(1,p2),nrow=1)
  clrX3_n <- log(x3_n)-matrix(rowMeans(log(x3_n)),ncol=1)%*%matrix(rep(1,p3),nrow=1)
  
  sigma1_clr_n <- cor(clrX1_n)
  sigma2_clr_n <- cor(clrX2_n)
  sigma3_clr_n <- cor(clrX3_n)
  
  clrX1_g <- log(x1_g)-matrix(rowMeans(log(x1_g)),ncol=1)%*%matrix(rep(1,p1),nrow=1)
  clrX2_g <- log(x2_g)-matrix(rowMeans(log(x2_g)),ncol=1)%*%matrix(rep(1,p2),nrow=1)
  clrX3_g <- log(x3_g)-matrix(rowMeans(log(x3_g)),ncol=1)%*%matrix(rep(1,p3),nrow=1)
  
  sigma1_clr_g <- cor(clrX1_g)
  sigma2_clr_g <- cor(clrX2_g)
  sigma3_clr_g <- cor(clrX3_g)
  
  corr1_clr_n <- c(corr1_clr_n,c(sigma1_clr_n[lower.tri(sigma1_clr_n)]))
  corr1_clr_g <- c(corr1_clr_g,c(sigma1_clr_g[lower.tri(sigma1_clr_g)]))
  
  corr2_clr_n <- c(corr2_clr_n,c(sigma2_clr_n[lower.tri(sigma2_clr_n)]))
  corr2_clr_g <- c(corr2_clr_g,c(sigma2_clr_g[lower.tri(sigma2_clr_g)]))
  
  corr3_clr_n <- c(corr3_clr_n,c(sigma3_clr_n[lower.tri(sigma3_clr_n)]))
  corr3_clr_g <- c(corr3_clr_g,c(sigma3_clr_g[lower.tri(sigma3_clr_g)]))
  
  ## Correlation on log X
  logX1_n <- log(x1_n)
  logX2_n <- log(x2_n)
  logX3_n <- log(x3_n)
  
  sigma1_logx_n <- cor(logX1_n)
  sigma2_logx_n <- cor(logX2_n)
  sigma3_logx_n <- cor(logX3_n)
  
  logX1_g <- log(x1_g)
  logX2_g <- log(x2_g)
  logX3_g <- log(x3_g)
  
  sigma1_logx_g <- cor(logX1_g)
  sigma2_logx_g <- cor(logX2_g)
  sigma3_logx_g <- cor(logX3_g)
  
  corr1_logx_n <- c(corr1_logx_n,c(sigma1_logx_n[lower.tri(sigma1_logx_n)]))
  corr1_logx_g <- c(corr1_logx_g,c(sigma1_logx_g[lower.tri(sigma1_logx_g)]))
  
  corr2_logx_n <- c(corr2_logx_n,c(sigma2_logx_n[lower.tri(sigma2_logx_n)]))
  corr2_logx_g <- c(corr2_logx_g,c(sigma2_logx_g[lower.tri(sigma2_logx_g)]))
  
  corr3_logx_n <- c(corr3_logx_n,c(sigma3_logx_n[lower.tri(sigma3_logx_n)]))
  corr3_logx_g <- c(corr3_logx_g,c(sigma3_logx_g[lower.tri(sigma3_logx_g)]))
  
  ## Correlation on log W
  sigma1_logw_n <- cor(log_w1_n)
  sigma2_logw_n <- cor(log_w2_n)
  sigma3_logw_n <- cor(log_w3_n)
  
  sigma1_logw_g <- cor(log_w1_g)
  sigma2_logw_g <- cor(log_w2_g)
  sigma3_logw_g <- cor(log_w3_g)
  
  corr1_logw_n <- c(corr1_logw_n,c(sigma1_logw_n[lower.tri(sigma1_logw_n)]))
  corr1_logw_g <- c(corr1_logw_g,c(sigma1_logw_g[lower.tri(sigma1_logw_g)]))
  
  corr2_logw_n <- c(corr2_logw_n,c(sigma2_logw_n[lower.tri(sigma2_logw_n)]))
  corr2_logw_g <- c(corr2_logw_g,c(sigma2_logw_g[lower.tri(sigma2_logw_g)]))
  
  corr3_logw_n <- c(corr3_logw_n,c(sigma3_logw_n[lower.tri(sigma3_logw_n)]))
  corr3_logw_g <- c(corr3_logw_g,c(sigma3_logw_g[lower.tri(sigma3_logw_g)]))
}

## Generate the boxplot
dat <- data.frame(Parameter=factor(rep(c(c(rep(50,length(corr1_n))),c(rep(100,length(corr2_n))),c(rep(200,length(corr3_n)))),8)), 
                  Correlation=c(corr1_n,corr2_n,corr3_n,corr1_clr_n,corr2_clr_n,corr3_clr_n,corr1_logx_n,corr2_logx_n,corr3_logx_n,corr1_logw_n,corr2_logw_n,corr3_logw_n,corr1_g,corr2_g,corr3_g,corr1_clr_g,corr2_clr_g,corr3_clr_g,corr1_logx_g,corr2_logx_g,corr3_logx_g,corr1_logw_g,corr2_logw_g,corr3_logw_g), 
                  Group=factor(rep(c(rep(c("X","clr","log X","log W"),each=length(corr1_n)+length(corr2_n)+length(corr3_n))),2)),
                  Model=factor(rep(c("Normal","Gamma"), each=4*(length(corr1_n)+length(corr2_n)+length(corr3_n)))))

dat$Model <- factor(dat$Model, levels=c("Normal","Gamma"))
dat$Parameter <- factor(dat$Parameter, levels=c(50,100,200))

require(ggplot2)
setEPS(width=8.4, height=4.2)
postscript(file="plot/boxplot.eps")

bp <- ggplot(dat, aes(x=Parameter, y=Correlation))
bp + geom_boxplot(aes(fill=Group), coef=10) +
  ylab("Sample correlation") + xlab(expression(italic("p"))) +
  facet_wrap(~ Model) +
  scale_fill_brewer(palette="RdBu", name="", breaks=c("clr","log W","log X","X"),
                    labels=c("clr",expression(bold("Y")),expression(paste("log ",bold("X"))),expression(bold("X")))) +
  theme_classic() +
  theme(axis.title=element_text(size=rel(1.2)),legend.text.align=0,
        strip.text=element_text(size=rel(1.2)),strip.background=element_blank())
dev.off()
