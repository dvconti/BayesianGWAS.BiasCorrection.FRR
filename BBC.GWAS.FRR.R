# Implements a Bayesian Weighted Correction
# This correction is similar to Zhong and Prentice's Weighted Average

# This code takes MLE for SNPs (betas and se.betas)
# and performs a hiearchical model for bias reduction with expected mean equal to a corrected beta (Zhong and Prentice)
# the code also calcualtes the proportion of FRR explained by a set of SNPS
# incorporating:
# the uncertainty of the SNP estimates
# bias reduction
# uncertainty in the pre-specified FRR
# David Conti

# set library and location
library(R2jags)
setwd("/Users/davidconti/Google Drive/BayesianBiasCorrection/BayesianGWAS.BiasCorrection.FRR")

# input data
beta.hat <- seq(from=-0.1, to=0.1, by=0.005)
M <- length(beta.hat)
sigma.hat <- rep(0.02, M)
p <- rep(0.15, M)
novel <- rbinom(M, 1, .5) + 1  # novel indexed with 2, known with 1
#novel <- rep(2, M)  # all treated novel
alpha.level <- 5e-08  # alpha level for determining significance - thus introducing bias

# include SNPs in FRR calculation
include.FRR <- rep(1, M)
include.FRR.known <- ifelse(novel==1,1,0)

# specify FRR and uncertainty on known FRR
lambda0.m <- 2.0
se0 <- 1/7
#c(lambda0.m-(1.96*se0),(lambda0.m+ 1.96*se0))


#Calculate needed quantities:
prec.beta <- 1/(sigma.hat^2)
q <- 1-p
zeros <- rep(0, M)
lambda0.prec <- 1/se0^2

# Estimate corrected beta via optimization of conditional expectation from Eq. 2.2 Zhong and Prentice 
c.stat <- qnorm(1-alpha.level/2)
phi <- function(x) { ((2*pi)^-.5)*exp(-(x^2)/2) }
conditionalExp <- function(theta, beta.hat.c, sigma.c)  {
  z <- theta/sigma.c
  phi.plus <- phi(z-c.stat)
  phi.neg <- phi(-z-c.stat)
  beta.e <- theta + sigma.c*((phi.plus-phi.neg)/(pnorm((z-c.stat))+pnorm((-z-c.stat))))
  (beta.e-beta.hat.c)^2
}
beta.cor <- unlist(lapply(1:length(beta.hat), FUN=function(m) { 
  r <- optim(beta.hat[m], conditionalExp, beta.hat.c=beta.hat[m], sigma.c=sigma.hat[m], method="Brent", upper=log(2.0), lower=-1*log(2.0))
  r$par
}))
# bias plot: Fig. 1 from Zhong and Prentice if across a range off beta.hat
#bias <- beta.hat-beta.cor
#plot(beta.hat, bias, type="l", lwd=2, ylim=c(-0.3,0.3))
beta.cor <- ifelse(novel==1, beta.hat, beta.cor) # correct only for novel SNPs

# calculate second-stage variance that corresponds to 95% prior certainty interval
sigma2 <- (abs(beta.hat-beta.cor))^2
sigma2 <- ifelse(sigma2==0, 1e-10, sigma2)
# 95% prior certainty interval =
#c(exp(0-1.96*sqrt(sigma2)),exp(0+1.96*sqrt(sigma2)))
#plot(exp(beta.hat), exp(beta.hat), type="l", lwd=2, xlab="OR.hat", ylab="OR")
#points(exp(beta.hat), exp(beta.cor), pch=16)
#lines(exp(beta.hat), exp(beta.cor-1.96*sqrt(sigma2)), lty=3)
#lines(exp(beta.hat), exp(beta.cor+1.96*sqrt(sigma2)), lty=3)

###### Run JAGs hierarchical model
# define JAGS model
model.string <-
  "model {
    C <- 10000 #this just has to be large enough to ensure all phi[m]'s > 0
    for(m in 1:M) {
      beta[m] ~ dnorm(beta.hat[m], prec.beta[m])  # generate the MLE effect estimates
      OR[m] <- exp(beta[m])
      r[m] <- exp(abs(beta[m])) # use for calculation in %FRR assuming p=RAF

      # normal prior on beta using the zeroes trick
      phi[m] <- C-l[m]
      l[m] <- -0.5*log(2*3.14) - 0.5*log(sigma2[m]) - 0.5 * pow((beta[m]-beta.cor[m]),2)/sigma2[m]
      zeros[m] ~ dpois(phi[m])

      # calculate components for FRR
      lambda[m] <- (p[m]*r[m]*r[m] + q[m])/((p[m]*r[m] + q[m])*(p[m]*r[m] + q[m]))
      log.lambda[m] <- log(lambda[m])
    }

    # calculate proportion of FRR
    FRR.all <- inprod(log.lambda[], include.FRR[])/log(lambda0)
    FRR.known <- inprod(log.lambda[], include.FRR.known[])/log(lambda0)

    # prior on lambda0 (known familial relative risk
    lambda0 ~ dnorm(lambda0.m, lambda0.prec)
  }"

jdata <- list(beta.hat=beta.hat,  prec.beta=prec.beta,
              beta.cor=beta.cor, sigma2=sigma2,
              p=p, q=q, M=M, 
              lambda0.m=lambda0.m, lambda0.prec=lambda0.prec,
              include.FRR=include.FRR, include.FRR.known=include.FRR.known,
              zeros=zeros)
var.s <- c("OR", "FRR.all", "FRR.known", "beta")
model.fit <- jags.model(file=textConnection(model.string), data=jdata, n.chains=1, n.adapt=5000)
update(model.fit, n.iter=10000, progress.bar="text")
model.fit <- coda.samples(model=model.fit, variable.names=var.s, n.iter=20000, thin=2, quiet=F)
model.coda <- as.data.frame(model.fit[[1]])
est <- apply(model.coda, 2, mean)

lambda.est <- est[grep("lambda", names(est))]
beta.est <- est[grep("beta", names(est))]
OR.est <- est[grep("OR", names(est))]
FRR.all.est <- est[grep("FRR.all", names(est))]
FRR.known.est <- est[grep("FRR.known", names(est))]

# output results
# plot of SNP bias reduction
pdf("BCC.Estimates.pdf")
plot(exp(beta.hat), exp(beta.est), pch=16, col=novel, xlab="MLE estimate", ylab="Biased Reduced Estimate")
abline(a=0, b=1)
legend(.9, 1.15, legend=c("Known", "Novel"), col=c(1,2), pch=16)
dev.off()

r <- summary(model.fit)
write.table(r$statistics, file="Summary.BBC.FRR.Estimates.txt", sep="\t")
write.table(r$quantiles, file="Summary.BBC.FRR.Estimates.txt", sep="\t", append=T)

