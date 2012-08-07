# Function to compute Bandi, Johannes etc estimators for time series x for Drift Diffusion Jump Nonparametrics Early Warning Signals
# Author: Stephen R Carpenter, 15 December 2011
# Modified by: Vasilis Dakos, January 4, 2012

Bandi5 <- function(x0,dx,nx,DT,bw,na,avec)  {
# Set up constants and useful preliminaries
SF <- 1/(bw*sqrt(2*pi))  # scale factor for kernel calculation
x02 <- x0*x0 # second power of x
dx2 <- dx*dx # second power of dx
dx4 <- dx2*dx2  # fourth power of dx
dx6 <- dx2*dx4  # sixth power of dx
# Compute matrix of kernel values
Kmat <- matrix(0,nrow=na,ncol=nx)
for(i in 1:(nx)) {  # loop over columns (x0 values)
  Kmat[,i] <- SF*exp(-0.5*(x0[i]-avec)*(x0[i]-avec)/(bw*bw))
  }
# Compute M1, M2, M4, moment ratio and components of variance for each value of a
M1.a <- rep(0,na)
M2.a <- rep(0,na)
M4.a <- rep(0,na)
M6M4r <- rep(0,na)  # vector to hold column kernel-weighted moment ratio
mean.a <- rep(0,na) # centering of conditional variance
SS.a <- rep(0,na)  # sum of squares
for(i in 1:na) {  # loop over rows (a values)
  Ksum <- sum(Kmat[i,])  # sum of weights
  M1.a[i] <- (1/DT)*sum(Kmat[i,]*dx)/Ksum
  M2.a[i] <- (1/DT)*sum(Kmat[i,]*dx2)/Ksum
  M4.a[i] <- (1/DT)*sum(Kmat[i,]*dx4)/Ksum
  M6.c <- (1/DT)*sum(Kmat[i,]*dx6)/Ksum
  M6M4r[i] <- M6.c/M4.a[i]
  mean.a[i] <- sum(Kmat[i,]*x0[2:(nx+1)])/Ksum 
  SS.a[i] <- sum(Kmat[i,]*x02[2:(nx+1)])/Ksum 
  }
# Compute conditional variance
S2.x <- SS.a - (mean.a*mean.a) # sum of squares minus squared mean
# Compute jump frequency, diffusion and drift
sigma2.Z <- mean(M6M4r)/(5) # average the column moment ratios
lamda.Z <- M4.a/(3*sigma2.Z*sigma2.Z)
sigma2.dx <- M2.a - (lamda.Z*sigma2.Z)
# set negative diffusion estimates to zero
diff.a <- ifelse(sigma2.dx>0,sigma2.dx,0)
sigma2.dx <- M2.a     # total variance of dx
mu.a <- M1.a
outlist <- list(mu.a,sigma2.dx,diff.a,sigma2.Z,lamda.Z,S2.x)
# outputs of function:
# mu.a is drift
# sigma2.dx is total variance of dx
# diff.a is diffusion
# sigma2.Z is jump magnitude
# lamda.Z is jump frequency
# S2.x is conditional variance
return(outlist)
} # end Bandi function