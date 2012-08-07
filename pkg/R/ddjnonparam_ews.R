# Drift Diffusion Jump Nonparametrics Early Warning Signals
# Author: Stephen R Carpenter, 15 December 2011
# Modified by: Vasilis Dakos, January 4, 2012


# Inputs:
#  x0 is the regressor
#  dx is the first difference of x0
#  nx is number of first differences
#  DT is time step
#  bw is the bandwidth for the kernel
#  na is number of a values for computing the kernel
#  avec is the mesh for the kernel

require(KernSmooth)

# MAIN FUNCTION
ddjnonparam_ews<-function(timeseries,bandwidth=0.6,na=500,logtransform=TRUE,interpolate=FALSE){
	
	timeseries<-ts(timeseries) #strict data-types the input data as tseries object for use in later steps
	if (dim(timeseries)[2]==1){
		Y=timeseries
		timeindex=1:dim(timeseries)[1]
		}else if(dim(timeseries)[2]==2){
		Y<-timeseries[,2]
		timeindex<-timeseries[,1]
		}else{
		warning("not right format of timeseries input")
		}
		
	# Interpolation
	if (interpolate){
		YY<-approx(timeindex,Y,n=length(Y),method="linear")
		Y<-YY$y
		}else{
		Y<-Y}
			
	# Log-transformation
	if (logtransform){
		Y<-log(Y+1)}
		
	# Preliminaries
	Xvec1<-Y
	Tvec1<-timeindex
	dXvec1 <- diff(Y)

	DT <- Tvec1[2]-Tvec1[1]
	bw <- bandwidth*sd(Xvec1) # bandwidth 
	alow <- min(Xvec1)
	ahigh <- max(Xvec1)
	na <- na
	avec <- seq(alow,ahigh,length.out=na)
	nx <- length(dXvec1)

	# Bandi-type estimates
	ParEst <- Bandi5(Xvec1,dXvec1,nx,DT,bw,na,avec)
	Drift.vec <- ParEst[[1]]
	TotVar.dx.vec <- ParEst[[2]]	
	Diff2.vec <- ParEst[[3]]
	Sigma2Z <- ParEst[[4]]
	LamdaZ.vec <- ParEst[[5]]
	S2.vec <- ParEst[[6]]

	# Interpolate time courses of indicators
	TotVar.i <- approx(x=avec,y=TotVar.dx.vec,xout=Xvec1)
	TotVar.t <- TotVar.i$y
	Diff2.i <- approx(x=avec,y=Diff2.vec,xout=Xvec1)
	Diff2.t <- Diff2.i$y
	Lamda.i <- approx(x=avec,y=LamdaZ.vec,xout=Xvec1)
	Lamda.t <- Lamda.i$y
	S2.i <- approx(x=avec,y=S2.vec,xout=Xvec1)
	S2.t <- S2.i$y

	# Plot the data
	dev.new()
	par(mfrow=c(2,1),mar=c(3, 3, 2, 2),mgp=c(1.5,0.5,0),oma=c(1,1,1,1))
	plot(Tvec1,Xvec1,type='l',col='black',lwd=2,xlab='',ylab='original data')
	grid()
	plot(Tvec1[1:length(Tvec1)-1],dXvec1,type='l',col='black',lwd=2,xlab='time',ylab='first-diff data')
	grid()
	
	# Plot indicators versus a
	dev.new()
	par(mfrow=c(2,2),mar=c(3, 3, 2, 2) ,cex.axis=1,cex.lab=1,mgp=c(2,1,0),oma=c(1,1,2,1))
	plot(avec,S2.vec,type='l',lwd=1,col='black',xlab='a',ylab='conditional variance')
	plot(avec,TotVar.dx.vec,type='l',lwd=1,col='blue',xlab='a',ylab='total variance of dx')
	plot(avec,Diff2.vec,type='l',lwd=1,col='green',xlab='a',ylab='diffusion')
	plot(avec,LamdaZ.vec,type='l',lwd=1,col='red',xlab='a',ylab='jump intensity')
	mtext("DDJ nonparametrics versus a",side=3,line=0.1,outer=TRUE)
	
	# Plot indicators versus time
	dev.new()
	par(mfrow=c(2,2),mar=c(3, 3, 2, 2),cex.axis=1,cex.lab=1,mgp=c(1.5,0.5,0),oma=c(1,1,2,1))
	plot(Tvec1,S2.t,type='l',lwd=1,col='black',xlab='time',ylab='conditional variance')
	plot(Tvec1,TotVar.t,type='l',lwd=1,col='blue',xlab='time',ylab='total variance of dx')
	plot(Tvec1,Diff2.t,type='l',lwd=1,col='green',xlab='time',ylab='diffusion')
	plot(Tvec1,Lamda.t,type='l',lwd=1,col='red',xlab='time',ylab='jump intensity')
	mtext("DDJ nonparametrics versus time",side=3,line=0.1,outer=TRUE)

	# Output
	nonpar_x<-data.frame(avec,S2.vec,TotVar.dx.vec,Diff2.vec,LamdaZ.vec)
	nonpar_t<-data.frame(Tvec1,S2.t,TotVar.t,Diff2.t,Lamda.t)
	return(c(nonpar_x,nonpar_t))
	}