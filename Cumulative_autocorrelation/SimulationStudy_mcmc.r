### study to investigate the accuracy of the assumption of asymptotic normality, 
### i.e. for one generated data set we will compute confidence regions based on 
### the observed Fisher information, but we also apply MCMC to sample from the 
### likelihood function. Then MCMC-sample and assumed normal distrbution will be 
### compared; cumulative autocorrelation is not considered here

### import required function and packages
source(file="effectivesamplesize.R");
source(file="sirseasonal.R");
source(file="datasimulationarma22.R");
source(file="datasimulationBB.R");
source(file="loglikelihood.ind.R");
source(file="loglikelihood.ar.R");
source(file="mcmc.R");

require('deSolve');
require('optimx');
require('mnormt');
require('numDeriv')
require('ks');

### set parameters for data simulation
# initial condition for sir-model (roughly 82,000,000 people) and time horizon
initial <- c(600500,190,399310,0)*82;
times <- 0:(8*52);
infodata <- list('initial'=initial,'times'=times)
# parameters for the sir-model
par.alpha <-5/3;
par.beta  <-0.3;
# set method for simulation of observation residuals, BrownianBridge or Arma(2,2)
simmethod <- 'BB'  ### either 'BB' or 'arma'
# parameters for the observation residuals process in case of ARMA(2,2) process
par.corr  <- c(0.4,0.4);
par.theta <- c(1,-0.6,0);
sigma     <- 0.05;


########### Start of the data simulation and parameter inference ##############
###############################################################################
# set the number of data samples K
K <- 1;
# define storevariables: for MLE and covariance of point estimated
mode_estimates <- matrix(0,nrow=K,ncol=2);
mode_covar     <- array(0,dim=c(K,2,2))                             
PV_mode        <- rep(0,K);

for (i in 1:K){
# simulation of data using the functions 'datasim.arma22' or 'datasim.BB', both
# functions utilize the model parameters and parameters for the correlation structure
# both functions return the logs of the weekly number of cases for eight years
if (simmethod=='arma'){
dat <- datasim.arma22(par.alpha=par.alpha,par.beta=par.beta,par.corr=par.corr,par.theta=par.theta,sigma=sigma,infolist=infodata)
}
if (simmethod=='BB'){
dat <- datasim.BB(par.alpha=par.alpha,par.beta=par.beta,sigma=sigma,infolist=infodata)
}

### Estimation of parameters alpha and beta based on the simulated data (assuming
### the initial condition of the SIR-model is known) ###########################
# define available information
infoest <- list('initial'=initial,'times'=times,'data'=dat)


### Estimation of the likelihood function by computing the MLE and corresponding
### Fisher information using a model assuming AR(1)-Residuals ##################
para <- c('par.alpha'=par.alpha,'par.beta'=par.beta,'par.corr'=par.corr[1],'par.sd'=sigma);
optpara_mode <- optim(para,loglikelihood.ar,method='Nelder-Mead',infolist=infoest,
                    control=list(maxit=10000,'trace'=0,reltol=1e-11),hessian=FALSE);      
covar_mode <- solve(hessian(func=loglikelihood.ar,x=optpara_mode[[1]],
                    method="Richardson",method.args=list(r=2),infolist=infoest));
# store MLE and covariance matrix
mode_estimates[i,] <- optpara_mode[[1]][1:2];
mode_covar[i,,]     <- covar_mode[1:2,1:2];
# compute PIT
d<-(c(par.alpha,par.beta)-optpara_mode[[1]][1:2])%*%solve(covar_mode[1:2,1:2])%*%(c(par.alpha,par.beta)-optpara_mode[[1]][1:2])
PV_mode[i] <- 1-exp(-d/2);

        

### Computation of a MCMC sample using the same likelihood function as target 
### distribution  ##############################################################
# transition kernel is based on the previously estimated covariance matrix
tr.kernel <- covar_mode;
# computation of the sample using the funciton 'mcmc'
parasample <- mcmc(fun=-loglikelihood.ar,initialpara=optpara_mode[[1]],
                   transitionkernel=tr.kernel,infolist=infoest,iter=100000,
                   burnin=10000,thinning=20);

### Checking for similarity between the MCMC sample and the normal distribution 
### first for each point of the MCMC sample we compute the PIT                     
# compute the effective distances according to the normal distribution
d <- rep(0,dim(parasample[[1]])[1]);
for (k in 1:dim(parasample[[1]])[1]){
d[k]<-(mode_estimates[i,]-parasample[[1]][k,1:2])%*%solve(mode_covar[i,,])%*%(mode_estimates[i,]-parasample[[1]][k,1:2])
}
# compute the PIT sample
PV_sample <- 1-exp(-d/2);

### compute a kernel density estimate based on the MCMC sample
# bandwith selection and density estimation
H.pi <- Hpi(parasample[[1]][,1:2], pilot="dscalar")
fhat <- kde(parasample[[1]][,1:2], H=H.pi, compute.cont=TRUE)

# evaluate the estimated normal density on the same grid
densapprox <- matrix(0,nrow=length(fhat$eval.points[[1]]),ncol=length(fhat$eval.points[[2]]));
for (k in 1:length(fhat$eval.points[[1]])){
  for (j in 1:length(fhat$eval.points[[2]])){
densapprox[k,j] <- dmnorm(x=c(fhat$eval.points[[1]][k],fhat$eval.points[[2]][j]), mean = mode_estimates[i,], varcov=mode_covar[i,,], log = FALSE);
  }
}
# compute the contour lines corresponding to the 10%, 50%, 75%, and 90% conficence
# interval of the normal density
confilevels <-  dmnorm(x=matrix(c(sqrt(-2*log(1-c(0.1,0.5,0.75,0.9))/solve(mode_covar[i,,])[1,1]),rep(0,4)),nrow=4),mean=c(0,0),varcov=mode_covar[i,,]);

### plotting of results
layout(matrix(1:3,nrow=1))
par(mar=c(5, 5.1, 4, 2) + 0.1)
# first: the MCMC sample and corresponding kernel estimate confidence intervals (10%, 50%, 75%, 90%)
contour(x=fhat$eval.points[[1]],y=fhat$eval.points[[2]],z=fhat$estimate,
       levels=c(fhat$cont[10],fhat$cont[50],fhat$cont[75],fhat$cont[90]),
       labels=c('10%','50%','75%','90%'),main='MCMC-Sample of loglikelihood function',
       labcex=0.8,ylab=expression(R[0]),xlab=expression(beta),cex.lab=2.6,cex.axis=2.2,cex.main=2.6);
points(parasample[[1]][,1:2],type='p',col='gray')
contour(x=fhat$eval.points[[1]],y=fhat$eval.points[[2]],z=fhat$estimate,
       levels=c(fhat$cont[10],fhat$cont[50],fhat$cont[75],fhat$cont[90]),
       labels=c('10%','50%','75%','90%'),labcex=0.8,add=TRUE);
# second: the confidence intervals of the approximate normal distribution (10%, 50%, 75%, 90%)
contour(x=fhat$eval.points[[1]],y=fhat$eval.points[[2]],z=densapprox,levels=confilevels,
        labels=c('10%','50%','75%','90%'),main='Approximated loglikelihood function',
        labcex=0.8,ylab=expression(R[0]),xlab=expression(beta),cex.lab=2.6,cex.axis=2.2,cex.main=2.6);
# third: the PIT histogram of the MCMC sample subject to the normal distribution        
hist(PV_sample,main='PIT-Histogram MCMC-Sample',xlab=expression(1-alpha),cex.lab=2.6,cex.axis=2.2,cex.main=2.6)
}

### save results
save(parasample,mode_estimates,mode_covar,PV_mode,file=paste('mcmc_simresults_',simmethod,'.Rdata');

 