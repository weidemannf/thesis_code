### Simulation study to investigate the performance of the inference appraoch
### incoorparating the cumulative autocorrelation to adjust credibility regions 
### around the posterior mode
### Here, incidence data is in the form of count data, which is simulated by a 
### seasonal sir model with added residuals allowing a possible autoregressive type
### autocorrelation
### investigated inference approaches are MLE assuming indipendent residuals,
### autregressive residuals, and indep. resid. accounting for cumul. autocorr.

### import required function and packages
source(file="effectivesamplesize.R");
source(file="anscombe.R");
source(file="sirseasonal.R");
source(file="discr.datasimulation.R");
source(file="discr.loglikelihood.ind.R");
source(file="discr.loglikelihood.ar.R");

require('deSolve');
require('optimx');
#require('mnormt');

### set parameters for data simulation
# initial condition for sir-model (roughly 82,000,000 people) and time horizon
initial <- c(600500,190,399310,0)*82;
times <- 0:(8*52);
infodata <- list('initial'=initial,'times'=times)
# parameters for the sir-model
par.alpha <-5/3;
par.beta  <-0.3;
# parameters for the observation residuals process
par.corr  <- 0.6;
par.size  <- 100;

########### Start of the data simulation and parameter inference ##############
###############################################################################
# set the number of samples K
K <- 1000;
# define storevariables: each method (ar,ind,ess) estimates two parameters 
# (alpha,beta) for all K data samples
ar_estimates  <- matrix(0,nrow=K,ncol=2);
ind_estimates <- matrix(0,nrow=K,ncol=2);
ess_estimates <- matrix(0,nrow=K,ncol=2);
# storevariables for the 2x2 covariance matrices of the K estimates according to 
# (possibly adjusted) Fisher-information (minus the hessian of loglikelihood)
ar_covar   <- array(0,dim=c(K,2,2))               
ind_covar  <- array(0,dim=c(K,2,2))
ess_covar  <- array(0,dim=c(K,2,2))
# storevariables for the results from the probability integral transform (PIT), 
# intuitively the p-value or quantile of the true parameters with respect to the  
# estimated distribution
PV_ar     <- rep(0,K);
PV_ind    <- rep(0,K);
PV_ess    <- rep(0,K);
# also stored are the results from the effective sample size (or cumulative 
# autocorrelation) estimates, i.e. K time the order of the fitted ARMA process 
# and the CA itself
effss <- matrix(0, nrow=K, ncol=3);


### start of the simulation loop #################
for (i in 1:K){
# simulation of data using the functions 'discr.datasim.ar', which utilize the model 
# and observation parameters to simulate count data for the weekly number of cases
dat <- discr.datasim.ar(par.alpha=par.alpha,par.beta=par.beta,par.corr=par.corr,par.size=par.size,infolist=infodata)

### Estimation of parameters alpha and beta based on the simulated data (assuming
### the initial condition of the SIR-model is known) ###########################
# define available information
infoest <- list('initial'=initial,'times'=times,'data'=dat)

### first: Max-Likelihood-Estimation assuming autoregressive observational residuals 
### for simplicity the optimisation procedure (Nelder-Mead) is initialized at the true value
# compute MLE for the parameters alpha, beta, corr and size
para <- c('par.alpha'=par.alpha,'par.beta'=par.beta,'par.corr'=par.corr,'par.size'=par.size);
optpara_ar <- optim(para,discr.loglikelihood.ar,method='Nelder-Mead',infolist=infoest,
           control=list(maxit=10000,'trace'=0,reltol=1e-10),hessian=TRUE);
# compute the estimated covariance matrix based on Fisher information
covar_ar <- solve(optpara_ar$hessian);
# store estimates and covariance
ar_estimates[i,] <- optpara_ar[[1]][1:2];
ar_covar[i,,]     <- covar_ar[1:2,1:2];
# compute the effective distance d of the true vector to the estimated vector subject to covariance 
# and compute its probability integral transform (PIT)
d<-(c(par.alpha,par.beta)-optpara_ar[[1]][1:2])%*%solve(covar_ar[1:2,1:2])%*%(c(par.alpha,par.beta)-optpara_ar[[1]][1:2])
PV_ar[i] <- 1-exp(-d/2);


### Second: MLE assuming independent residual ##################################
# same steps as above, but with a different loglikelihood function
para <- c('par.alpha'=par.alpha,'par.beta'=par.beta,'par.size'=par.size)
optpara_ind <- optim(para,discr.loglikelihood.ind,method='Nelder-Mead',infolist=infoest,
               control=list(maxit=10000,'trace'=0,reltol=1e-10),hessian=TRUE);
covar_ind <- solve(optpara_ind$hessian);
ind_estimates[i,] <- optpara_ind[[1]][1:2];
ind_covar[i,,]  <- covar_ind[1:2,1:2];
d<-(c(par.alpha,par.beta)-optpara_ind[[1]][1:2])%*%solve(covar_ind[1:2,1:2])%*%(c(par.alpha,par.beta)-optpara_ind[[1]][1:2])
PV_ind[i] <- 1-exp(-d/2);

### Third: CA-Method: accounting for autocorrelation ###########################
# point estimates are the same when assuming indipendent residuals
# observation residuals are computed from the expected observations according to the model
sol <- rk4(initial,times,sir,list('par.alpha'=optpara_ind[[1]][1],'par.beta'=optpara_ind[[1]][2]));
sol <- sol[,2:5];
obsmean <- sol[2:length(times),4]-sol[(1:length(times)-1),4]; 
resids <- anscombe.resi(dat,obsmean,disp=optpara_ind$par['par.size']);

effsamplesize <- effectivesamplesize(resids);
effss[i,] <- c(effsamplesize[[1]],effsamplesize[[2]],effsamplesize[[3]])
covar_ess <- covar_ind/effsamplesize[['ESS']];
ess_estimates[i,] <- optpara_ind[[1]][1:2];
ess_covar[i,,]     <- covar_ess[1:2,1:2];
d<-(c(par.alpha,par.beta)-optpara_ind[[1]][1:2])%*%solve(covar_ess[1:2,1:2])%*%(c(par.alpha,par.beta)-optpara_ind[[1]][1:2])
PV_ess[i] <- 1-exp(-d/2);


### Plot intermediate PIT histogram
layout(matrix(1:3,nrow=1))
hist(PV_ar[1:i],breaks=(0:10)/10,main='Normal Approximation: AR Residuals',xlab='PIT')
hist(PV_ind[1:i],breaks=(0:10)/10,main='Normal Approximation: Ind Residuals',xlab='PIT')
hist(PV_ess[1:i],breaks=(0:10)/10,main='Normal Approximation: CA Method',xlab='PIT')
}
### Plot and save final results
#pdf(file='PIT_comparison_discr.obs.pdf',width=15,height=5)
layout(matrix(1:3,nrow=1))
hist(PV_ar,breaks=(0:10)/10,main='Normal Approximation: AR Residuals',xlab='PIT')
hist(PV_ind,breaks=(0:10)/10,main='Normal Approximation: Ind Residuals',xlab='PIT')
hist(PV_ess,breaks=(0:10)/10,main='Normal Approximation: CA Method',xlab='PIT')
#dev.off()

save(par.corr,ar_estimates,ind_estimates,ess_estimates,ar_covar,
     ind_covar,ess_covar,PV_ar,PV_ind,PV_ess,file='mode_simresults_discr.obs.Rdata');

 
      
     