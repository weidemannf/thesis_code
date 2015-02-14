### simulation of count data for the weekly number of new cases based on a seasonal
### SIR-model and correlated residual
discr.datasim.ar <- function(par.alpha=5/3,par.beta=0.3,par.corr=0,par.size=100,
              infolist=list('initial'=c(600500,190,399310,0),'times'=0:(8*52)))
{
## arguments:
## par.alpha,par.beta: parameters of the SIR-model
## par.corr,par.size: parameters for the correlation and variance of the NegBinObservations
## infolist: initial conditions and time horizon

# compute SIR-model solutions and expected number of weekly cases
sol <- rk(infolist[['initial']],infolist[['times']],sir,list('par.alpha'=par.alpha,'par.beta'=par.beta),method='rk4');
sol <- sol[,2:5];
obsmean <- sol[2:length(times),4]-sol[(1:length(times)-1),4];                                   

# iteratively sample an observation using the mean and the last relative deviation 
# (i.e. the residual) from a negative binomial distribution 
obs <- rep(0,length(times)-1);
lastresi <- 0;
for (i in 1:(length(times)-1)){                                                                        
obs[i] <- rnbinom(n=1,mu=obsmean[i]*(1+par.corr*lastresi),size=par.size);
# the relative deviation from the mean
lastresi <- obs[i]/obsmean[i]-1;
}

# return the sampled count data series
return(obs)
}