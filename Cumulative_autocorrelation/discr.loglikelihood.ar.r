### loglikehood function of count data subject to the SIR-model parameters assuming 
### an autoregressive-type correlation structure within the count data
## arguments:
## para: parameter vector including alpha, beta, standard deviation, correlation coeff
## infolist: avaibable information, i.e. the initial condition of the system, time
##           horizon and data
discr.loglikelihood.ar <- function(para,infolist=list('initial'=c(600500,190,399310,0),'times'=0:(8*52),'data'=numeric(0)))
{
# compute solution to the SIR-ODE model
times <- infolist[['times']];
sol <- rk(infolist[['initial']],times,sir,list('par.alpha'=para['par.alpha'],'par.beta'=para['par.beta']),method='rk4');
# compute vector for expected number of new weekly cases
sol <- sol[,2:5];
obsmean <- sol[2:length(times),4]-sol[(1:length(times)-1),4];    
# compute relative residuals of each data point with respect ot predicted expectation 
# and shift by one time unit
lastresi <- c(0,infolist[['data']][1:length(obsmean)]/obsmean-1)
# compute log likelihood 
LL <- sum(dnbinom(x=infolist[['data']],mu=obsmean*(1+para['par.corr']*lastresi[1:length(obsmean)]),size=para['par.size'],log=TRUE))
# return negative loglikelihood
return(-LL)
}                                                                                     