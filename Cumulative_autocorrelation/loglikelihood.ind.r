### loglikehood function of the data subject to the SIR-model parameters assuming 
### independent observational resdiuals 
loglikelihood.ind <- function(para,infolist=list('initial'=c(600500,190,399310,0),'times'=0:(8*52),'data'=numeric(0)))
{
## arguments:
## para: parameter vector including alpha, beta, standard deviation
## infolist: avaibable information, i.e. the initial condition of the system, time
##           horizon and data
times <- infolist[['times']];
# solve the model with the given parameters
sol <- rk(infolist[['initial']],times,sir,list('par.alpha'=para['par.alpha'],'par.beta'=para['par.beta']),method='rk4');
# compute the model predicitons on the log scale
sol <- sol[,2:5];
log.obsmean <- log(sol[2:length(times),4]-sol[(1:length(times)-1),4]);
# compute the loglikelihood of the observed residuals to model predicted mean
LL <- sum(dnorm(x=(infolist[['data']]-log.obsmean),mean=0,sd=para['par.sd'],log=TRUE))
# return negative loglikelihood
return(-LL)
}