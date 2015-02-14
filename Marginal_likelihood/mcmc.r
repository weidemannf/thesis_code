### mcmc generates a Metropolis-Hastings sample from an (unnormalized)
### probability log density
### arguments: fun = unnormalized log density function
###            initialpara = initial parameter vector
###            transitionkernel = covariance matrix of the multivariate Gaussian 
###                               transition density
###            infolist = further arguments passed to fun
###            iter, burnin, thinning = number of iterations, burn in period, thinning                   
### output: list containing 1) a matrix for the posterior sample and 2) the acceptance ratio
mcmc <- function(fun,initialpara,transitionkernel,infolist=ls(),iter=10000,burnin=1000,thinning=10)
{
# setting up the storage matrix in momentaneous parameter vector
paramatrix <- initialpara;
steppara <- initialpara;
# compute density function at current parameter
stepLL <- fun(steppara,infolist=infolist)
logprobdens.sample <- stepLL;
jumpcounter <- 0;

### start the MH-chain  ##############################################
for (i in 2:iter){
# propose a parameter vector
  canpara <- steppara + rmnorm(n = 1, mean = rep(0, length(steppara)), varcov=transitionkernel);
  names(canpara) <- names(steppara);
# compute density function at parameter candidat and corresponding acceptance ratio
  canLL <- fun(canpara,infolist=infolist);
  acceptprob <- exp(canLL-stepLL);
# sample if transition happens
  u = runif(1);
  if(u<acceptprob){
    steppara <- canpara;
    stepLL <- canLL;
    jumpcounter <- jumpcounter+1;
  }
# store current parameter vector if necessary
  if((i>=burnin) && ((i-burnin)%%thinning==0)){      
      paramatrix <- rbind(paramatrix,steppara);
      logprobdens.sample <- c(logprobdens.sample,stepLL)
  }
}
# return final sample, the corresponding probability log density sample and the acceptance ratio
return(list('paramatrix'=paramatrix[2:dim(paramatrix)[1],],'probdens.sample'=logprobdens.sample[2:dim(paramatrix)[1]],'acceptance'=jumpcounter/iter))
}