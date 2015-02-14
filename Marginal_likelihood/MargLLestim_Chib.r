### Simulation study to investigate the performance of the marginal likelihood 
### estimation algortihm proposed by Chib et al. and modified version of their appraoch.
### Marginal likelihood is estimated for three model (posterior) types and variable
### dimension d

## import necessary functions and package
  source(file="mcmc.R");
  require('mnormt');

### define the model setting to be investigated
# set dimension 
  d <- 10;
# set type of posterior function, either 'normal', 'mixed' oder 'exponential'
  model <- 'exponential'; 

### define corresponding prior and likeliood function for each model ###
### note, that likelihood functions are defined without simulating actual data, and thus
### are subject to virtual data

### normal posterior with mean zero:
  if (model=='normal'){
## prior log density with variance sigma*diag(1,d)
    sigma <- 1;
    prior <- function(x){
      return(dmnorm(x, mean=rep(0,d), varcov=diag(sigma,d),log=TRUE))}
## log likelihood function with variance mu*diag(1,d)
    nu <- 1/100;
    LL <- function(x){
      return(dmnorm(x, mean=rep(0,d), varcov=diag(nu,d),log=TRUE))}
## transition kernel covariance for MCMC sampling
    transkern <- diag((nu*sigma)/(nu+sigma),d)*2.4^2/d;
    if(d==1){transkern <- (nu*sigma)/(nu+sigma)*2.4^2;}
## exact log marginal likelihood
    trueLML <- log(1/(2*pi*(sigma+nu))^(d/2))
  }

### mixed normal posterior:
  if (model=='mixed'){
## prior log density with variance sigma*diag(1,d) and mean zero
    sigma <- 1;
    prior <- function(x){
      return(dmnorm(x, mean=rep(0,d), varcov=diag(sigma,d),log=TRUE))}
## log likelihood function consisting of two multivariate normla densities with variance mu*diag(1,d)
    nu <- 1/100;
    LL <- function(x){
      return(log(0.6*dmnorm(x, mean=c(1,rep(0,d-1)), varcov=diag(nu,d),log=FALSE)+0.4*dmnorm(x, mean=c(-1,rep(0,d-1)), varcov=diag(nu,d),log=FALSE)))}
## transition kernel covariance for MCMC sampling
    transkern <- diag(c(2,rep((nu*sigma)/(nu+sigma),d-1)))*2.4^2/d;
    if(d==1){transkern <- 2*2.4^2;}
## exact log marginal likelihood
    trueLML <- log(1/(2*pi*(sigma+nu))^(d/2)*exp(-(sigma+nu)/(2*sigma*nu)*(nu/(sigma+nu)-(nu/(sigma+nu))^2)))
  }

### exponential posterior:
  if (model=='exponential'){
## prior log denisiy with positive support and mean sigma
    sigma <- 1;
    prior <- function(lambda){
      return(sum(dexp(x=lambda, rate = sigma, log = TRUE)))}
## log likelihood function corresponding to nsize exponetially distributed data points with 
## empirical mean nu/size, i.e. nu is sum of data values
    nu <- 3; nsize <- 3;
    LL <- function(lambda){
      if (min(lambda)<0){return(-Inf)}
      else{return(nsize*sum(dexp(x=rep(nu/nsize,length(lambda)), rate = lambda, log = TRUE)))}}
## transition kernel covariance for MCMC sampling
    transkern <- diag((nsize+1)/(sigma+nu)^2,d)*2.4^2/d;
    if(d==1){transkern <- (nsize+1)/(sigma+nu)^2*2.4^2;}
## exact log marginal likelihood
    trueLML <- d*log(gamma(nsize+1)*sigma/(sigma+nu)^(nsize+1));
}

### Fetch a posterior sample from the considered mode or, if none is available, 
### generate a sample by a MH-algorithm with Gaussian proposals
### using the above defined covariance matrix
## check is a sample exists
  sample.ia <- try(load(file=paste('mcmcsample',model,d,'.Rdata')));
  if (class(sample.ia)=="try-error"){
  ## initialze the MCMC algorithm
    initialvec <- rep(0,d);
    if (model=='exponential'){initialvec <- rep(1,d);} 
    logposterior <- function(x,...){
      return(prior(x)+LL(x))}
  ## generate the MCMC - Sample ####
    posteriorsample <- mcmc(fun=logposterior,initialpara=initialvec,transitionkernel=transkern,iter=1000000,burnin=0,thinning=10)  
    samplematrix  <- t(posteriorsample$paramatrix);
    logposteriorsample <- posteriorsample$probdens.sample;
  ## save the posterior sample 
    save(samplematrix, logposteriorsample,file=paste('mcmcsample',model,d,'.Rdata'));
  }


### Estimation of the log marginal likelihood (LML) according to a modified approach #
### proposed by Chib et al. ####################################################
## print true LML
print(c('true log marginal likelihood = ',trueLML));

### configurate the Estimation algorithm #############
# set choice of multivariate proposal kernel for LML estimation
  kern <- 'variance';
  if (kern == 'chib'){transkern <- cov(t(samplematrix))*2.4^2/d;}
  if (kern == 'variance'){transkern <- cov(t(samplematrix));}
## set number of estimator evaluations K (to be averaged later) and number of 
## respective point proposals N
  K <- 10;
  N <- 1000;
## set number M of evaluations of the overall estimator (to investigate its distribution)
## and storage vector
  M <- 1000;
  estMLL <- rep(0,M);
## loop over the number of overall evaluations M
  for (m in 1:M){
  # set storage vectors for estimator's numerators and denumerators
    fixLL <- rep(0,K);
    numerator <- rep(0,K);
    denumerator <- rep(0,K);
    margLL  <- rep(0,K);
  # loop over the number of single estimator evaluations K
    for (k in 1:K){
    # sample a parameter vector
      repeat{
        fixpara  <-  c(rmnorm(n = 1, mean = rowMeans(samplematrix), varcov=cov(t(samplematrix))));
        fixLL <- prior(fixpara)+LL(fixpara);
        if (fixLL>-Inf){break}
      }                       
    # compute the numerator    
      acceptsample    <- ((fixLL-logposteriorsample)-abs(fixLL-logposteriorsample))/2;
      transdenssample <- dmnorm(x=t(samplematrix), mean = fixpara, varcov=transkern);
      numerator <- mean(exp(acceptsample)*transdenssample);
    # compute the demunerator  ########
      acceptprob <- rep(0,N);
      canpara <-  rmnorm(n = N, mean = fixpara, varcov=transkern);
      if (model=='exponential'){
         canLL <- rep(0,N);
         for (n in 1:N)
         canLL[n]  <-  prior(canpara[n,])+LL(canpara[n,]);
      }else{
         canLL  <-  prior(canpara)+LL(canpara);}
      acceptprob <- 1-((1-exp(canLL-fixLL))+abs(1-exp(canLL-fixLL)))/2;
      denumerator <- mean(acceptprob);
    # compute the resulting log marginal likelihood estimate
      margLL[k] <- fixLL - log(numerator/denumerator);
    }
  # compute and print the averaged marginal likelihood over all K estimates
    estMLL[m] <- mean(margLL)
    print(c(m,'th estimated log marginal likelihood = ',estMLL[m])) 
  }
# save the overall LML estimator sample
  save(estMLL,file=paste('simresults_d',model,d,'K',K,'J',J,'kern',kern,'.Rdata'));
