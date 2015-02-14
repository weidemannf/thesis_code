##### Berechnung der Normalisierungskonstante #####
#### see Chip et al.
############# import necessary functions and packages ##########################
# functions which provide incidence data, demographic data or starting conditions 
source(file="demographicdataEFS.R")
# further functions providing, e.g., the likelihood or ODE system
source(file="functions.R")
source(file="loglikelihoodESS.R")
source(file="rotavacc.R")
source(file="likelihoodEFS.R")

# packages 
require('deSolve')
require('mnormt')


#### Setting of initial conditions and fixed parameters #######################
###############################################################################
      N <- initialEFS()        # initial condition of ODE-system 
      alpha <- logit(c(0.6, 0.667))
      theta <- logit(c(0.5, 0.25, 0))
      contacttype <- 6
  # demographic parameters
      birthcounts <- birthsEFS()
      death <- deathrate()
      agerates <-ageing()
      migra <- migrationEFS()
    
  # define time horizon   
      Stoptime  <- 23*52; times     <- 1:Stoptime;
  # import incidence data
      datEFS <- dataEFS(n=13)

### import inference results from previous modelling without vaccination as well 
### as the estimated cumulative autocorrelation of the time series (meaness.Rdata)
  # results.Rdata contains a vector para which is the posterior mode of the former modelling
    load(file="results.Rdata")
    load(file="meaness.Rdata")

### define the list with necessary information to compute the likelihood, this contains
### parameters which are not subject of estimation   
    infolist=list('N'=N,'times'=times,'births'=birthcounts,'death'=death,'migra'=migra,'ageing'=agerates,'dataEFS'=datEFS,'alpha'=alpha,'theta'=theta,'contacttype'=contacttype,'essEFS'=meanessEFS)

### define the prior distribution for each of the four parameter vectors #######
# para1 = modelpara, para2 = vaccpara, para3 = coverpara, para4 = immunpara as in likelihoodEFS.Rdata
    para1  <- para[c(1:13,15,17)]        
    para2  <- c(-0.147,-0.147,log(1/52))
    covar2 <- diag(c(0.871,0.871,0.5))             
    para3 <-  rep(0,8)
    covar3 <- diag(,nrow=8)     
    para4 <- c(0,0,0)
    covar4 <-  diag(c(1,0.01,0.0001)) 
    paradims <- c(length(para1),length(para2),length(para3),length(para4))
    block1 <- 1:paradims[1]
    block2 <- (cumsum(paradims)[1]+1):(cumsum(paradims)[2])
    block3 <- (cumsum(paradims)[2]+1):(cumsum(paradims)[3])
    block4 <- (cumsum(paradims)[3]+1):(cumsum(paradims)[4])   
    dimension <- sum(paradims)
# import of prior information on the modelparameters (para1)           
  priorparas <- read.table('priordensityparas.txt')           

## actual definition of log prior function 
## arguments: four parameter vectors as also defined in likelihoodEFS.R
## logprior is defined on the transformed space, i.e. R^d
  logprior <- function(modelpara, vaccpara, coverpara, immupara){
      logprio   <- 0;
    # transmission model parameters prior is defined through the information in 'priordensityparas.txt'
    # using a skewed normal distribution 
      for (i in 1:length(modelpara)){
           logprio <- logprio + log(dsnorm(modelpara[names(modelpara[i])],
                                      mean=priorparas$mu[priorparas$parameter==names(modelpara[i])],
                                      sd=priorparas$sigma[priorparas$parameter==names(modelpara[i])],
                                      a=priorparas$a[priorparas$parameter==names(modelpara[i])]))
            }
    # the three remaining parameter vectors are defined as multivariate normal using the above
    # defined means and covariance matrices, i.e. paraX and covarX                  
      logprio <- logprio + dmnorm(vaccpara, mean = para2, varcov=covar2, log = TRUE) 
                         + dmnorm(coverpara, mean = para3, varcov=covar3, log = TRUE) 
                         + dmnorm(immupara, mean = para4, varcov=covar4, log = TRUE);
    # additional priors to stabilize the transformed effectiveness parameters within [-5;5]
    # and also to penalize large relative differences of the contact parameters, i.e. paravec1[10:12]
      logprio <- logprio - exp(sum((vaccpara[1:2]/5)^10)) 
                         - exp(sum((modelpara[c(10,11,12)]-modelpara[c(12,10,11)])^2));
      return(logprio[[1]]);
  }       

## fetch the posterior sample for the model under consideration
  load(file='adaptivemcmcresults.Rdata')
                                                                                                                                                                
### Initialization of the Marginal Likelihood Estimation #######################
################################################################################
#### define transition kernel
  paramatrix  <- rbind(paramatrix1,paramatrix2,paramatrix3,paramatrix4)
  transkernel  <-  cov(t(paramatrix))
### set algorithm specifications: number of marginal likelihood estimations J, and 
### number of point proposals K
  J <- 100 
  K <- 100
  margLL <- rep(0,J)

### loop over marginal likelihood estimations
  for (j in 1:J){
### sample a parameter vector for point posterior density estimation
    fixLP <- NaN
    while (is.nan(fixLP)){
      fixpara <-  rmnorm(n = 1, mean = rowMeans(paramatrix), varcov=transkernel)
      fixpara1 <- fixpara[block1]; names(fixpara1) <- names(para1);
      fixpara2 <- fixpara[block2]
      fixpara3 <- fixpara[block3]
      fixpara4 <- fixpara[block4]
      fixLL <- likelihoodEFS(modelpara=fixpara1,vaccpara=fixpara2,coverpara=fixpara3,immupara=fixpara4,info=infolist)
      fixLP <- fixLL + logprior(modelpara=fixpara1, vaccpara=fixpara2, coverpara=fixpara3, immupara=fixpara4)
    }
    fixpara  <- c(fixpara1,fixpara2,fixpara3,fixpara4)

### compute the estimator's numerator  
    acceptsample    <- ((fixLP-LPsample)-abs(fixLP-LPsample))/2
    transdenssample <- dmnorm(x=t(paramatrix), mean = fixpara, varcov=transkernel)
    numerator <- mean(exp(acceptsample)*transdenssample)

### compute the estimator's denumerator 
    acceptprob <- rep(0,K)
### loop over the proposed points
    for (i in 1:K){
      canpara <-  rmnorm(n = 1, mean = fixpara, varcov=transkernel)
      canpara1 <- canpara[block1]; names(canpara1) <- names(para1);
      canpara2 <- canpara[block2]
      canpara3 <- canpara[block3]
      canpara4 <- canpara[block4]
      canLL <- likelihoodEFS(modelpara=canpara1,vaccpara=canpara2,coverpara=canpara3,immupara=canpara4,info=infolist)
      canLP <- canLL + logprior(modelpara=canpara1, vaccpara=canpara2, coverpara=canpara3, immupara=canpara4)
      if (is.nan(canLP)) {acceptprob[i] <- 0}
      else {acceptprob[i] <- min(exp(canLP-fixLP),1)}
    }
    denumerator <- mean(acceptprob)

### compute the marginal loglikelihood for the j-th point
    margLL[j] <- fixLP - log(numerator/denumerator)
    print(c('marginal logLikelihood = ',margLL[j]))
  }
### compute the geometric mean of the J marginal likelihood estimates
  meanmargLL <- mean(margLL)
### save results
  save(margLL,file='margLL.Rdata')


