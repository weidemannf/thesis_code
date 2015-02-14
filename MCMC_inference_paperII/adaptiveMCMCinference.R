#### adaptive Metropolis-Hastings algorithm for Bayesian parameter inference ######
###################################################################################

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
require('optimx')

#### Setting of initial conditions and fixed parameters #######################
###############################################################################
      begTime <- Sys.time()   # time measuring
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
                         + dmnorm(immupara, mean = para4, varcov=covar4, log = TRUE)
    # additional priors to stabilize the transformed effectiveness parameters within [-5;5]
    # and also to penalize large relative differences of the contact parameters, i.e. paravec1[10:12]
      logprio <- logprio - exp(sum((vaccpara[1:2]/5)^10)) 
                         - exp(sum((modelpara[c(10,11,12)]-modelpara[c(12,10,11)])^2))
      return(logprio[[1]])
}       
   
###### Initialisation of the Sample Chain using the adaptive Metropolis-Hastings
###### Algorithm ###############################################################
### Step 1: Definition of the initial transitions kernel, i.e. the covariance 
###         matrix of the centered Gaussian proposal
### Step 2: define initial parameter of the chain

## Checking if there is already a sample from the posterior available   
sample.ia <- try(load(file='adaptivemcmcresults.Rdata'))

## Case 1: there is an existing sample (e.g. from a short prerun), then define 
##         initial configuration based on that sample
if (class(sample.ia)=="character"){
    paramatrix <- rbind(paramatrix1,paramatrix2,paramatrix3,paramatrix4);
  # transition kernel
    transcovar <- cov(t(paramatrix));   
  # initial parameter vector
    steppara1 <- paramatrix1[,length(paramatrix1[1,])]; names(steppara1)<-names(para1);
    steppara2 <- paramatrix2[,length(paramatrix1[1,])]
    steppara3 <- paramatrix3[,length(paramatrix1[1,])]                              
    steppara4 <- paramatrix4[,length(paramatrix1[1,])]
}else{
## Case 2: there is no posterior sample so far, then define the caracteristics based
##         on the prior mode and curvature 
  # transition kernel
    optpara1 <- optim(para1,logprior,vaccpara=para2,coverpara=para3,immupara=para4,method="BFGS",control=list(fnscale=-1))$par
    hess     <- hessian(logprior,modelpara=optpara1,vaccpara=para2,coverpara=para3,immupara=para4)
    transcovar <- matrix(0,nrow=dimension,ncol=dimension)
    transcovar[block1,block1] <- -solve(hess)
    transcovar[block2,block2] <- covar2
    transcovar[block3,block3] <- covar3
    transcovar[block4,block4] <- covar4
  # downscale prior covariance as it is often to wide
    transcovar <- transcovar/10000
  # initial parameter vector
    steppara1 <- optpara1;  names(steppara1)<-names(para1);
    steppara2 <- para2
    steppara3 <- para3
    steppara4 <- para4
}   

### Step 3: variance mulitplicator of the transitionkernel (s_d)
    varscale <- 2.4^2/dimension
    
### Step 4: define length of sample chain  
  # total chain lenght
    K <- 20000
  # length of initial period  
    J <- 3000

### Step 5: define storage matrices    
    paramatrix1 <- matrix(0,nrow=length(para1),ncol=K)
    paramatrix2 <- matrix(0,nrow=length(para2),ncol=K)
    paramatrix3 <- matrix(0,nrow=length(para3),ncol=K)
    paramatrix4 <- matrix(0,nrow=length(para4),ncol=K)
    LPsample    <- rep(0,K)
    jumpcounter <- 0    
  # storage of initial parameter vectors             
    paramatrix1[,1] <- steppara1
    paramatrix2[,1] <- steppara2
    paramatrix3[,1] <- steppara3
    paramatrix4[,1] <- steppara4
       
    
### Step 6: compute log posterior denstiy of initial parameter vector      
    stepLL <- likelihoodEFS(modelpara=steppara1,vaccpara=steppara2,coverpara=steppara3,immupara=steppara4,info=infolist)
    stepLP <- stepLL + logprior(modelpara=steppara1, vaccpara=steppara2, coverpara=steppara3, immupara=steppara4)
    print(stepLP)
    LPsample[1] <- stepLP                    
   
############## Start of the adaptive Metropolis-Hastings Chain #################    
for (i in 2:K){
  # computation of the adaptive transition kernel after the initial period     
    if (i>J){
         interval <- 1:(i-1)
         paramatrix <-  rbind(paramatrix1[,interval],paramatrix2[,interval],paramatrix3[,interval],paramatrix4[,interval])
         transcovar <-  cov(t(paramatrix[,interval]))
    }
    
              
  # propose new parameter candidate 
      increment <-  rmnorm(n = 1, mean = rep(0,dimension), varcov=varscale*transcovar)
      canpara1 <- steppara1 + increment[block1]
      names(canpara1) <- names(para1)  
      canpara2 <- steppara2 + increment[block2]
      canpara3 <- steppara3 + increment[block3]
      canpara4 <- steppara4 + increment[block4]
      
  # compute log posterior density for candidate
      canLL <- likelihoodEFS(modelpara=canpara1,vaccpara=canpara2,coverpara=canpara3,immupara=canpara4,info=infolist)
      canLP <- canLL + logprior(modelpara=canpara1, vaccpara=canpara2, coverpara=canpara3, immupara=canpara4)  

  # compute acceptance probability
      acceptprob <- exp(canLP-stepLP)
  # if candidate is accepted set chain vector (steppara) to candidate
      if (!is.nan(canLP)){
         if(runif(1)<acceptprob){
            steppara1 <- canpara1
            steppara2 <- canpara2
            steppara3 <- canpara3
            steppara4 <- canpara4
            stepLP    <- canLP
            jumpcounter <- jumpcounter+1
        }
      }

  # store the chain parameter vector and its log posterior value
      paramatrix1[,i] <- steppara1
      paramatrix2[,i] <- steppara2
      paramatrix3[,i] <- steppara3
      paramatrix4[,i] <- steppara4
      LPsample[i]     <- stepLP
    
  # print intermediate step
      print(c(round(i,d=0),canLP,round(min(1,acceptprob),d=3),round(jumpcounter/(i-1),d=3)))
}

### save the posterior sample and the corresponding logposterior values and print runtime   
    save(paramatrix1,paramatrix2,paramatrix3,paramatrix4,LPsample,file='adaptivemcmcresults.Rdata')
    runTime <- Sys.time()-begTime                                       
    print(runTime)
