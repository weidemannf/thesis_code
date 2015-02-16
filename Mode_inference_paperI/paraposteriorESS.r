#### unnormalized logposterior function of model for variable parameter vectors
paraposteriorESS <- function(para=numeric(0),info=list(N,times,births,death,dataWFS,dataEFS,
                 alpha,theta,mu,omega,omega0,beta,p,sais,contactpara,h,fact,disp,
                 contacttype,migra,ageing,essEFS,essWFS))
### arguments:
### para = vector of posterior function arguments
###  , i.e. parameters which are subject of estimation
### info = list of all required posterior arguments (including those listed in para)
{    

# extract parameters and data from the information list
    # initial decompostion of population
      N <- info[['N']] 
    # time horizon of the model
      times <- info[['times']]; stoptime <- length(times);      
    # incidence data from EFS and WFS
      datEFS <- info[['dataEFS']]
      datWFS <- info[['dataWFS']]
    # Cumulative Autocorrelation or 'effective sample size' (ess) of all time series
      essEFS <- info[['essEFS']]
      essWFS <- info[['essWFS']]      
    # further arguments as listed in rota.R
    # parameters are also retransformed on (0,1) or (0,Inf) where necessary      
      birthcounts <- info[['births']]
      death  <- info[['death']]
      agerates <- info[['ageing']]
      migra  <- info[['migra']]    
      contacttype <- info[['contacttype']]
      alpha <- expit(info[['alpha']])                
      theta <- c(expit(info[['theta']]),0)       
      mu <- exp(info[['mu']])                     
      omega <- exp(info[['omega']])      
      omega0 <- exp(info[['omega0']])            
      beta <- exp(info[['beta']])                      
      p <- exp(info[['p']])                     
      sais <- c(info[['sais']][1],expit(info[['sais']][2]),
                info[['sais']][3],expit(info[['sais']][4]))  
      contactpara <- exp(info[['contactpara']])       
    # reporting rates and their relative increases  
      h <- expit(info[['h']])                    
      fact <- exp(info[['fact']])
    # dispersion parameter of observational distribution      
      disp <- exp(info[['disp']])                 

# now if parameters are included in 'para', these values are taken instead  
      if (!is.na(para['alpha1']))      {alpha[1]<-expit(para['alpha1'])}
      if (!is.na(para['alpha2']))      {alpha[2]<-expit(para['alpha2'])}
      if (!is.na(para['theta1']))      {theta[1]<-expit(para['theta1'])}
      if (!is.na(para['theta2']))      {theta[2]<-expit(para['theta2'])}
      if (!is.na(para['mu']))          {mu<-exp(para['mu'])}
      if (!is.na(para['omega']))       {omega<-exp(para['omega'])}
      if (!is.na(para['omega0']))      {omega0<-exp(para['omega0'])}
      if (!is.na(para['beta']))        {beta<-exp(para['beta'])}
      if (!is.na(para['p']))           {p<-exp(para['p'])}
      if (!is.na(para['a1']))          {sais[1]<-para['a1']}
      if (!is.na(para['b1']))          {sais[2]<-expit(para['b1'])}
      if (!is.na(para['a2']))          {sais[3]<-para['a2']}
      if (!is.na(para['b2']))          {sais[4]<-expit(para['b2'])}
      if (!is.na(para['contactpara1'])){contactpara[1]<-exp(para['contactpara1'])}
      if (!is.na(para['contactpara2'])){contactpara[2]<-exp(para['contactpara2'])}
      if (!is.na(para['contactpara3'])){contactpara[3]<-exp(para['contactpara3'])}
      if (!is.na(para['h1']))          {h[1]<-expit(para['h1'])}
      if (!is.na(para['h2']))          {h[2]<-expit(para['h2'])}
      if (!is.na(para['fact1']))       {fact[1]<-exp(para['fact1'])}
      if (!is.na(para['fact2']))       {fact[2]<-exp(para['fact2'])}
      if (!is.na(para['disp']))        {disp<-exp(para['disp'])}
      

### Computation of ODE solution for likelihood computation ##########
#####################################################################     
# next generation matrix for contact parameters and matrixtype  
    contact <- contactmatrix(cpara=contactpara,pattern=contacttype)
# computation of ODE solution for above parameters and demographic information
# using Runge-Kutta 4 scheme (package 'deSolve')  
    sol<-rk(N,times,rota,list('alpha'=alpha,'theta'=theta,'mu'=mu,'omega'=omega,
         'omega0'=omega0,'beta'=beta,'p'=p,'sais'=sais,'contact'=contact,
         'births'=birthcounts,'death'=death,'mig'=migra,'ageing'=agerates),method = "rk4")
# formatting solution array: time point, class, age group
    sol <- sol[times,2:(length(N)+1)]; dim(sol)<- c(length(times),dim(N)[1],dim(N)[2])
# calculate weekly age stratified number of new cases from solution: 
# (15th class counts number of occured cases)
    inci  <- sol[2:length(times),dim(N)[1],1:dim(N)[2]]-sol[1:(length(times)-1),dim(N)[1],1:dim(N)[2]]
# given the weekly case numbers, calculate loglikelihood of the data from WFS and EFS
# discard the first ten years from the solution (burn in)
    logLLEFS <- loglikelihoodESS(inci[(11*52):(stoptime-1),1:dim(N)[2]],datEFS[1:10,1:(stoptime-11*52)],h[1],fact[1],disp,essEFS,region='EFS')
    logLLWFS <- loglikelihoodESS(inci[(11*52):(stoptime-1),1:dim(N)[2]],datWFS[1:10,1:(stoptime-11*52)],h[2],fact[2],disp,essWFS,region='WFS')
   
#### Computation of log priorprobability of parameters to estimate  ###########
###############################################################################
# get prior information, i.e. parameters of skew normal distribution transformed scale
    priorparas <- read.table('priordensityparas.txt')    
# calculate logprior: 'dsnorm(x,mean,ds,a)' is density function of skew normal distrib.    
    logprior   <- 0
    if (length(para)>0){
      for (i in 1:length(para)){logprior<-logprior+log(dsnorm(para[names(para[i])],mean=priorparas$mu[priorparas$parameter==names(para[i])],sd=priorparas$sigma[priorparas$parameter==names(para[i])],a=priorparas$a[priorparas$parameter==names(para[i])]))}
    }
# unnormalized logposterior is sum of prior and likelihood
# negative sign as optimizatoin will look for minimum
  logposterior <- -logLLWFS-logLLEFS-logprior; names(logposterior) <- NULL;
  return(logposterior)
}