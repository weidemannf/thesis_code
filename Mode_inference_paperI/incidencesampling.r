### compute an incidence sample from the mixture of considered models for the EFS 
### and WFS from 2001 to 2008

############# import necessary functions and packages ##########################
# functions which provide incidence data, demographic data or starting conditions   
    source(file="demographicdata.R")    
# further functions providing, e.g., the posterior, likelihood, or ODE system     
    source(file="rota.R")           
    source(file="functions.R")    
# packages    
    require('deSolve')
    require('mnormt')
    
#### Setting of initial conditions and parameters #######################
#########################################################################
    begTime <- Sys.time()   # time measuring        
    N <- initial()          # initial condition of dynamic system 
    alpha <- c(0.6, 0.667)                         
    theta <- c(0.5, 0.25, 0)                           
    mu <- log(7/3)                                        
    omega <- log(7/8)                                     
    omega0 <- log(1/26)                      
    beta <- log(1/52)                        
    p <- log(1)                         
    sais <- c(0, logit(0.5),0 ,logit(0.5))     
    contactpara <- log(c(1,1,1))   
    contacttype <- 4  
    h <- logit(c(0.1, 0.1))          
    fact <- log(c(1, 1))             
    disp <- log(1)             

# demographic parameters    
    birthcounts <- births()
    death <-deathrate()
    agerates <-ageing()
    migra <- migration() 

# model time horizon    
    Stoptime  <- 19*52
    times     <- 1:Stoptime
# model time horizon for which data is available (exclude first eleven years)         
    lastyears <- (11*52):(Stoptime-1)

# import incidence data and population size data
    popEFS  <- read.table(file="PopulationEFS.txt")
    poptotal <- read.table(file="Populationtotal.txt")    
  
### fetch the considered ensemble of models  ##################################
# the model ensemble consists of models including or excluding parameters for 
# the duration of infection (mu) or (mu,omega) (three possibilities) as well as 
# accounting for six different contact structures. thus, overall possible 18 models 
    paralist <- list();  covarlist <-list();  

    for (contacttype in 1:6){
# results for models using contacttype one use different names
      if (contacttype == 1){
         direct1 <- "_omega0_beta_p_a1_b1_a2_b2_contactpara1_h1_h2_fact1_fact2_disp_ctype_1_ESS.Rdata"
         direct2 <- "_omega_omega0_beta_p_a1_b1_a2_b2_contactpara1_h1_h2_fact1_fact2_disp_ctype_1_ESS.Rdata"
         direct3 <- "_mu_omega_omega0_beta_p_a1_b1_a2_b2_contactpara1_h1_h2_fact1_fact2_disp_ctype_1_ESS.Rdata"
      }else{
         direct1 <- paste("_omega0_beta_p_a1_b1_a2_b2_contactpara1_contactpara2_contactpara3_h1_h2_fact1_fact2_disp_ctype",as.character(contacttype),"ESS.Rdata",sep = "_")
         direct2 <- paste("_omega_omega0_beta_p_a1_b1_a2_b2_contactpara1_contactpara2_contactpara3_h1_h2_fact1_fact2_disp_ctype",as.character(contacttype),"ESS.Rdata",sep = "_")
         direct3 <- paste("_mu_omega_omega0_beta_p_a1_b1_a2_b2_contactpara1_contactpara2_contactpara3_h1_h2_fact1_fact2_disp_ctype",as.character(contacttype),"ESS.Rdata",sep = "_")
      }
      try(load(file = direct1))
      paralist <- c(paralist,list(para)); covarlist <- c(covarlist,list(solve(hess)));
      try(load(file = direct2))
      paralist <- c(paralist,list(para)); covarlist <- c(covarlist,list(solve(hess)));
      try(load(file = direct3))
      paralist <- c(paralist,list(para)); covarlist <- c(covarlist,list(solve(hess)));
    }
      
# fetch model weights computed in advance
    try(load(file = 'weights.Rdata'))
   
### Start sampling the incidences using the sample of models and the corresponding 
### posterior distributions (each model has its own posterior) ##################
# set the number K of incidence time series to be drawn, i.e. for each age group 
# and region 
    K=400
# define the arrays of sampled incidences, once for the expected number of cases 
# from the model (aggreginci) and once for the reported incidence subject to the NegBin
# distribution (aggregsample), stratified by EFS (NB) and WFS (AB); 
# dimension = (number of age groups (10) , number of data weeks (8*52), K)   
    aggreginciWFS   <- array(0,dim=c(10,length(lastyears),K))
    aggreginciEFS   <- array(0,dim=c(10,length(lastyears),K))
    aggregsampleWFS <- array(0,dim=c(10,length(lastyears),K))
    aggregsampleEFS <- array(0,dim=c(10,length(lastyears),K))

# start sampling loop
for (k in 1:K){
# first sample the model according to model probabilities
# model order: first all models with contacttype 1, then contacttype 2, ...
# within these: first excluding mu and omega, then excluding only mu, then including both
   model <- sample(1:18,size=1,prob=totalweights)
   print(c(k,model))
   
# sample a parameter vector from the posterior distribution (corresponding to the sampled model)
   para  <- c(rmnorm(n = 1, mean = paralist[[model]],	varcov=covarlist[[model]]))
   names(para) <- names(paralist[[model]]) 

# retransform the sample parameter vector 
   mu <- 7/3      # these get overwritten if 'mu' or 'omega' are contained in 
   omega <- 7/8   # the sampled parameter vector
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

# define next generation matrix accorrding to sampled contacttype
   contacttype <- ((model-1)%/%3)+1
   contact <- contactmatrix(cpara=contactpara,pattern=contacttype)

### Compute solution of the ODE system using Runge-Kutta 4 scheme #########    
   sol<-rk(N,times,rota,list('alpha'=alpha,'theta'=theta,'mu'=mu,'omega'=omega,
         'omega0'=omega0,'beta'=beta,'p'=p,'sais'=sais,'contact'=contact,
         'births'=birthcounts,'death'=death,'mig'=migra,'ageing'=agerates),method = "rk4")     

# formatting of ODE solution array: time point, class, age group
    sol <- sol[times,2:(length(N)+1)]; dim(sol)<- c(length(times),dim(N)[1],dim(N)[2]);
# calculate weekly age stratified number of new cases from solution: 
# (15th class counts number of occured cases)
    inci  <- sol[2:length(times),dim(N)[1],1:dim(N)[2]]-sol[1:(length(times)-1),dim(N)[1],1:dim(N)[2]]

# accumulate weekly case expectation samples from model to match the data age resolution
    incistrata <- numeric(0);
    agestrata  <- list(1:6,7:8,9:10,11:12,13:14,15,16,17,18,19)
    for (ag in agestrata){incistrata <- rbind(incistrata,.rowSums(inci[,ag],m=dim(inci)[1],n=length(ag)))}
   
### calculate the decomposition into expected cases in the EFS and WFS
# define matrices for the expected number of cases (AB = WFS, NB = EFS)
   incistrataWFS <- matrix(0, nrow = 10, ncol=length(lastyears))
   incistrataEFS <- matrix(0, nrow = 10, ncol=length(lastyears))
# define the time dependent vector for the reporting rates
   vec_hospWFS <- numeric(0); vec_hospEFS <- numeric(0);
# calculate the expected number of reported cases in each (age group, region)
   for (i in 1:length(lastyears)){
# reporting rates depend on the baseline reporting rate (h) and the rel. incr. (fact)
     vec_hospWFS[i]<- h[2]+(h[2]*(fact[2]-1))*min(max(c(((i-(52*4))/26),0)),1)
     vec_hospEFS[i]<- h[1]+(h[1]*(fact[1]-1))*min(max(c(((i-(52*4))/26),0)),1)
# expected incidence is decomposed according to population size in EFS and WFS
     incistrataWFS[1:10,i]<-vec_hospWFS[i]*incistrata[1:10,i]*(1-(popEFS[[(1+floor((i-1)/52))]]/poptotal[[(1+floor((i-1)/52))]]))
     incistrataEFS[1:10,i]<-vec_hospEFS[i]*incistrata[1:10,i]*(popEFS[[(1+floor((i-1)/52))]]/poptotal[[(1+floor((i-1)/52))]])
   }

### sample reported incidences according to NegBin observational distribution    
   sampleinciWFS   <- rnbinom(n=10*length(lastyears),mu=incistrataWFS,size=disp)
   dim(sampleinciWFS) <- dim(incistrataWFS)
   sampleinciEFS   <- rnbinom(n=10*length(lastyears),mu=incistrataEFS,size=disp)
   dim(sampleinciEFS) <- dim(incistrataEFS)
### store the sample results   
   aggreginciWFS[,,k]   <- incistrataWFS
   aggreginciEFS[,,k]   <- incistrataEFS
   aggregsampleWFS[,,k] <- sampleinciWFS
   aggregsampleEFS[,,k] <- sampleinciEFS 
}

### save the samples and print runtime
   save(aggreginciWFS,aggreginciEFS,aggregsampleWFS,aggregsampleEFS,file=paste("incidencesample.Rdata"))
   runTime <- Sys.time()-begTime
   print(runTime)

