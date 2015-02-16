### Computation of the effective sample sizes of each time series for each model
### averaging of the age and region specific effective sample sizes over all considered
### models

# import necessary data, functions and packages
    source(file="demographicdata.R");
    source(file="functions.R");
    source(file="anscomberesiduals.R");
    source(file="effectivesamplesize.R");
# packages
    require('optimx');
    require('deSolve'); 
    require('polynom');

#### Setting of initial conditions and fixed parameters #######################
###############################################################################
    begTime <- Sys.time()   # time measuring
    N <- initial()          # initial condition of dynamic system
    alpha <- logit(c(0.6, 0.667))                          
    theta <- logit(c(0.5, 0.25))                           
    mu <- log(7/3)                                        
    omega <- log(7/8)                                     
    omega0 <- log(1/26)                      
    beta <- log(1/52)                        
    p <- log(1)                         
    sais <- c(0, logit(0.5),0 ,logit(0.5))     
    contactpara <- log(c(1,1,1))      
    h <- logit(c(0.1, 0.1))          
    fact <- log(c(1, 1))             
    disp <- log(1)                                

# demographic parameters    
    birthcounts <- births() 
    death <-deathrate()
    agerates <-ageing()
    migra <- migration()   

# import incidence data
    datWFS <- WFSdata(n=8)
    datEFS <- EFSdata(n=8)
      
# time horizon    
    Stoptime  <- 19*52
    times     <- 1:Stoptime

### now compute stratified effective sample sizes for all models based on the 
### respective posterior modes using a preceding Anscombe transformation of the residuals 
allmodels_essWFS <- array(0,dim=c(3,6,10)); allmodels_essEFS <- array(0,dim=c(3,6,10));
       
for (contacttype in 1:6){
    for (paraspace in c('full','no_mu','no_mu_omega')){
# define the full information list required by the ODE-model  
      infolist <- list('N'=N,'times'=times,'births'=birthcounts,'death'=death,'dataEFS'=datEFS,'dataWFS'=datWFS,'mu'=mu,'omega'=omega,
                   'alpha'=alpha,'theta'=theta,'omega0'=omega0,'beta'=beta,'p'=p,'sais'=sais,'contactpara'=contactpara,
                   'h'=h,'fact'=fact,'disp'=disp,'contacttype'=contacttype,'migra'=migra,'ageing'=agerates)      
# fetch optimal parameter vector of the model      
      if (contacttype == 1){
         if (paraspace=='no_mu_omega')   direct <- "_omega0_beta_p_a1_b1_a2_b2_contactpara1_h1_h2_fact1_fact2_disp_ctype_1_.Rdata"
         if (paraspace=='no_mu')   direct <- "_omega_omega0_beta_p_a1_b1_a2_b2_contactpara1_h1_h2_fact1_fact2_disp_ctype_1_.Rdata"
         if (paraspace=='full') direct <- "_mu_omega_omega0_beta_p_a1_b1_a2_b2_contactpara1_h1_h2_fact1_fact2_disp_ctype_1_.Rdata"
      }else{
         if (paraspace=='no_mu_omega')   direct <- paste("_omega0_beta_p_a1_b1_a2_b2_contactpara1_contactpara2_contactpara3_h1_h2_fact1_fact2_disp_ctype",as.character(contacttype),".Rdata",sep="_")
         if (paraspace=='no_mu')   direct <- paste("_omega_omega0_beta_p_a1_b1_a2_b2_contactpara1_contactpara2_contactpara3_h1_h2_fact1_fact2_disp_ctype",as.character(contacttype),".Rdata",sep="_")
         if (paraspace=='full') direct <- paste("_mu_omega_omega0_beta_p_a1_b1_a2_b2_contactpara1_contactpara2_contactpara3_h1_h2_fact1_fact2_disp_ctype",as.character(contacttype),".Rdata",sep="_")
      }
      try(load(file=direct))
# compute effective sample size (cumulative autocorrelation) of the model
      resids <- anscomberesiduals(para,info=infolist)
      for (i in 1:10){
        allmodels_essWFS[(1:3)[paraspace==c('full','no_mu','no_mu_omega')],(1:6)[contacttype==1:6],i] <- min(1,effectivesamplesize(resids[[1]][i,])[['ESS']])
        allmodels_essEFS[(1:3)[paraspace==c('full','no_mu','no_mu_omega')],(1:6)[contacttype==1:6],i] <- min(1,effectivesamplesize(resids[[2]][i,])[['ESS']])}
    }
}
### compute and store mean results for the time series specific effective sample sizes ###
    meanessWFS <- rep(0,10); meanessEFS <- rep(0,10)
    for (i in 1:10){
      meanessWFS[i] <- mean(allmodels_essWFS[,,i]) 
      meanessEFS[i] <- mean(allmodels_essEFS[,,i])
    }
    save(meanessWFS,meanessEFS,file='meaness.Rdata')

      
