### Computation of the posterior model probabilities (weights) for each model  ###
##################################################################################

############# import necessary functions and packages ##########################
# functions which provide incidence data, demographic data or starting conditions   
    source(file="demographicdata.R")    
# further functions providing, e.g., the posterior, likelihood, or ODE system    
    source(file="loglikelihoodESS.R")
    source(file="rota.R")    
    source(file="paraposteriorESS.R")    
    source(file="functions.R")    
# packages    
    require('optimx')
    require('deSolve')

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
    contacttype <- 1   
    h <- logit(c(0.1, 0.1))          
    fact <- log(c(1, 1))             
    disp <- log(1)                                

# demographic parameters    
    geburt <- geburten() 
    death <-deathrate()
    agerates <-ageing()
    migra <- migration()   

# import incidence data
    datAB <- ABDaten(n=8)
    datNB <- NBDaten(n=8)
      
# time horizon    
    Stoptime  <- 19*52
    times     <- 1:Stoptime
      
# import cumulative autocorrelation measures for each of the 20 time series
    essdata.ia <- try(load(file="meaness.Rdata"))
    if (class(essdata.ia)=='character'){
    essNB<-meanessNB
    essAB<-meanessAB
    }else{
    essNB<-rep(1,10)
    essAB<-rep(1,10)
    }

# define the full information list required by the posterior function
    infolist <- list('N'=N,'times'=times,'geburt'=geburt,'death'=death,'datenNB'=datNB,
             'datenAB'=datAB,'mu'=mu,'omega'=omega,'alpha'=alpha,'theta'=theta,
             'omega0'=omega0,'beta'=beta,'p'=p,'sais'=sais,'contactpara'=contactpara,
             'h'=h,'fact'=fact,'disp'=disp,'contacttype'=contacttype,'migra'=migra,
             'ageing'=agerates,'essNB'=essNB,'essAB'=essAB)

### now compute marginal likelihoods for all models based on the 
### assumed normality of the posterior distribution 
marginallikelihoods <- matrix(0,nrow=3,ncol=6); 
       
for (contacttype in 1:6){
    for (paraspace in c('full','no_mu','no_mu_omega')){
# fetch posterior mode, respective value and hessian of the model      
      if (contacttype == 1){
         if (paraspace=='no_mu_omega')   direct <- "_omega0_beta_p_a1_b1_a2_b2_contactpara1_h1_h2_fact1_fact2_disp_ctype_1_ESS.Rdata"
         if (paraspace=='no_mu')   direct <- "_omega_omega0_beta_p_a1_b1_a2_b2_contactpara1_h1_h2_fact1_fact2_disp_ctype_1_ESS.Rdata"
         if (paraspace=='full') direct <- "_mu_omega_omega0_beta_p_a1_b1_a2_b2_contactpara1_h1_h2_fact1_fact2_disp_ctype_1_ESS.Rdata"
      }else{
         if (paraspace=='no_mu_omega')   direct <- paste("_omega0_beta_p_a1_b1_a2_b2_contactpara1_contactpara2_contactpara3_h1_h2_fact1_fact2_disp_ctype",as.character(contacttype),"ESS.Rdata",sep="_")
         if (paraspace=='no_mu')   direct <- paste("_omega_omega0_beta_p_a1_b1_a2_b2_contactpara1_contactpara2_contactpara3_h1_h2_fact1_fact2_disp_ctype",as.character(contacttype),"ESS.Rdata",sep="_")
         if (paraspace=='full') direct <- paste("_mu_omega_omega0_beta_p_a1_b1_a2_b2_contactpara1_contactpara2_contactpara3_h1_h2_fact1_fact2_disp_ctype",as.character(contacttype),"ESS.Rdata",sep="_")
      }
      try(load(file=direct))
# compute log marginal likelihood of the model
      infolist$contacttype <- contacttype    
      pointpost <- paraposteriorESS(para,info=infolist)      
      mineiwert <- min(eigen(hess)[[1]])   
      if (mineiwert>0){
      A <- solve(hess)                                                                       
      marginallikelihoods[(1:3)[paraspace==c('full','no_mu','no_mu_omega')],contacttype] <- length(para)/2*log(2*pi)+log(det(A))/2-pointpost
      }else{
      print('Error: Hessian not positive definite')
      marginallikelihoods[(1:3)[paraspace==c('full','no_mu','no_mu_omega')],contacttype] <- -Inf}       
    }
}
### compute and store probability weights for each model ###
   marginallikelihoods <- exp(marginallikelihoods-max(marginallikelihoods))
   totalweights <- marginallikelihoods/sum(marginallikelihoods)
   save(totalweights,file='weights.Rdata')
   

      
