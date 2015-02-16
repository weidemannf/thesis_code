### functions for computing the Anscombe residuals of NegBin distributed data
### for more details see Hardin and Hilbe: Generalized Linear Models and Extensions 
V <- function(mu,r){
return(((mu^2)/r+mu))
}
integrand <- function(mu,r){
return(((mu^2)/r+mu)^(-1/3))
}
A <- function(y,mu=1,r=1){
A <- integrate(integrand,lower=mu,upper=y,r=r)
return(A[[1]])
}
a.resi <- function(y,mu=1,r=1){
return(A(y,mu,r)/(V(mu,r)^(1/6)))
}


### anscomberesiduals computes the Anscombe resiuduals of the incidence data with respect
### to the ODE model output at the posterior mode.
### arguments:
### para = vector of posterior function arguments
###  , i.e. parameters which are subject of estimation
### info = list of all required arguments for solving the model ODE and data (including those listed in para)
anscomberesiduals <- function(para,info=list(N,times,birthcounts,death,dataWFS,dataEFS,alpha,theta,mu,omega,
                                          omega0,beta,p,sais,contactpara,h,fact,disp,contacttype,migra,agerates))
{
# extract parameters and data from the information list
    # initial decompostion of population
      N <- info[['N']] 
    # time horizon of the model
      times <- info[['times']]; stoptime <- length(times);  
    # incidence data from EFS (NB) and WFS (AB)
      datWFS <- info[['dataWFS']]
      datEFS <- info[['dataEFS']]  
      datayears <- (11*52-1)+(1:dim(datWFS)[2])
    # further arguments as listed in rota.R
    # parameters are also retransformed on (0,1) or (0,Inf) where necessary      
      birthcounts <- info[['births']]
      death  <- info[['death']]
      ageing <- info[['ageing']]
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
         'births'=birthcounts,'death'=death,'mig'=migra,'ageing'=ageing),method = "rk4")
# formatting solution array: time point, class, age group
    sol <- sol[times,2:(length(N)+1)]; dim(sol)<- c(length(times),dim(N)[1],dim(N)[2]);
# calculate weekly age stratified number of new cases from solution: 
# (15th class counts number of occured cases)
    inci  <- sol[2:length(times),dim(N)[1],1:dim(N)[2]]-sol[1:(length(times)-1),dim(N)[1],1:dim(N)[2]]

#### cumulate incidences according to data and decompose by region #####
    incistrata <- numeric(0)
    agestrata  <- list(1:6,7:8,9:10,11:12,13:14,15,16,17,18,19)
    for (i in 1:length(agestrata)){
      if (length(agestrata[[i]])>1){
         incistrata <- rbind(incistrata,rowSums(inci[datayears,agestrata[[i]]]))
      }else{
         incistrata <- rbind(incistrata,inci[datayears,agestrata[[i]]])
      }
    }

# get population sizes from EFS and all German states (year specific since 2001)
    popEFS   <- read.table(file="PopulationEFS.txt")
    poptotal <- read.table(file="Populationtotal.txt")
 
# compute age and region specific model output  (see loglikelihoodESS.R)
    incistrataWFS <- matrix(0,nrow=10,ncol=length(datayears))
    incistrataEFS <- matrix(0,nrow=10,ncol=length(datayears))
    vec_hospWFS <- rep(0,length(datayears)); vec_hospEFS <- rep(0,length(datayears))
    for (i in 1:length(datayears)){
       vec_hospWFS[i]<- h[2]+(h[2]*(fact[2]-1))*min(max(c(((i-(52*4))/26),0)),1)
       vec_hospEFS[i]<- h[1]+(h[1]*(fact[1]-1))*min(max(c(((i-(52*4))/26),0)),1)

       incistrataWFS[,i]<-vec_hospWFS[i]*incistrata[,i]*(1-(popEFS[[(1+floor((i-1)/52))]]/poptotal[[(1+floor((i-1)/52))]]))
       incistrataEFS[,i]<-vec_hospEFS[i]*incistrata[,i]*(popEFS[[(1+floor((i-1)/52))]]/poptotal[[(1+floor((i-1)/52))]])
     }

# compute and center Anscombe residuals
    ResidualsWFS <- datWFS*0;
    ResidualsEFS <- datEFS*0;     
    for (i in 1:length(datayears)){
      for(j in 1:10){
         ResidualsWFS[j,i] <- a.resi(datWFS[j,i],incistrataWFS[j,i],disp)
         ResidualsEFS[j,i] <- a.resi(datEFS[j,i],incistrataEFS[j,i],disp)
      }
    }
    ResidualsWFS <- ResidualsWFS - rowMeans(ResidualsWFS);
    ResidualsEFS <- ResidualsEFS - rowMeans(ResidualsEFS);    
# return the matrices of residuals   
    return(list(ResidualsWFS,ResidualsEFS))
}
     
