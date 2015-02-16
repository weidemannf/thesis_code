### Bayesian parameter inference for for variable parameter spaces
### inference is based on finding posterior mode and the corresponding fisher information
### posterior distribution is then defined as multivariate normal using these quantities
    
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
    contacttype <- 4  
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
      
# import cumulative autocorrelation measures for each of the 20 time series
# for optimization not accounting for autocorrelation, the vectors must be set to one
    essdata.ia <- try(load(file="meaness.Rdata"))
    if (class(essdata.ia)=='character'){
      essEFS<-meanessEFS
      essWFS<-meanessWFS
    }else{
      essEFS<-rep(1,10)
      essWFS<-rep(1,10)
    }

# define the full information list required by the posterior function
    infolist <- list('N'=N,'times'=times,'births'=birthcounts,'death'=death,'dataEFS'=datEFS,
             'dataWFS'=datWFS,'mu'=mu,'omega'=omega,'alpha'=alpha,'theta'=theta,
             'omega0'=omega0,'beta'=beta,'p'=p,'sais'=sais,'contactpara'=contactpara,
             'h'=h,'fact'=fact,'disp'=disp,'contacttype'=contacttype,'migra'=migra,
             'ageing'=agerates,'essEFS'=essEFS,'essWFS'=essWFS)

################# start the inference procedure####################  ###########
################################################################################
# define parameters which are subject of estimation 
# parameters included in this vector will be the variables of the posterior function
# other parameters contained in the information list will be considered fixed values
# when using contactpattern 1: Use 'contactpara1'= contactpara[1] in the parameter vector      
    para <- c('omega0'=omega0,'beta'=beta,'p'=p,'a1'=sais[1],
            'b1'=sais[2],'a2'=sais[3],'b2'=sais[4],'contactpara'=contactpara,'h'=h,
            'fact'=fact,'disp'=disp)
    
# compute posterior value at the initial parameter vector  
    pointpost <- paraposteriorESS(para,info=infolist)
    print(c('initialposterior: ',pointpost))
                                                                                             
    
### if inference (posterior optimisation) for this parameter without considering 
### autocorrelation was already conducted, fetch these results
### if this vector is available, it will be our new starting point for maximisation 
    direct<-character(0)
    for (i in 1:length(names(para))){direct<-paste(direct,names(para)[i],sep = "_")}
    direct<-paste(direct,"ctype",as.character(contacttype),".Rdata",sep = "_")
    try(load(file=direct))
# compute the posterior value for this vector by also considering autocorrelation     
    pointpost <- paraposteriorESS(para,info=infolist)
    print(c('initialposterior: ',pointpost))

### if inference (posterior optimisation) for this parameter including  
### autocorrelation was also already conducted, fetch these results
### if this vector is available, instead this will be our new starting point for maximisation       
    direct<-character(0)
    for (i in 1:length(names(para))){direct<-paste(direct,names(para)[i],sep = "_")}
    direct<-paste(direct,"ctype",as.character(contacttype),"ESS.Rdata",sep = "_")
    try(load(file=direct))
    
# compute posterior for this parameter vector, accounting for autocorrelation
    pointpost <- paraposteriorESS(para,info=infolist) 
    print(c('initialposterior: ',pointpost))

### start the posterior maximisation procedure #################################
# the hybrid scheme consists of a Nelder-Mead optimisation (capped at a maximum 
# number of iterations) and afterwards optimisation along the gradient of the 
# posterior at the reached parameter vector
# set the number K of iteration steps for the hybrid optimisation scheme
    K <- 20
    for (i in 1:K){     
# use Nelder-Mead to optimise posterior function, only varying parameters contained 
# in 'para'; the algorithm is capped at 2000 iteration;
      optpara <- optim(para,paraposteriorESS,method='Nelder-Mead',info=infolist,control=list('trace'=1,maxit=2000,reltol=1e-15))
      para  <- optpara[[1]]
# evaluate the posterior at the 'new' best vector
      pointpost <- paraposteriorESS(para,info=infolist)
# compute and normalize the posterior gradient at the new vector
      vect <- grad(paraposteriorESS,para,method='Richardson',method.args=list(r=6),info=infolist)
      vect <- vect/(10000*mean(abs(vect)))
# compute the optimum of the one-dimesional function g(x)=posterior(point-x*vect),
# where point is the so far optimal value and vect the normalized gradient
# the function 'optimgrad' computes the optimal value for 'x', using the posterior 
# function, the best known parameter vector (para) and its gradient (vect) as args.
# usage of vect/10000 yields faster convergence of 'optimgrad' 
      optval <- optimgrad(paraposteriorESS,para,vect,info=infolist)
      print(c('multiplicator: ',optval[[1]]))
      newpara <- para-(optval[[1]]*vect)
      # Auswertung der Posterior an diesem neuen Optimum
      newpointpost <- paraposteriorESS(newpara,info=infolist)
# if 'para' was already near to the optimum, 'optimgrad' does not always improve the 
# parameter vector, thus actual improvement is checked
      if (newpointpost < pointpost){
        para <- newpara
        pointpost <- newpointpost
      }
      print(c(i,'. finalposterior: ',pointpost))
    }

### compute the Fisher information, i.e. minus the hessian of the logposterior, at
### the optimal parameter vector    
    hess <- hessian(paraposteriorESS, para, method="Richardson",method.args=list(r=6),info=infolist)
      
### save results, filename determined by parameter subject of estimation  
    direct<-character(0)
    for (i in 1:length(names(para))){direct<-paste(direct,names(para)[i],sep = "_")}
    if (class(essdata.ia)=='character'){
    direct<-paste(direct,"ctype",as.character(contacttype),"ESS.Rdata",sep = "_")
    }else{
    direct<-paste(direct,"ctype",as.character(contacttype),".Rdata",sep = "_")
    }
    save(para,hess,file=direct)

### print runtime 
    runTime <- Sys.time()-begTime
    print(runTime)                            