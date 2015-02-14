### compute an incidence sample from model for the EFS incidence following introduction
### of routine vaccination from 1999 up to a prespecified time horizon

############# import necessary functions and packages ##########################
# functions which provide incidence data, demographic data or starting conditions 
source(file="demographicdataEFS.R")
source(file="projectdemographics.R")
# further functions providing, e.g., the likelihood or ODE system
source(file="functions.R")
source(file="rotavacc.R")
# packages 
require('deSolve')

#### Setting of initial conditions and fixed parameters #######################
###############################################################################
      N <- initialEFS()        # initial condition of ODE-system 
      alpha <- c(0.6, 0.667)
      theta <- c(0.5, 0.25, 0)
      contacttype <- 6
  # demographic parameters
      birthcounts <- birthsEFS()
      death <- deathrate()
      agerates <-ageing()
      migra <- migrationEFS()    
  # define prediction time horizon   
      periods   <- 40             
      stoptime  <- periods*52;     times <- 1:stoptime;
      plotperiods <- periods-10;   lastyears <- (stoptime-plotperiods*52+1):stoptime;
  # set contactpattern
      contacttype <- 6
  # fetch results from the posteior sampling 
      load(file='adaptivemcmcresults.Rdata')
    
### Initilizing the predictive incidence sampling  ###########################  
  # set the vector of coverage levels in [0,1] which should be investigated
      coverage.levels <- c(0.5,0.9)
  # set sample size for each investigated long term coverage level
      K <- 100
  # define index vector for drawing from the posterior sample    
      index <- round(dim(paramatrix1)[2]*(1:min(K,dim(paramatrix1)[2]))/min(K,dim(paramatrix1)[2]))
    
  # sampling future death, migration and fertitlity rates
      prodeathsample  <- array(0,dim=c(dim(N)[2],periods,length(index)))                          
      promigrasample  <- array(0,dim=c(dim(N)[2],periods,length(index)))
      fertratessample <- array(0,dim=c(stoptime,length(index))) 
      for (i in 1:length(index)){
          prodeathsample[,,i] <- projectdeath(olddeaths=death,timehorizon=periods)                            
          promigrasample[,,i] <- projectmigra(oldmigras=migra,timehorizon=periods)
          fertratessample[,i] <- weeklyfertility(timehorizon=periods)
      }

### Start of the samping procedure ############################################
### Step1: Loop over the vector of coverage levels
for (cr in coverage.levels){
# define storage matrices
  storeexpcases    <- array(0,dim=c(10,52*plotperiods,length(index)))
  storereportcases <- array(0,dim=c(10,52*plotperiods,length(index)))
  storepops        <- array(0,dim=c(10,52*plotperiods,length(index)))

### Loop over the sample size for each coverage level
  for (i in 1:length(index)){
    print(c(cr,i))
# fetch parameter from the posterior sample
    modelpara <- paramatrix1[,index[i]] 
    names(modelpara) <- c('mu','omega','omega0','beta','p','a1','b1','a2','b2',
              'contactpara1','contactpara2','contactpara3','h1','fact1','disp')
    vaccpara  <- paramatrix2[,index[i]]
    coverpara <- paramatrix3[,index[i]]
    immupara  <- paramatrix4[,index[i]]

# compute age specific immunities 
    ages <- c(seq(0,1,by=1/6),seq(1.5,4.5,by=0.5),10,25,45,65,85)
    phirates <- 2*expit(immupara[1]+immupara[2]*ages+immupara[3]*ages^2)     

# compute coverage process, i.e. weekly coverage rates
    coverage <- vacccoverage(coveragerates=coverpara)
    lastcover <- coverage[length(coverage)]
    if (stoptime > 52*24){coverage <- c(coverage, lastcover+0.5*(1-cos((1:104)*pi/104))*(cr-lastcover))}
    if (stoptime > 52*25){coverage <- c(coverage, rep(cr,stoptime-52*25))}     

### compute solution of the ODE-system according to the parameters #############
### from the posterior sample                                      #############
    sol<-rk(N,times,rotavacc,list('alpha'=alpha,'theta'=theta,'phi'=phirates,
                                'mu'=exp(modelpara['mu']),'omega'=exp(modelpara['omega']),
                                'omega0'=exp(modelpara['omega0']),'beta'=exp(modelpara['beta']),
                                'p'=exp(modelpara['p']),
                                'sais'=c(modelpara['a1'],expit(modelpara['b1'])-0.5,modelpara['a2'],expit(modelpara['b2'])-0.25),
                                'contact'=contactmatrix(cpara=exp(modelpara[c('contactpara1','contactpara2','contactpara3')]),pattern=contacttype),
                                'births'=birthcounts,'birthrates'=fertratessample[,i],
                                'death'=prodeathsample[,,i],'mig'=promigrasample[,,i],
                                'ageing'=agerates,'coverage'=coverage,
                                'vacceff'=c(expit(vaccpara[1:2]),exp(vaccpara[3]))),
             method='rk4')

# formatting of solution array
    sol <- sol[times,-1];  dim(sol)<- c(length(times),dim(N)[1],dim(N)[2]);
                    
### compute the simulated population counts and the overall number of cases 
# from the solution array and store the population results
    solperm <- aperm(sol[lastyears,1:(dim(N)[1]-1),], perm=c(1,3,2))
    modelpopulation <- rowSums(solperm,dims=2)
    inci  <- sol[lastyears,dim(N)[1],]-sol[lastyears-1,dim(N)[1],]    
    population <- numeric(0) 
    incistrata <- numeric(0)
    aggregagegroups <- list(1:6,7:8,9:10,11:12,13:14,15,16,17,18,19)
    for (ag in aggregagegroups){   
      population <- rbind(population,.rowSums(modelpopulation[,ag],m=dim(modelpopulation)[1],n=length(ag)))
      incistrata <- rbind(incistrata,.rowSums(inci[,ag],m=dim(inci)[1],n=length(ag)))
    }    
    storepops[,,i] <- population

### compute the simulated expected and observed numbers of reported case 
# set time dependent detection rate
    hosp <- expit(modelpara['h1']); names(hosp) <- NULL;
    fact <- exp(modelpara['fact1']); names(fact) <- NULL;
    vec.rep <- c(rep(hosp,(4*52)-1),
                 seq(hosp,hosp*fact,by=(hosp*(fact-1))/26),
                 rep(hosp*fact,dim(inci)[1]-4.5*52))

# compute expectation and observation sample                  
    expcases <- t(vec.rep*t(incistrata[,]))
    storeexpcases[,,i] <- expcases
    reportcases <- rnbinom(n=length(expcases),mu=expcases,size=exp(modelpara['disp']))
    dim(reportcases) <- dim(expcases)
    storereportcases[,,i] <- reportcases
}
### store results for each coverage level
save(storereportcases,file=paste('incidencesampling',cr,'.Rdata'))
save(storeexpcases,file=paste('expincisampling',cr,'.Rdata'))
save(storepops,file=paste('populationsample',cr,'.Rdata'))
}