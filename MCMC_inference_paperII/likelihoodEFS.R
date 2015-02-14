#### unnormalized loglikelihood function with respect to incidence data from EFS
#### loglikelihood also takes ESS adjustment into account

likelihoodEFS <- function(modelpara,vaccpara,coverpara,immupara,info=list('N'=N,
                 'times'=times,'births'=birthcounts,'death'=death,'migra'=migra,'ageing'=agerates,
                 'dataEFS'=datEFS,'alpha'=alpha,'theta'=theta,'contacttype'=contacttype,'essEFS'=essEFS))
# arguments:
# modelpara: includes transmission specific parameters  (see rotavacc.R)
# vaccpara : vaccine effectiveness parameters
# coverpara: coveragerate parameters
# immupara : parameters for age-specific immunity
# info     : list with further fixed model parameters not subject of estimation,
#            i.e. incidence data, demographic parameters, model state N and time horizon
{    
# extract time horizon of the model
      times <- info[['times']]; stoptime <- length(times);

# age specific immunity parameters and computation
      ages <- c(0,1/6,2/6,3/6,4/6,5/6,1,1.5,2,2.5,3,3.5,4,4.5,10,25,45,65,85)
      phirates <- 2*expit(immupara[1]+immupara[2]*ages+immupara[3]*ages^2)
    
### Computation of ODE solution for likelihood computation ##########
#####################################################################    
  # ODE solution using Runge-Kutta 4 scheme (package 'deSolve') 
    sol<-rk(info[['N']],times,rotavacc,list('alpha'=expit(info[['alpha']]),'theta'=expit(info[['theta']]),'phi'=phirates,
                                'mu'=exp(modelpara['mu']),'omega'=exp(modelpara['omega']),
                                'omega0'=exp(modelpara['omega0']),'beta'=exp(modelpara['beta']),
                                'p'=exp(modelpara['p']),
                                'sais'=c(modelpara['a1'],expit(modelpara['b1'])-0.5,modelpara['a2'],expit(modelpara['b2'])-0.25),
                                'contact'=contactmatrix(cpara=exp(modelpara[c('contactpara1','contactpara2','contactpara3')]),pattern=info[['contacttype']]),
                                'births'=info[['births']],
                                'death'=info[['death']],'mig'=info[['migra']],
                                'ageing'=info[['ageing']],'coverage'=vacccoverage(coveragerates=coverpara),
                                'vacceff'=c(expit(vaccpara[1:2]),exp(vaccpara[3]))),
             method='rk4')
      
  # formatting solution array: time point, class, age group
    sol <- sol[times,-1]                  #2:362
    dim(sol)<- c(length(times),dim(info[['N']])[1],dim(info[['N']])[2])                
  # calculate weekly age stratified number of new cases from solution: 
  # (19th col (2nd dim) counts number of occured cases)    
    inci  <- sol[2:length(times),dim(info[['N']])[1],]-sol[1:(length(times)-1),dim(info[['N']])[1],]
       
  # given the weekly case numbers, calculate loglikelihood of the EFS data
  # discard the first ten years from the solution (burn in)
    logLLEFS <- loglikelihoodESS(inci[(10*52):(stoptime-1),],dat=info[['dataEFS']],hosp=expit(modelpara['h1']),fact=exp(modelpara['fact1']),disp=exp(modelpara['disp']),ESS=info[['essEFS']])
  return(logLLEFS)
}