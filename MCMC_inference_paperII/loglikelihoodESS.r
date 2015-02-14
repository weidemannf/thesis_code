# this function computes the likelihood of observed incidence data in the EFS
# given the expected case numbers computed the model
# assumes negative binomial distribution of the data 

loglikelihoodESS <- function(inci,dat,hosp,fact,disp,ESS)
{
# inci: model predicted age stratified weekly cases (since 2001, 19 age groups according)
# dat:  data on weekly case series (since 2001), 10 age groups due to data resolution
# hosp: baseline reproting rate in EFS
# fact: relative reporting rate increase in 2005 in EFS
# disp: overdispersion of observational ditribution
# ESS:  vector of the estimated inverse of the cumulative autocorrelation of each data series
#       (i.e. ten components for ten age groups in the EFS data), (measure of effective sample size)

# compute minimal time horizon of available case and model predictions
    datsize <- min(dim(inci)[1],dim(dat)[2])
# accumulate weekly case predicton from model to match the data age resolution
    incistrata <- numeric(0)
    agestrata  <- list(1:6,7:8,9:10,11:12,13:14,15,16,17,18,19)
    for (ag in agestrata){incistrata <- rbind(incistrata,.rowSums(inci[,ag],m=dim(inci)[1],n=length(ag)))} 
### calculate the expected number of reported cases per week and age group ##
#############################################################################   
# week specific reporting rate in the EFS (i=1 is first week of 2001)
    vec.rep <- c(rep(hosp,(4*52)-1),
                 seq(hosp,hosp*fact,by=(hosp*(fact-1))/26),
                 rep(hosp*fact,datsize-4.5*52))
# calculate model-predicted reported number of cases in EFS 
    repinci <- vec.rep*t(incistrata[1:10,1:datsize])

#### Calculate loglikelihood of the data utilizing the cumulative autocorrelations
# 0 < ESS < 1; ESS = (cumulative autocorrelation)^(-1); specified for each
# assuming NegBin distribution of the data with mean coming from the model  
    logLL <- sum(ESS*dnbinom(x=dat[1:10,1:datsize],mu=t(repinci),size=disp,log = TRUE))
    return(logLL)
}