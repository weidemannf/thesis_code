### these functions compute the loglikelihood of observed incidence data in the EFS
### and WFS given the expected case numbers computed the model
### assumes negative binomial distribution of the data 

loglikelihoodESS <- function(inci,dat,hosp,fact,disp,ESS,region)
{
# inci: model predicted age stratified weekly cases (since 2001, 19 age groups according)
# dat:  data on weekly case series (since 2001), 10 age groups due to data resolution
# hosp: baseline reproting rate
# fact: relative reporting rate increase in 2005
# disp: overdispersion of observational ditribution
# ESS:  vector of the estimated inverse of the cumulative autocorrelation of each data series
#       (i.e. ten components for ten age groups in the data), (measure of effective sample size)
# region: character determining the corresponding region of the data, either 'AB' or 'NB'

# compute minimal time horizon of available case and model predictions
    datsize <- min(dim(inci)[1],dim(dat)[2])
# accumulate weekly case predicton from model to match the data age resolution
    incistrata <- numeric(0)
    agestrata  <- list(1:6,7:8,9:10,11:12,13:14,15,16,17,18,19)
    for (ag in agestrata){incistrata <- rbind(incistrata,.rowSums(inci[,ag],m=dim(inci)[1],n=length(ag)))}

# get population sizes from EFS and all German states (year specific since 2001)
    popEFS   <- read.table(file="PopulationEFS.txt")
    poptotal <- read.table(file="Populationtotal.txt")
 
### calculate the expected number of reported cases per week and age group ##
#############################################################################
    vec_hosp <- rep(0,datsize)
    reg.incistrata <- matrix(0,nrow=10,ncol=datsize)
# distiguishing between the EFS (NB) and WFS (AB)  
    if (region == 'EFS'){
       for (i in 1:datsize){
# week specific reporting rate (i=1 is first week of 2001)
          vec_hosp[i]<- hosp+(hosp*(fact-1))*min(max(c(((i-(52*4))/26),0)),1)   
# age, and time specific number of predicted cases 
# (assuming homogeneous distrbution over both regions, EFS and WFS)
          reg.incistrata[,i]<-incistrata[,i]*popEFS[[(1+floor((i-1)/52))]]/poptotal[[(1+floor((i-1)/52))]]
    }}else{
        for (i in 1:datsize){
          vec_hosp[i]<- hosp+(hosp*(fact-1))*min(max(c(((i-(52*4))/26),0)),1)   
          reg.incistrata[,i]<-incistrata[,i]*(1-(popEFS[[(1+floor((i-1)/52))]]/poptotal[[(1+floor((i-1)/52))]]))
        }
    }
            
# calculate model-predicted reported number of cases in EFS 
    repinci <- vec_hosp*t(reg.incistrata[1:10,1:datsize])

#### Calculate loglikelihood of the data utilizing the cumulative autocorrelations
# 0 < ESS <= 1; ESS = (cumulative autocorrelation)^(-1); specified for each time series
# assuming NegBin distribution of the data with mean coming from the model  
    logLL <- sum(ESS*dnbinom(x=dat[1:10,1:datsize],mu=t(repinci),size=disp,log = TRUE))
    return(logLL)
}