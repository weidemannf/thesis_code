### initial provides the initial decomposition of the population into the model states,
### the decomposition is based on model steady state reached after running the model for
### 200 years and reweighting along the age groups according to the 1990 demography
### arguments: none
### output: a n x m matrix with n = 1 + # disease states (for cumulative number of cases in the last row),
###         and m = # age groups
initial <- function(){
  demo1990 <- 79753000*c(0.217,0.316,0.263,0.166,0.038)
  dat <- as.matrix(read.table(file="unnormedInitial.txt"))
  dat <- (dat+abs(dat))/2
  aggreggroup <- list(1:15,16,17,18,19)
  for (i in 1:length(aggreggroup)){dat[,aggreggroup[[i]]] <- dat[,aggreggroup[[i]]]*demo1990[i]/sum(dat[,aggreggroup[[i]]])}
  dat <- rbind(dat,rep(0,19))
  return(dat)
}



### births provides the weekly number of births for the years 1990-2009
### arguments: none
### output: vector of length 20*52
births <- function(){
birthdata <- c(t(as.matrix(read.table(file="births.txt",sep = ""))))
birthcounts   <- rep(0,20*52)
for (i in 1:length(birthcounts)){birthcounts[i] <- 3/13*(birthdata[ceiling(3*i/13)]/2+ birthdata[min(ceiling(3*(i+1)/13),240)]/2)}
return(birthcounts)
}



### deathrate provides the year-specific death rates (with respect to deaths per week)
### for each age group (1990-2012)
### arguments: none
### value : a n x m matrix, where n is #age groups and m is #years
deathrate <- function(){
survmale <- read.table(file='survivalcounts1990-2012male.txt')
survfema  <- read.table(file='survivalcounts1990-2012female.txt')
survprob <- log((survmale[[1]]+survfema[[1]])/(2*100000))
for (i in 2:10){ survprob <- cbind(survprob,log((survmale[[i]]+survfema[[i]])/(survmale[[i-1]]+survfema[[i-1]])))}
groupsize <- c(rep(52,5),(15*52),rep(20*52,3),10*52)
deathrate <- -t(survprob)/groupsize
modeldeathrate <- t(matrix(0,nrow=19,ncol=length(deathrate[1,])))
aggreggroup <- list(1:6,7:8,9:10,11:12,13:14,15,16,17,18,19)
for (i in 1:length(aggreggroup)){modeldeathrate[,aggreggroup[[i]]]  <- t(deathrate[i,])}
return(t(modeldeathrate))
}



### ageing provides the ageing rates for each age group,
### i.e. 1/(length of each age group in weeks)
### arguments: none
### value: a vector with length equal to the number of age groups
ageing <- function(){
x <- c(rep(6/52,6),rep(2/52,8),1/(15*52),rep(1/(20*52),3),0);
return(x)}



### migration provides weekly migration rates in absolute counts (1990-2009)
### stratified by age groups and years
### arguments: none
### output: a n x m matrix with n = #age groups and m = #years
migration <- function() {
migdata    <- as.matrix(read.table(file='migration.txt'))
groupsizes <- c(rep(1/6,6),rep(0.5,8),15,rep(20,3))
aggreggroup <- list(1:14,15,16,17,18)
weeklymigmatrix <- t(matrix(0,nrow=19,ncol=20))
for (i in 1:length(aggreggroup)){weeklymigmatrix[,aggreggroup[[i]]] <- (migdata[i,]/52)%*%t(groupsizes[aggreggroup[[i]]]/sum(groupsizes[aggreggroup[[i]]]))}
return(t(weeklymigmatrix))
}



### WFSData and EFSData provide the rotavirus weekly case data from 2001 until 2009 stratified
### by 10 age groups for the WFS and EFS, respectively
### arguments: n = the number of years from 2001 on, for which data should be provided
### output: a m x r matrix, with m = #agegroups and r = #dataweeks (52*#years except for n=13)
EFSdata <- function(n=8){
EFSdat<-numeric(0)
for (i in 1:n){
   direct<-paste("RotadataEFS200",i,".txt")
   rotadata <- read.table(file=direct)
   for (j in 1:min(52,length(rotadata))){EFSdat<-cbind(EFSdat,rotadata[[j]])}
}
return(EFSdat)
}

WFSdata <- function(n=8){
WFSdat<-numeric(0)
for (i in 1:n){
   direct<-paste("RotadataWFS200",i,".txt")
   rotadata <- read.table(file=direct)
   for (j in 1:min(52,length(rotadata))){WFSdat<-cbind(WFSdat,rotadata[[j]])}
}
return(WFSdat)
}
