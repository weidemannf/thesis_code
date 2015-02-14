### initialNB provides the initial decomposition of the EFS population into the model states
### the decomposition is based on model steady state reached after running the model for
### 200 years and reweighting along the age groups according to the 1990 demography
### arguments: none
### output: a n x m matrix with n = 1 + # disease states (for cumulative number of cases in the last row),
###         and m = # age groups
initialEFS <- function(){
  demo1990 <- as.matrix(read.table(file="popEFS.txt"))[1,]
  dat <- as.matrix(read.table(file="unnormedInitialNB.txt"))
  dat <- (dat+abs(dat))/2
  aggreggroup <- list(1:6,7:8,9:10,11:12,13:14,15,16,17,18,19)
  for (i in 1:length(aggreggroup)){dat[,aggreggroup[[i]]] <- dat[,aggreggroup[[i]]]*demo1990[i]/sum(dat[,aggreggroup[[i]]])}
  dat <- rbind(dat,matrix(0,nrow=5,ncol=19))
  return(dat)
}



### birthsEFS provides the weekly number of births in the EFS for the years 1990-2012
### arguments: none
### output: vector of length 23*52
birthsEFS <- function(){
birthdata <- read.table(file="birthsEFS.txt")[[1]]
birthcounts   <- rep(0,23*52)
for (i in 1:length(birthcounts)){
birthcounts[i] <- 3/13*(birthdata[ceiling(3*i/13)]/2+ birthdata[min(ceiling(3*(i+1)/13),264)]/2)}
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
for (i in 2:10){
  survprob <- cbind(survprob,log((survmale[[i]]+survfema[[i]])/(survmale[[i-1]]+survfema[[i-1]])))
}
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
x <- c(rep(6/52,6),rep(2/52,8),1/(15*52),rep(1/(20*52),3),0)
return(x)}



### migrationEFS provides weekly migration rates (in absolute counts) for the EFS (1990-2013)
### stratified by age groups and years, recalculated from changing population sizes
### arguments: none
### output: a n x m matrix with n = #age groups and m = #years
migrationEFS <- function() {
pop <- read.table(file='popEFS.txt')
deathrates <- (52*deathrate())[c(1,7,9,11,13,15:19),]
agerates   <- c(rep(1,5),1/(15),rep(1/(20),3),0)
birth      <- birthsEFS(); dim(birth) <- c(52,23); birth <- colSums(birth)
migmatrix <- matrix(0,nrow=10,ncol=20)
for (j in 1:dim(migmatrix)[2]){
    migmatrix[1,j] <- (deathrates[1,j+1]+agerates[1])*(pop[[1]][j+1]-pop[[1]][j]*exp(-deathrates[1,j+1]-agerates[1]))/(1-exp(-deathrates[1,j+1]-agerates[1]))-birth[j]
    for (i in 2:dim(migmatrix)[1]){
        migmatrix[i,j] <- (deathrates[i,j+1]+agerates[i])*(pop[[i]][j+1]-pop[[i]][j]*exp(-deathrates[i,j+1]-agerates[i]))/(1-exp(-deathrates[i,j+1]-agerates[i])) -agerates[i-1]*pop[[i-1]][j]
}}
weeklymigmatrix <- t(matrix(0,nrow=19,ncol=20))
aggreggroup <- list(1:6,7:8,9:10,11:12,13:14,15,16,17,18,19)
for (i in 1:length(aggreggroup)){weeklymigmatrix[,aggreggroup[[i]]] <- migmatrix[i,]/(52*length(aggreggroup[[i]]))}
weeklymigmatrix <- cbind(t(weeklymigmatrix),matrix(rep(colMeans(weeklymigmatrix[17:19,]),4),ncol=4))
return(weeklymigmatrix)
}



### dataEFS provides the rotavirus weekly case data from 2001 until 2013 stratified
### by 10 age groups
### arguments: n = the number of years from 2001 on, for which data should be provided
### output: a m x r matrix, with m = #agegroups and r = #dataweeks (52*#years except for n=13)
dataEFS <- function(n=8){
datEFS<-numeric(0);
for (i in 1:n){
   direct<-paste("RotadataEFS200",i,".txt")
   rotadata <- read.table(file=direct);
   for (j in 1:min(52,length(rotadata))){
      datEFS<-cbind(datEFS,rotadata[[j]]);
   }
}
return(datEFS)
}
