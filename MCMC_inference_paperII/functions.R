### expit function
expit <- function(x){
return(exp(x)/(1+exp(x)))
}



### logit function
logit <- function(p){
return(log(p/(1-p)))
}



### dsnorm is density function of skew normal distribution
dsnorm <- function(x,mean=0,sd=1,a=0){
return(2/sd*dnorm((x-mean)/sd)*pnorm(a*(x-mean)/sd))
}



### vacccoverage computes the weekly vaccination coverage rates from 2006 until 2013
### by applying Spline-Interpolation on the transformed coveragefunction on the real axis
### Retransformation (expit) projects into the predefined intervals:
### [0.007,0.0068]; [0.027,0.118]; [0.31,0.48]; [0.51,0.67]; [0.45,0.68]; [0.45,0.68]; [0.45,0.68]; [0.45,0.68]
### for vaccination coverage at the end of year 2006, ..., 2013 , respectively
### arguments: coveragerates: vector in R^8
### output:    vector of length 25*23 (momentaneous coverage for each week in 1991-2013)
vacccoverage <- function(coveragerates=rep(0,8))
{
cover     <-rep(0,(52*23))
crates    <- c(0.0001,0.007,0.027,0.31,0.51,0.45,0.45,0.45,0.45)+c(0,0.061,0.091,0.17,0.16,0.23,0.23,0.23,0.23)*c(0,expit(coveragerates))
vacccover <-spline(52*(0:8),y=logit(crates),xout=1:416)
cover[(52*15+1):(52*23)]<-expit(vacccover$y)
return(cover)
}



### contactmatrix computes the 19x19 WAIFW matrix for disease transmission
### within the age structured population, six distinct contact patterns are possible
### arguments: c1,c2,c3 are contact frequencies and fill components of the matrix
### output: a 19x19 matrix
contactmatrix <- function(cpara=rep(1,3),pattern=1)
{                                    
contact <- matrix(cpara[1],nrow=19,ncol=19);                   
if (pattern==2){
contact[,] <- 0
contact[1:14,1:14]  <-cpara[1]
contact[15:17,15:17]<-cpara[2]
contact[18:19,18:19]<-cpara[3]}
if (pattern==3){
contact[,15:17] <-cpara[2]
contact[,18:19] <-cpara[3]}                                            
if (pattern==4){
contact[15:17,15:17] <-cpara[2]
contact[1:19,18:19]  <-cpara[3]
contact[18:19,1:17]  <-cpara[3]}
if (pattern==5){
contact[,]         <-cpara[3]
contact[1:17,1:17] <-cpara[2]
contact[1:14,1:14] <-cpara[1]}                                         
if (pattern==6){
contact[,]          <-cpara[3]
contact[15:17,15:17]<-cpara[2]
contact[1:14,1:14]  <-cpara[1]}
return(contact)
}