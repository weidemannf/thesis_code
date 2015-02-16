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



### optimgrad optimizes a function on a multidimnensional space along a predefined
### direction given through a vector. 
### arguments: fun is a function, initial is the point from which optimization begins
###            grad is the vector defining the (two-sided) direction in which to optimize
###            info is further information passed to fun
### output:    a new parameter vector containing a new optimum
optimgrad <- function(fun,initial,grad,info){
newfun <- function(x){
value <- fun(initial-(x*grad),info)
return(value)}
newpara <- optimise(newfun,interval=c(-1,1))
return(newpara);
}



### contactmatrix computes the 19x19 WAIFW matrix for disease transmission
### within the age structured population, six distinct contact patterns are possible
### arguments: c1,c2,c3 are contact frequencies and fill components of the matrix
### output: a 19x19 matrix
contactmatrix <- function(cpara=rep(1,3),pattern=1)
{                                    
contact <- matrix(cpara[1],nrow=19,ncol=19);                   
if (pattern==2){
contact[,] <- 0;
contact[1:14,1:14]  <-cpara[1];
contact[15:17,15:17]<-cpara[2];
contact[18:19,18:19]<-cpara[3];}
if (pattern==3){
contact[,1:14]  <-cpara[1]
contact[,15:17] <-cpara[2];
contact[,18:19] <-cpara[3];}                                            
if (pattern==4){
contact[15:17,15:17] <-cpara[2];
contact[1:19,18:19]  <-cpara[3];
contact[18:19,1:17]  <-cpara[3];}
if (pattern==5){
contact[,]         <-cpara[3];
contact[1:17,1:17] <-cpara[2];
contact[1:14,1:14] <-cpara[1];}                                         
if (pattern==6){
contact[,]          <-cpara[3];
contact[15:17,15:17]<-cpara[2];
contact[1:14,1:14]  <-cpara[1];}
return(contact);
}