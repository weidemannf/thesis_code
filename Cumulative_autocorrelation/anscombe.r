### functions for computing the Anscombe residuals of NegBin distributed data
### for more details see Hardin and Hilbe: Generalized Linear Models and Extensions
V <- function(mu,r){
return((mu^2)/r+mu)}

integrand <- function(mu,r){
return(((mu^2)/r+mu)^(-1/3))}

A <- function(y,mu=1,r=1){
return(integrate(integrand,lower=mu,upper=y,r=r)[[1]])}

resi <- function(y,mu=1,r=1){
return(A(y,mu,r)/(V(mu,r)^(1/6)))}

### main function to compute the Anscombe residuals of NegBin distributed data
### arguments: obs = vector of observations 
###            obsmean = vector of expected values
###            disp = dispersion parameter
anscombe.resi <- function(obs,obsmean,disp)
{
anscomberesi <- rep(0,length(obs));
for (i in 1:length(obs)){
anscomberesi[i] <- resi(obs[i],obsmean[i],disp)
}
return(anscomberesi-mean(anscomberesi))
}