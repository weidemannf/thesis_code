### to a given univariate time series 'x' this function copmutes the cumulative 
### autocorrelation by fitting the time series to an ARMA process and calculation
### all autocorrelations (for all lags) based on the fitted coefficients 

effectivesamplesize <- function(x,max.timelag=length(x))
{
# arguments
# x: a real vector (an univariate time series)
# max.timelag: the maximum time lag T for later summation of the lag-t autocorr.

### first find the best fitting ARMA model to match the time series (using AIC)
### the search algorthm searches within models with max(p,q)<=n.max, where p and q 
### are the order of the AR and MA component, respectively
### the algorithm then increases n.max until no further model improvement (smaller AIC) 
### is achieved
# inital optimal p and q  equal zero
p.order <- 0; q.order <- 0; 
# compute aic of the fitted ARMA-model
model <-arima(x,order=c(p.order,0,q.order),include.mean=FALSE)
aic <- model$aic;

# the search begin
notstop <- TRUE;
n.max <- 0;
# loop over n.max
while (notstop){
# assume that this is the last n.max to be incestigated
  notstop <- FALSE;
# increase nmax  
  n.max <- n.max+1;
# first keep q=n.max and search along p in 0:n.max, computing the corresponding
# model fits, i.e. AICs  
  for(p in 0:n.max){
      testmodel <-try(arima(x,order=c(p,0,n.max),include.mean=FALSE));
      if(!(attr(testmodel,'class')=="try-error")){
# if one model yields smaller AIC, define this as the new best model and tell 
# the algorithm to also investigate the next n.max (notstop <- TRUE)
        if (testmodel$aic < aic){
         aic <- testmodel$aic;
          p.order <- p;
          q.order <- n.max;
          notstop <- TRUE;
      }}
  }
# do the same for p=n.max and q in 0:n.max  
  for(q in 0:n.max){
      testmodel <-try(arima(x,order=c(n.max,0,q),include.mean=FALSE));
      if(!(attr(testmodel,'class')=="try-error")){
      if (testmodel$aic < aic){
        aic <- testmodel$aic;
        p.order <- n.max;
        q.order <- q;
        notstop <- TRUE;
      }}
  }
}
# eventually p.order and q.order store the best fitting model orders
# compute the best fitting model, i.e. its coefficients
model <- arima(x,order=c(p.order,0,q.order),include.mean=FALSE);

### Based on the ARMA-coefficients we compute the cumulative autocorrelation, i.e.
### the sum of all autocorrelations (all lags) of the model
### the single autocorrelation are computed according to P.Brockwell/R.Davis:
### Introduction to Time Series and Forecasting, Second Edition page 90

m <- max(p.order,(q.order+1));

# get phi (AR-coeffs), theta (MA-coeffs) und sigma (standard deviation of gaussian incr.)
phi   <- model[[1]][1:p.order]; theta <- numeric(0);
if (q.order >0){theta <- model[[1]][(p.order+1):(p.order+q.order)]}
# include theta_0=1 and the missing components theta_k=0 for q.order < k <=m
theta <- c(1,theta); theta <- c(theta,rep(0,m-length(theta)));
sigma2 <- model[[2]];

### compute the zeroth to m-th autocovariance  #####################
####################################################################
### vector of autocovariances gamma will be computed by A*gamma=b
### with suitably defined A and b

### compute matrix A, where A = A1 + A2  , both (m+1)^2 matrices with i,j in 0:m
### A1[i,j] = -phi_(i+j), with -phi_0=1 and phi_k=0 for k>p.order
### A2[i,j] = -phi_(i-j), with -phi_0=1 and phi_k=0 for k>p.order or k<0
# add zeroth component to phi for further use: phi_0=-1
phi <- c(-1,phi)
A1 <- matrix(0,nrow=m+1,ncol=m+1);
A2 <- matrix(0,nrow=m+1,ncol=m+1);
for (j in 1:(m+1)){
  for (i in 1:(m+1)){
    if(!is.na(phi[i+j-1])){A1[i,j] <- -phi[i+j-1]}
    if(j>1 & i>=j){if(!is.na(phi[1+i-j])){A2[i,j] <- -phi[1+i-j]}}
  }
}
A <- A1+A2;

### compute b, where b is vector of length m+1 (i in 0:m)
### with b[i]= sigma^2* sum_{j=(0:Inf)} theta_i+j * psi_j,
### where psi_j= theta_j + sum_{k=(1:min(p,j)} theta_k * psi_(j-k), and psi_0=1
### see page 85 

psi <- (1:m)*0;
psi[1] <- 1;
if (m>1){
for (i in 2:m){psi[i] <- theta[i] +phi[2:min(length(phi),i)]%*%psi[(i-1):(1+i-min(length(phi),i))];}
}
b <- (1:m)*0;
for (i in 1:m){b[i] <- sigma2*(theta[i:length(theta)])%*%psi[1:(length(theta)-i+1)];}

### now solve for the autocovariances gamma (here acvs)
acvs <- solve(A,c(b,0));          

### now compute the autocovariances for lag > m #############################
### these take the form gamma(k)= alpha_1*xi_1^(-k)+  ... + alpha_p*xi_p^(-k)
### where the xi_j are the roots (Phi(z)==0) of Phi(z)=1- phi_1*z- ... - phi_p*z^p
### the alphas are determined as the above equation must also hold for gamma_k with k<=m

if (p.order > 0){
# compute xi
require('polynom');
xi <-  solve(polynomial(phi))
# define matrix XI with XI[k,i] = xi_i^(m-p.order+k) 
XI <- matrix(0,nrow = p.order, ncol = p.order);
for (i in 1:p.order){XI[i,]<-xi^-(m-p.order-1+i)}
# solve for alpha: XI %*% alpha = gamma((m-p.order):(m-1))  
alpha <- solve(XI,acvs[(m-p.order+1):m]);

### compute the autocorrelations 
# first up to (m-1)th autocorrelation
acr <- acvs[1:m]/acvs[1]
# then from m up to max.timelag
for (i in m:max.timelag){
acr <- c(acr,alpha%*%(xi^-(i))/acvs[1]);
}
### compute the inverse cumulative autocorrelation
CA <- 1/(1+2*sum(acr[2:max.timelag]));  
### return the best fitting model and the inverse cumulative autocorrelation (ESS)
return(list('p.order'=p.order,'q.order'=q.order,'ESS'=Re(CA)));
}

#### for the case p.order=0, the autocorrelations acr_k vanish for k>q.order 
if (p.order == 0){
acr <- c(acvs[1:m]/acvs[1],rep(0,max.timelag-(m-1)))
CA <- 1/(1+2*sum(acr[2:max.timelag])); 
return(list('p.order'=p.order,'q.order'=q.order,'ESS'=Re(CA)));
}
}