### function to simulate incidence data from the SIR-model using a Brownian Bridge
### process (which returns to zeor after each year) for the process of observation 
### residuals; 
datasim.BB <- function(par.alpha=5/3,par.beta=0.3,sigma=0.1,
              infolist=list('initial'=c(600500,190,399310,0),'times'=0:(8*52)))
{
## arguments:
## par.alpha, par.beta, sigma: SIR-model parameters and volatility of the BrownBridge
## infolist: initial conditions (S,I,R)_0, time horizon

# compute solution of the SIR-model
sol <- rk4(infolist[['initial']],infolist[['times']],sir,list('par.alpha'=par.alpha,'par.beta'=par.beta));
sol <- sol[,2:5];
# Compute logs of the expected weekly case numbers
log.obsmean <- log(sol[2:length(times),4]-sol[(1:length(times)-1),4]); 

### compute added noise (the Browian Bridge process)
# number of seasons
m <- ceiling((length(times)-1)/52);
# compute gaussian increments and the corresponding random walk
errors <- rnorm(n=m*52,sd=sigma);
cumerrors <- cumsum(errors);
# compute order-1-spline to connect the randomwalk points after each year
interpol <- approx(x=(((0:m)*52)+1),y=c(0,cumerrors[52*(1:m)]),xout=1:(52*m+1));
# difference between the random walk and the spline yields a brownian bridge process
# which returns to zero after each year, these are the residuals
resids <- cumerrors - interpol$y[2:length(interpol$y)];

# add the residuals to the mean process and return the log observation process
log.obs <- log.obsmean+resids[1:length(log.obsmean)];
return(log.obs)
}