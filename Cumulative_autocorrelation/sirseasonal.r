### ODE system governing the simple seasonal SIR-model
### the system includes a seasonal term: (1+0.2*sin(t*2*pi/52))
### and births and deaths at 1/4000 per week
sir <- function(t,x,para=list('par.alpha'=0.4,'par.beta'=0.4))
{
## arguments: 
## t is time
## x is state of the system (4-dim)
## para is list of alpha (R0) and beta (recovery rate)

# define left hand side of the ODE system
dx <- numeric(0);
# extract information from the para list
par.alpha <- para[['par.alpha']];
par.beta  <- para[['par.beta']];
# calculate the population size N
N<-sum(x[1:3]);
# compute left hand side
dx[1] <- -par.alpha*par.beta*(1+0.2*sin(t*2*pi/52))*x[1]*x[2]/N + (N-x[1])/4000;
dx[2] <- par.alpha*par.beta*(1+0.2*sin(t*2*pi/52))*x[1]*x[2]/N - par.beta*x[2] - x[2]/4000;
dx[3] <- par.beta*x[2] - x[3]/4000;                            
# last component governs the infinitesimal increase in the number of occured infections
dx[4] <- par.alpha*par.beta*(1+0.2*sin(t*2*pi/52))*x[1]*x[2]/N;
# return results
return(list(dx));
}