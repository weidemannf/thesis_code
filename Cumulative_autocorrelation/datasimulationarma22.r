datasim.arma22 <- function(par.alpha=0.5,par.beta=0.3,par.corr=c(0.0,0),par.theta=c(1,0,0),sigma=0.1,infolist=list('initial'=c(990,10,0),'times'=0:52),savedata=FALSE)
{
x <- infolist[['initial']];
times <- infolist[['times']];


sol <- rk4(x,times,sir,list('par.alpha'=par.alpha,'par.beta'=par.beta));
sol <- sol[,2:4];
# Compute observations and noise

log.obsmean <- log(par.alpha*sol[1:(length(times)-1),1]*sol[1:(length(times)-1),2]/sum(x));

lastresi <- c(0,0);
lasterror <- c(0,0,0);

log.obs <- numeric(0);

for (i in 1:(length(times)-1)){
lasterror[2:3] <- lasterror[1:2];
lasterror[1] <- rnorm(n=1,sd=sigma);
log.obs[i] <- log.obsmean[i]+par.corr%*%lastresi+par.theta%*%lasterror;
lastresi[2] <- lastresi[1];
lastresi[1] <-  log.obs[i]-log.obsmean[i];
}

if (savedata){
save(log.obs,file=paste('simdataarma22',as.character(par.alpha),as.character(par.beta),'.Rdata'));
plot(log.obs-log.obsmean,type='l');
}

return(log.obs)
}