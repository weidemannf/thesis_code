### projectdeath provides projected yearly age stratified deathrates for the years
### 2011 and following. deathrates are predicted by assuming the increments in the
### age specific log-deathrates to be normally distributed with mean and variance estimated
### from the data
### arguments: olddeaths = a matrix containing the past age stratified deathrates
###            timehorizon = maximum time for simulation of new death rates
### output: a matrix containing input matrix olddeaths and the simulated future death rates
###         for each age group
projectdeath <- function(olddeaths=deathrate(),timehorizon=40){
deathrates        <- matrix(0,nrow=dim(olddeaths)[1],ncol=timehorizon)
deathrates[,1:20] <- olddeaths[,1:20]
for (i in 1:dim(olddeaths)[1]){
loginc <- log(olddeaths[i,2:20])-log(olddeaths[i,1:19])
project.loginc <- rnorm(n=timehorizon-20,mean=mean(loginc[9:19]),sd=sqrt(var(loginc)))
logrates <- log(olddeaths[i,20])+cumsum(project.loginc)
deathrates[i,21:timehorizon] <- exp(logrates)
}
return(deathrates);
}



### projectmigra provides projected yearly age stratified migration rates for the years
### 2012 and following. Migration rates are simulated by assuming the age specific
### rates to follow an AR(1)-process fitted to the past migration rates up to 2011.
### arguments: oldmigras = a matrix containing the past age stratified migration rates
###            timehorizon = maximum time for simulation of new migration rates
### output: a matrix containing input matrix oldmigras and the simulated future migration rates
###         for each age group
projectmigra <- function(oldmigras=migrationNB(),timehorizon=40){
migras        <- matrix(0,nrow=dim(oldmigras)[1],ncol=timehorizon)
migras[,1:21] <- oldmigras[,1:21]
for (i in 1:dim(oldmigras)[1]){
ar.model <- ar(migras[i,10:21]-mean(migras[i,10:21]),order.max=1,AIC=FALSE)
promigs <- mean(migras[i,10:21])+arima.sim(model=list(ar=ar.model$ar), n=timehorizon-21, sd=sqrt(ar.model$var),n.start=1,start.innov=migras[i,21]-mean(migras[i,1:21]))
migras[i,22:timehorizon] <- c(promigs)
}
return(migras)
}



### weeklyfertility provides projected weekly fertility rates for the years 2013
### and following. Fertility rate is defined as the number of birth per week per person in
### the age group 20-39 years. Log fertility rates are simulated from a fitted time series
### including autogression (order 1), linear trend, and seasonality. A
### random walk was added to the simulated rates to represent future uncertainty.
### arguments: timehorizon = maximum time for simulation of new migration rates
### output: a matrix containing simulated weekly fertility rates from 1990 until end of timehorizon
# Model for ARIMA-fitting of weekly birthrates:
# y is vector of log(weekly births devided by weekly population (20-39years)) 2006-2012:
# y <- log(geburt[833:1196]/population[7,833:1196]);
# times <- 833:1196;
# coeffs <- arima(y,order = c(1, 0, 0),xreg=cbind(times,cos((times)*2*pi/52),sin((times)*2*pi/52)));
# coeffs$coef
#           ar1     intercept         times
# 0.8297456990 -7.5311359991  0.0002007732 -0.0703888503 -0.0534658917
# sigma^2 estimated as 0.0005973:  log likelihood = 833.93,  aic = -1655.87
weeklyfertility <- function(timehorizon=40){
times <- 1:(timehorizon*52)
logfert <- -7.5311 + 0.0002007*times - 0.070388*cos((times)*2*pi/52) - 0.053465*sin((times)*2*pi/52)
resids <- c(rep(0,23*52), arima.sim(model=list(ar=0.8237398),n=(timehorizon-23)*52,sd=sqrt(0.0005973),n.start=1,start.innov=-7.2837-logfert[23*52]))
deviation <- c(rep(0,23*52),cumsum(rnorm(n=(timehorizon-23)*52,mean=0,sd=0.003)))
fert.predict <- logfert+resids+deviation
return(exp(fert.predict))
}
