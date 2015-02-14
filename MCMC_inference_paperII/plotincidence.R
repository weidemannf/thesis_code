#### Generate plots of the incidence predictons from the posterior distribution  ##
###  and various parameter distributions ##########################################
###################################################################################

############ import necessary functions and packages and data ##################
## functions which provide incidence data, demographic data or other functions
source(file="demographicdataEFS.R")
source(file="functions.R")

## fetch results from the posterior sampling used for incidence prediction ####
    load(file='adaptivemcmcresults.Rdata')
## fetch incidence data
    datEFS <- dataEFS(n=13)
## fetch incidence prediction sample for a specific coverage
    coveragelevel <- 0.9
load(file=paste('populationsample',coveragelevel,'.Rdata'))      
load(file=paste('incidencesampling',coveragelevel,'.Rdata'))     
load(file=paste('expincisampling',coveragelevel,'.Rdata'))

### compute pointwise prediction bands for the expected and reported incidence 
### based on the available predictive sample for aggregated age groups
# define, which age groups should be pooled together for the prupose of illustration
agestrata <- list(1:5,6:8,9:10)
# set up storage matrices
expquantile  <- numeric(0); reportquantile  <- numeric(0);
popmatrix  <- numeric(0);   datmatrix  <- numeric(0);
# loop over the pooled age groups
for (i in agestrata){
if (length(i)>1){
expmatrix    <- colSums(storeexpcases[i,,])
reportmatrix <- colSums(storereportcases[i,,])
popmatrix   <-  colSums(storepops[i,,])  
# observed incidence until 2013 not subject to stochastic variation of the population,
# i.e. equal for all population samples
datmatrix   <-  rbind(datmatrix,100000*colSums(datEFS[i,])/popmatrix[1:dim(datEFS)[2],1])
}else{
expmatrix    <- storeexpcases[i,,]
reportmatrix <- storereportcases[i,,]
popmatrix   <-  storepops[i,,]  
datmatrix   <-  rbind(datmatrix,100000*datEFS[i,]/popmatrix[1:dim(datEFS)[2],1])
}
## loop over the available samples
for (j in 1:dim(storeexpcases)[2]){
expquantile     <- rbind(expquantile,100000*quantile(expmatrix[j,]/popmatrix[j,],p=c(0.025,0.5,0.975)))
reportquantile  <- rbind(reportquantile,100000*quantile(reportmatrix[j,]/popmatrix[j,],p=c(0.025,0.975)))
}}
dim(expquantile) <- c(dim(storeexpcases)[2],length(agestrata),3)
dim(reportquantile) <- c(dim(storereportcases)[2],length(agestrata),2)
datmatrix  <- t(datmatrix)



### Plotting of the incidence curve according to the data and the prediction bands
### from the incidences sampled by the model 
# define plotting years, i.e. from plotstart to plotend
    plotstart <- 2004
    plotend   <- 2021
    plotyears <- (1+52*(plotstart-2001)):(52*(plotend-2001))
    if (plotstart < 2014){
    datyears  <-  (1+52*(plotstart-2001)):min(650,52*(plotend-2000))
    }

# define some vectors for flexible axis and graphic notation 
    axis.spots <- (0:(periods-1))*52+1;  axis.labels <- as.character((2001):(2001+periods-1));  
    l.agebound <- c(0,1,2,3,4,5,20,40,60,80)
    u.agebound <- c(0,1,2,3,4,19,39,59,79,100)
       
    #pdf(file='vaccincidencelong.pdf',width=18,height=10.5)
    layout(matrix(1:length(agestrata),nrow=length(agestrata)),width=c(1),height=c(2,rep(1,length(agestrata)-1)))
    par(mar=c(5, 5, 4, 0.5) + 0.1)   
    for (j in 1:length(agestrata)){
    if (j==1){plot(plotyears,expquantile[plotyears,j,2],main=paste(l.agebound[min(agestrata[[j]])],'-',u.agebound[max(agestrata[[j]])],'years of age'),type='l',col='black',
              ylim=c(0,1.16*max(c(datmatrix[datyears,j],reportquantile[plotyears,j,2]))),ylab='weekly incidence (per 100,000)',xlab='time (years)',xaxt='n',cex.main=1.9,cex.lab=1.9,cex.axis=1.2)
    }else{
    plot(plotyears,expquantile[plotyears,j,2],main=paste(l.agebound[min(agestrata[[j]])],'-',u.agebound[max(agestrata[[j]])],'years of age'),type='l',col='black',
         ylim=c(0,1.4*max(c(reportquantile[plotyears,2:length(agestrata),2]))),ylab='weekly incidence',xlab='time (years)',xaxt='n',cex.main=1.9,cex.lab=1.9,cex.axis=1.2)}
    axis(1,at=axis.spots,labels=axis.labels,las=0,cex.axis=1.5)
    polygon(c(plotyears,rev(plotyears)),c(reportquantile[plotyears,j,1],rev(reportquantile[plotyears,j,2])),col='grey85',border = NA)
    polygon(c(plotyears,rev(plotyears)),c(expquantile[plotyears,j,1],rev(expquantile[plotyears,j,3])),col='grey55',border = NA)
    lines(plotyears,expquantile[plotyears,j,2],col='green',type='l',lwd=2)
    if (plotstart < 2014){lines(datyears,datmatrix[datyears,j],col='red',lty=1,lwd=2)}
    if (j==1){
      leg1<-legend("topright",legend=c('95% prediction interval for reported incidence','95% prediction interval for expected incidence','mean incidence prediction','reported incidence data'),
                  fill=c('grey85','grey55','white','white'),border=c('grey85','grey55','white','white'),horiz=FALSE,cex=2.3)
      legend(x=leg1$rect$left-7,y=leg1$rect$top,legend=c('','','',''),col=c(NA,NA,'green','red'),lty=c(1,1),lwd=c(0,0,2,2),box.col = NA, , bg=NA, horiz=FALSE,cex=2.3)
      }
    }
    #dev.off()
    
### plotting of the momenaneous vaccination coverage rates according to ########
#### posterior distribution ####################################################
# set size for coverage parameter sample from the posterior
  K <- 1000
# define index vector for drawing from the posterior sample    
  index <- round(dim(paramatrix1)[2]*(1:min(K,dim(paramatrix1)[2]))/min(K,dim(paramatrix1)[2]))
# calculate samples for the pointwise (wrt. time) coverage rate 
  coveragesample <- numeric(0)
  for (i in index){
      coverage <- vacccoverage(coveragerates=paramatrix3[,i])  
      coveragesample <- rbind(coveragesample,coverage[521:length(coverage)])
  }
  covquantile <- numeric(0)
  for (j in 1:dim(coveragesample)[2]){
      covquantile    <- rbind(covquantile,quantile(coveragesample[,j],p=c(0.025,0.5,0.975)))
  }
# prior data
  priorups <- c(0.07,0.115,0.485,0.67,0.675,0.675,0.675,0.675)
  priorbottoms <- c(0.01,0.025,0.315,0.51,0.445,0.445,0.445,0.445)
  priormeans <- logit((expit(priorups)+expit(priorbottoms))/2)

# plotting of the prediction bands and priors 
  plotyears <- (1+52*3):(52*13)    
  #pdf(file='vacccoverage.pdf',width=12,height=5)
  plot(plotyears,covquantile[plotyears,2],main='Vaccination coverage',type='l',col='blue',ylim=c(0,1),ylab='coverage rate',xlab='time (years)',yaxt='n',xaxt='n',cex.main=1.4,cex.lab=1.6,cex.axis=1.3)
  axis(1,at=(3:12)*52+1,labels=as.character(2004:2013),las=0,cex=1.2)
  axis(2,at=c(0,0.2,0.4,0.6,0.8,1),labels=c('0%','20%','40%','60%','80%','100%'),las=0,cex=1.2)
  polygon(c(plotyears,rev(plotyears)),c(covquantile[plotyears,1],rev(covquantile[plotyears,3])),col='lightblue',border = NA)
  lines(plotyears,covquantile[plotyears,2],col='blue',type='l',lwd=2)
  for(i in 1:8){
    lines(c(1+52*(5+i),1+52*(5+i)),c(priorups[i],priorbottoms[i]),col='grey59',lwd=2)
    lines(c(-3+52*(5+i),5+52*(5+i)),c(priorups[i],priorups[i]),col='grey59',lwd=2)
    lines(c(-3+52*(5+i),5+52*(5+i)),c(priorbottoms[i],priorbottoms[i]),col='grey59',lwd=2)
    points(1+52*(5+i),priormeans[i],pch=20,col='grey59',lwd=1)
  }
  legend("topleft",legend=c(expression(paste('95% credibility interval for ',phi[t])),'mean and range of prior distribution'),
            col=c('white','white'),fill=c(color[1],'grey59'),border=c(color[1],'white'),horiz=FALSE,cex=1.4)
  #dev.off()
    
    
### Plotting of histograms of the vaccine effectiveness parameters #############
################################################################################
# set size for effectiveness parameter sample from the posterior
    K <- 10000
# define index vector for drawing from the posterior sample    
    index <- round(dim(paramatrix1)[2]*(1:min(K,dim(paramatrix1)[2]))/min(K,dim(paramatrix1)[2]))    
    etaI <- expit(paramatrix2[1,index]); Iquant <- quantile(etaI,p=c(0.025,0.5,0.975))
    etaS <- expit(paramatrix2[2,index]); Squant <- quantile(etaS,p=c(0.025,0.5,0.975))
    etaP <- etaI*etaS;                   Pquant <- quantile(etaP,p=c(0.025,0.5,0.975))
    etaW <- exp(paramatrix2[3,index]);        Wquant <- quantile(1/(52*etaW),p=c(0.025,0.5,0.975))

    #pdf(file='vacceffectiveness.pdf',width=13,height=13)
    layout(t(matrix(1:4,nrow=2)),width=c(1,1),height=c(1,1))
    par(mar=c(5, 5, 4, 0.5) + 0.1)
    plotdata <- hist(etaI,main=expression(atop('(a)','Risk ratio for acquiring infection')),xlim=c(0,1),breaks=40,xlab=expression(eta[I]),freq=FALSE,ylab='density',
                cex.main=1.8,cex.sub=1.9,cex.lab=2.2,cex.axis=1.8)
    lines(rep(Iquant[1],2),c(-max(plotdata$density)/20,max(plotdata$density)/20),col='lightblue',lwd=4)
    lines(rep(Iquant[2],2),c(-max(plotdata$density)/20,max(plotdata$density)/20),col='blue',lwd=4)
    lines(rep(Iquant[3],2),c(-max(plotdata$density)/20,max(plotdata$density)/20),col='lightblue',lwd=4)
    plotdata <- hist(etaS,main=expression(atop('(b)','Risk ratio for developing symptoms')),xlim=c(0,1),breaks=40,xlab=expression(eta[S]),freq=FALSE,ylab='density',
                cex.main=1.8,cex.sub=1.9,cex.lab=2.2,cex.axis=1.8)
    lines(rep(Squant[1],2),c(-max(plotdata$density)/20,max(plotdata$density)/20),col='lightblue',lwd=4)
    lines(rep(Squant[2],2),c(-max(plotdata$density)/20,max(plotdata$density)/20),col='blue',lwd=4)
    lines(rep(Squant[3],2),c(-max(plotdata$density)/20,max(plotdata$density)/20),col='lightblue',lwd=4)
    plotdata <- hist(etaP,main=expression(atop('(c)','Risk ratio for acquiring symptomatic infection')),xlim=c(0,1),breaks=20,xlab=expression(eta[I]*eta[S]),freq=FALSE,ylab='density',
                cex.main=1.8,cex.sub=1.9,cex.lab=2.2,cex.axis=1.8)
    lines(rep(Pquant[1],2),c(-max(plotdata$density)/20,max(plotdata$density)/20),col='lightblue',lwd=4)
    lines(rep(Pquant[2],2),c(-max(plotdata$density)/20,max(plotdata$density)/20),col='blue',lwd=4)
    lines(rep(Pquant[3],2),c(-max(plotdata$density)/20,max(plotdata$density)/20),col='lightblue',lwd=4)
    plotdata <- hist((1/(52*etaW))[(1/(52*etaW))<10],main=expression(atop('(d)','Mean duration of immunity loss')),breaks=seq(0,10,by=0.25),xlab=expression(1/(52*eta[W])),
                cex.main=1.8,cex.sub=1.9,cex.lab=2.2,cex.axis=1.8,freq=FALSE,ylab='density')
    lines(rep(Wquant[1],2),c(-max(plotdata$density)/20,max(plotdata$density)/20),col='lightblue',lwd=4)
    lines(rep(Wquant[2],2),c(-max(plotdata$density)/20,max(plotdata$density)/20),col='blue',lwd=4)
    lines(rep(Wquant[3],2),c(-max(plotdata$density)/20,max(plotdata$density)/20),col='lightblue',lwd=4)
    #dev.off()
  
  
  
      
#### Plot the shift in the seasonality of incidence ############################
### this requires the computation of the incidence predictions at the beginning
### of this script #############################################################
# define, which seasons to plot
    refyear <- 2006;  refseason     <- (-11:40)+(52*(refyear-2001))
    datyear <- 2013;  lastdatseason <- (-11:40)+(52*(datyear-2001))
    futyear <- 2026;  futureseason  <- (-11:40)+(52*(futyear-2001))
# plot graphics    
    axen <- -11:40 
    #pdf(file='seasincidence.pdf',width=9*length(agestrata),height=9) 
    par(mar=c(6, 5, 4, 0.5) + 0.1)   
    plot(axen,reportquantile[refseason,1,1],main=paste(l.agebound[min(agestrata[[1]])],'-',u.agebound[max(agestrata[[1]])],'years of age'),type='l',col='white',
         ylim=c(0,1.4*max(c(reportquantile[refseason,1,2]))),ylab='weekly incidence (per 100,000)',xlab='Week of Year from January 1',xaxt='n',cex.main=1.8,cex.lab=1.7,cex.axis=1.5)
    axis(1,at=seq(-10,40,by=5),labels=as.character(seq(-10,40,by=5)),las=0,cex.axis=1.4)
    lines(axen,expquantile[refseason,1,2],col=rgb(0, 1, 0,1),type='l',lwd=2)
    lines(axen,expquantile[lastdatseason,1,2],col=rgb(0, 0, 1,1),type='l',lwd=2)
    lines(axen,expquantile[futureseason,1,2],col=rgb(1, 0, 0,1),type='l',lwd=2)
    polygon(c(axen,rev(axen)),c(reportquantile[refseason,1,1],rev(reportquantile[refseason,1,2])),col=rgb(0, 1, 0,0.25),border = NA)   
    polygon(c(axen,rev(axen)),c(reportquantile[lastdatseason,1,1],rev(reportquantile[lastdatseason,1,2])),col=rgb(0, 0, 1,0.25),border = NA)  
    polygon(c(axen,rev(axen)),c(reportquantile[futureseason,1,1],rev(reportquantile[futureseason,1,2])),col=rgb(1, 0, 0,0.25),border = NA)   
    legend("topright",legend=c(paste(refyear,'model prediction (Exp. and 95% PI)'),paste(datyear,'model prediction (Exp. and 95% PI)'),paste(futyear,'model prediction (Exp. and 95% PI)')),
            col=c('green','blue','red'),fill=c(rgb(0, 1, 0,0.25),rgb(0, 0, 1,0.25),rgb(1, 0, 0,0.25)),border=c(rgb(0, 1, 0,0.25),col=rgb(0, 0, 1,0.25),rgb(1, 0, 0,0.25)),horiz=FALSE,lwd=2,cex=1.3)  
    #dev.off()