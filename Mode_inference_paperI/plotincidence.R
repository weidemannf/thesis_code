### Plotting the results from the predictive incidence sampling #######

## import of data function
source(file="demographicdata.R")
    datWFS <- WFSdata(n=8)
    datEFS <- EFSdata(n=8)
## fetch the generated incidence sample 
    load(file="incidencesample.Rdata")
    
## set agestrate to define which of the ten data agegroups are pooled and plotted together
   agestrata=list(1:2,3:5)
    
### compute prediction bands and sample means for expected and observed number 
### of cases stratifed by week and pooled age group
# plot time horizon
    plotyears <- 1:length(aggreginciWFS[1,,1])
    periods <- length(plotyears)/52
# define storage arrays    
    sampleWFS.quantiles <- array(0,dim=c(length(agestrata),length(plotyears),2))
    sampleEFS.quantiles <- array(0,dim=c(length(agestrata),length(plotyears),2))   
    meanWFS.quantiles   <- array(0,dim=c(length(agestrata),length(plotyears),2))
    meanEFS.quantiles   <- array(0,dim=c(length(agestrata),length(plotyears),2))  
    meanWFS  <- array(0,dim=c(length(agestrata),length(plotyears)))
    meanEFS  <- array(0,dim=c(length(agestrata),length(plotyears)))
    dataWFS  <- array(0,dim=c(length(agestrata),length(plotyears)))
    dataEFS  <- array(0,dim=c(length(agestrata),length(plotyears)))
# compute prediction bands and sample mean by the empirical quantiles of the sampled incidences    
    for (j in 1:length(agestrata)){
        for (i in 1:length(plotyears)){
            sampleWFS.quantiles[j,i,]  <- quantile(colSums(aggregsampleWFS[agestrata[[j]],i,]), probs=c(0.975,0.025),na.rm=TRUE)
            sampleEFS.quantiles[j,i,]  <- quantile(colSums(aggregsampleEFS[agestrata[[j]],i,]), probs=c(0.975,0.025),na.rm=TRUE)   
            meanWFS.quantiles[j,i,]  <- quantile(colSums(aggreginciWFS[agestrata[[j]],i,]), probs=c(0.975,0.025),na.rm=TRUE)
            meanEFS.quantiles[j,i,]  <- quantile(colSums(aggreginciEFS[agestrata[[j]],i,]), probs=c(0.975,0.025),na.rm=TRUE) 
            meanWFS[j,i]    <- mean(colSums(aggreginciWFS[agestrata[[j]],i,]),na.rm=TRUE)
            meanEFS[j,i]    <- mean(colSums(aggreginciEFS[agestrata[[j]],i,]),na.rm=TRUE) 
        }
        dataWFS[j,] <- colSums(datWFS[agestrata[[j]],])
        dataEFS[j,] <- colSums(datEFS[agestrata[[j]],])
    }
    
### Plotting the graphical output  #########
############################################
    
# define some vectors for flexible axis and graphic notation    
    axis.spots <- (0:(periods-1))*52+1;  axis.labels <- as.character((2001):(2001+periods-1));  
    l.agebound <- c(0,1,2,3,4,5,20,40,60,80)
    u.agebound <- c(0,1,2,3,4,19,39,59,79,100)

## generate the plot: set plot window; then for each age group and region:
## 1) draw prediction bands for observations and expectations;
## 2) draw prediction mean; 
## 3) draw data; 
             
    #pdf(file='incidence.pdf',width=18,height=10)
    par(mar=c(5.5, 5, 4, 0.5) + 0.1)
    layout(matrix(1:(2*length(agestrata)),nrow=length(agestrata)),width=c(1,1),height=rep(0.5,length(agestrata)))
    leg <- TRUE
    for (j in 1:length(agestrata)){
    plot(plotyears,meanWFS[j,],main=paste(l.agebound[min(agestrata[[j]])],'-',u.agebound[max(agestrata[[j]])],'years of age, WFS'),type='l',col='green',ylim=c(0,1.4*max(c(sampleWFS.quantiles[j,,1],sampleEFS.quantiles[j,,1]))),ylab='weekly reported cases',xlab='time (weeks)',xaxt='n',cex.main=2.0,cex.lab=2.0,cex.axis=1.6)
    axis(1,at=axis.spots,labels=axis.labels,las=0,cex.axis=1.8)
    polygon(c(plotyears,rev(plotyears)),c(sampleWFS.quantiles[j,,2],rev(sampleWFS.quantiles[j,,1])),col='grey85',border = NA)
    polygon(c(plotyears,rev(plotyears)),c(meanWFS.quantiles[j,,2],rev(meanWFS.quantiles[j,,1])),col='grey55',border = NA)
    lines(meanWFS[j,],col='green',type='l')
    lines(dataWFS[j,],col='red',type='l')
    if (leg){
      legend("topleft",legend=c(expression(paste('95% prediction interval for ',X[t])),expression(paste('95% prediction interval for ',h %.% Y[t](vartheta)))),
            col=c('white','white'),fill=c('grey85','grey55'),border=c('grey85','grey55'),horiz=FALSE,cex=2.2)
      legend("topright",legend=c(expression(paste('E(',X[t],')')),expression(D[t])),
            col=c('green','red'),lty=c(1,1),horiz=FALSE,cex=2.0)
      leg <- FALSE      
      }
    }
    for (j in 1:length(agestrata)){
    plot(plotyears,meanEFS[j,],main=paste(l.agebound[min(agestrata[[j]])],'-',u.agebound[max(agestrata[[j]])],'years of age, EFS'),type='l',col='green',ylim=c(0,1.4*max(c(sampleWFS.quantiles[j,,1],sampleEFS.quantiles[j,,1]))),ylab='weekly reported cases',xlab='time (weeks)',xaxt='n',cex.main=2.0,cex.lab=2.0,cex.axis=1.6)
    axis(1,at=axis.spots,labels=axis.labels,las=0,cex.axis=1.8)
    polygon(c(plotyears,rev(plotyears)),c(sampleEFS.quantiles[j,,2],rev(sampleEFS.quantiles[j,,1])),col='grey85',border = NA)
    polygon(c(plotyears,rev(plotyears)),c(meanEFS.quantiles[j,,2],rev(meanEFS.quantiles[j,,1])),col='grey55',border = NA)
    lines(meanEFS[j,],col='green',type='l')
    lines(dataEFS[j,],col='red',type='l')
    }
    #dev.off()

