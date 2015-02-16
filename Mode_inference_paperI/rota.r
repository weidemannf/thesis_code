#### ODE system for compartmentmodel with 14 classes und 19 agegroups  ###
#### this function computes the left hand side of the ODE system for given parameters

# t=time in weeks
# N= current decomposition of population (matrix) + vector of symptomatic cases so far (last row)
##### parameters in the list function arguments:
# alpha  = sequential susceptibility decrease after one and two infections
# theta  = chance of developing symptoms within first, second and third infection
# mu     = rate of revovery from symptomatic infection
# omega  = rate of revovery from asymptomatic infection
# omega0 = rate of waning maternal protection
# p      = relative infectiousness of symptomatically infected
# beta   = rate of waning natural immunity
# sais   = transmission seasonality parameters
# contact= next generation matrix
# births = vector of weekly number of births
# death  = age and time(year) dependent death rates
# mig    = age and time(year) dependent migration rates
# ageing = ageing rates

rota <- function(t,N,params=list('alpha'=c(0.6,0.667),'theta'=c(0.5,0.25,0),
      'mu'=7/3,'omega'=7/9,'omega0'=1/8,'beta'=1/50,'p'=0.5,'sais'=c(0.2,0.2,0,0),
      'contact'=contactmatrix(cpara=rep(1,3),pattern=1),'births'=births(),
      'death'=deathrate(),'mig'=migration(),'ageing'=ageing()))        
{   
# set matrix dimensions and overall population size            
    dim(N) <- c(15,length(N)/15); Pop <- sum(N[1:(dim(N)[1]-1),1:dim(N)[2]]);  
# compute modelyear
    year <- 1+floor((t-1)/52)                     
    
# extract parameters from the parameters list
    delta <- params[['ageing']]
    deathrate <- params[['death']]
    migra <- params[['mig']]
    alpha <- params[['alpha']]
    theta <- params[['theta']]
    mu <- params[['mu']]
    omega <- params[['omega']]
    omega0 <- params[['omega0']]
    beta <- params[['beta']]  
    p <- params[['p']]        
    sais <- params[['sais']]; dim(sais)<-c(2,2);
    contact <- params[['contact']]
    birthcounts <- params[['births']]    
    
# define vector of new births per age group; youngest age group has index dim(N)[2]; 
# births obviously only apply to youngest age group;                         
    birth <- c(rep(0,dim(N)[2]-1),birthcounts[round(t)])   
    
# compute the age specific force of infection using seasonality, 
# number of infectives per age group, and next generation matrix
    FOI <- exp(sais[1,1]*cos(2*pi*(t/52-sais[2,1]))+sais[1,2]*cos(2*pi*(t/26-sais[2,2])))/Pop
    Infect <- p*N[9,]+N[10,]+(p*N[11,]+N[12,])+(p*N[13,]+N[14,])
    FOI <- FOI*contact%*%Infect
        
# define left hand side of ODE system (sol)
# same size as population matrix N including number of new symptomatic infections
    sol    <- (1:length(N))*0; dim(sol) <- dim(N)
# define age-specific proportions of each model class 
    props <- t(t(N[1:(dim(N)[1]-1),1:dim(N)[2]])/colSums(N[1:(dim(N)[1]-1),1:dim(N)[2]]))
    
# M                
    sol[1,dim(N)[2]:1]  <- delta[c((dim(N)[2]-1):1,dim(N)[2])]*N[1,c((dim(N)[2]-1):1,dim(N)[2])] + props[1,dim(N)[2]:1]*migra[dim(N)[2]:1,year] + 
    birth - (delta[dim(N)[2]:1]+deathrate[dim(N)[2]:1,year]+omega0)*N[1,dim(N)[2]:1];                                                                                     
# S1
    sol[2,dim(N)[2]:1]  <- delta[c((dim(N)[2]-1):1,dim(N)[2])]*N[2,c((dim(N)[2]-1):1,dim(N)[2])] + props[2,dim(N)[2]:1]*migra[dim(N)[2]:1,year] + 
    omega0*N[1,dim(N)[2]:1] + beta*N[5,dim(N)[2]:1] - (FOI[dim(N)[2]:1]+delta[dim(N)[2]:1]+deathrate[dim(N)[2]:1,year])*N[2,dim(N)[2]:1];                                                      
# S2    
    sol[3,dim(N)[2]:1]  <- delta[c((dim(N)[2]-1):1,dim(N)[2])]*N[3,c((dim(N)[2]-1):1,dim(N)[2])] + props[3,dim(N)[2]:1]*migra[dim(N)[2]:1,year] + 
    omega*N[10,dim(N)[2]:1] + beta*N[8,dim(N)[2]:1] - (alpha[1]*FOI[dim(N)[2]:1]+beta+delta[dim(N)[2]:1]+deathrate[dim(N)[2]:1,year])*N[3,dim(N)[2]:1];                                        
# S2a    
    sol[4,dim(N)[2]:1]  <- delta[c((dim(N)[2]-1):1,dim(N)[2])]*N[4,c((dim(N)[2]-1):1,dim(N)[2])] + props[4,dim(N)[2]:1]*migra[dim(N)[2]:1,year] +  
    beta*N[3,dim(N)[2]:1] - (alpha[1]*FOI[dim(N)[2]:1]+beta+delta[dim(N)[2]:1]+deathrate[dim(N)[2]:1,year])*N[4,dim(N)[2]:1];                                                           
# S2b    
    sol[5,dim(N)[2]:1]  <- delta[c((dim(N)[2]-1):1,dim(N)[2])]*N[5,c((dim(N)[2]-1):1,dim(N)[2])] + props[5,dim(N)[2]:1]*migra[dim(N)[2]:1,year] + 
    beta*N[4,dim(N)[2]:1] - (alpha[1]*FOI[dim(N)[2]:1]+beta+delta[dim(N)[2]:1]+deathrate[dim(N)[2]:1,year])*N[5,dim(N)[2]:1];                                                           
# S3
    sol[6,dim(N)[2]:1]  <- delta[c((dim(N)[2]-1):1,dim(N)[2])]*N[6,c((dim(N)[2]-1):1,dim(N)[2])] + props[6,dim(N)[2]:1]*migra[dim(N)[2]:1,year] + 
    omega*(N[12,dim(N)[2]:1] + N[14,dim(N)[2]:1]) - (alpha[1]*alpha[2]*FOI[dim(N)[2]:1]+beta+delta[dim(N)[2]:1]+deathrate[dim(N)[2]:1,year])*N[6,dim(N)[2]:1];                                 
# S3a    
    sol[7,dim(N)[2]:1]  <- delta[c((dim(N)[2]-1):1,dim(N)[2])]*N[7,c((dim(N)[2]-1):1,dim(N)[2])] + props[7,dim(N)[2]:1]*migra[dim(N)[2]:1,year] + 
    beta*N[6,dim(N)[2]:1] - (alpha[1]*alpha[2]*FOI[dim(N)[2]:1]+beta+delta[dim(N)[2]:1]+deathrate[dim(N)[2]:1,year])*N[7,dim(N)[2]:1];                                                  
# S3b    
    sol[8,dim(N)[2]:1]  <- delta[c((dim(N)[2]-1):1,dim(N)[2])]*N[8,c((dim(N)[2]-1):1,dim(N)[2])] + props[8,dim(N)[2]:1]*migra[dim(N)[2]:1,year] + 
    beta*N[7,dim(N)[2]:1] - (alpha[1]*alpha[2]*FOI[dim(N)[2]:1]+beta+delta[dim(N)[2]:1]+deathrate[dim(N)[2]:1,year])*N[8,dim(N)[2]:1];
# I1                                                      
    sol[9,dim(N)[2]:1]  <- delta[c((dim(N)[2]-1):1,dim(N)[2])]*N[9,c((dim(N)[2]-1):1,dim(N)[2])] + props[9,dim(N)[2]:1]*migra[dim(N)[2]:1,year] + 
    theta[1]*FOI[dim(N)[2]:1]*N[2,dim(N)[2]:1] - (mu+delta[dim(N)[2]:1]+deathrate[dim(N)[2]:1,year])*N[9,dim(N)[2]:1];                                                                  
# A1        
    sol[10,dim(N)[2]:1] <- delta[c((dim(N)[2]-1):1,dim(N)[2])]*N[10,c((dim(N)[2]-1):1,dim(N)[2])]+ props[10,dim(N)[2]:1]*migra[dim(N)[2]:1,year] + 
    (1-theta[1])*FOI[dim(N)[2]:1]*N[2,dim(N)[2]:1] + mu*N[9,dim(N)[2]:1] - (omega+delta[dim(N)[2]:1]+deathrate[dim(N)[2]:1,year])*N[10,dim(N)[2]:1];                                           
# I2    
    sol[11,dim(N)[2]:1] <- delta[c((dim(N)[2]-1):1,dim(N)[2])]*N[11,c((dim(N)[2]-1):1,dim(N)[2])]+ props[11,dim(N)[2]:1]*migra[dim(N)[2]:1,year] + 
    alpha[1]*theta[2]*FOI[dim(N)[2]:1]*(N[3,dim(N)[2]:1]+N[4,dim(N)[2]:1]+N[5,dim(N)[2]:1]) - 
    (mu+delta[dim(N)[2]:1]+deathrate[dim(N)[2]:1,year])*N[11,dim(N)[2]:1];
# A2                                      
    sol[12,dim(N)[2]:1] <- delta[c((dim(N)[2]-1):1,dim(N)[2])]*N[12,c((dim(N)[2]-1):1,dim(N)[2])]+ props[12,dim(N)[2]:1]*migra[dim(N)[2]:1,year] + 
    alpha[1]*(1-theta[2])*FOI[dim(N)[2]:1]*(N[3,dim(N)[2]:1]+N[4,dim(N)[2]:1]+N[5,dim(N)[2]:1]) + mu*N[11,dim(N)[2]:1] - 
    (omega+delta[dim(N)[2]:1]+deathrate[dim(N)[2]:1,year])*N[12,dim(N)[2]:1];     
# I3          
    sol[13,dim(N)[2]:1] <- delta[c((dim(N)[2]-1):1,dim(N)[2])]*N[13,c((dim(N)[2]-1):1,dim(N)[2])]+ props[13,dim(N)[2]:1]*migra[dim(N)[2]:1,year] + 
    alpha[1]*alpha[2]*theta[3]*FOI[dim(N)[2]:1]*(N[6,dim(N)[2]:1]+N[7,dim(N)[2]:1]+N[8,dim(N)[2]:1]) - 
    (mu+delta[dim(N)[2]:1]+deathrate[dim(N)[2]:1,year])*N[13,dim(N)[2]:1];
# A3                             
    sol[14,dim(N)[2]:1] <- delta[c((dim(N)[2]-1):1,dim(N)[2])]*N[14,c((dim(N)[2]-1):1,dim(N)[2])]+ props[14,dim(N)[2]:1]*migra[dim(N)[2]:1,year] + 
    alpha[1]*alpha[2]*(1-theta[3])*FOI[dim(N)[2]:1]*(N[6,dim(N)[2]:1]+N[7,dim(N)[2]:1]+N[8,dim(N)[2]:1]) + mu*N[13,dim(N)[2]:1] - 
    (omega+delta[dim(N)[2]:1]+deathrate[dim(N)[2]:1,year])*N[14,dim(N)[2]:1];  
# number of new infections    
    sol[15,dim(N)[2]:1] <- theta[1]*FOI[dim(N)[2]:1]*N[2,dim(N)[2]:1] + alpha[1]*theta[2]*FOI[dim(N)[2]:1]*(N[3,dim(N)[2]:1]+
    N[4,dim(N)[2]:1]+N[5,dim(N)[2]:1]) + alpha[1]*alpha[2]*theta[3]*FOI[dim(N)[2]:1]*(N[6,dim(N)[2]:1]+N[7,dim(N)[2]:1]+N[8,dim(N)[2]:1]);
                              
return(list(c(sol))) # return left hand side of ODE system
}                                                   