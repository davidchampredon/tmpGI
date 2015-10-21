library("deSolve")


SEmInR.beta.gamma <- function(t, x, parms)
{
  ### (Code inspired from Ben Bolker's)
  
  with(c(as.list(parms),as.list(x)), {
    
    S = x[1]
    E = x[2:(nE+1)]
    I = x[(nE+2):(nE+nI+1)]
    R = x[nE+nI+2]
    C = x[(nE+2)]
    D = x[nE+nI+2]
    K = x[(nE+2)] # cumul symptoms onset
    
    sigma2=sigma*nE ## transition rate through E boxes
    gamma2=gamma*nI ## transition rate through I boxes
    
    infrate = beta*S*sum(I)
    
    dS = mu * (1  - S)
    dS = dS - infrate
    
    ## Exposed
    if (nE>1) {
      dE = sigma2*(c(0,E[1:(nE-1)])-E)
    } else dE=-E*sigma
    dE[1] = dE[1] + infrate
    dE = dE - mu*E
    
    ## Infected
    if (nI>1) {
      dI = gamma2*(c(0,I[1:(nI-1)])-I)
    } else dI=-I*gamma
    onset = sigma2*E[nE]
    
    dI[1] = dI[1]+onset
    dI = dI - mu*I
    recovery = gamma2*I[nI]
    
    ## Recovered 
    dR = (1-f)*recovery - mu*R
    
    ## Cumulative Incidence
    dC = infrate
    
    ## Cumulative symptoms onset
    dK = onset
    
    ## Cumulative RECOVERY
    dD = recovery
    
    res=c(dS, dE, dI, dR, dC, dD, dK)
    list(res)
  })
}



calc.Iall <- function(dat)
{
  ### Sum number of all individuals in any "I" compartment
  
  col.I <- which(grepl("I",names(dat)))
  if(length(col.I)>1) dat$Iall <- rowSums(dat[,col.I])
  if(length(col.I)==1) dat$Iall <- dat[,col.I]
  return(dat)
}

calc.Eall <- function(dat)
{
  ### Sum number of all individuals in any "I" compartment
  col.I <- which(grepl("E",names(dat)))
  if(length(col.I)>1) dat$Eall <- rowSums(dat[,col.I])
  if(length(col.I)==1) dat$Eall <- dat[,col.I]
  return(dat)
}


calc.Jall <- function(dat)
{
  col.I <- which(grepl("J",names(dat)))
  if(length(col.I)>1) dat$Jall <- rowSums(dat[,col.I])
  if(length(col.I)==1) dat$Jall <- dat[,col.I]
  return(dat)
}

calc.Yall <- function(dat)
{
  col.I <- which(grepl("Y",names(dat)))
  if(length(col.I)>1) dat$Yall <- rowSums(dat[,col.I])
  if(length(col.I)==1) dat$Yall <- dat[,col.I]
  return(dat)
}



SEmInR.FX <- function(t, x, parms)
{
  ### CALCULATE PROBABILITY TO BE INFECTIOUS
  
  with(c(as.list(parms),as.list(x)), {
    
    S = x[1]
    E = x[2:(nE+1)]
    I = x[(nE+2):(nE+nI+1)]  # prevalence
    R = x[nE+nI+2]
    Y = x[(nE+nI+3):(2*nE+nI+2)]  # proba to be in Ek
    J = x[(2*nE+nI+3):(2*nE+2*nI+2)] # proba to be in Ik
    Z = x[2*nE+2*nI+3]  # cumulative incidence
    
    sigma2=sigma*nE ## transition rate through E boxes
    gamma2=gamma*nI ## transition rate through I boxes
    
    infrate = beta*S*sum(I)
    
    dS = mu * (1  - S)
    dS = dS - infrate
    
    ## Exposed
    if (nE>1) {
      dE = sigma2*(c(0,E[1:(nE-1)])-E)
    } else dE=-E*sigma
    dE[1] = dE[1] + infrate
    dE = dE - mu*E
    
    ## Infected
    if (nI>1) {
      dI = gamma2*(c(0,I[1:(nI-1)])-I)
    } else dI=-I*gamma
    onset = sigma2*E[nE]
    
    dI[1] = dI[1]+onset
    dI = dI - mu*I
    recovery = gamma2*I[nI]
    
    ## Recovered 
    dR = (1-f)*recovery - mu*R
    
    
    ## Proba being in Ek
    if (nE>1) {
      dY = sigma2*(c(0,Y[1:(nE-1)])-Y)
    } else dY=-Y*sigma
    dY[1] = dY[1] + 0*infrate
    
    ## Proba being in Ik
    if (nI>1) {
      dJ = gamma2*(c(0,J[1:(nI-1)])-J)
    } else dJ=-J*gamma
    onset.proba = sigma2*Y[nE]
    
    dJ[1] = dJ[1]+onset.proba
    
    dZ = infrate
    
    res=c(dS, dE, dI, dR, dY, dJ, dZ)
    list(res)
  })
}

solve.SEmInR <- function(R0,horizon,
                         nE,nI,
                         latent_mean,
                         infectious_mean,
                         N,
                         I.init){
  ### Solve the SEmInR
  ### ODEs equations
  
  # Parameters
  dt <- seq(0,horizon,1)
  sigma <- 1/latent_mean
  gamma <- 1/infectious_mean
  beta <- R0*gamma
  f <- 0.0
  mu <- 0.00
  
  # Initial conditions
  S.init <- 1 - I.init
  E.init <- 0
  
  params.SEmInR <- c(mu=mu, 
                     beta=beta,
                     sigma=sigma,
                     gamma=gamma,
                     f=f,
                     nE=nE, nI=nI)
  
  ### Inital conditions
  EIvec <- c(E=rep(I.init, # <-- WARNING - FOR PLOT ESTHETIC ONLY! should be 0, but would introduce an artificial spike in incidence. 
                   ifelse(nE==Inf,1,nE)),
             I=c(I.init,rep(0,ifelse(nI==Inf,0,nI-1))))
  inits.SEmInR <- c(S=1-I.init,EIvec,R=0,C=I.init,D=0,K=1)
  
  #### Solutions
  SEmInR <- as.data.frame(lsoda(inits.SEmInR, dt, 
                                SEmInR.beta.gamma, 
                                parms=params.SEmInR))
  SEmInR <- calc.Iall(SEmInR)  # <- global prevalence
  SEmInR$inc <- c(I.init,diff(SEmInR$C))  # <- global incidence
  return(SEmInR)
}

solve.SEmInR.FX <- function(R0,horizon,
                            nE,nI,
                            latent_mean,
                            infectious_mean,
                            N,
                            I.init){
  ### Solve the SEmInR
  ### ODEs equations
  ### with probability to be infectious
  
  # Parameters
  dt <- seq(0,horizon,0.2) #0.2
  sigma <- 1/latent_mean
  gamma <- 1/infectious_mean
  beta <- R0*gamma
  f <- 0.0
  mu <- 0.00
  
  params.SEmInR <- c(mu=mu, 
                     beta=beta,
                     sigma=sigma,
                     gamma=gamma,
                     f=f,
                     nE=nE, nI=nI)
  
  # Inital conditions
  # W A R N I N G : 
  # initial infectious in I[1] compartment (not E[1])
  
  S.init <- 1 - I.init
  E.init <- 0
  
  EIvec <- c(E= rep(0,ifelse(nE==Inf,1,nE)),
             I= c(I.init,rep(0,ifelse(nI==Inf,0,nI-1))) )
  
  YJvec <- c(Y=c(1,rep(0,ifelse(nE==Inf,1,nE-1))),
             J=rep(0,ifelse(nI==Inf,0,nI) ))
  
  inits.SEmInR <- c(S=1-I.init, EIvec, R=0, YJvec, Z=I.init)
  
  # Solutions 
  SEmInR <- as.data.frame(lsoda(inits.SEmInR, dt,
                                SEmInR.FX, 
                                parms=params.SEmInR))
  SEmInR <- calc.Iall(dat = SEmInR)
  SEmInR <- calc.Jall(dat = SEmInR)
  SEmInR <- calc.Yall(dat = SEmInR)
  
  # incidence (from solved cumulative incidence)
  SEmInR$inc <- c(SEmInR$Z[1],diff(SEmInR$Z))
  
  return(SEmInR)
}





