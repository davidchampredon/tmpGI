source("SEmInR_deterministic.R")

require("deSolve")
require("ggplot2")
require("gridExtra")


GI.intrinsic <- function(tau, df){
  ### CALCULATE INTRINSIC GENERATION INTERVAL DISTRIBUTION
  ### AT A GIVEN TIME SINCE INFECTION
  ### USE THE EXPRESSION: g(tau) = K(tau)/R0
  ### WITH ADDITIONAL ASSUMPTION MEAN INFECTIOUSNESS (LAMBDA)
  ### IS CONSTANT:
  ### g(tau) = F(tau)/integral(F(x)dx) 
  ### WITH F(x) PROBABILITY TO BE INFECTIOUS x TIME UNITS AFTER INFECTION
  res = 0
  N = length(df$time)
  dt = df$time[2]-df$time[1]
  JJ = df$Jall[1:N]
  integ  = sum(JJ)*dt
  res = df$Jall[tau]/integ
  return(res)
}

GI.fwd <- function(tau, s, df){
  ### CALCULATE THEORETICAL FORWARD GENERATION INTERVAL
  res = 0
  N = length(df$time)
  dt = df$time[2]-df$time[1]
  
  JJ = df$Jall[1:(N-s)]
  SS = df$S[(s+1):N]
  integ  = sum(JJ*SS)*dt
  res = df$Jall[tau]*df$S[tau+s]/integ
  return(res)
}


GI.bck <- function(tau, t, df, correct.first.inc=TRUE){
  ### CALCULATE THEORETICAL BACKWARD GENERATION INTERVAL
  res = 0
  
  # Incidence solved from the ODE (=I[0])
  # has an artificial spike at time t=0
  # b/c the incidence at the very first time step
  # (t=1*dt) is beta*S[0]*I[0] << I[0]
  # The integration is very sensitive to the shape
  # of incidence. If the incidence 'initial spike' is kept
  # then an artificial 'bump' for backward GI appears
  # especially for low R0 (for high R0, the explosive
  # incidence dwarfs the initial spike).
  
  # To solve this issue, get rid of the spike:
  if(correct.first.inc) df$inc[1]<-df$inc[2]
  #-----
  
  if(tau<t){
    N = length(df$time)
    dt = df$time[2]-df$time[1]
    
    JJ = df$Jall[1:(t-1)]
    II = df$inc[(t-1):1]
    integ = sum(JJ*II)*dt
    res = df$Jall[tau]*df$inc[t-tau]/integ
  }
  return(res)
}




calc.theoretical.GI <- function(file.prmset, 
                                n.points.GI.crv,
                                horizon,
                                R0.override = NULL,
                                do.plot = FALSE)
{
  ### Parameters
  prm <- read.csv(file.prmset,header=F)
  N <- prm[prm$V1=="popSize",2]
  dt <- seq(0,horizon,0.2) #0.2
  # Time unit = DAYS
  latent_mean = prm[prm$V1=="latent_mean",2]
  infectious_mean = prm[prm$V1=="infectious_mean",2]
  sigma <- 1/latent_mean
  gamma <- 1/infectious_mean
  if(is.null(R0.override))  R0 <- prm[prm$V1=="R0",2]
  if(!is.null(R0.override))  R0 <- R0.override
  beta <- R0*gamma
  f <- 0.0
  mu <- 0.00
  
  # Initial conditions
  I.init <- prm[prm$V1=="init_I1",2]/N
  print(paste("init_I1 =",prm[prm$V1=="init_I1",2]))
  print(paste("N =",N))
  print(paste("I.init =",I.init))
  S.init <- 1 - I.init
  E.init <- 0
  
  nE <- prm[prm$V1=="nE",2]
  nI <- prm[prm$V1=="nI",2]
  
  params.SEmInR <- c(mu=mu, 
                     beta=beta,
                     sigma=sigma,
                     gamma=gamma,
                     f=f,
                     nE=nE, nI=nI)
  
  ### Inital conditions
  # W A R N I N G : 
  # initial infectious in I[1] compartment (not E[1])
  
  EIvec <- c(E= rep(0,ifelse(nE==Inf,1,nE)),
             I= c(I.init,rep(0,ifelse(nI==Inf,0,nI-1))) )
  
  YJvec <- c(Y=c(1,rep(0,ifelse(nE==Inf,1,nE-1))),
             J=rep(0,ifelse(nI==Inf,0,nI) ))
  
  inits.SEmInR <- c(S=1-I.init, EIvec, R=0, YJvec, Z=I.init)
  
  #### Solutions of the SEmInR model ####
  SEmInR <- as.data.frame(lsoda(inits.SEmInR, dt, SEmInR.FX, 
                                parms=params.SEmInR))
  SEmInR <- calc.Iall(dat = SEmInR)
  SEmInR <- calc.Jall(dat = SEmInR)
  SEmInR <- calc.Yall(dat = SEmInR)
  
  # incidence (from solved cumulative incidence)
  SEmInR$inc <- c(SEmInR$Z[1],diff(SEmInR$Z))
  
  
  ##############################
  #### Generation intervals ####
  ##############################
  
  NT = nrow(SEmInR)
  
  # -- mean forward & backward generation interval
  tt = vector()
  f.bar = vector()
  g.bar = vector()
  
  # Sequence of indices used in the loop that goes through
  # all duration-since-infection points:
  loop.idx = round(seq(1,NT-1,length.out=n.points.GI.crv))
  dt = SEmInR$time[2]-SEmInR$time[1]
  
  ### --- Calculate gi.intrinsic ---
  gi.intrinsic = vector()
  for(tau in 1:NT) gi.intrinsic[tau] <- GI.intrinsic(tau,SEmInR)
  # theoretical mean: 
  theo.mean.gii <- latent_mean + infectious_mean*(nI+1)/2/nI
  # numerical mean:
  mean.gii <- sum(SEmInR$time*gi.intrinsic*dt)
  # variance:
  var.gii <- sum(SEmInR$time^2*gi.intrinsic*dt) - mean.gii^2
  
  gi.fwd.lst = list()
  gi.bck.lst = list()
  gi.fwd.time.lst = list()
  gi.bck.time.lst = list()
  
  cnt = 1
  print(paste("Integrating bck & fwd GI with",length(loop.idx),"steps"))
  
  for(s in loop.idx){
    cat(paste0(cnt,"."))
    
    # --- calculate gi.fwd:
    gi.fwd = vector()
    for(tau in 1:(NT-s))  gi.fwd[tau] <- GI.fwd(tau, s, SEmInR)
    
    # calculate expectation of gi.fwd:
    tt[cnt] = SEmInR$time[s]
    
    tvec = SEmInR$time[c(1:(NT-s))]
    f.bar[cnt] = sum(gi.fwd*tvec*dt)
    
    # --- calculate gi.bck:
    gi.bck = vector()
    for(tau in 1:(s-1)){
      gi.bck[tau] <- GI.bck(tau, s, SEmInR)
    }  
    
    # calculate expectation of gi.bck:
    tvec.bck = SEmInR$time[c(1:(s-1))]
    g.bar[cnt] = sum(gi.bck*tvec.bck*dt)
    
    gi.fwd.lst[[cnt]] <- gi.fwd
    gi.bck.lst[[cnt]] <- gi.bck
    gi.fwd.time.lst[[cnt]] <- tvec
    gi.bck.time.lst[[cnt]] <- tvec.bck
    cnt = cnt+1
  }
  
  GI.ODE <- data.frame(time=tt,
                       GI.fwd.mean = f.bar,
                       GI.bck.mean = g.bar)
  
  #############################
  ###### PLOTS TO CHECK #######
  #############################
  
  TIME = SEmInR$time
  # time since infection
  t.inf <- TIME[which(TIME<3*(latent_mean+infectious_mean))]
  n.t.inf <- length(t.inf)
  
  if(do.plot){
    par(mfrow=c(2,1))
    
    ### PREVALENCE & INCIDENCE
    plot(x=TIME, y=SEmInR$Iall, typ="l", col="red",lwd=6, main="Prevalence")
    plot(x=TIME, y=SEmInR$inc, typ="l", col="red",lwd=6, main="Incidence")
    
    ### PROBABILITIES
    plot(x=t.inf, y=SEmInR$Yall[1:n.t.inf], 
         typ="l",lwd=6, 
         xlab="Time since infection", ylab="",las=1,
         main="Probability of being in E[k]")
    col.Y <- which(grepl("Y",names(SEmInR)))
    
    for(i in 1:length(col.Y)){
      col.i = 1-i/length(col.Y)
      lines(x=t.inf, y=SEmInR[1:n.t.inf,col.Y[i]], 
            col=rgb(col.i,col.i,col.i))
    }
    abline(v=latent_mean,lty=3)
    
    plot(x=t.inf, y=SEmInR$Jall[1:n.t.inf], 
         typ="l",lwd=6, 
         xlab="Time since infection", ylab="",las=1,
         main="Probability of being in I[k]")
    col.J <- which(grepl("J",names(SEmInR)))
    
    for(i in 1:length(col.J)){
      col.i = 1-i/length(col.J)
      lines(x=t.inf, y=SEmInR[1:n.t.inf,col.J[i]], 
            col=rgb(col.i,col.i,col.i))  
    }
    abline(v=latent_mean+infectious_mean*(nI+1)/2/nI,lty=2)
    
    ### GENERATION INTERVALS
    
    # Forward:
    plot(x=tt,y=f.bar,pch=16,typ="o",
         main="Mean Forward GI",
         xlab = "Calendar time",
         ylim=c(0,max(f.bar,latent_mean+infectious_mean,na.rm = T)),
         las=1, lwd=6)
    abline(h=mean.gii,lty=2)
    
    # Backward:
    plot(x=tt,y=g.bar,pch=16,typ="o",
         main="Mean Backward GI",
         xlab = "Calendar time",
         ylim=c(0,max(g.bar,latent_mean+infectious_mean,na.rm = T)),
         las=1, lwd=6)
    abline(h=mean.gii,lty=2)
    
    # Intrinsic:
    # equivalent gamma distribution
    shape <- mean.gii^2/var.gii
    rate <- mean.gii/var.gii
    equiv.gamma <- dgamma(x=t.inf, shape=shape, rate=rate)
    
    title=paste0("Intrinsic GI \n",
                 "lat.mean=",latent_mean,
                 " ; infec.mean=",infectious_mean,
                 " ; nE=",nE, " ; nI=",nI)
    
    par(mfrow=c(1,1))
    plot(x=t.inf, y=gi.intrinsic[1:n.t.inf], 
         typ="l", 
         main = title,
         xlab="Time Since Infection",
         ylab="",las=1,
         ylim=range(gi.intrinsic[1:n.t.inf],equiv.gamma),
         lwd=6)
    lines(x=t.inf,y=equiv.gamma,
          col=rgb(1,0,0,0.6),lwd=6,lty=3)
    abline(v=mean.gii,lty=2,lwd=3)
    abline(v=theo.mean.gii,lty=1,col="green",lwd=3)
    legend(x = "topright", legend = c("GI integrated","Gamma(m,v)"),
           lwd=6, col=c("black","red"),lty=c(1,3))
  }
  return(list(GI.ODE=GI.ODE, 
              SEmInR=SEmInR, 
              GI.intrinsic=gi.intrinsic,
              GI.fwd.theo = gi.fwd.lst,
              GI.bck.theo = gi.bck.lst,
              GI.fwd.theo.time = gi.fwd.time.lst,
              GI.bck.theo.time = gi.bck.time.lst,
              time.vec = tt)
  )
}


calc.theoretical.GI.base <- function(N,
                                     latent_mean,
                                     infectious_mean ,
                                     R0,
                                     init_I1,
                                     nE,
                                     nI,
                                     n.points.GI.crv,
                                     horizon,
                                     do.plot = FALSE)
{
  ### Parameters
  dt <- seq(0,horizon,0.2) #0.2
  # Time unit = DAYS
  sigma <- 1/latent_mean
  gamma <- 1/infectious_mean
  beta <- R0*gamma
  f <- 0.0
  mu <- 0.00
  
  # Initial conditions
  I.init = init_I1/N
  S.init <- 1 - I.init
  E.init <- 0
  print(paste("N =",N))
  print(paste("init_I1 =",init_I1))
  print(paste("I.init =",I.init))
  
  params.SEmInR <- c(mu=mu, 
                     beta=beta,
                     sigma=sigma,
                     gamma=gamma,
                     f=f,
                     nE=nE, nI=nI)
  
  ### Inital conditions
  # W A R N I N G : 
  # initial infectious in I[1] compartment (not E[1])
  
  EIvec <- c(E= rep(0,ifelse(nE==Inf,1,nE)),
             I= c(I.init,rep(0,ifelse(nI==Inf,0,nI-1))) )
  
  YJvec <- c(Y=c(1,rep(0,ifelse(nE==Inf,1,nE-1))),
             J=rep(0,ifelse(nI==Inf,0,nI) ))
  
  inits.SEmInR <- c(S=1-I.init, EIvec, R=0, YJvec, Z=I.init)
  
  #### Solutions of the SEmInR model ####
  SEmInR <- as.data.frame(lsoda(inits.SEmInR, dt, SEmInR.FX, 
                                parms=params.SEmInR))
  SEmInR <- calc.Iall(dat = SEmInR)
  SEmInR <- calc.Jall(dat = SEmInR)
  SEmInR <- calc.Yall(dat = SEmInR)
  
  # incidence (from solved cumulative incidence)
  SEmInR$inc <- c(SEmInR$Z[1],diff(SEmInR$Z))
  
  
  ##############################
  #### Generation intervals ####
  ##############################
  
  NT = nrow(SEmInR)
  
  # -- mean forward & backward generation interval
  tt = vector()
  f.bar = vector()
  g.bar = vector()
  f.var = vector()
  g.var = vector()
  
  
  # Sequence of indices used in the loop that goes through
  # all duration-since-infection points:
  loop.idx = round(seq(1,NT-1,length.out=n.points.GI.crv))
  dt = SEmInR$time[2]-SEmInR$time[1]
  
  ### --- Calculate gi.intrinsic ---
  gi.intrinsic = vector()
  for(tau in 1:NT) gi.intrinsic[tau] <- GI.intrinsic(tau,SEmInR)
  # theoretical mean: 
  theo.mean.gii <- latent_mean + infectious_mean*(nI+1)/2/nI
  # numerical mean:
  mean.gii <- sum(SEmInR$time*gi.intrinsic*dt)
  # variance:
  var.gii <- sum(SEmInR$time^2*gi.intrinsic*dt) - mean.gii^2
  
  gi.fwd.lst = list()
  gi.bck.lst = list()
  gi.fwd.time.lst = list()
  gi.bck.time.lst = list()
  
  cnt = 1
  print(paste("Integrating bck & fwd GI with",length(loop.idx),"steps"))
  
  for(s in loop.idx){
    cat(paste0(cnt,"."))
    
    # --- calculate gi.fwd:
    gi.fwd = vector()
    for(tau in 1:(NT-s))  gi.fwd[tau] <- GI.fwd(tau, s, SEmInR)
    
    # calculate expectation of gi.fwd:
    tt[cnt] = SEmInR$time[s]
    
    tvec = SEmInR$time[c(1:(NT-s))]
    f.bar[cnt] <- sum(gi.fwd*tvec*dt)
    f.var[cnt] <- sum(gi.fwd*tvec^2*dt)-f.bar[cnt]^2
    
    # --- calculate gi.bck:
    gi.bck = vector()
    for(tau in 1:(s-1)){
      gi.bck[tau] <- GI.bck(tau, s, SEmInR)
    }  
    
    # calculate expectation of gi.bck:
    tvec.bck = SEmInR$time[c(1:(s-1))]
    g.bar[cnt] <- sum(gi.bck*tvec.bck*dt)
    g.var[cnt] <- sum(gi.bck*tvec.bck^2*dt)-g.bar[cnt]^2
    
    gi.fwd.lst[[cnt]] <- gi.fwd
    gi.bck.lst[[cnt]] <- gi.bck
    gi.fwd.time.lst[[cnt]] <- tvec
    gi.bck.time.lst[[cnt]] <- tvec.bck
    cnt = cnt+1
  }
  
  GI.ODE <- data.frame(time=tt,
                       GI.fwd.mean = f.bar,
                       GI.bck.mean = g.bar,
                       GI.fwd.var = f.var,
                       GI.bck.var = g.var)
  
  #############################
  ###### PLOTS TO CHECK #######
  #############################
  
  TIME = SEmInR$time
  # time since infection
  t.inf <- TIME[which(TIME<3*(latent_mean+infectious_mean))]
  n.t.inf <- length(t.inf)
  
  if(do.plot){
    par(mfrow=c(2,1))
    
    ### PREVALENCE & INCIDENCE
    plot(x=TIME, y=SEmInR$Iall, typ="l", col="red",lwd=6, main="Prevalence")
    plot(x=TIME, y=SEmInR$inc, typ="l", col="red",lwd=6, main="Incidence")
    
    ### PROBABILITIES
    plot(x=t.inf, y=SEmInR$Yall[1:n.t.inf], 
         typ="l",lwd=6, 
         xlab="Time since infection", ylab="",las=1,
         main="Probability of being in E[k]")
    col.Y <- which(grepl("Y",names(SEmInR)))
    
    for(i in 1:length(col.Y)){
      col.i = 1-i/length(col.Y)
      lines(x=t.inf, y=SEmInR[1:n.t.inf,col.Y[i]], 
            col=rgb(col.i,col.i,col.i))
    }
    abline(v=latent_mean,lty=3)
    
    plot(x=t.inf, y=SEmInR$Jall[1:n.t.inf], 
         typ="l",lwd=6, 
         xlab="Time since infection", ylab="",las=1,
         main="Probability of being in I[k]")
    col.J <- which(grepl("J",names(SEmInR)))
    
    for(i in 1:length(col.J)){
      col.i = 1-i/length(col.J)
      lines(x=t.inf, y=SEmInR[1:n.t.inf,col.J[i]], 
            col=rgb(col.i,col.i,col.i))  
    }
    abline(v=latent_mean+infectious_mean*(nI+1)/2/nI,lty=2)
    
    ### GENERATION INTERVALS
    
    # Forward:
    plot(x=tt,y=f.bar,pch=16,typ="o",
         main="Mean Forward GI",
         xlab = "Calendar time",
         ylim=c(0,max(f.bar,latent_mean+infectious_mean,na.rm = T)),
         las=1, lwd=6)
    abline(h=mean.gii,lty=2)
    
    # Backward:
    plot(x=tt,y=g.bar,pch=16,typ="o",
         main="Mean Backward GI",
         xlab = "Calendar time",
         ylim=c(0,max(g.bar,latent_mean+infectious_mean,na.rm = T)),
         las=1, lwd=6)
    abline(h=mean.gii,lty=2)
    
    # Intrinsic:
    # equivalent gamma distribution
    shape <- mean.gii^2/var.gii
    rate <- mean.gii/var.gii
    equiv.gamma <- dgamma(x=t.inf, shape=shape, rate=rate)
    
    title=paste0("Intrinsic GI \n",
                 "lat.mean=",latent_mean,
                 " ; infec.mean=",infectious_mean,
                 " ; nE=",nE, " ; nI=",nI)
    
    par(mfrow=c(1,1))
    plot(x=t.inf, y=gi.intrinsic[1:n.t.inf], 
         typ="l", 
         main = title,
         xlab="Time Since Infection",
         ylab="",las=1,
         ylim=range(gi.intrinsic[1:n.t.inf],equiv.gamma),
         lwd=6)
    lines(x=t.inf,y=equiv.gamma,
          col=rgb(1,0,0,0.6),lwd=6,lty=3)
    abline(v=mean.gii,lty=2,lwd=3)
    abline(v=theo.mean.gii,lty=1,col="green",lwd=3)
    legend(x = "topright", legend = c("GI integrated","Gamma(m,v)"),
           lwd=6, col=c("black","red"),lty=c(1,3))
  }
  return(list(GI.ODE=GI.ODE, 
              SEmInR=SEmInR, 
              GI.intrinsic=gi.intrinsic,
              GI.fwd.theo = gi.fwd.lst,
              GI.bck.theo = gi.bck.lst,
              GI.fwd.theo.time = gi.fwd.time.lst,
              GI.bck.theo.time = gi.bck.time.lst,
              time.vec = tt)
  )
}