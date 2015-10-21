

calc.mean <- function(f,x){
  dx = x[2]-x[1]
  return(sum(f*x)*dx)
}


calc.GI.fwdbck.theo <- function(R0, max.horizon, nE, nI, 
                                latent_mean, infectious_mean, 
                                N, I.init){ 
  ### CALCULATE BITH FWD & BCK GENERATION INTERVALS
  ### WITHIN A SEmInR FRAMEWORK
  
  tt <- c(1:max.horizon)
  
  ### SEMINR time series 
  # This one is used to calculate S and I:
  seminr <- solve.SEmInR(R0,max.horizon,
                         nE,nI,
                         latent_mean,
                         infectious_mean,
                         N,
                         I.init)
  # And this one to calculate the proba to be infectious F:
  seminr.proba <- solve.SEmInR.FX(R0,max.horizon,
                                  nE,nI,
                                  latent_mean,
                                  infectious_mean,
                                  N,
                                  I.init)
  
  S <- seminr$S
  I <- seminr$inc
  NT = nrow(seminr.proba)
  gi.time <- seminr.proba$time
  
  horizon <- max(seminr$time[I>1E-7])
  tth <- tt[1:horizon]
  mult.horiz <- 1.6  # <-- make sure this is >>1 in orderto have clean integration for large calendar times
  theo.GI <- calc.theoretical.GI(paste0(path.model,file.prm), 
                                 n.points.GI.crv = min(200,mult.horiz*horizon),
                                 horizon = mult.horiz*horizon,
                                 R0.override = R0,
                                 do.plot = FALSE)
  GI.ODE <- theo.GI[["GI.ODE"]]
  theo.gi.fwd <- theo.GI[["GI.fwd.theo"]]
  theo.gi.bck <- theo.GI[["GI.bck.theo"]]
  theo.gi.fwd.time <- theo.GI[["GI.fwd.theo.time"]]
  theo.gi.bck.time <- theo.GI[["GI.bck.theo.time"]]
  
  
  ### --- Calculate gi.intrinsic ---
  
  # fine time scale intrinsic GI
  gi.intrinsic = vector()
  for(tau in 1:NT) gi.intrinsic[tau] <- GI.intrinsic(tau,seminr.proba)
  mean.g <- calc.mean(gi.intrinsic,gi.time)
  
  ### plot esthetic features
  ylim.g.mult <- 1.35
  ylim.g.mult.bck <- 1.35
  lwd.gi = 5
  xlim.g <- 2*mean.g
  
  # coarser time scale intrinsic GI
  # where g[i] means g at i time units
  idx = which(seminr.proba$time==round(seminr.proba$time))
  g <- gi.intrinsic[idx]
  g <- g[tt] # <- makes sure same size (else larger by 1)
  
  ### ----- Backward GI -----
  b.bar <- GI.ODE$GI.bck.mean
  b.bar.time <- GI.ODE$time
  # Ends when incidence is tiny ~ end of epidemic
  b.bar <- b.bar[b.bar.time<horizon]
  b.bar.time <- b.bar.time[b.bar.time<horizon]
  
  ### ----- Forward GI -----
  f.bar <- GI.ODE$GI.fwd.mean
  f.bar.time <- GI.ODE$time
  # Ends when incidence is tiny ~ end of epidemic
  f.bar <- f.bar[f.bar.time<horizon]
  f.bar.time <- f.bar.time[f.bar.time<horizon]
  
  return(list(b.bar.time=b.bar.time,
              b.bar=b.bar,
              f.bar.time=f.bar.time,
              f.bar=f.bar,
              mean.g = mean.g,
              horizon=horizon))
}

