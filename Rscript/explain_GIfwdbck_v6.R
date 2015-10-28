############################################################
###   
###   CREATE FIGURES EXPLAING TEMPORAL
###   VARIATION OF BACKWARD AND 
###   FORWARD GENERATION INTERVAL DISTRIBUTIONS
###
###   Created 2015-03-16 by David Champredon
###   Modified 2015-07-04 by David Champredon
###
############################################################

source("SEmInR_deterministic.R")
source("calc_theoretical_GI.R")

### Global Variables
save.to.file <- TRUE
do.bck <- TRUE
do.fwd <- TRUE
do.vert.density <- FALSE

# Path to the C++ model generating simulations
path.model <- "../Gillespie_SEmInR/"
simprm.list <- as.character(read.csv(paste0(path.model,"param_all_list.csv"),header=F)[,1])  
file.prm <- simprm.list[1]

# Retrieve SEmInR model parameters
prm <- read.csv(paste0(path.model,file.prm),header=FALSE)
max.horizon <- prm[prm$V1=="horizon",2]
nE <- prm[prm$V1=="nE",2]
nI <- prm[prm$V1=="nI",2]
latent_mean <- prm[prm$V1=="latent_mean",2]
infectious_mean <- prm[prm$V1=="infectious_mean",2]
N <- prm[prm$V1=="popSize",2] ; popSize <- N
R0 <- prm[prm$V1=="R0",2]
I.init <- prm[prm$V1=="init_I1",2]/N

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

horizon <- max(seminr$time[seminr$inc>1E-7])
tth <- tt[1:horizon]
mult.horiz <- 1.6  # <-- make sure this is >>1 in orderto have clean integration for large calendar times
theo.GI <- calc.theoretical.GI(paste0(path.model,file.prm), 
                               n.points.GI.crv = min(200,mult.horiz*horizon),
                               horizon = mult.horiz*horizon,
                               do.plot = FALSE)
GI.ODE <- theo.GI[["GI.ODE"]]
theo.gi.fwd <- theo.GI[["GI.fwd.theo"]]
theo.gi.bck <- theo.GI[["GI.bck.theo"]]
theo.gi.fwd.time <- theo.GI[["GI.fwd.theo.time"]]
theo.gi.bck.time <- theo.GI[["GI.bck.theo.time"]]

calc.mean <- function(f,x){
  dx = x[2]-x[1]
  return(sum(f*x)*dx)
}


### --- Calculate gi.intrinsic ---

# fine time scale intrinsic GI
gi.intrinsic = vector()
for(tau in 1:NT) gi.intrinsic[tau] <- GI.intrinsic(tau,seminr.proba)
mean.g <- calc.mean(gi.intrinsic,gi.time)

# coarser time scale intrinsic GI
# where g[i] means g at i time units
idx = which(seminr.proba$time==round(seminr.proba$time))
g <- gi.intrinsic[idx]
g <- g[tt] # <- makes sure same size (else larger by 1)

### plot esthetic features
ylim.g.mult <- 1.35
ylim.g.mult.bck <- 1.35
lwd.gi = 5
xlim.g <- 2*mean.g




plot.g <- function(caltime){
  ### Plot the intrinsic generation interval distribution
  plot(x= seminr.proba$time, y=gi.intrinsic,
       typ="l", 
       main=paste0("Intrinsic GI Distribution\n (at calendar time t=",caltime,")"), 
       lwd=4,
       xlim=c(0,xlim.g),ylim=c(0,max(g)*ylim.g.mult),
       xlab="time since infection", ylab="",yaxt="n")  
  abline(v=mean.g,lty=2,lwd=1)
}

plot.S <- function(slope,inflection,time.s){
  plot(tt[1:horizon],S[1:horizon],typ="l", 
       main=paste0("Susceptible Proportion\n(at calendar time s=",time.s,")"), 
       lwd=1, lty=3,
       xlim=c(0,horizon), ylim=c(0,1),
       xlab="calendar time", ylab="",yaxt="n")
  yval <- c(0,0.5,1)
  axis(side = 2,at = yval,labels = yval)
  
  # Time direction arrows
  y.arrow = max(S)/3
  arrows(x0 = time.s,y0 = y.arrow, 
         x1=time.s+xlim.g, y1=y.arrow,
         length=0.1,col="darkgray", lwd=2)
  
  # bold parts
  # lines(tt[time.s:horizon],S[time.s:horizon],typ="l", lwd=2 )
  lines(tt[1:time.s],S[1:time.s],typ="l", lwd=2 )
  lines(tt[time.s:(time.s+xlim.g)],S[time.s:(time.s+xlim.g)],typ="l", lwd=6 )
  points(tt[time.s],S[time.s],cex=3)
  
  # Gray-shaded area
  bx1 = c(time.s,time.s+xlim.g)
  by1 = c(0,0)
  bx2 = rev(bx1)
  by2 = c(1,1)  
  polygon(x = c(bx1,bx2),y = c(by1,by2), col=rgb(0,0,0,0.15),border = NA)
}


plot.I <- function(time.t)
{
  plot(tt[1:horizon],I[1:horizon],typ="l", 
       main=paste0("Incidence\n(at calendar time t=",time.t,")"), 
       lwd=1,lty=3,
       xlim=c(0,horizon),
       xlab="calendar time", yaxt="n",ylab="")
  ymax <- round(max(I[1:horizon])*0.97,digits = 2)
  yval <- c(0,ymax/2,ymax)
  axis(side = 2,at = yval,labels = yval)
  
  # bold part
  bp = max(0,time.t-xlim.g):time.t
  lines(tt[1:time.t],I[1:time.t],lwd=2)
  lines(tt[bp],I[bp],lwd=6)
  points(tt[time.t],I[time.t],cex=3)
  # Gray-shaded area
  bx1 = c(time.t-xlim.g,time.t)
  by1 = c(0,0)
  bx2 = rev(bx1)
  by2 = rep(max(I),2) 
  polygon(x = c(bx1,bx2),y = c(by1,by2), col=rgb(0,0,0,0.15),border = NA)
  # Direction of time arrow
  y.arrow = max(I)/5
  arrows(x0 = time.t,y0 = y.arrow, x1=time.t-xlim.g, y1=y.arrow,
         length=0.1,col="darkgray", lwd=2)
}



plot.GI.fwd <- function(time.s, g, S, horizon,tt)
{
  
  fwd = theo.gi.fwd[[time.s]]
  fwd.time = theo.gi.fwd.time[[time.s]]
  m.GI.fwd = calc.mean(fwd, fwd.time)
  
  plot(x=fwd.time, y=fwd,
       typ="l", lwd=8,
       xlim=c(0,xlim.g),ylim=c(0,max(g)*ylim.g.mult ),
       xlab="time since infection", ylab="",yaxt="n",
       main = paste0("Forward GI Distribution\n(at calendar time s=",time.s,")") )
  
  lines(gi.time,gi.intrinsic,col="gray")
  
  abline(v=mean.g, lty=2,lwd=1)
  abline(v=m.GI.fwd, lty=1,lwd=3)
}



plot.GI.bck <- function(time.t,tt,g, I, correct.first.inc)
{
  bck = theo.gi.bck[[time.t]]
  bck.time = theo.gi.bck.time[[time.t]]
  m.GI.bck = calc.mean(bck, bck.time)
  
  plot(x=bck.time, y=bck, 
       typ="l", lwd=8,
       xlim=c(0,xlim.g),ylim=c(0,max(g)*ylim.g.mult.bck ),
       xlab="time since infection", ylab="",yaxt="n",
       main = paste0("Backward GI Distribution\n(at calendar time s=",time.t,")") )
  
  lines(gi.time,gi.intrinsic,col="gray")
  abline(v=mean.g, lty=2,lwd=1)
  abline(v=m.GI.bck, lty=1,lwd=3)
}


##################
###    PLOTS   ###
##################

if(do.bck){
  if(save.to.file) pdf("explain_GI_bck.pdf", height=15,width=10)
  
  layout(matrix(c(1,2,3,4,5,6,7,8,9,10,10,10),
                nrow = 4, ncol=3, byrow = TRUE), 
         widths=c(1,1,1), heights=c(1,1,1,1.5))
  
  par(cex.main=1.7, cex.lab=1.5, 
      cex.axis=1.5, las = 1)
  
  time.t.vec <- c(20,48,70)
  
  correct.first.inc <- FALSE
  
  ### ----- First Row -----
  for (i in 1:3) plot.g(caltime=time.t.vec[i])
  
  ### ----- Second Row -----
  for (i in 1:3) plot.I(time.t=time.t.vec[i])
  
  ### ----- Third Row -----
  for (i in 1:3) plot.GI.bck(time.t = time.t.vec[i],tt,g,I,correct.first.inc)
  
  ### ----- Fourth Row -----
  
  b.bar.plot <- GI.ODE$GI.bck.mean
  b.bar.plot.time <- GI.ODE$time
  
  plot(x=b.bar.plot.time, y=b.bar.plot, 
       main = "Backward GI Mean",
       xlim=c(0,horizon), xlab="calendar time",
       ylim=c(0,max(b.bar.plot,
                    xlim.g*do.vert.density,na.rm = T)), 
       ylab="",
       lwd=lwd.gi,
       typ="l",
       las=1)
  
  idx <- which(floor(b.bar.plot.time)%in%time.t.vec)
  points(x=b.bar.plot.time[idx], y=b.bar.plot[idx],cex=3,pch=1)
  points(x=b.bar.plot.time[idx], y=b.bar.plot[idx],cex=2,pch=16)
  
  if(do.vert.density){
    for(i in 1:length(time.t.vec)){
      gib <- GI.bck.g(time.t.vec[i],g, I)[1:xlim.g]
      xi <- rep(time.t.vec[i],length(gib))
      lines(x=xi+gib*30, y=1:xlim.g,lwd=1,col="grey")
    }
  }
  abline(h=mean.g,lty=2)
  
  if(save.to.file) dev.off()
}

if(do.fwd){
  if(save.to.file) pdf("explain_GI_fwd.pdf", height=15,width=10)
  
  layout(matrix(c(1,2,3,4,5,6,7,8,9,10,10,10),
                nrow = 4, ncol=3, byrow = TRUE), 
         widths=c(1,1,1), heights=c(1,1,1,1.5))
  
  par(cex.main=1.7, cex.lab=1.5, 
      cex.axis=1.5, las = 1)
  
  time.s.vec <- c(10,38,50)
  
  ### ----- First Row -----
  for (i in 1:3) plot.g(time.s.vec[i])	
  
  ### ----- Second Row -----
  for (i in 1:3) plot.S(slope = slope,
                        inflection = inflection,
                        time.s=time.s.vec[i])
  
  ### ----- Third Row -----
  for (i in 1:3) plot.GI.fwd(time.s = time.s.vec[i],g, S, horizon,tt)
  
  ### ----- Fourth Row -----
  
  f.bar.plot <- GI.ODE$GI.fwd.mean
  f.bar.plot.time <- GI.ODE$time
  t.idx = which(round(f.bar.plot.time)%in%time.s.vec)
  f.bar <- f.bar.plot[t.idx]
  
  plot(x=f.bar.plot.time, y=f.bar.plot, 
       main = "Forward GI Mean",
       xlim=c(0,horizon*0.7), xlab="calendar time",
       ylim=c(min(f.bar)*0.9,mean.g*1.1), ylab="",
       lwd=lwd.gi,
       typ="l",
       las=1)
  
    points(x=f.bar.plot.time[t.idx], y=f.bar.plot[t.idx],cex=3,pch=1)
    points(x=f.bar.plot.time[t.idx], y=f.bar.plot[t.idx],cex=2,pch=16)
  
  abline(h=mean.g,lty=2)
  
  if(save.to.file) dev.off()
}
