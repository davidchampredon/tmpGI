#####################################################################
###
###   READ SIMULATIONS RUN IN C++
###
###   PLOT GENERATION INTERVAL DISTRIBUTIONS OF
###   SIMULATIONS AND ANALYTICAL FORMULA (SOLVED NUMERICALLY)
###
###   Created: 2015-07-03 by David Champredon
###
#####################################################################


library(ggplot2); theme_set(theme_bw())
library(plyr)
library("gridExtra")

source("calc_theoretical_GI.R")
source("figures_ms.R")
source("read_simul_GI_FCT.R")

saveplot = TRUE
doParallel = FALSE  # if multiple CPUs used (WARNING: not really faster!)

### Check flags
chk.prev <- T#FALSE
chk.cumInc <- T#FALSE


### Retrieve the number of MC iterations and parameter sets
simprm.list <- as.character(read.csv("param_all_list.csv",header=F)[,1])  
n.prm.set <- length(simprm.list)


### Loop on all parameter sets ###

for(p in 1:n.prm.set)
{
  print(paste0("Processing parameter set #",p," (",simprm.list[p],")...  "))
  
  ### Retrieve simulation parameters
  file.param <- simprm.list[p]

  simprm <- read.csv(file.param,header=F)
  R0 <- simprm[simprm$V1=="R0",2]
  nE <- simprm[simprm$V1=="nE",2]
  nI <- simprm[simprm$V1=="nI",2]
  mc <- simprm[simprm$V1=="mc_iter",2]
  popSize <- simprm[simprm$V1=="popSize",2]
  latent_mean <- simprm[simprm$V1=="latent_mean",2]
  infectious_mean <- simprm[simprm$V1=="infectious_mean",2]
  horiz <- simprm[simprm$V1=="horizon",2]
  fileplot = paste0("_R0_",R0,"_nE_",nE,"_nI_",nI,
                    "_lat_",latent_mean,"_inf_",infectious_mean,
                    "_pop_",popSize/1000,"k_MC_",mc,
                    ".pdf")
  
  ########################
  ### Retrieve simulations
  ########################
  
  # Slim data frames
  t.bucket = 0.02
  
  # -- Prevalence --
  
  prev <- read.mc.files(doParallel,"prev",file.paramset=simprm.list[p])[["dat"]]
  names(prev)=c("time","n","mc")
  
  prev$time2 = round(prev$time/t.bucket)*t.bucket
  d.prev = ddply(prev,c("time2"),summarize,
                 min=min(n),
                 m=mean(n),
                 max=max(n))
  
  max.horizon = max(d.prev$time2)
  
  # -- Cumulative incidence --
  if(chk.cumInc){
    cumInc <- read.mc.files(doParallel,"cumInc",file.paramset=simprm.list[p])[["dat"]]
    names(cumInc)=c("time","n","mc")
    
    cumInc$time2 = round(cumInc$time/t.bucket)*t.bucket
    d.cumInc = ddply(cumInc,c("time2"),summarize,
                     min=min(n),
                     m=mean(n),
                     max=max(n))  
  }
  
  myq <- 0.50
  
  # -- Backward GI --
  GIbck <- get.GI.bck.sim(doParallel,file.param,t.bucket)
  
  # -- Forward GI --
  GIfwd <- get.GI.fwd.sim(doParallel,file.param,t.bucket)
  GIfwd2 = melt(GIfwd,id.vars = c("time.infector","mc"),na.rm = T)
  GIfwd2$time.infector2 = round(GIfwd2$time.infector/t.bucket)*t.bucket
  names(GIfwd2)[names(GIfwd2)=="value"] <- "GI"
  
  ### Calculate theoretical 
  ### forward & backward generation intervals
  theo.GI <- calc.theoretical.GI(simprm.list[p], 
                                 n.points.GI.crv = min(200,max.horizon),
                                 horizon = 1.02*max.horizon,
                                 do.plot = FALSE)
  
  GI.ODE <- theo.GI[["GI.ODE"]]
  SEmInR <- theo.GI[["SEmInR"]]
  
  # Merge all results in one single data frame
  GIfwd2$param.set <- p
  GIbck$param.set <- p
  GI.ODE$param.set <- p
  
  if(p==1){
    GIfwd2.all <- GIfwd2
    GIbck.all <- GIbck
    GI.ODE.all <- GI.ODE
  }
  
  if(p>1){
    GIfwd2.all <- rbind(GIfwd2.all, GIfwd2)
    GIbck.all <- rbind(GIbck.all, GIbck)
    GI.ODE.all <- rbind(GI.ODE.all, GI.ODE)
  }
  
  # Mean intrinsic GI 
  mean.gi.intrinsic = latent_mean + infectious_mean*(nI+1)/nI/2
  
  
  ####################################################################
  ###    P L O T S
  ####################################################################
  
  info.prm <- paste0("R0=",R0," ; nE=",nE," ; nI=",nI,
                     " ; lat=",latent_mean," ; infect=",infectious_mean,
                     " ; PopSize=",popSize/1000,"k ; MC=",mc)
  
  if(saveplot) {
    
    # Forward generation interval
    
    pdf(file=paste0("plot_fwd",fileplot),width=15,height=10)
    g <- plot.gilles.vs.analytic(dat.gil = GIfwd,
                            fwdOrBck = "fwd",
                            dat.ode = GI.ODE,
                            n.times = 50, 
                            myq = myq,
                            title = "Mean forward GI: theory vs. simulations",
                            infoprm = info.prm,
                            mean.gi.intrinsic = mean.gi.intrinsic,
                            show.quantile = FALSE,
                            min.mc = mc/4)
    plot(g)
    dev.off()
    
    # Backward generation interval
    
    pdf(file=paste0("plot_bck",fileplot),width=15,height=10)
    g <- plot.gilles.vs.analytic(dat.gil = GIbck,
                            fwdOrBck = "bck",
                            dat.ode = GI.ODE,
                            n.times = 50, myq,
                            title = "Mean forward GI: theory vs. simulations",
                            infoprm=info.prm,
                            mean.gi.intrinsic=mean.gi.intrinsic,
                            show.quantile=FALSE,
                            min.mc = mc/4)
    plot(g)
    dev.off()
    
    # Checks:
    if(chk.prev) check.prevalence(d.prev)
    if(chk.cumInc) check.cumInc(d.cumInc)
  }
  print(" - DONE! -")
}

### PLOT DIFFERENCES B/W PARAMETER SETS

if(saveplot & n.prm.set>1){
  pdf("plot_all_GI_prmset.pdf")
  ggplot(GI.ODE.all)+geom_line(aes(x=time,y=GI.fwd.mean,colour=factor(param.set)),size=3)
  ggplot(GI.ODE.all)+geom_line(aes(x=time,y=GI.bck.mean,colour=factor(param.set)),size=3)
  dev.off()
}