require(snowfall)
require(data.table)
library(reshape2)

read.mc.files.parallel <- function(filerootname,file.paramset)
{
  ### Function to merge all MC files into one data frame
  
  ### Find the successfull executions only
  lsname  = paste0(path.model,"OUT/",file.paramset,"__",filerootname,"*.out")
  fname <- system(paste0("ls ", lsname ),intern = T)
  mc.success <- length(fname)
  print(paste("Reading",mc.success,"MC",filerootname, "output files (parallel)..."))
  
  ### Parallel execution:
  n.cpus <- sfCpus()
  sfInit(parallel = T, cpu = n.cpus)
  
  read.one.file <- function(dummy_index,flist) {
    message(paste("DEBUG",dummy_index))
    tmp <- read.csv(file = flist[dummy_index], header = F)
    tmp$mc <- dummy_index
    return(tmp)
  }
  sfExportAll()
  idx <- c(1:mc.success)
  res <- sfLapply(idx,fun=read.one.file, flist=fname)
  sfStop()
  
  d <- as.data.frame(rbindlist(res))
  return(list(dat=d, mc.success=mc.success))
}

read.mc.files.serial <- function(filerootname,file.paramset,path.model)
{
  ### Function to merge all MC files into one data frame
  
  ### Find the successfull executions only
  lsname  = paste0(path.model,"OUT/",file.paramset,"__",filerootname,"*.out")
  fname <- system(paste0("ls ", lsname ),intern = T)
  mc.success <- length(fname)
  print("---")
  print(paste("Reading",mc.success,"MC",filerootname, "output files (serial)..."))
  
  read.one.file <- function(i,flist) {
    cat(paste0(i,"."))
    tmp <- read.csv(file = flist[i], header = F)
    tmp$mc <- i
    return(tmp)
  }
  idx <- c(1:mc.success)
  res <- lapply(idx,FUN=read.one.file, flist=fname)
  
  d <- as.data.frame(rbindlist(res))
  return(list(dat=d, mc.success=mc.success))
}

read.mc.files <- function(doParallel=FALSE,filerootname,file.paramset,path.model){
  if(doParallel)  res<-read.mc.files.parallel(filerootname,file.paramset,path.model)
  if(!doParallel)  res<-read.mc.files.serial(filerootname,file.paramset,path.model)
  return(res)
}


# Check prevalence
check.prevalence <- function(d.prev){
  g.chk = ggplot(d.prev)+geom_ribbon(aes(x=time2,y=m,ymin=min,ymax=max),
                                     fill="lightgray")
  g.chk = g.chk + geom_line(aes(x=time2,y=m))
  g.chk = g.chk + geom_line(data=SEmInR,
                            aes(x=time,y=popSize*Iall),
                            colour="red",lwd=2)
  g.chk = g.chk + ggtitle(paste0("Check Prevalence Gillespie vs Deterministic ODE\n",info.prm))
  g.chk = g.chk + xlab("Time")+ylab("")
  
  pdf(file=paste0("plot_prev",fileplot),width=15,height=10)
  plot(g.chk)
  dev.off()
}


# Check cumulative incidence
check.cumInc <- function(d.cumInc){
  g.cuminc = ggplot(d.cumInc)+geom_ribbon(aes(x=time2,y=m,ymin=min,ymax=max,dummy=factor(mc)),
                                          fill="lightgray")
  g.cuminc = g.cuminc + geom_line(aes(x=time2,y=m))
  g.cuminc = g.cuminc + geom_line(data=SEmInR,aes(x=time,y=Z*popSize),colour="red",size=2)
  g.cuminc = g.cuminc + ggtitle(paste0("Check Cum Incidence Gillespie vs Deterministic ODE\n",info.prm))
  
  pdf(file=paste0("plot_cumInc",fileplot),width=15,height=10)
  plot(g.cuminc)
  dev.off()
}



get.GI.bck.sim <- function(doParallel,file.paramset,t.bucket,path.model){
  ### RETRIEVE BACKWARD GI FROM SIMULATIONS
  
  GIbck <- read.mc.files(doParallel,"GIbck",file.paramset,path.model)[["dat"]]
  names(GIbck)=c("time.infectee","GI","mc")
  GIbck <- GIbck[GIbck$GI>0,]
  GIbck$time.infectee2 = round(GIbck$time.infectee/t.bucket)*t.bucket
  return(GIbck)
}


get.GI.fwd.sim <- function(doParallel,file.paramset,t.bucket,path.model){
  ### RETRIEVE FORWARD GI FROM SIMULATIONS
  
  GIfwd <- read.mc.files(doParallel,"GIfwd",file.paramset,path.model)[["dat"]]
  names(GIfwd)[1]="time.infector"
  # Remove the initial infectious individuals (those infected at time t=0)
  # b/c they are They are put in the I[1] compartment directly at time 0 (not in E[1]), 
  # and the clock starts counting from that moment, so we miss the exposed duration! 
  GIfwd <- subset(GIfwd,time.infector>0)
  
  # DELETE?
  # GIfwd2 = melt(GIfwd,id.vars = c("time.infector","mc"),na.rm = T)
  # GIfwd2$time.infector2 = round(GIfwd2$time.infector/t.bucket)*t.bucket
  # GIfwd2$GIf <- GIfwd2$value
  return(GIfwd)
}

get.GI.fwd.sim.melt <- function(doParallel,file.paramset,t.bucket,path.model){
  ### RETRIEVE FORWARD GI FROM SIMULATIONS
  ### AND REFORMAT IN TALL SLIM FORMAT
  
  GIfwd <-get.GI.fwd.sim(doParallel,file.paramset,t.bucket,path.model)
  GIfwd2 = melt(GIfwd,id.vars = c("time.infector","mc"),na.rm = T)
  GIfwd2$time.infector2 = round(GIfwd2$time.infector/t.bucket)*t.bucket
  GIfwd2$GI <- GIfwd2$value
  return(GIfwd2)
}


get.simulation.param <- function(file.param){
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
  return(c(R0=R0,nE=nE,nI=nI))
}
