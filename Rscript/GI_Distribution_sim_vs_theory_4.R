#################################################################
###   
###   COMPARE EMPIRICAL & THEORETICAL GI DISTRIBUTIONS
###
###   Created 2015-07-03 by David Champredon
###
#################################################################

source("calc_theoretical_GI.R")
source("read_simul_GI_FCT.R")
source("figures_ms.R")

save.to.file <- TRUE
add.info.filename <- F #TRUE
info.in.title <- F #TRUE

# Path to the C++ model generating simulations
path.model <- "../Gillespie_SEmInR/"

# Read the simulation parameter values:
simprm.list <- as.character(read.csv(paste0(path.model,"param_all_list.csv"),header=F)[,1])  
file.prm <- simprm.list[1]
simprm <- read.csv(paste0(path.model,file.prm),header=F)
R0 <- simprm[simprm$V1=="R0",2]
nE <- simprm[simprm$V1=="nE",2]
nI <- simprm[simprm$V1=="nI",2]
mc <- simprm[simprm$V1=="mc_iter",2]
popSize <- simprm[simprm$V1=="popSize",2]
latent_mean <- simprm[simprm$V1=="latent_mean",2]
infectious_mean <- simprm[simprm$V1=="infectious_mean",2]
horiz <- simprm[simprm$V1=="horizon",2]
prm.info = paste0("_R0_",R0,"_nE_",nE,"_nI_",nI,
                  "_lat_",latent_mean,"_inf_",infectious_mean,
                  "_pop_",popSize/1000,"k_MC_",mc)

# File name for output plots
fname.fwd <- ifelse(add.info.filename,paste0("plot_fwd_dist",prm.info,".pdf"),"plot_fwd_dist.pdf")
fname.bck <- ifelse(add.info.filename,paste0("plot_bck_dist",prm.info,".pdf"),"plot_bck_dist.pdf")

# Mean intrinsic GI 
mean.gi.intrinsic = latent_mean + infectious_mean*(nI+1)/nI/2

# Slim data frames
t.bucket = 0.002

### Retrieve generation intervals data from simulations:
GIbck.sim <-  get.GI.bck.sim(doParallel=FALSE, file.prm,t.bucket,path.model)
GIbck.sim$t <- ceiling(GIbck.sim$time.infectee)
GIfwd.sim <-  get.GI.fwd.sim.melt(doParallel=FALSE,file.prm,t.bucket,path.model)
GIfwd.sim$t <- round(GIfwd.sim$time.infector)
max.horizon <- round(max(GIbck.sim$time.infectee))+1

### Calculate theoretical (using a SEmInR model)
### forward & backward generation intervals
theo.GI <- calc.theoretical.GI(file.prmset = paste0(path.model,file.prm), 
                               n.points.GI.crv = min(200,max.horizon),
                               horizon = 1.02*max.horizon,
                               do.plot = FALSE)
GI.ODE <- theo.GI[["GI.ODE"]]
theo.gi.fwd <- theo.GI[["GI.fwd.theo"]]
theo.gi.bck <- theo.GI[["GI.bck.theo"]]
theo.gi.fwd.time <- theo.GI[["GI.fwd.theo.time"]]
theo.gi.bck.time <- theo.GI[["GI.bck.theo.time"]]
theo.time <- theo.GI[["time.vec"]]



###############
###  PLOTS  ###
###############

# Calendar times where we look:
tsvec.fwd <- c(5,40,60)
tsvec.bck <- c(5,40,60)
plot.w <- 10
plot.h <- 5

if(save.to.file) pdf(fname.fwd,width = plot.w, height = plot.h)
layout(matrix(c(1,1,1,2,3,4),
              nrow = 2, ncol=3, byrow = TRUE), 
       widths=c(1,1,1), heights=c(1,1))

thetitle = "Mean forward GI: theory vs. simulations"
if(info.in.title) thetitle = paste(thetitle,"\n",prm.info)

plot.theo.vs.sim(dat.gil = GIfwd.sim,
                 fwdOrBck = "fwd",
                 dat.ode = GI.ODE,
                 n.times = 50, 
                 title = thetitle,
                 mean.gi.intrinsic = mean.gi.intrinsic,
                 min.mc = mc/4,
                 tsvec = tsvec.fwd)

sapply(tsvec.fwd,FUN=compare.sim.theo.distrib,
       GIfwd.sim,g,I,
       theo.gi.fwd,theo.gi.fwd.time,"Forward")
if(save.to.file) dev.off()


if(save.to.file) pdf(fname.bck,width = plot.w, height = plot.h)
layout(matrix(c(1,1,1,2,3,4),
              nrow = 2, ncol=3, byrow = TRUE), 
       widths=c(1,1,1), heights=c(1,1))

thetitle = "Mean backward GI: theory vs. simulations"
if(info.in.title) thetitle = paste(thetitle,"\n",prm.info)

plot.theo.vs.sim(dat.gil = GIbck.sim,
                 fwdOrBck = "bck",
                 dat.ode = GI.ODE,
                 n.times = 50, 
                 title = thetitle,
                 mean.gi.intrinsic = mean.gi.intrinsic,
                 min.mc = mc/4,
                 tsvec = tsvec.bck)

sapply(tsvec.bck,FUN=compare.sim.theo.distrib,
       GIbck.sim,g,I,
       theo.gi.bck,theo.gi.bck.time,"Backward")

if(save.to.file) dev.off()
