#################################################################
#####     
#####     FIT SEmInR MODEL PARAMETERS (R0,DOL,DOI,nE,nI)
#####     TO BCKWD GENERATION INTERVAL DATA
#####
#####     WARNING: This code is supposed to use simulations
#####     that are typically different from other R codes of 
#####     this project! 
#####     Here, simulation with one MC iteration is needed 
#####     (i.e., fit to one realization of data, as we'd do in real life)
#####     
#####     The simulation must have been run beforehand!	
#####
#####     Created: 2015-08-02 by David Champredon
#####     Modified: 
#####     v2: 2015-10-21 by David Champredon to include reviewers' comments
#####
#################################################################


source("calc_theoretical_GI.R")
source("read_simul_GI_FCT.R")
source("figures_ms.R")
source("fit_model_bckGI_FCT_2.R")

fitrdata <- "fit.RData"
system(paste("rm -rf",fitrdata),intern=FALSE)

### Retrieve fit parameters:
###
fitprm <- read.csv("fit_model_bckGI_param.csv")
# ABC iterations:
iter.abc <- fitprm$value[fitprm$name=="iter.abc"]
# Number of CPUs used in parallel:
ncpus <- fitprm$value[fitprm$name=="ncpus"]
# Number of dates to fit to
max.x <- fitprm$value[fitprm$name=="max.x"]

save.to.file <- TRUE


### ==== Retrieve simulation data (target to fit) ====
  
# Read the simulation parameter values:
path.model <- as.character(read.csv(file="path_model.txt",header = F)[1,1])
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

# Slim data frames
t.bucket = 0.002

### Retrieve generation intervals data 
### from (previously run) simulations:
###
GIbck.sim <-  get.GI.bck.sim(doParallel=FALSE, file.prm,t.bucket,path.model)

# Force to take only one MC realization
# (in the case the model previously run had more than one MC)
GIbck.sim <- subset(GIbck.sim, mc==1)

GIbck.sim$t <- ceiling(GIbck.sim$time.infectee)
max.horizon <- round(max(GIbck.sim$time.infectee))+1
tlab <- "time.infectee"  # <-- legacy from previous code where could switch b/w fwd andbck GI (keep it in case)

# Define calendar time vector
tmax <- max(GIbck.sim[tlab])
n.times <- 100
tt <- seq(0,tmax, length.out=n.times+2)	
dt <- tt[2]-tt[1]


# merge time points into time buckets for plots 
#(designed when simulation's MC>1, but here MC should be =1, 
# so kind of useless but easier to keep it - it's harmless)
GIbck.sim$timeround = round(GIbck.sim[tlab][,1]/dt)*dt
GIbck.sim.melt = ddply(GIbck.sim,c("timeround"),summarize,
                       m = mean(GI),
                       sd = sd(GI),
                       n = length(GI))
# 'True' model bck GI implied by 
# 'true' value parameters (the one used in simulations)
theo.GI.truth <- calc.theoretical.GI(file.prmset = paste0(path.model,file.prm), 
                               n.points.GI.crv = min(200,max.horizon),
                               horizon = 1.02*max.horizon,
                               do.plot = FALSE)

GI.ODE.truth <- theo.GI.truth[["GI.ODE"]]


### ==== PARAMETER FIT FOR BACKWARD GI ====
### 
###  Using Approximate Bayesian Computation 
###

# ABC priors:
latent_mean <- runif(n=iter.abc,min = 1, max=15)  
infectious_mean <- runif(n=iter.abc,min = 1, max=15)  
R0 <- pmax(rep(x=0.1,times=iter.abc),rnorm(n=iter.abc, mean = 2, sd=0.6))
nE <- round(runif(n=iter.abc, min = 1, max=12))
nI <- round(runif(n=iter.abc, min = 1, max=12))

### ABC trials:
###
require(snowfall)
sfInit(parallel = TRUE, cpu = ncpus)
sfLibrary(deSolve)
idx.apply <- 1:iter.abc
sfExportAll()
system.time(res <- sfSapply(idx.apply, ABC.trials.unit, 
							GIbck.sim.melt,
							latent_mean,
							infectious_mean,R0,nE,nI,max.x))
sfStop()

# Store all trials in a data frame:
abc.trials <- data.frame(t(res))

# Posterior
distance.threshold <- 2*min(abc.trials$dist) 
abc.post <- ABC.posterior(abc.trials, thresh=distance.threshold)


### NAIVE FIT OF INTERINSIC GI
###
require(snowfall)
sfInit(parallel = TRUE, cpu = ncpus)
sfLibrary(deSolve)
idx.apply <- 1:iter.abc
sfExportAll()
system.time(res.intrinsic <- sfSapply(idx.apply, 
							ABC.trials.unit.intrinsic, 
							GIbck.sim.melt,
							latent_mean,
							infectious_mean,R0,nE,nI,max.x))
sfStop()

# Store all trials in a data frame:
abc.trials.intrinsic <- data.frame(t(res.intrinsic))

distance.threshold.intrinsc <- 2*min(abc.trials.intrinsic$dist) 
abc.post.intrinsic <- ABC.posterior(abc.trials.intrinsic, thresh=distance.threshold)


######################
### PLOT THE FITS  ###
######################

system("rm -rf fit_model_bckGI*.pdf",intern = FALSE)

# Fitting on BACKWARD GI,
# Mean values from ABC posterior distributions:
latent_mean.mean = abc.post$m[abc.post$variable=="latent_mean"]
infectious_mean.mean =  abc.post$m[abc.post$variable=="infectious_mean"]
R0.mean = abc.post$m[abc.post$variable=="R0"]
nE.mean = round(abc.post$m[abc.post$variable=="nE"])
nI.mean = round(abc.post$m[abc.post$variable=="nI"])

# Fitting on BACKWARD GI,
# Value from best (=shortest distance) ABC trial
idx.best <- which.min(abc.trials$dist)
latent_mean.best = abc.trials$latent_mean[idx.best]
infectious_mean.best =  abc.trials$infectious_mean[idx.best]
R0.best = abc.trials$R0[idx.best]
nE.best = abc.trials$nE[idx.best]
nI.best = abc.trials$nI[idx.best]
print(abc.trials[idx.best,])


# Fitting on INTRINSIC GI,
# Mean values from ABC posterior distributions:
latent_mean.mean.int = abc.post.intrinsic$m[abc.post.intrinsic$variable=="latent_mean"]
infectious_mean.mean.int =  abc.post.intrinsic$m[abc.post.intrinsic$variable=="infectious_mean"]
R0.mean.int = abc.post.intrinsic$m[abc.post.intrinsic$variable=="R0"]
nE.mean.int = round(abc.post.intrinsic$m[abc.post.intrinsic$variable=="nE"])
nI.mean.int = round(abc.post.intrinsic$m[abc.post.intrinsic$variable=="nI"])


# plotting fit using mean of posterior for fitted param values:

# DELETE:
# plot.fit.results(GI.ODE.truth,
#                  GIbck.sim.melt,
#                  deg.poly=9,
#                  popSize, 
#                  init_I1=7,
#                  latent_mean = latent_mean.mean,
#                  infectious_mean = infectious_mean.mean,
#                  R0 = R0.mean,
#                  nE = nE.mean,
#                  nI = nI.mean,
#                  max.x = max.x
# )

prm.fit = list( latent_mean = latent_mean.mean,
				infectious_mean = infectious_mean.mean,
				R0 = R0.mean,
				nE = nE.mean,
				nI = nI.mean)
prm.fit.naive = list( latent_mean = latent_mean.mean.int,
				infectious_mean = infectious_mean.mean.int,
				R0 = R0.mean.int,
				nE = nE.mean.int,
				nI = nI.mean.int)

if(save.to.file) pdf("fit_model_bckGI.pdf",width=15, height=10)

plot.fit.results.new(GI.ODE.truth,
				 GIbck.sim.melt,
				 deg.poly=9,
				 popSize, 
				 init_I1=7,
				prm.fit=prm.fit,
				prm.fit.naive = prm.fit.naive,
				 max.x = max.x
)


# plotting fit using fitted param values that minimizes the error:
# plot.fit.results(GI.ODE.truth,
#                  GIbck.sim.melt,
#                  deg.poly=9,
#                  popSize, 
#                  init_I1=7,
#                  latent_mean = latent_mean.best,
#                  infectious_mean = infectious_mean.best,
#                  R0 = R0.best,
#                  nE = nE.best,
#                  nI = nI.best,
#                  max.x = max.x
# )

if(save.to.file) dev.off()



# Prior coverage:
if(save.to.file) pdf("fit_model_bckGI_info.pdf",width=15, height=10)
plot(abc.trials)
if(save.to.file) dev.off()

save.image(fitrdata)

