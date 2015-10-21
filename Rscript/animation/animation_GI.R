library("animation")

source("calc_theoretical_GI.R")
source("SEmInR_deterministic.R")
source("animation_GI_FCT.R")

t0<-Sys.time()

moviename <- "movie_GI.gif"
system(paste0("rm -f ",moviename))

# Path to the C++ model generating simulations
path.model <- "./"
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
R0 <-  prm[prm$V1=="R0",2]
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

#Set delay between frames when replaying
ani.options(interval=0.15)
#Generate files under the current working directory.
ani.options(outdir = getwd())
ani.options(ani.width=1000, ani.height=1100)
t.min=10
t.max=horizon - 5
nbFrames = (t.max-t.min)

times.movie <- round(seq(t.min,t.max,length.out=nbFrames))

saveGIF(expr = {
	for(t in times.movie) plot.one.frame(t)
},movie.name = moviename)

message(paste("time elapsed:",Sys.time()-t0))