####################################################
###   PLOT BACKWARD AND FORWARD 
###   INTERVAL MEAN AS A FUNCTION OF R0
###
###   Created 2015-05-19 by David Champredon
###
####################################################


source("SEmInR_deterministic.R")
source("calc_theoretical_GI.R")
source("GIfwdbck_multi_R0_FCT.R")

### Global Variables
save.to.file <- TRUE

# Path to the C++ model generating simulations
path.model <- "../Gillespie_SEmInR/"
simprm.list <- as.character(read.csv(paste0(path.model,"param_all_list.csv"),header=F)[,1])  
file.prm <- simprm.list[1]

# Retrieve SEmInR model parameters
# (but not R0!!!!)
prm <- read.csv(paste0(path.model,file.prm),header=FALSE)
max.horizon <- prm[prm$V1=="horizon",2]
nE <- prm[prm$V1=="nE",2]
nI <- prm[prm$V1=="nI",2]
latent_mean <- prm[prm$V1=="latent_mean",2]
infectious_mean <- prm[prm$V1=="infectious_mean",2]
N <- prm[prm$V1=="popSize",2] ; popSize <- N
I.init <- prm[prm$V1=="init_I1",2]/N

### R0 values that will be on the plot
R0.values <- c(1.5, 2.0, 4.0, 8.0, 16.0)
n.R0 <- length(R0.values)

### Calculate fwd & bck GI for all R0 values:
GI <- list()
for(i in 1:n.R0){
  print(paste("R0 =",R0.values[i]))
  GI[[i]] <- calc.GI.fwdbck.theo(R0=R0.values[i], max.horizon, nE, nI, 
                                latent_mean, infectious_mean, 
                                N, I.init)
}
# Mean intrinsic GI
mean.igi <- GI[[1]][["mean.g"]]


###  ==== PLOTS ====

if(save.to.file) pdf("GIfwdbck_multi_R0.pdf",width=12,height=6)
par(cex.main=1.7, cex.lab=1.5, 
    cex.axis=1.5, las = 1)

par(mfrow=c(1,2))
lwd.gi = 4
pch.ep <- 16
cex.ep <- 1.5
pos.lab.bck <- c(3,4,4,4,3)
pos.lab.fwd <- 3 
col.line <- rgb(0,0,0,0.7)

### Backward GI:
plot(x=NULL,y=NULL,
     lwd=lwd.gi,
     main = "Mean Backward Generation Interval",
     xlim=c(0,GI[[1]][["horizon"]]), 
     ylim= c(0,3.1*mean.igi),
     xlab="calendar time",
     ylab="Days",
     typ="l", 
     las=1)
grid()
abline(h=mean.igi,lty=5,lwd=lwd.gi/2)

for(i in 1:n.R0){
  lines(x=GI[[i]][["b.bar.time"]], y=GI[[i]][["b.bar"]],
        lwd=lwd.gi, col=col.line)
  nn = length(GI[[i]][["b.bar.time"]])
  xi <- GI[[i]][["b.bar.time"]][[nn]]
  yi <- GI[[i]][["b.bar"]][[nn]]
  points(x=xi, y=yi, pch=pch.ep,cex=cex.ep)
  text(x=xi,y=yi,labels = R0.values[i],col="grey",pos = pos.lab.bck)
}


### Forward GI:
plot(x=NULL,y=NULL,
     lwd=lwd.gi ,
     main = "Mean Forward Generation Interval",
     xlim=c(0,GI[[1]][["horizon"]]), 
     ylim = c(0,1.1*mean.igi),
     xlab="calendar time",
     ylab="Days",
     typ="l")
grid()
abline(h=mean.igi,lty=5,lwd=lwd.gi/2)

for(i in 1:n.R0){
  nn = length(GI[[i]][["f.bar.time"]])
  lines(x=GI[[i]][["f.bar.time"]], y=GI[[i]][["f.bar"]],
        lwd=lwd.gi, col=col.line)
  xi <- GI[[i]][["f.bar.time"]][[nn]]
    yi <- GI[[i]][["f.bar"]][[nn]]
  points(x=xi, y=yi, pch=pch.ep,cex=cex.ep)
  text(x=xi,y=yi,labels = R0.values[i],col="grey",pos = pos.lab.fwd)
}

if(save.to.file) dev.off()
