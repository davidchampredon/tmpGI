#################################################################################
###
###   GENERATE ALL FILES THAT SPECIFY
###   PARAMETERS FOR SEmInR Gillespie MODEL
###
###   Created 2015-03-16 by David Champredon
###
#################################################################################

### CLEAN-UP PREVIOUS FILES
system("rm param_all_list.csv")
system("rm param_*.csv")

### =======================
### Non-variable parameters:
### =======================

prm <- read.csv("simul_param.csv",header = T)

popSize <- prm$value[prm$parameter=="popSize"]
horizon <-  prm$value[prm$parameter=="horizon"] # <-- should be highest possible horizon across all param combinations
mc_iter <- prm$value[prm$parameter=="mc_iter"]
njobs   <- prm$value[prm$parameter=="njobs"]
nE <- prm$value[prm$parameter=="nE"]
nI <- prm$value[prm$parameter=="nI"]
init_I1 <- prm$value[prm$parameter=="init_I1"]
calc_WIW_Re <- prm$value[prm$parameter=="calc_WIW_Re"]
doExact <- prm$value[prm$parameter=="doExact"]
timeStepTauLeap <- prm$value[prm$parameter=="timeStepTauLeap"]

nv.prm <- data.frame(name=c("popSize","horizon","mc_iter","njobs",
                            "nE","nI","init_I1","calc_WIW_Re",
                            "doExact","timeStepTauLeap"),
                     value=c(popSize,horizon,mc_iter,njobs,
                             nE,nI,init_I1,calc_WIW_Re,
                             doExact,timeStepTauLeap)
                     )
nv.prm$name <- as.character(nv.prm$name)

### =======================
### Variable parameters:
### =======================

R0.vec <- read.csv("simul_param_R0.csv",header = F,comment.char = "#")[,1]
latent.mean.vec <- read.csv("simul_param_DOL.csv",header = F,comment.char = "#")[,1]
infect.mean.vec <- read.csv("simul_param_DOI.csv",header = F,comment.char = "#")[,1]
n.R0<-length(R0.vec)
n.GI<-length(latent.mean.vec)
stopifnot(length(latent.mean.vec)==length(infect.mean.vec))

message(paste0("Generating ",n.R0*n.GI," parameter files..."))

filelist <- "param_all_list.csv"
filename <- vector()

cnt <- 1
for(r in 1:n.R0){
  for(g in 1:n.GI){
    filename[cnt] = paste0("param_",cnt,".csv")
    
    df.prm <- rbind(nv.prm,c("R0",R0.vec[r]))
    df.prm <- rbind(df.prm,c("latent_mean",latent.mean.vec[g]))
    df.prm <- rbind(df.prm,c("infectious_mean",infect.mean.vec[g]))
    
    write.table(x=df.prm,file=filename[cnt],sep=",",
                row.names=F,col.names=F,quote=F)
    message(paste0("File ",filename[cnt]," written."))
    cnt <- cnt+1
  }
}

write.table(x=filename,file=filelist,row.names=F,col.names=F,quote=F)


