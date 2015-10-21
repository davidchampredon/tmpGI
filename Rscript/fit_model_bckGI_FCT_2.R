


distance.interpol <- function(x1,y1, 
							  x2,y2, 
							  deg.poly,
							  max.x,
							  do.plot = F){
	### CALCULATE DISTANCE OF POLYNOMIAL INTERPOLATION 
	### OF 2 FUNCTIONS
	### THAT DO NOT HAVE SAME Xs
	### (can't easily do "f(x)-g(x)" because x do not match)
	
	fit1 <- lm(y1~poly(x1,deg.poly))
	interpol1 <- predict(fit1, data.frame(x1=x1))
	
	fit2 <- lm(y2~poly(x2,deg.poly))
	interpol2 <- predict(fit2, data.frame(x2=x1))
	
	if(do.plot){ # plots to debug
		plot(x1,y1)
		lines(x1,interpol1,typ="l")
		points(x2,y2,col="red",pch=4)
		lines(x1,interpol2,typ="l",col="red")
		abline(v=max.x,lty=2,lwd=3)
	}
	idx.max = which(x1 <= max.x)
	dist <- sqrt(sum((interpol1[idx.max]-interpol2[idx.max])^2))
	return(dist)
}



fct.to.minimize <- function(GIbck.sim.melt,
							latent_mean ,
							infectious_mean,
							R0,
							nE,
							nI,
							max.x){
	### RETURNS THE DISTANCE B/W THE OBSERVED 
	### BACKWARD GI AND THE ONE IMPLIED BY THE MODEL
	### DISTANCE = SUM OF DISTANCE OF INTERPOLATED CURVES OF MEAN AND VARIANCE 
	### BCK GI AT ALL CALENDAR TIMES < max.x
	
	theo.GI.new <- calc.theoretical.GI.base(N=popSize,
											latent_mean = latent_mean,
											infectious_mean = infectious_mean,
											R0 = R0,
											init_I1 = 7,
											nE = nE,
											nI = nI,
											n.points.GI.crv = min(200,max.horizon),
											horizon = 1.02*max.horizon,
											do.plot = FALSE)
	GI.ODE.new <- theo.GI.new[["GI.ODE"]]

	do.plot <- F
	par(mfrow=c(1,2))
	err.mean <- distance.interpol(x1 = GIbck.sim.melt$timeround,
								  y1 = GIbck.sim.melt$m, 
								  x2 = GI.ODE.new$time,
								  y2 = GI.ODE.new$GI.bck.mean,
								  deg.poly=9,
								  max.x=max.x,
								  do.plot=do.plot)
	err.var <- distance.interpol(x1 = GIbck.sim.melt$timeround,
								 y1 = GIbck.sim.melt$sd^2, 
								 x2 = GI.ODE.new$time,
								 y2 = GI.ODE.new$GI.bck.var,
								 deg.poly=9,
								 max.x=max.x,
								 do.plot=do.plot)
	return(err.mean+err.var)
}





fct.to.minimize.intrinsic <- function(GIbck.sim.melt,
							latent_mean ,
							infectious_mean,
							R0,
							nE,
							nI,
							max.x){
	### RETURNS THE DISTANCE B/W THE OBSERVED 
	### BACKWARD GI AND THE ONE IMPLIED BY THE MODEL
	### DISTANCE = SUM OF DISTANCE OF INTERPOLATED CURVES OF MEAN AND VARIANCE 
	### BCK GI AT ALL CALENDAR TIMES < max.x
	
	theo.GI.new <- calc.theoretical.GI.base(N=popSize,
											latent_mean = latent_mean,
											infectious_mean = infectious_mean,
											R0 = R0,
											init_I1 = 7,
											nE = nE,
											nI = nI,
											n.points.GI.crv = min(200,max.horizon),
											horizon = 1.02*max.horizon,
											do.plot = FALSE)
	#GI.ODE.new <- theo.GI.new[["GI.ODE"]]
	
	# Retrieve intrinsic GI (g) from this simulation:
	g.sim  <- theo.GI.new[["GI.intrinsic"]]
	g.m  <- mean(g.sim)
	g.sd  <- sd(g.sim)
	
	# Calculate mean and sd from the target data:
	target.m <- mean(GIbck.sim.melt$m)
	target.sd <- sd(GIbck.sim.melt$m)
	
	
	err.mean <- sqrt((g.m-target.m)^2)
	err.var <- sqrt((g.sd-target.sd)^2)
	
	return(err.mean+err.var)
}




ABC.trials.unit <- function(i, GIbck.sim.melt,latent_mean,infectious_mean,
							R0,nE,nI, max.x){
	### APPROXIMATE BAYESIAN CALCULATION
	### (UNIT FUNCTION TO BE RUN IN PARALLEL)
	dist = fct.to.minimize(GIbck.sim.melt,
						   latent_mean[i] ,
						   infectious_mean[i],
						   R0[i],
						   nE[i],
						   nI[i], 
						   max.x)
	return(c(dist=dist,
			 latent_mean=latent_mean[i], 
			 infectious_mean=infectious_mean[i], 
			 R0=R0[i], nE=nE[i],nI=nI[i]) )
}


ABC.trials.unit.intrinsic <- function(i, GIbck.sim.melt,latent_mean,infectious_mean,
							R0,nE,nI, max.x){
	### APPROXIMATE BAYESIAN CALCULATION
	### (UNIT FUNCTION TO BE RUN IN PARALLEL)
	dist = fct.to.minimize.intrinsic(GIbck.sim.melt,
									 latent_mean[i] ,
									 infectious_mean[i],
									 R0[i],
									 nE[i],
									 nI[i], 
									 max.x)
	return(c(dist=dist,
			 latent_mean=latent_mean[i], 
			 infectious_mean=infectious_mean[i], 
			 R0=R0[i], nE=nE[i],nI=nI[i]) )
}



ABC.posterior <- function(abc, thresh){
	### ANALYSIS OF ABC TRIALS BELOW ACCEPTANCE THRESHOLD
	
	abc.keep <- subset(abc,dist<thresh)
	plot(sort(abc.keep$dist),typ="s")
	
	print("Minimum distance parameters:")
	print(abc.keep[which.min(abc.keep$dist),])
	
	library(reshape2)
	library(ggplot2)
	library(plyr)
	df <- melt(abc.keep,id.vars = "dist")
	
	g<- ggplot(df)+geom_histogram(aes(x=value),binwidth=0.5)+facet_wrap(~variable,scales="free")
	plot(g)
	
	df2 <- ddply(df,.variables = "variable",summarize,m=mean(value),s=sd(value),n=length(value))
	return(df2)
}


ABC.trials<- function(GIbck.sim.melt, iter = 100){
	### APPROXIMATE BAYESIAN CALCULATION
	### (NOT PARALLELIZED)
	for(i in 1:iter){
		print(paste("\n ABC iter ", i))
		
		latent_mean <- max(0.1, rnorm(n = 1, mean = 5, sd=2))
		infectious_mean <- max(0.1, rnorm(n = 1, mean = 5, sd=2))
		R0 <- max(0.1,rnorm(n=1, mean = 2, sd=0.5))
		nE <- round(runif(n=1, min = 1, max=12))
		nI <- round(runif(n=1, min = 1, max=12))
		
		dist = fct.to.minimize(GIbck.sim.melt,
							   latent_mean ,
							   infectious_mean,
							   R0,
							   nE,
							   nI)
		
		if(i==1) df = data.frame(dist=dist,latent_mean=latent_mean, infectious_mean=infectious_mean, R0=R0, nE=nE,nI=nI)
		if(i>1) df = rbind(df,c(dist=dist,latent_mean=latent_mean, infectious_mean=infectious_mean, R0=R0, nE=nE,nI=nI))
	}
	return(df)
}


plot.fit.results <- function(GI.ODE.truth,
							 GIbck.sim.melt,
							 deg.poly,
							 popSize,
							 init_I1,
							 latent_mean ,
							 infectious_mean ,
							 R0 ,
							 nE,
							 nI,
							 max.x
){
	
	# Calculate theoretical generation interval:
	
	theo.GI.fit <- calc.theoretical.GI.base(N=popSize,
											latent_mean=latent_mean,
											infectious_mean=infectious_mean ,
											R0=R0,
											init_I1=init_I1 ,
											nE=nE ,
											nI=nI ,
											n.points.GI.crv = min(200,max.horizon),
											horizon = 1.02*max.horizon,
											do.plot = FALSE)
	GI.ODE.fit <- theo.GI.fit[["GI.ODE"]]
	
	do.plot <- FALSE
	err.mean <- distance.interpol(x1 = GIbck.sim.melt$timeround,
								  y1 = GIbck.sim.melt$m, 
								  x2 = GI.ODE.fit$time,
								  y2 = GI.ODE.fit$GI.bck.mean,
								  deg.poly=9,
								  max.x=max.x,
								  do.plot=do.plot)
	err.var <- distance.interpol(x1 = GIbck.sim.melt$timeround,
								 y1 = GIbck.sim.melt$sd^2, 
								 x2 = GI.ODE.fit$time,
								 y2 = GI.ODE.fit$GI.bck.var,
								 deg.poly=9,
								 max.x=max.x,
								 do.plot=do.plot)
	
	par(mfrow=c(1,2))
	yrng = c(0,1.1*max(GIbck.sim.melt$m,GI.ODE.fit$GI.bck.mean,na.rm = T))
	plot(x=GI.ODE.truth$time, y=GI.ODE.truth$GI.bck.mean, 
		 typ="l",lwd=1,
		 col="grey",
		 ylim=yrng,
		 main="Backward GI Mean",
		 xlab = "Calendar time",
		 ylab = "Mean of backward GI")
	# abline(a = 0,b=1, lty=2, col="grey")
	lines(x=GI.ODE.fit$time, GI.ODE.fit$GI.bck.mean,col="red",lwd=6)
	points(x=GIbck.sim.melt$timeround, y = GIbck.sim.melt$m,pch=16)
	abline(v=max.x,lty=2)
	grid()
	
	plot(x=GIbck.sim.melt$timeround, y = GIbck.sim.melt$sd^2, 
		 main="Backward GI Variance",
		 xlab = "Calendar time",
		 ylab = "Variance of backward GI")
	lines(x=GI.ODE.fit$time,y=GI.ODE.fit$GI.bck.var,col="red",lwd=6)
	abline(v=max.x,lty=2)
	grid()
}




plot.fit.results.new <- function(GI.ODE.truth,
							 GIbck.sim.melt,
							 deg.poly,
							 popSize,
							 init_I1,
							 prm.fit,
							 prm.fit.naive,
# 							 latent_mean ,
# 							 infectious_mean ,
# 							 R0 ,
# 							 nE,
# 							 nI,
							 max.x
){
	
	# unpack parameters:
	latent_mean = prm.fit[["latent_mean"]]
	infectious_mean = prm.fit[["infectious_mean"]]
	R0 = prm.fit[["R0"]]
	nE = prm.fit[["nE"]]
	nI = prm.fit[["nI"]]
	
	latent_mean.naive = prm.fit.naive[["latent_mean"]]
	infectious_mean.naive = prm.fit.naive[["infectious_mean"]]
	R0.naive = prm.fit.naive[["R0"]]
	nE.naive = prm.fit.naive[["nE"]]
	nI.naive = prm.fit.naive[["nI"]]
	
	
	# Calculate theoretical generation interval
	# with fitted parameters
	theo.GI.fit <- calc.theoretical.GI.base(N=popSize,
											latent_mean=latent_mean,
											infectious_mean=infectious_mean ,
											R0=R0,
											init_I1=init_I1 ,
											nE=nE ,
											nI=nI ,
											n.points.GI.crv = min(200,max.horizon),
											horizon = 1.02*max.horizon,
											do.plot = FALSE)
	GI.ODE.fit <- theo.GI.fit[["GI.ODE"]]
	
	# Calculate theoretical generation interval
	# with NAIVELY fitted parameters
	theo.GI.fit.naive <- calc.theoretical.GI.base(N=popSize,
											latent_mean=latent_mean.naive,
											infectious_mean=infectious_mean.naive ,
											R0=R0.naive,
											init_I1=init_I1 ,
											nE=nE.naive ,
											nI=nI.naive ,
											n.points.GI.crv = min(200,max.horizon),
											horizon = 1.02*max.horizon,
											do.plot = FALSE)
	GI.ODE.fit.naive <- theo.GI.fit.naive[["GI.ODE"]]
	
	
	par(mfrow=c(1,2))
	yrng = c(0,1.1*max(GIbck.sim.melt$m,GI.ODE.fit$GI.bck.mean,na.rm = T))
	
	### Mean GI
	###
	plot(x=GI.ODE.truth$time, 
		 y=GI.ODE.truth$GI.bck.mean, 
		 typ="l",lwd=1,
		 col="grey",
		 ylim=yrng,
		 main="Backward GI Mean",
		 xlab = "Calendar time",
		 ylab = "Mean of backward GI")
	
	lines(x=GI.ODE.fit$time, GI.ODE.fit$GI.bck.mean,col="red",lwd=6)
	lines(x=GI.ODE.fit.naive$time, GI.ODE.fit.naive$GI.bck.mean,col="red",lwd=6, lty=2)
	
	points(x=GIbck.sim.melt$timeround, y = GIbck.sim.melt$m,pch=16)
	abline(v=max.x,lty=2)
	grid()
	
	### Std dev of GI
	###
	plot(x=GIbck.sim.melt$timeround, 
		 y = GIbck.sim.melt$sd^2, 
		 main="Backward GI Variance",
		 xlab = "Calendar time",
		 ylab = "Variance of backward GI")
	
	lines(x=GI.ODE.fit$time,y=GI.ODE.fit$GI.bck.var,col="red",lwd=6)
	lines(x=GI.ODE.fit.naive$time,y=GI.ODE.fit.naive$GI.bck.var,col="red",lwd=6, lty=2)
	
	abline(v=max.x,lty=2)
	grid()
}



fct.to.minimize.wrap <- function(x,GIbck.sim.melt){
	### Wrap for 'optim()'
	res<-fct.to.minimize(GIbck.sim.melt,
						 latent_mean=x[1] ,
						 infectious_mean=x[2],
						 R0=x[3],
						 nE=x[4],
						 nI=x[5])
	write.table(x = t(c(x,res)),file="optim.csv",append = TRUE, 
				quote=F,sep = ",",col.names = F,row.names = F)
	return(res)
}