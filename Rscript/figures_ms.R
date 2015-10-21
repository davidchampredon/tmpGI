###
### PRODUCES FIGURES FOR MANUSCRIPT
###

require(ggplot2)
require(reshape2)
require(plyr)

plot.theo.vs.sim <- function(dat.gil,
							 fwdOrBck,
							 dat.ode,
							 n.times,
							 title,
							 mean.gi.intrinsic,
							 min.mc,
							 tsvec)
{
	if(fwdOrBck=="fwd") timelabel <- "time.infector"
	if(fwdOrBck=="bck") timelabel <- "time.infectee"
	
	# Define calendar time vector
	tmax <- max(dat.gil[timelabel])
	tt <- seq(0,tmax, length.out=n.times+2)	
	dt <- tt[2]-tt[1]
	
	# merge time points into time buckets for plots
	dat.gil$timeround = round(dat.gil[timelabel][,1]/dt)*dt
	
	# Summary statistics (=mean, for now) for each time buckets
	d.dat.gil.melt = ddply(dat.gil,c("timeround"),summarize,
						   m = mean(GI),
						   n = length(GI))
	
	# filter out the time points that do not have
	# enough data points (in term of MC realizations)
	df.tmp <- ddply(dat.gil,c("timeround"),summarize,
					n.mc = length(unique(mc)))
	time.enough.mc <- df.tmp$timeround[df.tmp$n.mc>=min.mc] 
	d.dat.gil.melt <- d.dat.gil.melt[d.dat.gil.melt$timeround%in%time.enough.mc,]
	dat.ode <- dat.ode[dat.ode$time<=max(time.enough.mc),]
	
	# Plots
	
	if(fwdOrBck=="fwd") gi.ode = dat.ode$GI.fwd.mean
	if(fwdOrBck=="bck") gi.ode = dat.ode$GI.bck.mean
	
	plot.rng <- range(mean.gi.intrinsic, gi.ode,d.dat.gil.melt$m)
	
	plot(x=dat.ode$time, y=gi.ode, 
		 typ = "l", lwd =3,
		 ylab="", xlab="Calendar time",
		 ylim = c(0.9*plot.rng[1],plot.rng[2]),
		 las = 1, 
		 main = title)
	grid()
	points(x=d.dat.gil.melt$timeround, y=d.dat.gil.melt$m)
	
	abline(h=mean.gi.intrinsic, lty=2)
	
	tmp = floor(dat.ode$time)
	
	# DEBUG
	print(dat.ode$time)
	
	idx.dist <- which(tmp%in%tsvec)
	points(x=tsvec, y=gi.ode[idx.dist],pch=15, cex=1.5)
	points(x=tsvec, y=gi.ode[idx.dist],pch=0, cex=2.5)
}


compare.sim.theo.distrib <- function(ts,sim.gi,g,I,
									 theo.gi, theo.gi.time, 
									 gitype){
	
	# Simulated GI distribution at given calendar time
	sgi = sim.gi$GI[sim.gi$t==ts]
	# Theoretical GI distribution at given calendar time
	ts2 <- ifelse(gitype=="Forward",ts,ts+1)  # <-- the way the backward gi time is calculated is shifted by 1, so re-align it...
	tgi <- theo.gi[[ts2]]   
	tt <- theo.gi.time[[ts2]]   
	dt = tt[2]-tt[1]
	# means
	m.tgi = sum(tgi*tt)*dt  #m.tgi = sum(tgi*theo.gi.time[[ts]])*dt
	m.sgi = mean(sgi)
	
	# Compare historam of simualted GI
	# with theoretical curves
	tmp <- hist(sgi,plot=F)
	hist(sgi,
		 breaks = c(0:(max(tmp$breaks)+1)),
		 probability = T,
		 main=paste0(gitype," GI empirical and \n theoretical distributions \n at calendar time t=",ts),
		 xlab = "Time since infection (days)", ylab="",
		 ylim = c(0,max(tgi,tmp$density)*1.1),
		 yaxt="n",
		 col="lightgrey",border = "grey")
	lines(x = tt, y=tgi,
		  typ="l", lwd=4)
	abline(v=c(m.tgi,m.sgi),lwd=1,
		   lty=c(1,2),col=c("black","black"))
}



#####################################################################
#####################################################################

plot.gilles.vs.analytic <- function(dat.gil,
									fwdOrBck,
									dat.ode,
									n.times,
									myq,
									title,
									infoprm,
									mean.gi.intrinsic,
									show.quantile,
									min.mc)
{
	if(fwdOrBck=="fwd") timelabel <- "time.infector"
	if(fwdOrBck=="bck") timelabel <- "time.infectee"
	
	# Define calendar time vector
	tmax <- max(dat.gil[timelabel])
	tt <- seq(0,tmax, length.out=n.times+2)	
	dt <- tt[2]-tt[1]
	
	# reformat in one tall, simple dataframe when fwd case
	dat.gil.melt = dat.gil
	if(fwdOrBck=="fwd"){
		dat.gil.melt = melt(dat.gil,
							id.vars = c(timelabel,"mc"),
							na.rm = T)
		names(dat.gil.melt)[names(dat.gil.melt)=="value"]<-"GI"
	}
	
	# merge time points into time buckets for plots
	dat.gil.melt$timeround = round(dat.gil.melt[timelabel][,1]/dt)*dt
	
	# Summary statistics for each time buckets
	d.dat.gil.melt = ddply(dat.gil.melt,c("timeround"),summarize,
						   min = quantile(GI,probs = (1-myq)/2),
						   m = mean(GI),
						   max = quantile(GI,probs=myq+(1-myq)/2),
						   n = length(GI))
	
	# filter out the time points that do not have
	# enough data points (in term of MC realizations)
	df.tmp <- ddply(dat.gil.melt,c("timeround"),summarize,
					n.mc = length(unique(mc)))
	time.enough.mc <- df.tmp$timeround[df.tmp$n.mc>=min.mc] 
	d.dat.gil.melt <- d.dat.gil.melt[d.dat.gil.melt$timeround%in%time.enough.mc,]
	dat.ode <- dat.ode[dat.ode$time<=max(time.enough.mc),]
	
	# Plots
	
	if(fwdOrBck=="fwd") gi.ode = "GI.fwd.mean"
	if(fwdOrBck=="bck") gi.ode = "GI.bck.mean"
	
	g = ggplot() + geom_line(data = dat.ode, 
							 aes_string(x="time",y=gi.ode),
							 size=4,
							 alpha=0.5)
	if(show.quantile){
		g = g + geom_pointrange(data=d.dat.gil.melt,
								aes(x=timeround,
									y = m,
									ymin=min,
									ymax=max),
								colour="gray",
								shape=22,
								fill="black",
								size=1)
	}
	if(!show.quantile){
		g = g + geom_point(data=d.dat.gil.melt,
						   aes(x = timeround,
						   	y = m),
						   colour="gray",
						   shape=22,
						   fill="black",
						   size=5)
	}
	
	g = g + geom_hline(yintercept=mean.gi.intrinsic, linetype=2)
	g = g + ylab("")+ggtitle(paste(title,"\n",infoprm))
	g = g + theme(axis.text=element_text(size=20))
	g = g + theme(axis.title=element_text(size=20))
	g = g + theme(title=element_text(size=20))
	
	g = g + geom_text(data=d.dat.gil.melt,
					  aes(x=timeround,y = m,label=n,vjust=2))
	
	#g = g + geom_point(data=dat.gil.melt,aes(x=timeround,y=GI),colour="grey")
	return(g)
}
