cex.point <- 5


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
		 main=paste0("Susceptible Proportion\n(at calendar time s=",round(time.s),")"), 
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
	points(tt[time.s],S[time.s],cex=cex.point, pch=1)
	points(tt[time.s],S[time.s],cex=cex.point-1, pch=16)
	
	
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
		 main=paste0("Incidence\n(at calendar time t=",round(time.t),")"), 
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
	points(tt[time.t],I[time.t],cex=cex.point,pch=1)
	points(tt[time.t],I[time.t],cex=cex.point-1,pch=16)
	
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
		 main = paste0("Forward GI Distribution\n(at calendar time s=",round(time.s),")") )
	
	lines(gi.time,gi.intrinsic,col="gray",lwd=2)
	
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
		 main = paste0("Backward GI Distribution\n(at calendar time s=",round(time.t),")") )
	
	lines(gi.time,gi.intrinsic,col="gray",lwd=2)
	abline(v=mean.g, lty=2,lwd=1)
	abline(v=m.GI.bck, lty=1,lwd=3)
}


####################################################
###    SINGLE FRAME    ###
####################################################




plot.one.frame <-function(time.t.one)
{
	par(mfrow=c(3,2))
	par(cex.main=2, cex.lab=1.6, 
		cex.axis=1.6, las = 1)
	
	i <- 1
	
	correct.first.inc <- FALSE
	
	# LEFT COLUMN = backward GI
	# RIGHT COLUMN = forward GI
	
	### ----- first Row -----
	plot.I(time.t=time.t.one[i])
	plot.S(slope = slope,
		   inflection = inflection,
		   time.s=time.t.one[i])
	
	### ----- second Row -----
	plot.GI.bck(time.t = time.t.one[i],tt,g,I,correct.first.inc)
	plot.GI.fwd(time.s = time.t.one[i],g, S, horizon,tt)
	
	
	### ----- third Row, LEFT COLUMN(backwd GI) -----
	
	b.bar.plot <- GI.ODE$GI.bck.mean
	b.bar.plot.time <- GI.ODE$time
	
	dt <- max(diff(b.bar.plot.time))
	idx = which(abs(b.bar.plot.time - time.t.one[i])<=dt/2)[1]
	
	plot(x=b.bar.plot.time, y=b.bar.plot, 
		 main = "Backward GI Mean",
		 xlim=c(0,horizon), xlab="calendar time",
		 ylim=c(0,max(b.bar.plot,na.rm = T)),		 #delete --> xlim.g*do.vert.density,na.rm = T)), 
		 ylab="",
		 lwd=1,col="lightgrey",
		 typ="l",
		 las=1)
	lines(x=b.bar.plot.time[1:idx], y=b.bar.plot[1:idx],
	      lwd = lwd.gi,
	      col = "black")
	points(x=b.bar.plot.time[idx], y=b.bar.plot[idx],cex=cex.point,pch=1)
	points(x=b.bar.plot.time[idx], y=b.bar.plot[idx],cex=cex.point-1,pch=16)

	abline(h=mean.g,lty=2)
	
	
	### ----- Third Row, RIGHT COLUMN (Fwd GI) -----
	
	f.bar.plot <- GI.ODE$GI.fwd.mean
	f.bar.plot.time <- GI.ODE$time
	dt <- max(diff(f.bar.plot.time))

	# t.idx = which(round(f.bar.plot.time)%in%time.t.one)
	t.idx = which(abs(f.bar.plot.time - time.t.one[i])<=dt/2)[1]
	f.bar <- f.bar.plot[t.idx]

	plot(x=f.bar.plot.time, y=f.bar.plot, 
		 main = "Forward GI Mean",
		 xlim=c(0,horizon*0.95),
		 xlab="calendar time",
		 ylim=c(mean.g*0.75, mean.g*1.05), 
		 ylab="",
		 lwd=1,col="lightgrey",
		 typ="l",
		 las=1)
	
	lines(x=f.bar.plot.time[1:t.idx], y=f.bar.plot[1:t.idx], 
	      lwd = lwd.gi,
	      col="black")
	
	points(x=f.bar.plot.time[t.idx], y=f.bar.plot[t.idx],cex=cex.point,pch=1)
	points(x=f.bar.plot.time[t.idx], y=f.bar.plot[t.idx],cex=cex.point-1,pch=16)
	abline(h=mean.g,lty=2)
}



