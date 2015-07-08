library(stringr)
library(ggplot2)

getDirs = function(name){
	dirs = list.dirs(recursive=FALSE)
	grep(name,dirs,value=TRUE)
}

params = function(folder){
  text <- readLines(paste(folder,'/lotus.f90',sep=''))
  amp = as.numeric(str_match(text[grep('amp =',text)],"[0?.0-9]+"))
  freq = as.numeric(str_match(text[grep('freq =',text)],"[0?.0-9]+"))
  lambda = round(4/freq)/4
  data.frame(folder=folder, amp=amp, freq=freq, lambda = lambda)
}

forces = function(p,first=1){
  data = read.table(paste(p$folder,'/fort.9',sep=''),
	col.names = c("time","CFL","drag","lift","Fz","dragf","liftf","Fzf"))
  data$time = data$time-0.25/p$freq
  data$cycle = data$time*p$freq
  data$phase = data$cycle%%1
  last = floor(max(data$cycle))
  first = ceiling(100*p$freq)
  return(subset(data,cycle>=first & cycle<last))
}

myfft <- function(x,t){
	len <- length(t)
	t2 = seq(t[1],t[len],len=len)
	x2 = spline(t,x,xout=t2)$y
	dt <- t2[len]-t2[(len-1)]
	coeffs <- fft(x2)/len*2
	mag <- Mod(coeffs)
	phase <- Arg(coeffs)
	freq <- 0:(len-1)/dt/len
	fourier <- data.frame(freq,mag,phase,coeffs)
	return(fourier)
}

zetaA <- function(tData){
	pres = subset(myfft(tData$lift,tData$time),freq<1)
	fric = subset(myfft(tData$liftf,tData$time),freq<1)
	return(data.frame(zeta=pres$freq,pres=pres$coeffs,fric=fric$coeffs))
}

peaks <- function(data,x){
	i = which(diff(sign(diff(x)))==-2)+1
	data[i,]
}

extrema <- function(data,x){
	i = which(abs(diff(sign(diff(x))))==2)+1
	data[i,]
}

freqA <- function(zData,p,fs=0.188){
	if(nrow(zData)==0) return(data.frame(pres=NA,fric=NA,fs=NA))
	f = zData[which.min(abs(zData$zeta-p$freq)),]
	zp = peaks(zData,Mod(zData$pres))
	zp$mag = Mod(zp$pres)*exp(-(zp$zeta-fs)^2*1e3)
	fs = zp$zeta[which.min(-zp$mag)]
	data.frame(pres = f$pres, fric = f$fric, fs=fs)
}
