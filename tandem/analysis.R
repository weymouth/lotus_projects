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
  if(first<0) first <- last+first
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

freqA <- function(tData,p,fs=0.188){
	if(nrow(tData)==0) return(data.frame(freq=p$freq, amp=p$amp, 
		Rp = NA, Ip = NA, Rf = NA, If = NA, fs=NA))
	fData = subset(myfft(tData$lift,tData$time),freq<1)
	ffData = subset(myfft(tData$liftf,tData$time),freq<1)
	f = fData[which.min(abs(fData$freq-p$freq)),]
	ff = ffData[which.min(abs(fData$freq-p$freq)),]
	fData$prob = -fData$mag*exp(-(fData$freq-fs)^2*1e3)
	fs = fData$freq[which.min(fData$prob)]
	data.frame(Rp = Re(f$coeffs),Ip = Im(f$coeffs), 
		   Rf = Re(ff$coeffs),If = Im(ff$coeffs),fs=fs)
}
