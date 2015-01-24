library(stringr)
library(ggplot2)

params = function(folder){
  text <- readLines(paste(folder,'/lotus.f90',sep=''))
  amp = as.numeric(str_match(text[grep('amp =',text)],"[0?.0-9]+"))
  freq = as.numeric(str_match(text[grep('freq =',text)],"[0?.0-9]+"))
  list(folder=folder,amp=amp,freq=freq)
}

forces = function(folder,amp,freq,first=1){
  omega = 2*pi*freq
  data = read.table(paste(folder,'/fort.9',sep=''),
	col.names = c("time","CFL","drag","lift","Fz","dragf","liftf","Fzf"))
#  data$y = amp*sin(omega*data$time)
#  data$v = omega*amp*cos(omega*data$time)
#  data$a = -omega^2*amp*sin(omega*data$time)
#  data$La = (data$lift+data$liftf)*data$a
#  data$Lv = (data$lift+data$liftf)*data$v
  data$cycle = data$time*freq
  data$phase = data$cycle%%1
  return(subset(data,cycle>=first))
}

percycle = function(data,amp,freq){
  omega = 2*pi*freq
  m = 1
  n = m*max(data$cycle)-m+1
  cycle=rep(0,n)
  cv   =rep(0,n)
  ca   =rep(0,n)
  for(i in 1:n){
    c0 = (i-1)/m
    per = subset(data,cycle>=c0 & cycle<c0+1)
    cycle[i] = c0
    ca[i] = mean(per$La)/(amp^2*omega^4*pi/4)
    cv[i] = mean(per$Lv)/(amp^2*omega^2)
  }
  data.frame(cycle,ca,cv,amp,freq)
}

stats = function(data,first=max(data$cycle-15)){
  late = subset(data,cycle>=first)
  list(	meanCa=mean(late$ca),ampCa=sqrt(2*var(late$ca)),
	meanCv=mean(late$cv),ampCv=sqrt(2*var(late$cv)),
	amp=data$amp[1],freq=data$freq[1])
}

myfft <- function(x,t){
	len <- length(t)
	t2 = seq(t[1],t[len],len=len)
	x2 = spline(t,x,xout=t2)$y
	dt <- t2[len]-t2[(len-1)]
	coeffs <- fft(x2)/len*2
	mag <- Mod(coeffs)^2
	phase <- Arg(coeffs)
	freq <- 0:(len-1)/dt/len
	fourier <- data.frame(freq,mag,phase,coeffs)
	return(fourier)
}

freqA <- function(tData,p){
	fData = subset(myfft(tData$lift+tData$liftf,tData$time),freq<1)
	f = fData[which.min(abs(fData$freq-p$freq)),]
	fs = fData[order(-fData$mag),]$freq
	fs = fs[fs>p$freq/1.25 & fs<p$freq*1.25]
	data.frame(freq=p$freq, amp=p$amp, Lv = f$mag*cos(f$phase),La = f$mag*sin(f$phase), fs=fs[1])
}

ppdf = function(plot,name,x=8,y=4){
     pdf(name,x,y)
     print(plot)
     dev.off()
}

ppng = function(plot,name){
     png(name,8,4,units='in',res=210)
     print(plot)
     dev.off()
}

