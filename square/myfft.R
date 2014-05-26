myfft <- function(x,t){
	len <- length(t)
	t2 = seq(t[1],t[len],len=len)
	x2 = spline(t,x,xout=t2)$y
	dt <- t2[len]-t2[(len-1)]
	coeffs <- fft(x2)/len*2
	mag <- Mod(coeffs)^2*dt^2; mag[1] = 0
	phase <- Arg(coeffs)
	freq <- 0:(len-1)/dt/len
	fourier <- data.frame(freq,mag,phase,coeffs)
	return(fourier)
	}

peakFreq <- function(data,d){
	f = subset(myfft(data[,d],data[,1]),freq<1)
	peak = 0
	for( i in 1:length(f$mag)){
		if(f$mag[i]>peak){
			peak = f$mag[i]
			j = i
		}
	}
	peak = sum(f$mag[(j-1):(j+1)])
	freq = f$freq[j]
	return(data.frame(freq,peak))
	}

dfft <- function(data,d,cut=50){
	d2 = subset(data,time>cut)
	f = subset(myfft(d2[,d],d2[,1]),freq<1 & freq>0)
	f$name = data$name[1]
	f$freq = round(f$freq,2)
	f$mag = f$mag/max(f$mag)
	return(f)
	}
