source('analysis.R')
library(scales) # for muted
################### functions ############################## 

getF <- function(folder){
  p = params(folder)
  f = freqA(forces(p,20),p) # drop first 20 cycles
  cbind(p,f)
}

getData <- function(fData=data.frame()){
	dirs <- getDirs('caseC')
	dirs = dirs[!(dirs %in% fData$folder)] # keep new folders
	for(folder in dirs){
	  f = getF(folder)
	  fData = rbind(fData,f)
	}
	return(fData)
}

qfill <- function(aesfill,data){
	ggplot(data=data,aes(x=lambda,y=amp))+geom_tile(aesfill)+
	scale_fill_gradient2()+labs(y='a',x=expression(lambda))+
	theme_minimal()+theme(title=element_text(size=15),
#	plot.background = element_rect(fill='#A4AEB5'),
	axis.title.y=element_text(angle=0))
}

interp = function(fData){
	newData = fData
	newlambda <- seq(4,12,0.25)
	newlambda = newlambda[!(newlambda %in% unique(fData$lambda))]
	for( a in unique(fData$amp) ){
	for( l in newlambda ){
		left = subset(fData,amp==a & lambda==l+0.25)
		right = subset(fData,amp==a & lambda==l-0.25)
		newData = rbind(newData,
			data.frame(folder="interp",amp=a,freq=1/l,lambda=l,
			Rp = 0.5*(left$Rp+right$Rp),
			Ip = 0.5*(left$Ip+right$Ip),
			Rf = 0.5*(left$Rf+right$Rf),
			If = 0.5*(left$If+right$If),
			fs = 0.5*(left$fs+right$fs)))
	}}
	return(newData)
}

################### runtime ############################## 

fData <- getData(if(exists('fData')){fData})

newData <- interp(fData)

qfill(aes(fill=Rp),newData)+labs(fill=expression('\U211C'(hat(F)[yp])))
ggsave('RFyp.png',width=6,height=3)
qfill(aes(fill=Rf),newData)+labs(fill=expression('\U211C'(hat(F)[yf])))
ggsave('RFyf.png',width=6,height=3)

qfill(aes(fill=Ip+If),newData)+labs(fill=expression('\U2111'(hat(F)[y])))
ggsave('IFy.png',width=6,height=3)
qfill(aes(fill=tanh(Ip+If)),newData)+labs(fill=expression('\U2111'(hat(F)[y])))
ggsave('tanhIFy.png',width=6,height=3)
qfill(aes(fill=Ip),newData)+labs(fill=expression('\U2111'(hat(F)[yp])))
ggsave('IFyp.png',width=6,height=3)
qfill(aes(fill=Ip),subset(newData,lambda>=7))+labs(fill=expression('\U2111'(hat(F)[yp])))
ggsave('IFyp_long.png',width=6,height=3)

fs = qfill(aes(fill=fs/mean(newData$fs)-1),newData)+labs(fill=expression(Delta~f[s]~'(%)'))
ggsave('fs.png',width=6,height=3)

