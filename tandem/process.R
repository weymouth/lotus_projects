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
	axis.title.y=element_text(angle=0))
}

interp = function(fData){
	newData = fData
	newlambda <- seq(4,15,0.25)
	newlambda = newlambda[!(newlambda %in% unique(fData$lambda))]
	for( l in newlambda ){
		leftR = subset(fData,lambda==l+0.25)
		rightR = subset(fData,lambda==l-0.25)
		for( a in intersect(leftR$amp,rightR$amp)){
			left = subset(leftR,amp==a)
			right = subset(rightR,amp==a)
			newData = rbind(newData,
				data.frame(folder="interp",amp=a,freq=1/l,lambda=l,
				pres = 0.5*(left$pres+right$pres),
				fric = 0.5*(left$fric+right$fric),
				fs = 0.5*(left$fs+right$fs)))
		}
	}
	return(newData)
}

################### runtime ############################## 

fData <- getData(if(exists('fData')){fData})

newData <- interp(fData)

qfill(aes(fill=Arg(pres)),subset(newData,amp<0.9 & lambda<=12))+labs(fill=expression('Arg'(hat(F)[yp])))
ggsave('AFyp.png',width=6,height=3)
qfill(aes(fill=Arg(pres)),subset(newData,amp<0.9 & lambda<=12))+labs(fill=expression('Arg'(hat(F)[yp])))
ggsave('AFyp.png',width=6,height=3)

qfill(aes(fill=Im(pres+fric)),subset(newData,amp<0.9 & lambda<=12))+labs(fill=expression('\U2111'(hat(F)[y])))
ggsave('IFy.png',width=6,height=3)
qfill(aes(fill=Im(pres)),subset(newData,amp<0.9 & lambda<=12))+labs(fill=expression('\U2111'(hat(F)[yp])))
ggsave('IFyp.png',width=6,height=3)
qfill(aes(fill=Im(pres)),subset(newData,amp<0.9 & lambda<=12 & lambda>=8))+labs(fill=expression('\U2111'(hat(F)[yp])))
ggsave('IFyp_long.png',width=6,height=3)

qfill(aes(fill=sin(Arg(pres))),subset(newData,amp<0.9 & lambda<=12))+labs(fill=expression(sin(phi[p])))
ggsave('sinp.png',width=6,height=3)
qfill(aes(fill=cos(Arg(pres))),subset(newData,amp<0.9 & lambda<=12))+labs(fill=expression(cos(phi[p])))
ggsave('cosp.png',width=6,height=3)

qfill(aes(fill=Re(pres+fric)),subset(newData,amp<0.9 & lambda<=12))+labs(fill=expression('\U211C'(hat(F)[y])))
ggsave('RFy.png',width=6,height=3)
qfill(aes(fill=Re(pres)),subset(newData,amp<0.9 & lambda<=12))+labs(fill=expression('\U211C'(hat(F)[yp])))
ggsave('RFyp.png',width=6,height=3)

fs = qfill(aes(fill=fs/mean(newData$fs)-1),newData)+labs(fill=expression(Delta~f[s]~'(%)'))
ggsave('fs.png',width=6,height=3)

