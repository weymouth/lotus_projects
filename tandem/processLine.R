source('analysis.R')
library(scales) # for muted
################### functions ############################## 

gettData <- function(tData=data.frame()){
	dirs <- getDirs('caseCamp0.4')
	dirs = dirs[!(dirs %in% tData$folder)] # keep new folders
	for(folder in dirs){
		p = params(folder)
		f = forces(p,20) # drop first 20 cycles
		tData = rbind(tData,cbind(p,f))
	}
	return(tData)
}

getzData <- function(tData,zData = data.frame()){
	dirs <- unique(tData$folder)
	dirs = dirs[!(dirs %in% zData$folder)] # keep new folders
	for(fold in dirs){
		p = params(fold)
		t = subset(tData,folder==fold)
		z = zetaA(t)
		zData = rbind(zData,cbind(p,z))
	}
	return(zData)
}

getfData <- function(zData){
	fData = data.frame()
	for(fold in unique(zData$folder)){
		p = params(fold)
		z = subset(zData,folder==fold)
		f = freqA(z,p,fs=0.172)
		fData = rbind(fData,cbind(p,f))
	}
	return(fData)
}

################### runtime ############################## 

tData = gettData()
zData = getzData(tData)
fData = getfData(zData)

sub = c(4,5.75,6,10,12,16.5,30,62.5,125)
stData = subset(tData,lambda %in% sub)
szData = subset(zData,lambda %in% sub)
sfData = subset(fData,lambda %in% sub)

sstData = stData[seq(1,nrow(stData),len=2e4),]

x = seq(0,1,len=41)
xdat = data.frame(phase=x,lift=cos(2*pi*x))
qplot(phase,lift,data=sstData,geom='line',alpha=I(0.5),color=factor(floor(cycle)))+
	geom_line(data=xdat,color='blue',linetype=2)+facet_wrap(~lambda)+
	scale_color_discrete(c=0,guide=FALSE)+theme_bw()+
	theme(axis.title.y=element_text(angle=0))+ylab(expression(F[yp]))+
	scale_x_continuous('t/T',labels=c("0","1/4","1/2","3/4","1"))
ggsave('phasesAmp0.4.png',width=8,height=6)

qplot(freq,Mod(pres),data=sfData,facets=~lambda,shape=I(19),color=I('blue'))+
	geom_point(aes(x=2*fs-freq,y=-0.05),shape=17,color='grey')+
	geom_point(aes(x=2*fs+freq,y=-0.05),shape=17,color='grey')+
	geom_line(data=subset(szData,zeta<=0.65),aes(x=zeta),color='grey40')+
	theme_bw()+theme(axis.title.y=element_text(angle=0))+
	ylab(expression(group('|',hat(F)[y],'|')))+
	scale_x_continuous(expression(zeta/f[s]),breaks=mean(sfData$fs)*(0:3),labels=0:3)
ggsave('ModzetaAmp0.4.png',width=8,height=6)

ggplot(data=subset(fData,lambda<300),aes(x=lambda))+
	scale_x_log10(expression(lambda/lambda[s]),breaks=5.75*2^(0:5),labels=2^(0:5))+
	scale_y_log10(expression(group('|',hat(F)[y],'|')) )+
	geom_abline(intercept=2,slope=-3,linetype=2,color='blue')+
	geom_abline(intercept=2/3,slope=-2,linetype=2,color='red')+
	geom_line(data=data.frame(x=c(20,300)),aes(x=x,y=2.4/x),linetype=2,color='blue')+
	geom_point(aes(y=Mod(pres)),color='blue',shape=19)+
	geom_point(aes(y=Mod(fric)),color='red',shape=18)+
	theme_bw()+theme(axis.title.y=element_text(angle=0))+
	annotate("text",label="lambda^{-2}",x=5,y=.25,color='red',parse=TRUE)+
	annotate("text",label="lambda^{-1}",x=30,y=.1,color='blue',parse=TRUE)+
	annotate("text",label="lambda^{-3}",x=5,y=.5,color='blue',parse=TRUE)
ggsave('ModlambdaAmp0.4.png',width=8,height=5)

ggplot(data=subset(fData,lambda<300),aes(x=lambda))+
	scale_x_log10(expression(lambda/lambda[s]),breaks=5.75*2^(0:5),labels=2^(0:5))+
	scale_y_continuous(expression(phi),breaks=c(pi/2,0,-pi/2),labels=expression(pi/2,0,-pi/2))+
	theme_bw()+theme(axis.title.y=element_text(angle=0))+
	geom_point(aes(y=Arg(pres)),color='blue',shape=19)+
	geom_point(aes(y=Arg(fric)),color='red',shape=18)
ggsave('philambdaAmp0.4.png',width=8,height=5)

ggplot(data=subset(fData,lambda<300),aes(x=lambda))+
	scale_x_log10(expression(lambda/lambda[s]),breaks=5.75*2^(0:5),labels=2^(0:5))+
	theme_bw()+theme(axis.title.y=element_text(angle=0))+
	geom_point(aes(y=Im(pres)/Mod(pres)),color='blue',shape=19)+
	geom_point(aes(y=Im(fric)/Mod(fric)),color='red',shape=18)

