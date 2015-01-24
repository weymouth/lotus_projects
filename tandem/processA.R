source('analysis.R')
source('GLS.R')

dirs = list.dirs(recursive=FALSE)
dirs = grep("caseA",dirs,value=TRUE)

getPC = function(folder){
  p = params(folder)
  forces(folder,p$amp,p$freq)
  percycle(data,p$amp,p$freq)
}

getF = function(folder){
  p = params(folder)
  freqA(forces(folder,p$amp,p$freq,10),p)
}

fData = data.frame(freq=NULL,amp=NULL,Lv=NULL,La=NULL,fs=NULL)
for(folder in dirs){
  fData = rbind(fData,getF(folder))
}

#pcData = data.frame(cycle=NULL,ca=NULL,cv=NULL,amp=NULL,freq=NULL)
#sData = data.frame(meanCa=NULL,ampCa=NULL,meanCv=NULL,ampCv=NULL,amp=NULL,freq=NULL)
#for(folder in dirs){
#  p = getPC(folder)
#  s = as.data.frame(stats(p,first=15))
#  pcData = rbind(pcData,p)
#  sData = rbind(sData,s)
#}

#dat = subset(pcData,cycle<40 & freq %in% c(0.1,0.162,0.18,0.2,0.22,0.275))
#g = ggplot(dat,aes(cycle,cv*amp^2))+theme_bw()+geom_hline(yintercept=0,color='grey')+facet_grid(amp~freq,scales="free_x")+stat_smooth(method='loess')+ylab('Fluid excitation') 
#ppdf(g,'ClV_A_pc.pdf',x=12,y=8)

#GLSout = GLS(cbind(sData$freq,sData$amp), sData$meanCv*sData$amp, lambda=0.001)
#grid = expand.grid(freq=seq(3.2,8,length.out=100),amp=seq(0.05,0.8,length.out=100))
#grid$freq <- 1/grid$freq
#grid$yhat <- GLS.predict(GLSout,cbind(grid$freq,grid$amp))

#q = qplot(1/freq,amp,data=grid,fill=yhat,geom='tile')+geom_point(data=sData,aes(fill=NULL))
#q = q+scale_fill_gradient2('Excitation',limits=c(-0.25,0.25),na.value="red")
#q = q+stat_contour(aes(z=yhat),binwidth=0.01)+stat_contour(aes(z=yhat),bins=1,size=1)
#q = q+labs(x='U*',y='A/D')
#ppdf(q,'ClV_A_GLS.pdf',x=10,y=6)
