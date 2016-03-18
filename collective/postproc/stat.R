library(dplyr);library(ggplot2);library(zoo)
ppdf = function(plot,name){
     pdf(name,8,4)
     print(plot)
     dev.off()
}

data = read.table("fort.9",col.names = c("time","CFL","drag","side","lift"))
k = 4
data$drag = rollmean(data$drag,2*k+1,fill=TRUE)
data$side = rollmean(data$side,2*k+1,fill=TRUE)
data$lift = rollmean(data$lift,2*k+1,fill=TRUE)

l = length(data$time)
n = 2000
j = round(seq(1+k,l-k,len=min(l-k,n)))
data = data[j,]

CFL = qplot(time,CFL,data=data,geom="line")
ppdf(CFL,"11CFL.pdf")

easy = data %>% filter(time>mean(time)) %>%
  summarize(t = min(time),
			mdrag = mean(drag), adrag = mad(drag), ndrag = 1.05*min(drag),
      mside = mean(side), aside = mad(side), nside = 1.05*min(side),
      mlift = mean(lift), alift = mad(lift), nlift = 1.05*min(lift))
attach(easy)

drag = qplot(time,drag,data=data,geom="line")+ylab("Fx")
drag = drag+annotate("text",x=t,y=ndrag,label=paste("mean=",round(mdrag,3)," amp=",round(adrag,3)))
ppdf(drag,"01drag.pdf")

side = qplot(time,side,data=data,geom="line")+ylab("Fy")
side = side+annotate("text",x=t,y=nside,label=paste("mean=",round(mside,3)," amp=",round(aside,3)))
ppdf(side,"02side.pdf")

lift = qplot(time,lift,data=data,geom="line")+ylab("Fz")
lift = lift+annotate("text",x=t,y=nlift,label=paste("mean=",round(mlift,3)," amp=",round(alift,3)))
ppdf(lift,"03lift.pdf")

data = read.table("fort.13",col.names = c("time","position","velocity"))
data = data[c(FALSE,TRUE,FALSE,FALSE),]
l = length(data$time)
n = 2000
j = round(seq(1,l,len=min(l,n)))

data = data[j,]
s = qplot(time,position,data=data,geom="line")+ylab("length")
u = qplot(time,velocity,data=data,geom="line")+ylab("rate")
ppdf(s,"04length.pdf")
ppdf(u,"05rate.pdf")

data = read.table("mgs.txt", col.names = c("itr","res"))

l = length(data$itr)
n = 2000
data$i = seq(1,l)
j = round(seq(1,l,len=min(l,n)))
#j = unique(c(which(data$itr>1),j))
data = data[j,]

itr = qplot(i,log(itr,2),data=data)+ylab(expression(log[2](iteration)))
res = qplot(i,log(res,10),data=data)+ylab(expression(log[10](residual)))
ppdf(itr,"13itr.pdf")
ppdf(res,"12res.pdf")

try(source('analysis.R'))
