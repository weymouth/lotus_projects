library(dplyr);library(ggplot2)

data = read.table("fort.9",col.names = c("time","CFL","drag","lift","side"))
l = length(data$time)
n = 2000
j = round(seq(4,l,len=min(l,n)))
data = data[j,]
CFL = qplot(time,CFL,data=subset(data,CFL<2),geom="line")

stats = data %>% filter(time>mean(time)) %>%
	summarize(t = min(time),
			mdrag = mean(drag), mlift = mean(lift),
		  adrag = mad(drag), alift = mad(lift),
			ndrag = min(drag), nlift=min(lift))
attach(stats)

drag = qplot(time,drag,data=data,geom="line")
drag = drag+annotate("text",x=t,y=ndrag,label=paste("mean=",round(mdrag,3)," amp=",round(adrag,3)))
lift = qplot(time,lift,data=data,geom="line")
lift = lift+annotate("text",x=t,y=nlift,label=paste("mean=",round(mlift,3)," amp=",round(alift,3)))

ppdf = function(plot,name){
     pdf(name,8,4)
     print(plot)
     dev.off()
}

ppdf(CFL,"11CFL.pdf")
ppdf(drag,"01drag.pdf")
ppdf(lift,"02lift.pdf")

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
