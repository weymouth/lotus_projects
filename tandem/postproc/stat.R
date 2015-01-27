require(ggplot2)
data = read.table("fort.9",col.names = c("time","CFL","drag","lift","Fz","dragf","liftf","Fzf"))
l = length(data$time)
n = 2000
j = round(seq(30,l,len=min(l,n)))
data = data[j,]
CFL = qplot(time,CFL,data=data,geom="line")
t = 0.5*(min(data$time)+max(data$time))
late = subset(data,time>t)
mdrag = mean(late$drag)
adrag = sqrt(2*mean((late$drag-mdrag)^2))
mlift = mean(late$lift)
alift = sqrt(2*mean((late$lift-mlift)^2))
drag = qplot(time,drag,data=data,geom="line")
y = (mdrag+adrag)*1.2
drag = drag+annotate("text",x=t,y=y,label=paste("mean=",round(mdrag,3)," amp=",round(adrag,3)))
lift = qplot(time,lift,data=data,geom="line")
y = (mlift+alift)*1.2
lift = lift+annotate("text",x=t,y=y,label=paste("mean=",round(mlift,3)," amp=",round(alift,3)))

ppdf = function(plot,name){
     pdf(name,8,4)
     print(plot)
     dev.off()
}

ppdf(CFL,"CFL.pdf")
ppdf(drag,"drag.pdf")
ppdf(lift,"lift.pdf")

data = read.table("mgs.txt", col.names = c("itr","res"))

l = length(data$itr)
n = 2000
data$i = seq(1,l)
j = round(seq(1,l,len=min(l,n)))
#j = unique(c(which(data$itr>1),j))
data = data[j,]

itr = qplot(i,log(itr,2),data=data)+ylab(expression(log[2](iteration)))
res = qplot(i,log(res,10),data=data)+ylab(expression(log[10](residual)))
ppdf(itr,"itr.pdf")
ppdf(res,"res.pdf")

source('analysis.R')
