require(ggplot2)
data = read.table("fort.9",col.names = c("time","CFL","thrust","lift"))
data$drag = -data$thrust
l = length(data$time)
n = 2000
j = round(seq(1,l,len=min(l,n)))
data = data[j,]
CFL = qplot(time,CFL,data=data,geom="line")
data = subset(data,time>0.1)
t = 0.5*max(data$time)
late = subset(data,time>t)
mdrag = mean(late$drag)
plift = max(late$lift)
drag = qplot(time,drag,data=data,geom="line")
drag = drag+annotate("text",x=1,y=mdrag,label=paste("mean",round(mdrag,2)))
lift = qplot(time,lift,data=data,geom="line")
lift = lift+annotate("text",x=1,y=plift*0.95,label=paste("peak",round(plift,2)))

ppdf = function(plot,name){
     pdf(name,8,4)
     print(plot)
     dev.off()
}

ppdf(CFL,"CFL.pdf")
ppdf(drag,"drag.pdf")
ppdf(lift,"lift.pdf")
print(mdrag)
print(plift)

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
