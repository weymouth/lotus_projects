library(dplyr);library(ggplot2)

data = read.table("fort.9",col.names = c("time","CFL","Re","Ret"))
l = length(data$time)
n = 2000
j = round(seq(5,l,len=min(l,n)))
data = data[j,]

CFL = qplot(time,CFL,data=data,geom="line")
Re = qplot(time,Re,data=data,geom="line")
Ret = qplot(time,Ret,data=data,geom="line")

ppdf = function(plot,name){
     pdf(name,8,4)
     print(plot)
     dev.off()
}

ppdf(CFL,"11CFL.pdf")
ppdf(Re,"01Re.pdf")
ppdf(Ret,"02Ret.pdf")

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
