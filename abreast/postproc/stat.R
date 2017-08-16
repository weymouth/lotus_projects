library(dplyr);library(ggplot2)

runfilter <- function(x,k=21){
  n <- length(x); l=(k-1)/2
  y = vector(mode = "numeric", length = n)
  f = l+1-abs(l+1-seq(1,k))
  f = f/sum(f)
  for(i in (1+l):(n-l)){
    y[i] = mean(x[(i-l):(i+l)]*f)
  }
  y[1:l] = y[1+l]; y[(n-l+1):n] = y[n-l]
  y
}
peaks <- function(x,tol=1e-10){
  n <- length(x)
  j<-2:(n-1)
  x[j] <- (x[j]-x[(j-1)])*(x[(j+1)]-x[j])
  x[1] <- 1; x[n] <- 1
  x<0
}

ppdf = function(plot,name){
  pdf(name,8,4)
  print(plot)
  dev.off()
}

data = read.table("fort.9",col.names = c("time","CFL","drag","lift","Fz"))

easy_range = data %>% filter(peaks(runfilter(drag))) %>% 
  top_n(21,time) %>% summarise(i=min(time),e=max(time))

easy = data %>% filter(time>=easy_range$i,time<=easy_range$e) %>%
  summarize(t = min(time),
            mdrag = mean(drag), mlift = mean(lift),
            adrag = mad(drag), alift = mad(lift),
            ndrag = 1.05*min(drag), nlift=1.05*min(lift)
  )
attach(easy)

l = length(data$time)
n = 2000
j = round(seq(5,l,len=min(l,n)))
data = data[j,]

CFL = qplot(time,CFL,data=subset(data,CFL<1),geom="line")

data$drag[data$time<1] = mdrag
drag = qplot(time,drag,data=data,geom="line")
drag = drag+annotate("text",x=t,y=ndrag,label=paste("mean=",round(mdrag,3)," amp=",round(adrag,3)))
lift = qplot(time,lift,data=data,geom="line")
lift = lift+annotate("text",x=t,y=nlift,label=paste("mean=",round(mlift,3)," amp=",round(alift,3)))

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
