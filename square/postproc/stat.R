library(dplyr);library(ggplot2)

data = read.table("fort.9",col.names = c("time","CFL","drag","lift","Fz","y","v"))
l = length(data$time)
n = 2000
j = round(seq(5,l,len=min(l,n)))
data = data[j,]
CFL = qplot(time,CFL,data=subset(data,CFL<1),geom="line")

easy = data %>% filter(time>mean(time)) %>%
  summarize(t = min(time),
			mdrag = mean(drag), mlift = mean(lift), mpos = mean(y),
		  adrag = mad(drag), alift = mad(lift), apos = mad(y),
			ndrag = 1.05*min(drag), nlift=1.05*min(lift), npos = 1.05*min(y)
			)
attach(easy)

medium = data %>% filter( y>lead(y) & y>lag(y) ) %>%
    mutate(period = time-lag(time))

harder = data %>% filter( (y>lead(y) & y>lag(y)) | (y<lead(y) & y<lag(y)), time>mean(time)) %>%
    transmute(amp = abs(y-lag(y))/2) %>%
    filter(cume_dist(desc(amp))<0.1) %>%
    summarise(As=max(amp), As10=mean(amp))
attach(harder)

data$drag[data$time<1] = mdrag
data$lift[data$time<1] = mlift
drag = qplot(time,drag,data=data,geom="line")
drag = drag+annotate("text",x=t,y=ndrag,label=paste("mean=",round(mdrag,3)," amp=",round(adrag,3)))
lift = qplot(time,lift,data=data,geom="line")
lift = lift+annotate("text",x=t,y=nlift,label=paste("mean=",round(mlift,3)," amp=",round(alift,3)))
pos = qplot(time,y,data=data,geom="line")
pos = pos+annotate("text",x=t,y=npos,label=paste("mean=",round(mpos,3)," A*=",round(As,3),
                                                 " A*[10]=",round(As10,3)))

ppdf = function(plot,name){
     pdf(name,8,4)
     print(plot)
     dev.off()
}

ppdf(CFL,"11CFL.pdf")
ppdf(drag,"01drag.pdf")
ppdf(lift,"02lift.pdf")
ppdf(pos,"03pos.pdf")

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
