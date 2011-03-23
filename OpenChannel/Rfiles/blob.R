step <- function(t,dt=2) floor(t/dt)*dt
require(ggplot2);require(RColorBrewer)

blob <- read.csv(file="blob.csv",head=TRUE,sep=",")

blob$radius <- (3./4./pi*blob$vol)^(1./3.)
blob$max <- apply(blob[6:8],1,max)
blob$min <- apply(blob[6:8],1,min)
blob$ratio <- blob$max/blob$min
blob$bigtime <- apply(blob[1],1,step)

ggrt <- qplot(data=blob,radius,time,size=radius,alpha=0.01)+scale_size(to = c(0, 2))+scale_alpha(legend=FALSE)
ggxt <- qplot(data=blob,x,time,size=radius,alpha=0.01)+scale_size(to = c(0, 2.5))+scale_alpha(legend=FALSE)

crop <- subset(blob,time<42)
ggcrop <- qplot(data=crop,x=radius,weight=vol/sum(vol)*10,binwidth=0.05, facets=bigtime~.)

steady <- subset(blob,time>34)
ggsteady <- qplot(data=steady,x=radius,weight=vol/sum(vol),binwidth=0.05,xlim=c(0,3.5))

sized <- subset(blob,min>0&ratio<=10)
ggrtc <- qplot(data=sized,radius,time,size=radius,colour=ratio)+scale_size(to = c(0, 2))+scale_color_gradient(low="yellow",high="red")

