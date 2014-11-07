require(ggplot2)
require(dplyr)

read = function(name,f=1,Xi='1/4'){
	raw = read.table(paste(name,"/fort.9",sep=''),col.names = c("time","CFL","thrust","lift","Fz"))
	raw$drag = -raw$thrust
	raw = subset(raw,time>-0.3)
	time = seq(-0.2,1.2,0.01)
	drag = predict(loess(drag~time,raw,span=0.03),time)
	lift = predict(loess(lift~time,raw,span=0.03),time)
	smooth = data.frame(time = time, lift = lift, drag = drag, L = 200/f, Xi = Xi)
}

f0p5.Xii4 = read("f0p5_Xii4",0.5)
f0p7.Xii4 = read("f0p7_Xii4",sqrt(0.5))
f1.Xii4 = read("f1_Xii4",1)
f1p4.Xii4 = read("f1p4_Xii4",sqrt(2))
f2.Xii4 = read("f2_Xii4",2)
f2p8.Xii4 = read("f2p8_Xii4",sqrt(8))
f4.Xii4 = read("f4_Xii4",4)
f5p6.Xii4 = read("f5p6_Xii4",sqrt(32))
f8.Xii4 = read("f8_Xii4",8)

ref = subset(f0p5.Xii4,time>0 & time<1)
all = rbind(f0p5.Xii4,f0p7.Xii4,f1.Xii4,f1p4.Xii4,f2.Xii4,f2p8.Xii4,f4.Xii4,f5p6.Xii4,f8.Xii4)

plt = subset(all,L==400|L==200|L==100|L==50)
qdrag = qplot(time,drag,data=plt,linetype=factor(L),geom='path')+theme_bw()+scale_linetype_manual(values=c(3,4,2,1))+labs(x="tU/L",y=expression(C[D]),legend="h/L")+theme(legend.position=c(0.85,0.6))

pdf("drag_time_L.pdf",4,3)
plot(qdrag)
dev.off()

lengths <- group_by(subset(all,time>0 & time<1),L)
lengths <- mutate(lengths,Edrag = drag-ref$drag,Elift = lift-ref$lift)
sums <- summarize(lengths,ml=mean(lift),md=mean(drag),pl=max(lift),pd=max(drag),L2l = sqrt(mean(Elift^2)),L2d = sqrt(mean(Edrag^2)))

elift = ggplot(data=data.frame(x=c(1/250,1/25),y=c(1.2e-1,1.2),z=c(3e-2,3)))+geom_line(aes(x,y),color='grey',linetype=2)+geom_line(aes(x,z),color='grey')+geom_point(data=sums[1:7,],aes(1/L,max(ref$lift)-pl))+theme_bw()+scale_x_log10('h/L',breaks=c(0.005,0.01,0.02,0.04))+scale_y_log10(expression(peak~C[L]~error))

edrag = ggplot(data=data.frame(x=c(1/250,1/25),y=c(0.7e-1,0.7),z=c(1.5e-2,1.5)))+geom_line(aes(x,y),color='grey',linetype=2)+geom_line(aes(x,z),color='grey')+geom_point(data=sums[1:7,],aes(1/L,max(ref$drag)-pd))+theme_bw()+scale_x_log10('h/L',breaks=c(0.005,0.01,0.02,0.04))+scale_y_log10(expression(peak~C[D]~error))

rmselift = ggplot(data=data.frame(x=c(1/250,1/25),y=c(0.63e-1,0.63),z=c(1.8e-2,1.8)))+geom_line(aes(x,y),color='grey',linetype=2)+geom_line(aes(x,z),color='grey')+geom_point(data=sums[1:7,],aes(1/L,L2l))+theme_bw()+scale_x_log10('h/L',breaks=c(0.005,0.01,0.02,0.04))+scale_y_log10(expression(RMSE(C[L])))

pdf("peak_lift_error.pdf",3,3)
plot(elift)
dev.off()

pdf("peak_drag_error.pdf",3,3)
plot(edrag)
dev.off()

pdf("lift_rmse.pdf",3,3)
plot(rmselift)
dev.off()

f1.Xii2 = read("f1_Xii2",1,'1/2')
f1.Xii8 = read("f1_Xii8",1,'1/8')
f1.Xii16 = read("f1_Xii16",1,'1/16')
f1.Xii32 = read("f1_Xii32",1,'1/32')

all = rbind(f1.Xii2,f1.Xii4,f1.Xii8,f1.Xii16,f1.Xii32)
Xlift = qplot(time,lift,data=all,color=Xi,geom='path')+theme_bw()+theme(legend.position=c(0.85,0.6))
Xdrag = qplot(time,drag,data=all,color=Xi,geom='path')+theme_bw()+theme(legend.position=c(0.85,0.6))

pdf("Xi_lift.pdf",4,4)
plot(Xlift)
dev.off()

pdf("Xi_drag.pdf",4,4)
plot(Xdrag)
dev.off()


write.csv(all,"smooth.csv")
