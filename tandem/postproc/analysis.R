source('../analysis.R')
p = params('./')
tData = forces(p)

x = seq(0,1,0.01); y = max(tData$lift)
q = qplot(phase,lift,data=tData,alpha=I(0.05),color=I('blue'))+geom_point(aes(y=liftf),color='red',alpha=0.05)+geom_line(data=data.frame(phase=x,lift=y*cos(2*pi*x)),linetype=2)+ylab(expression(F[y]))+theme_bw()+theme(axis.title.y=element_text(angle=0))
ggsave("Phase.png",width=6,height=4)

fData = subset(myfft(tData$lift,tData$time),freq<1)
ffData = subset(myfft(tData$liftf,tData$time),freq<1)

fplot = function(aesy){
	ggplot(data=fData,aes(x=freq))+xlab(expression(zeta))+
	geom_line(aesy,color='blue')+
	geom_line(data=ffData,aesy,color='red')+
	theme_bw()+theme(axis.title.y=element_text(angle=0))
}
q = fplot(aes(y=Mod(coeffs)))+ylab(expression(group('|',hat(F)[y],'|')))
ggsave("Mod.png",width=6,height=4)
q = fplot(aes(y=Re(coeffs)))+ylab(expression("\U211c"(hat(F)[y])))
ggsave("Re.png",width=6,height=4)
q = fplot(aes(y=Im(coeffs)))+ylab(expression("\U2111"(hat(F)[y])))
ggsave("Im.png",width=6,height=4)
