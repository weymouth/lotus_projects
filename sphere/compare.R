library(dplyr);library(ggplot2)

ppdf = function(plot,name){
     pdf(name,8,4)
     print(plot)
     dev.off()
}

A = read.table("MoveTp2D100/fort.9",col.names = c("time","CFL","drag","lift","Fz"))
A$name = "move"
B = read.table("GustTp2D100/fort.9",col.names = c("time","CFL","drag","lift","Fz"))
B$name = "gust"

data = rbind(A,B)
data$bouy = 0
data$bouy[data$time<=2.2 & data$name=="gust"] = 2.5*4/3

drag = qplot(time,-4*drag,data=data,linetype=name,geom='line')+
  geom_line(aes(y=bouy),data=subset(data,bouy>0))+
  labs(y=expression(C[Dp]),x='tU/D',linetype="")
ppdf(drag,"drag_compare.pdf")
