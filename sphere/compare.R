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

drag = qplot(time,4*drag,data=data,linetype=name,geom='line')+labs(y=expression(C[Dp]),x='tU/D',linetype="")
ppdf(drag,"drag_compare.pdf")
