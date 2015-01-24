source('../analysis.R')
p = params('./')
omega = 2*pi*p$freq
tData = subset(forces(p$folder,p$amp,p$freq),time>1)
pcData = percycle(tData,p$amp,p$freq)
s = stats(pcData)

q = qplot(phase,lift,data=tData,color=cycle,alpha=I(0.05))+ylab(expression(C[L]))
ppng(q,"ClAll.png")

q = qplot(phase,La,data=tData,alpha=I(0.1))
q = q+labs(y=expression(C[L]*a),title = paste("fluid mass coefficient =",round(s$meanCa,3)))
ppng(q,"ClaLate.png")

q = qplot(phase,Lv,data=tData,alpha=I(0.1))
q = q+labs(y=expression(C[L]*v),title = paste("fluid excitation coefficient =",round(s$meanCv,3)))
ppng(q,"ClvLate.png")

q = qplot(cycle,ca,data=pcData,geom=c('line','point','smooth'),method='loess')
ppng(q,'ClaCycle.png')
q = qplot(cycle,cv,data=pcData,geom=c('line','point','smooth'),method='loess')
ppng(q,'ClvCycle.png')

