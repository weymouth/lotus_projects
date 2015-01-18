require(ggplot2)
source('GLS.R')

fout = function(x){x[1]*sin(2*pi*x[2])}

dat = expand.grid(x1=seq(0,1,1/3),x2=seq(0,1,1/3))
# dat$x1 = dat$x1+runif(nrow(dat),-0.1,0.1)
# dat$x2 = dat$x2+runif(nrow(dat),-0.1,0.1)
dat$y = apply(dat,1,fout)
GLSout = GLS(dat[,1:2], dat[,3])

grid = expand.grid(x1=seq(0,1,.05),x2=seq(0,1,.05))
grid$y = apply(grid,1,fout)

grid$yhat = GLS.predict(GLSout,grid[,1:2])
grid$s2 = GLS.se(GLSout,grid[,1:2])

grid$D2 = 0
dyhat = function(grid){
	for(i in 2:20){
	for(j in 2:20){
		k = i+(j-1)*21
		di2 = grid$yhat[k-1]-2*grid$yhat[k]+grid$yhat[k+1]
		dj2 = grid$yhat[k-21]-2*grid$yhat[k]+grid$yhat[k+21]
		dij = (grid$yhat[k+22]-grid$yhat[k+20]-grid$yhat[k-20]+grid$yhat[k-22])/4
		grid$D2[k] = di2^2+dj2^2+dij^2
	}}
	return(grid$D2)
}
grid$D2 <- dyhat(grid)


#sub = grid[sample(nrow(grid), 4), ]

#qplot(x1,x2,data=grid,fill=y,geom='tile')+geom_point(data=sub,aes(fill=NULL))+scale_fill_gradient2()
