i = which(grid$s2+grid$D2>max(grid$s2+grid$D2)-1e-10)
dat = rbind(dat,grid[i,1:3])
GLSout = GLS(dat[,1:2], dat[,3])
grid$yhat = GLS.predict(GLSout,grid[,1:2])
grid$s2 = GLS.se(GLSout,grid[,1:2])
grid$D2 <- dyhat(grid)

print(paste("real SE=",max((grid$yhat-grid$y)^2),", SE+D2=",max(grid$s2+grid$D2),sep=""))
