cholSolve = function(L,b){backsolve(L,forwardsolve(t(L),b))}

GLS = function(X,y,lambda=0){
	X <- scale(X)
	X.mu = attributes(X)$`scaled:center`
	X.sd = attributes(X)$`scaled:scale`
	L <- chol(kern(X,X)+lambda*diag(nrow(X)))
	sRi1 = sum(cholSolve(L,matrix(1,nrow(L))))
	mu = sum(cholSolve(L,y))/sRi1
	yp = y-mu
	Riyp = cholSolve(L,yp)
	sig2 = t(yp) %*% Riyp / nrow(L)
	list(X=X,X.mu=X.mu,X.sd=X.sd,L=L,sRi1=sRi1,mu=mu,Riyp=Riyp,sig2=sig2)
}

kern = function(X,new){
	apply(new,1,function(y)apply(X,1,function(x)exp(-sum((x-y)^2))))
}

GLS.predict = function(GLS,new){
	new <- scale(new,center=GLS$X.mu,scale=GLS$X.sd)
	GLS$mu+t(kern(GLS$X,new)) %*% GLS$Riyp
}

GLS.se = function(GLS,new){
	new <- scale(new,center=GLS$X.mu,scale=GLS$X.sd)
	r = kern(GLS$X,new)
	Rir = cholSolve(GLS$L,r)
	as.matrix(GLS$sig2*(1-colSums(r*Rir)+(1-colSums(Rir))^2/GLS$sRi1))
}

