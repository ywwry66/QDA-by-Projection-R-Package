require(MASS)
require(nloptr)
#require(NlcOptim)

MisRate <- function(a,mu0,mu1,sigma0,sigma1){
  m0=t(a) %*% mu0; m1=t(a) %*% mu1
  s0=sqrt(t(a) %*% sigma0 %*% a); s1=sqrt(t(a) %*% sigma1 %*% a)
  delta=sqrt((m0-m1)^2+(s0^2-s1^2)*log(s0^2/s1^2))
  if(s0==s1){
    s=s0
    return(pnorm(-abs(m0-m1)/(2*s)))
  }else 
    return((1/2)*(pnorm((s1*(m1-m0)-s0*delta)/(s0^2-s1^2))+1-pnorm((s1*(m1-m0)+s0*delta)/(s0^2-s1^2)))+
      (1/2)*(pnorm((s0*(m1-m0)+s1*delta)/(s0^2-s1^2))-pnorm((s0*(m1-m0)-s1*delta)/(s0^2-s1^2))))
}
drt_2 <- function(mu0,mu1,sigma0,sigma1,sigma,lambda,iter){
  p=length(mu0)
  a=matrix(0,p,iter)
  target <- function(a) MisRate(a,mu0,mu1,sigma0,sigma1)+lambda*sum(abs(a))
  const <- function(a) 1-sum(a^2)
  #const <- function(a) return(list(c=1-sum(a^2)))
  x0=ginv(sigma)%*%(mu0-mu1)
  x0=x0/sqrt(sum(x0^2))
  #print(x0)
  #d=optim(par=ginv(sigma)%*%(mu0-mu1),fn=MisRate,method="BFGS",mu0=mu0,mu1=mu1,sigma0=sigma0,sigma1=sigma1)$par
  #d=d/sqrt(sum(d^2))
  #d[which(abs(d)<lambda)]=0
  #a[,1]=d
  s=nloptr(x0=x0,eval_f=target,eval_g_ineq =const,opts = list("algorithm"="NLOPT_LN_COBYLA","xtol_rel"=1.0e-4,"maxeval"=30000),lb=rep(-1,p),ub=rep(1,p))
  #s=solnl(x0, objfun = target,confun = const,lb=rep(-1,p),ub=rep(1,p))
  #print(s)
  a[,1]=s$solution
  #a[,1]=s$par
  if(iter>1){
      for(j in 2:iter){
        Q=qr.Q(qr(a[,1:(j-1)]),complete = TRUE)[,j:p]
        d=Q%*%optim(par=ginv(t(Q)%*%sigma%*%Q)%*%t(Q)%*%(mu0-mu1),fn=MisRate,method="BFGS",mu0=t(Q)%*%mu0,mu1=t(Q)%*%mu1,sigma0=t(Q)%*%sigma0%*%Q,sigma1=t(Q)%*%sigma1%*%Q)$par
        d=d/sqrt(sum(d^2))
        d[which(abs(d)<lambda)]=0
        a[,j]=d
      }
  }
  return(a)
}

QDAp_2 <- function(x,y,xnew,lambda=0,iter=1){
  x=data.matrix(x)
  xnew=data.matrix(xnew)
  x0=x[which(y==0),]
  x1=x[which(y==1),]
  mu0=colMeans(x0)
  mu1=colMeans(x1)
  sigma0=cov(x0)
  sigma1=cov(x1)
  sigma=((nrow(x0)-1)*sigma0+(nrow(x1)-1)*sigma1)/(nrow(x0)+nrow(x1)-1)
  a=drt_2(mu0,mu1,sigma0,sigma1,sigma,lambda,iter)
  qda.fit=qda(x%*%a,y)
  ynew=predict(qda.fit,xnew%*%a)$class
  return(ynew)
}
