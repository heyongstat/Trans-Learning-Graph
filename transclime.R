library(mvtnorm)
library(fastclime)
library(Matrix)
library(lavaSearch2) #only to symmetrize matrices
source("ADMM_proj.R")#admm for spd projection
source("GaussianCdfTransform.R")

#generate Sigma, Sigma^{1},...,Sigma^{K} of the data for simulation
rdmatrix=function(K, A.size, h, p=100, type='Toep'){
  if(type=='Toep'){
    Theta<-toeplitz(0.6^(1:p)*2)
    Theta[which(abs(Theta)<=0.05, arr.ind=T)]<- 0
  }else if(type=='Bdiag'){
    Theta<-kronecker(diag(p/4), toeplitz(c(1.2,0.9,0.6,0.3)))
  }
  Sig<- solve(Theta)
  d=sqrt(diag(Sig))
  d=diag(d)
  Theta=d%*%Theta%*%d #inverse correlation matrix#
  Sig<- cov2cor(Sig) #correlation matrix#
  Omega.vec<-0
  Delta.k=array(0,dim=c(p,p,A.size))
  Theta.k=array(0,dim = c(p,p,(K-A.size)))
  Sig.k=array(0,dim = c(p,p,K))
  for(k in 1 : K){
    if(k<=A.size){
      Delta.k[,,k]<-matrix(rbinom(p^2,size=1,prob=0.1)*runif(p^2,-h/p,h/p),ncol=p)
      Sig.k[,,k]<-(diag(1,p)+Delta.k[,,k])%*%Sig 
      Sig.k[,,k]<-lavaSearch2:::symmetrize(Sig.k[,,k], update.upper = TRUE)
      if(min(eigen(Sig.k[,,k])$values)<0.05){
        Sig.k[,,k]<-Spd.proj(Sig.k[,,k])$mat
      }
      Sig.k[,,k] <- cov2cor(Sig.k[,,k])
      Omega.vec<-c(Omega.vec, max(colSums(abs(diag(1,p)-Sig.k[,,k]%*%Theta))))
    }
    else{
      Theta.k[,,(k-A.size)]<-diag(1.5,p)+matrix(rbinom(p^2,size=1,prob=0.1)*0.2,ncol=p)
      Theta.out<-lavaSearch2:::symmetrize(Theta.k[,,(k-A.size)], update.upper = TRUE)
      if(min(eigen(Theta.out)$values)<0.05){
        Theta.out<-Spd.proj(Theta.out)$mat
      }
      Sig.k[,,k]=solve(Theta.out)
      Sig.k[,,k] <- cov2cor(Sig.k[,,k])
    }
  }
  list(Sig.k=Sig.k,Sig=Sig,Omega.l1=max(Omega.vec),Theta0=Theta) 
}

#generate the data for simulation
DataGen<-function(K, n.vec, Sig.k, Sig, p=100, tf.type="cdf"){
  if (tf.type == "cdf"){
    X<-rmvnorm(n.vec[1],rep(0,p), sigma=Sig)
    X = GaussianCdfTransform(z=X,mu=rep(0,p), sigma=Sig)
    for(k in 1 : K){
      Y = rmvnorm(n.vec[k+1],rep(0,p), sigma=Sig.k[,,k])
      Y = GaussianCdfTransform(z=Y,mu=rep(0,p), sigma=Sig.k[,,k])
      X<- rbind(X,Y)
    }
  } else if (tf.type == "exp"){
    X<-rmvnorm(n.vec[1],rep(0,p), sigma=Sig)
    X<-exp(X)
    for(k in 1 : K){
      Y = rmvnorm(n.vec[k+1],rep(0,p), sigma=Sig.k[,,k])
      Y = exp(Y)
      X<- rbind(X,Y)
    }
  } else{
    X<-rmvnorm(n.vec[1],rep(0,p), sigma=Sig)
    for(k in 1 : K){
      Y = rmvnorm(n.vec[k+1],rep(0,p), sigma=Sig.k[,,k])
      X<- rbind(X,Y)
    }
  }
  list(X=X)
}

###algorithm based on fastclime package
### min \|Theta\|_1
###subject to \|X^TX%*%Theta-Bmat\|_max <=lambda (if X is the raw sample)
###subject to \|X%*%Theta-Bmat\|_max <=lambda (if X is the sample covariance matrix)
Myfastclime.s<-function(X,Y,Bmat,lambda=0.1, scale=T){
  p<-ncol(X)
  obj=rep(-1,2*p)
  obj_bar=rep(0,2*p)
  rhs_bar<-rep(1, 2*p)
  
  Sig.hat<-X
  
  feasible=T
  Theta.hat<-NULL
  mat=rbind(cbind(Sig.hat,-Sig.hat),
            cbind(-Sig.hat,Sig.hat))
  
  for(j in 1:p){
    rhs <- c(Bmat[,j],-Bmat[,j])
    out.txt<-capture.output(  fastlp.re<-fastlp(obj=obj, mat=mat, rhs=rhs+rhs_bar*lambda))
    if(!grepl("optimal", out.txt)){
      feasible=F
      break
    }
    Theta.hat<-cbind(Theta.hat,(fastlp.re[1:p]-fastlp.re[-(1:p)]))
    if(all(Theta.hat[,j]==0)) scale2=F else scale2=T
    if(scale & scale2){
      Theta.hat[,j]<-
        as.numeric(Theta.hat[j,j]/ (t(Theta.hat[,j])%*%Y%*%Theta.hat[,j]))*Theta.hat[,j]
    }
  }
  if(!feasible){
    cat('Theta.hat not found','\n')
    Theta.hat<-diag(1,p)
  }
  list(Theta.hat=Theta.hat, conv=feasible)
}

## compute the semi-positive definite projection of SigA.hat 
## with smallest eigenvalue lower bounded by eps
Spd.proj<-function(SigA.hat, eps=NULL){
  p=ncol(SigA.hat)
  if(is.null(eps)){
    eps<-5/200
  }
  feasible=1
  SigA.t<-SigA.hat
  if(min(eigen(SigA.t)$values) <=eps ){
    feasible=2
    SigA.t<-ADMM_proj(SigA.hat, epsilon=eps)$mat
  }
  SigA.t<-lavaSearch2:::symmetrize(SigA.t, update.upper = TRUE)
  list(mat=SigA.t, conv=feasible)
}

####the main Copula Trans-CLIME algorithm###
##corX: sample correlation matrix of 2/3 primary data;
##n0: 2n/3, n is the sample size of primary data;
##corX1: sample correlation matrix of primary data;
##corX.til: sample correlation matrix of the rest 1/3 primary data;
##corX.A: {corX^{(k)}, k in A}; lambda:lambda.Omega; 
##conv: whether calculation of CLIME estimator is feasible or not;
### agg: perform LS aggregation or not;  Theta.cl: CLIME estimator
Trans.CLIME<-function(corX, corX.A, corX1, corX.til, lambda, conv, Theta.cl=NULL, n0, agg=T){
  p<-ncol(corX) #X,Y,Bmat,lambda=0.1, scale=T
  if(conv){
    lam.delta<-2*sqrt(log(p)/n0) 
    Delta.re <- Myfastclime.s(X=corX, Y=NULL, Bmat=corX-corX.A, lambda=lam.delta, scale=F)
    if(Delta.re$conv){ 
      Delta.init<-Delta.re$Theta.hat
      Theta.db<-Delta.init+ Theta.cl%*%(corX-corX.A)
      Delta.re <- Myfastclime.s(X=diag(1,p), Y=NULL, Bmat=Theta.db, lambda=2*lam.delta, scale=F)
      if(Delta.re$conv){
        Delta.hat<-Delta.re$Theta.hat
        Theta.re <- Myfastclime.s(X=corX.A, Y=corX1,Bmat=diag(1,p)-t(Delta.hat), lambda=lambda)
        if(Theta.re$conv){
          if(agg){
            Omega.hat<-Agg(Theta.init=cbind(Theta.cl, Theta.re$Theta.hat), Sigtilde=corX.til)
          }
          else{
            Omega.hat<-Theta.tl.re$Theta.hat
          }
        }else{
          Omega.hat=diag(1,p);conv=F
        }
      }else{
        Omega.hat=diag(1,p);conv=F
      }
    }else{
      Omega.hat=diag(1,p);conv=F
    }
  }else{
    Omega.hat=diag(1,p);conv=F
  }
  list(Omega.hat=Omega.hat,conv=conv)
}

####LS aggregation function with 
### Theta.init=(Omega.clime, Omega.hat);
### Sigtilde: sample correlation matrix of some primary samples
Agg<-function(Theta.init, Sigtilde){
  p<-ncol(Sigtilde)
  v.mat<-sapply(1:p, function(j){
    W.j<-diag(0,2)
    W.j[1,1]=t(Theta.init[,j])%*%Sigtilde%*%Theta.init[,j]
    W.j[1,2]=t(Theta.init[,j])%*%Sigtilde%*%Theta.init[,p+j]
    W.j[2,1]=W.j[1,2]
    W.j[2,2]=t(Theta.init[,p+j])%*%Sigtilde%*%Theta.init[,p+j]
    if(eigen(W.j)$values[2]<=10^(-6)){
      matrix(c(0,1),nrow=2,ncol=1)
    }else{
      solve(W.j)%*%c(Theta.init[j,j], Theta.init[j,p+j])
    }
  })
  
  v.mat<-apply(v.mat,2, function(x) {
    if(any(x>=0))
      x*(x>=0)/(sum(x*(x>=0)))
    else
      x
  })
  Theta.hat<-sapply(1:p, function(j) cbind(Theta.init[,j], Theta.init[,p+j])%*% v.mat[,j])
  
  Theta.hat
}

###compute the estimation errors in Frobenius norm, max norm and spectral norm,
###false negative rate, and false positive rate.
Dist<- function(Theta.hat,Theta0){
  p<-ncol(Theta.hat)
  Theta.hat<-lavaSearch2:::symmetrize(Theta.hat, update.upper = TRUE)
  Frob=sum((Theta.hat-Theta0)^2)/p
  maxnorm= max(abs(Theta.hat-Theta0))
  S=max(abs(svd(Theta.hat-Theta0)$d))^2
  a <- array(dim=(p^2))
  for (b in 1:(p^2)) {
    if (Theta.hat[b]==0 & Theta0[b]==0) #TN
    {
      a[b]<-1
    } 
    else if(Theta.hat[b]!=0 & Theta0[b]!=0){  #TP
      a[b]<-2
    }
    else if(Theta.hat[b]==0 & Theta0[b]!=0){   #FN
      a[b]<-3
    }
    else if(Theta.hat[b]!=0 & Theta0[b]==0){  #FP
      a[b]<-4
    }
    else{
      a[b]<-0
    }
  }
  FNR = sum(a==3)/sum(Theta0!=0)
  FPR = sum(a==4)/sum(Theta0==0)
  list(Frob=Frob, maxnorm= maxnorm, S=S, FNR=FNR, FPR=FPR)
}

###cross validation for selecting tuning parameters
###X: primary data; cor.n: sample correlation matrix of primary data;
###copula: cv for Copula Trans-CLIME or Trans-CLIME
cv.clime<-function(p,X,cor.n, nfold=5,copula=T){
  library(caret)
  folds<-createFolds(1:nrow(X), k=nfold)
  te<-NULL
  lam.seq<-seq(0,3.5,length.out=176)*2*sqrt(log(p)/nrow(X)*nfold/(nfold-1))
  for(i in 1:nfold){
    if(copula){
      cor.cv1=sin(pi/2*cor(X[-folds[[i]],],method="kendall"))
      cor.cv2=sin(pi/2*cor(X[folds[[i]],],method="kendall"))
    }
    else{
      cor.cv1=cor(X[-folds[[i]],])
      cor.cv2=cor(X[folds[[i]],])
    }
    te<-rbind(te,sapply(lam.seq, function(lam){
      cur.clime<-Myfastclime.s(X=cor.cv1,Y=cor.n,Bmat=diag(1,p), lambda=lam)$Theta.hat
      cur.clime<-lavaSearch2:::symmetrize(cur.clime, update.upper = TRUE)
      Theta<-Spd.proj(SigA.hat=cur.clime, eps=0.001)$mat
      eigens<-eigen(Theta)$values
      sum(diag(cor.cv2%*%Theta))/2-sum(log(eigens[eigens>0]))/2})
    )
  }
  te.ave<-colMeans(te)
  te.min<-which.min(te.ave)
  cat(te.ave,'\n')
  if(te.ave[te.min]==te.ave[1] | te.ave[te.min]==te.ave[176]){
    lam<-seq(0,3.5,length.out=176)[which(diff(sign(diff(te)))==2)+1]
    if(length(lam)==0){
      lam<-seq(0,3.5,length.out=176)[te.min]
    }
  }else{
    lam<-seq(0,3.5,length.out=176)[te.min]
  }
  lam
}