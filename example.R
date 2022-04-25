source('copulatransclime.R')
###an example for ROC curve
set.seed(213)
p=100
K=5
A.size=3
DeThe=rdmatrix(K, A.size, h=10, p, type='Bdiag')
#list(Sig.k=Sig.k,Sig=Sig,Omega.l1=max(Omega.vec),Theta0=Theta) 
Theta0<-DeThe$Theta0 ##true Omega
n.vec<-c(200,rep(200,K))
n.cum=cumsum(n.vec[-1])
copindi=array(0,dim=c(528,5,100))
indi=array(0,dim=c(528,5,100))
constall=seq(0,3.5,length.out=176)
for(i in 1:100){
  dat.all<-DataGen(K=K, p=p, n.vec=n.vec, Sig.k=DeThe$Sig.k, Sig=DeThe$Sig, tf.type="cdf")
  X.all<-dat.all$X
  n0=round(n.vec[1]*2/3)
  X0=X.all[1:n0,]
  X.til=X.all[(n0+1):n.vec[1],]
  X.A=X.all[-(1:n.vec[1]),]
  
  ##calculate rank-based sample correlation matrix
  corX.n0_k = sin(pi/2*cor(X0,method="kendall"))
  corX.n1_k=sin(pi/2*cor(X.all[1:n.vec[1],],method="kendall"))
  corX.til_k=sin(pi/2*cor(X.til,method="kendall"))
  nA<-sum(n.vec[2:(A.size+1)])
  corX.A_k=n.vec[2]/nA*sin(pi/2*cor(X.A[1:n.vec[2],],method="kendall"))
  for (l in 2:A.size) {
    corX.A_k=corX.A_k+n.vec[l+1]/nA*sin(pi/2*cor(X.A[(n.cum[l-1]+1):n.cum[l],],method="kendall"))
  }
  nB<-sum(n.vec[2:(K+1)])
  corX.Apool_k=n.vec[2]/nB*sin(pi/2*cor(X.A[1:n.vec[2],],method="kendall"))
  for (l in 2:K) {
    corX.Apool_k=corX.Apool_k+n.vec[l+1]/nB*sin(pi/2*cor(X.A[(n.cum[l-1]+1):n.cum[l],],method="kendall"))
  }
  
  ##calculate the sample correlation matrix
  corX.n0 = cor(X0)
  corX.n1=cor(X.all[1:n.vec[1],])
  corX.til=cor(X.til)
  corX.A=n.vec[2]/nA*cor(X.A[1:n.vec[2],])
  for (l in 2:A.size) {
    corX.A=corX.A+n.vec[l+1]/nA*cor(X.A[(n.cum[l-1]+1):n.cum[l],])
  }
  corX.Apool=n.vec[2]/nB*cor(X.A[1:n.vec[2],])
  for (l in 2:K) {
    corX.Apool=corX.Apool+n.vec[l+1]/nB*cor(X.A[(n.cum[l-1]+1):n.cum[l],])
  }
  
  j=0
  for (const in constall){
    ###original copula clime####
    copTheta.clime<-Myfastclime.s(X=corX.n1_k,Y=corX.n1_k,Bmat=diag(1,p), lambda=const*2*sqrt(log(p)/n.vec[1]))
    copTheta.hatclime0<-copTheta.clime$Theta.hat
    if(copTheta.clime$conv){
      copindi[(3*j+1), ,i]=unlist(Dist(copTheta.hatclime0, Theta0))
    }else{copindi[(3*j+1), ,i]=NA}
    
    ###copula Trans-CLIME algorithm
    copTheta.re0<-Myfastclime.s(X=corX.n0_k,Y=corX.n1_k,Bmat=diag(1,p), lambda=const*2*sqrt(log(p)/n0))
    copTheta.hat0<-copTheta.re0$Theta.hat
    #Trans.CLIME<-function(corX, corX.A, corX1, corX.til, lambda, conv, Theta.cl=NULL, n0, agg=T)
    copOmega.re<-Trans.CLIME(corX=corX.n0_k, corX.A=corX.A_k, corX1=corX.n1_k,
                             corX.til=corX.til_k,n0=n0,lambda=const*2*sqrt(log(p)/nA), 
                             conv=copTheta.re0$conv, Theta.cl=copTheta.hat0)
    if(copOmega.re$conv){
      copindi[(3*j+2), ,i]=unlist(Dist(copOmega.re$Omega.hat, Theta0))
    }else{copindi[(3*j+2), ,i]=NA}
    #pooled copula Trans-CLIME
    copPooled.clime <- Trans.CLIME(corX=corX.n0_k, corX.A=corX.Apool_k, corX1=corX.n1_k,
                                   corX.til=corX.til_k,n0=n0,lambda=const*2*sqrt(log(p)/nB), 
                                   conv=copTheta.re0$conv, Theta.cl=copTheta.hat0)
    if(copPooled.clime$conv){
      copindi[(3*j+3), ,i]=unlist(Dist(copPooled.clime$Omega.hat, Theta0))
    }else{copindi[(3*j+3), ,i]=NA}
    
    ###original clime####
    Theta.clime<-Myfastclime.s(X=corX.n1,Y=corX.n1,Bmat=diag(1,p), lambda=const*2*sqrt(log(p)/n.vec[1]))
    Theta.hatclime0<-Theta.clime$Theta.hat
    if(Theta.clime$conv){
      indi[(3*j+1), ,i]=unlist(Dist(Theta.hatclime0, Theta0))
    }else{indi[(3*j+1), ,i]=NA}
    ###Trans-CLIME algorithm
    Theta.re0<-Myfastclime.s(X=corX.n0,Y=corX.n1,Bmat=diag(1,p), lambda=const*2*sqrt(log(p)/n0))
    Theta.hat0<-Theta.re0$Theta.hat
    Omega.re<-Trans.CLIME(corX=corX.n0, corX.A=corX.A, corX1=corX.n1,
                          corX.til=corX.til,n0=n0,lambda=const*2*sqrt(log(p)/nA), 
                          conv=Theta.re0$conv, Theta.cl=Theta.hat0)
    if(Omega.re$conv){
      indi[(3*j+2), ,i]=unlist(Dist(Omega.re$Omega.hat, Theta0))
    }else{indi[(3*j+2), ,i]=NA}
    #pooled Trans-CLIME
    Pooled.clime<-Trans.CLIME(corX=corX.n0, corX.A=corX.Apool, corX1=corX.n1,
                              corX.til=corX.til,n0=n0,lambda=const*2*sqrt(log(p)/nB), 
                              conv=Theta.re0$conv, Theta.cl=Theta.hat0)
    if(Pooled.clime$conv){
      indi[(3*j+3), ,i]=unlist(Dist(Pooled.clime$Omega.hat, Theta0))
    }else{indi[(3*j+3), ,i]=NA}
    j=j+1
  }
}

#save(list = ls(all.names = TRUE), file = "copulatransclime-cdfG-B-h10-cor.RData")

###an example using cv
set.seed(213)
p=100
K=5
A.size=3
DeThe=rdmatrix(K, A.size, h=10, p, type='Bdiag')
#list(Sig.k=Sig.k,Sig=Sig,Omega.l1=max(Omega.vec),Theta0=Theta) 
Theta0<-DeThe$Theta0 ##true Omega
n.vec<-c(200,rep(200,K))
n.cum=cumsum(n.vec[-1])

dat.all<-DataGen(K=K, p=p, n.vec=n.vec, Sig.k=DeThe$Sig.k, Sig=DeThe$Sig, tf.type="cdf")
X.all<-dat.all$X
n0=round(n.vec[1]*2/3)
X0=X.all[1:n0,]
X.til=X.all[(n0+1):n.vec[1],]
X.A=X.all[-(1:n.vec[1]),]

corX.n0_k = sin(pi/2*cor(X0,method="kendall"))
corX.n1_k=sin(pi/2*cor(X.all[1:n.vec[1],],method="kendall"))
corX.til_k=sin(pi/2*cor(X.til,method="kendall"))
nA<-sum(n.vec[2:(A.size+1)])
corX.A_k=n.vec[2]/nA*sin(pi/2*cor(X.A[1:n.vec[2],],method="kendall"))
for (l in 2:A.size) {
  corX.A_k=corX.A_k+n.vec[l+1]/nA*sin(pi/2*cor(X.A[(n.cum[l-1]+1):n.cum[l],],method="kendall"))
}
nB<-sum(n.vec[2:(K+1)])
corX.Apool_k=n.vec[2]/nB*sin(pi/2*cor(X.A[1:n.vec[2],],method="kendall"))
for (l in 2:K) {
  corX.Apool_k=corX.Apool_k+n.vec[l+1]/nB*sin(pi/2*cor(X.A[(n.cum[l-1]+1):n.cum[l],],method="kendall"))
}

corX.n0 = cor(X0)
corX.n1=cor(X.all[1:n.vec[1],])
corX.til=cor(X.til)
corX.A=n.vec[2]/nA*cor(X.A[1:n.vec[2],])
for (l in 2:A.size) {
  corX.A=corX.A+n.vec[l+1]/nA*cor(X.A[(n.cum[l-1]+1):n.cum[l],])
}
corX.Apool=n.vec[2]/nB*cor(X.A[1:n.vec[2],])
for (l in 2:K) {
  corX.Apool=corX.Apool+n.vec[l+1]/nB*cor(X.A[(n.cum[l-1]+1):n.cum[l],])
}

#cv.clime<-function(p,X,cor.n, nfold=5,copula=T)
copconst<-cv.clime(p=p,X=X.all[1:n.vec[1],],cor.n=corX.n1_k, nfold=5,copula=T)

###original copula clime####
copTheta.clime<-Myfastclime.s(X=corX.n1_k,Y=corX.n1_k,Bmat=diag(1,p), lambda=copconst*2*sqrt(log(p)/n.vec[1]))
copTheta.hatclime0<-copTheta.clime$Theta.hat
if(copTheta.clime$conv){
  unlist(Dist(copTheta.hatclime0, Theta0))
}else{print("NA")}

###copula Trans-CLIME algorithm
copTheta.re0<-Myfastclime.s(X=corX.n0_k,Y=corX.n1_k,Bmat=diag(1,p), lambda=copconst*2*sqrt(log(p)/n0))
copTheta.hat0<-copTheta.re0$Theta.hat
#Trans.CLIME<-function(corX, corX.A, corX1, lambda, conv, Theta.cl=NULL, n0, agg=T,corX.til){
copOmega.re<-Trans.CLIME(corX=corX.n0_k, corX.A=corX.A_k, corX1=corX.n1_k,
                         corX.til=corX.til_k,n0=n0,lambda=copconst*2*sqrt(log(p)/nA), 
                         conv=copTheta.re0$conv, Theta.cl=copTheta.hat0)
if(copOmega.re$conv){
  unlist(Dist(copOmega.re$Omega.hat, Theta0))
}else{print("NA")}
#pooled copula Trans-CLIME
copPooled.clime <- Trans.CLIME(corX=corX.n0_k, corX.A=corX.Apool_k, corX1=corX.n1_k,
                               corX.til=corX.til_k,n0=n0,lambda=copconst*2*sqrt(log(p)/nB), 
                               conv=copTheta.re0$conv, Theta.cl=copTheta.hat0)
if(copPooled.clime$conv){
  unlist(Dist(copPooled.clime$Omega.hat, Theta0))
}else {print("NA")}

const<-cv.clime(p=p,X=X.all[1:n.vec[1],],cor.n=corX.n1, nfold=5,copula=F)
##original clime
Theta.clime<-Myfastclime.s(X=corX.n1,Y=corX.n1,Bmat=diag(1,p), lambda=const*2*sqrt(log(p)/n.vec[1]))
Theta.hatclime0<-Theta.clime$Theta.hat
if(Theta.clime$conv){
  unlist(Dist(Theta.hatclime0, Theta0))
}else {print("NA")}

###Trans-CLIME algorithm
Theta.re0<-Myfastclime.s(X=corX.n0,Y=corX.n1,Bmat=diag(1,p), lambda=const*2*sqrt(log(p)/n0))
Theta.hat0<-Theta.re0$Theta.hat
Omega.re<-Trans.CLIME(corX=corX.n0, corX.A=corX.A, corX1=corX.n1,
                      corX.til=corX.til,n0=n0,lambda=const*2*sqrt(log(p)/nA), 
                      conv=Theta.re0$conv, Theta.cl=Theta.hat0)
if(Omega.re$conv){
  unlist(Dist(Omega.re$Omega.hat, Theta0))
}else {print("NA")}
#pooled Trans-CLIME
Pooled.clime<-Trans.CLIME(corX=corX.n0, corX.A=corX.Apool, corX1=corX.n1,
                          corX.til=corX.til,n0=n0,lambda=const*2*sqrt(log(p)/nB), 
                          conv=Theta.re0$conv, Theta.cl=Theta.hat0)
if(Pooled.clime$conv){
  unlist(Dist(Pooled.clime$Omega.hat, Theta0))
}else{print("NA")}