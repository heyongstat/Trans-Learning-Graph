### Functionality: Gaussian CDF transformation 
GaussianCdfTransform = function(z,mu, sigma){
  d = ncol(z)
  # define and standardize g = f^{-1}
  g = NULL
  bFunc = NULL
  for(k in c(1:d)){
    g = c(g, function(x) x)
    bFunc = c(bFunc, function(x) x)
  }
  gMean = rep(0,d)
  gRatio= rep(0,d)  
  for (i in 1:d) {
    grid = seq(mu[i]-5*sigma[i,i],mu[i]+5*sigma[i,i],length=1000)
    delta = grid[2] - grid[1]
    w = dnorm(grid,mu[i],sigma[i,i])
    bFunc[[i]] = function(x) {pnorm(x,0.05,0.4)}
    gMean[i] = delta*sum(w*bFunc[[i]](grid))
    gRatio[i]= sigma[i,i] /sqrt( delta*sum((w*(bFunc[[i]](grid) - gMean[i])^2)))
    g[[i]] = function(x) {(bFunc[[i]](x)- gMean[i] )*gRatio[i]   + mu[i]}      
  }
  # conduct the transformation
  x = 0*z
  for (i in 1:d) {x[,i] = g[[i]](z[,i])}
  return(x)
} 