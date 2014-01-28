# !! To do: Implement input checking (dimension/size, class), print warnings/errors if input has unexpected format

library(compiler)

##############################################
# crps.sim - compute CRPS of MCMC sample
# Input:
# - rlz, (scalar) realization
# - dat, MCMC forecast sample
# Output:
# - list of two elements: score - CRPS score, bw - bandwidth used for smoothing
# Notes:
# - Output is smoothed using a Gaussian Kernel and the Sheather and Jones (JRSSB, 1991) bandwidth choice as implemented in R
# - Sample size of dat should be at most (!) 10,000 to be computationally tractable
#   (larger samples can be thinned, which also reduces autocorrelation in the MCMC draws)
##############################################
crps.sim <- cmpfun(function(rlz,dat){
  auxcrps <- cmpfun(function(mu,s){
    return(2*s*dnorm(mu/s)+mu*(2*pnorm(mu/s)-1))
  })
  h <- bw.SJ(dat)
  mat1 <- 0.5*mean(outer(dat,dat,function(x,y) auxcrps(x-y,sqrt(2*h^2))))
  list(score = mat1 - mean(auxcrps(rlz-dat,h)), bw = h)
})

##############################################
# qs.sim - compute QS of MCMC sample
# same as crps.sim above, but for QS
##############################################
qs.sim <- cmpfun(function(rlz,dat){
  h <- bw.SJ(dat)
  mat1 <- outer(dat,dat,"-")
  mat1 <- mean(exp(-(mat1^2)/(4*h^2))*((4*pi*h^2)^(-0.5)))
  list(score=2*mean(dnorm((rlz-dat)/h))/h - mat1,bw=h)  
})

##############################################
# ls.sim - compute LS of MCMC sample
# same as crps.sim above, but for LS
##############################################
ls.sim <- cmpfun(function(rlz,dat){
  h <- bw.SJ(dat)
  list(score=log(mean(dnorm(rlz,dat,h))),bw=h)
})

##############################################
# main.sim - compute score of MCMC sample
# Input:
# - rlz, (scalar) realization
# - dat, MCMC forecast sample
# - sel, choice of scoring rule 1 = LS, 2 = QS, 3 = CRPS, 4 = squared error
# Output:
# - appropriate score
# Notes:
# - see crps.sim above for notes
##############################################
main.sim <- cmpfun(function(rlz,dat,sel){
if (sel == 1){
ls.sim(rlz,dat)
} else if (sel == 2){
qs.sim(rlz,dat)
} else if (sel == 3){
crps.sim(rlz,dat)
} else if (sel == 4){
list(score=-(rlz-mean(dat))^2,bw=NA)
}})

##############################################
# p2pnorm - compute cdf of 2-piece normal distribution
# Input:
# - x, real number/vector, position(s) at which to compute cdf
# - m, scalar, mode of distribution
# - s1, s2, positive scalars, standard deviation parameters
# Output:
# - real number/vector, cdf at positions given by x
##############################################
p2pnorm <- cmpfun(function(x,m,s1,s2){
  ifelse(x <= m,((2*s1)/(s1+s2))*pnorm((x-m)/s1),((s1-s2)/(s1+s2))+((2*s2)/(s1+s2))*pnorm((x-m)/s2))  
})

##############################################
# d2pnorm - compute density of 2-piece normal
# Input/output analogous to p2pnorm above, but for pdf
# Setting lg = TRUE yields log density
##############################################
d2pnorm <- cmpfun(function(x,m,s1,s2,lg=FALSE){
  a <- 2/(sqrt(2*pi)*(s1+s2))
  ret <- ifelse(x <= m,a*sqrt(2*pi)*s1*dnorm(x,m,s1),a*sqrt(2*pi)*s2*dnorm(x,m,s2))  
  if (lg==FALSE) ret else log(ret)
})

##############################################
# r2pnorm - Simulate from 2-piece normal
# Inputs m, s1, s2 analogous to p2pnorm above
# n gives the number of draws to be simulated
##############################################
r2pnorm <- cmpfun(function(n,m,s1,s2){
  u <- runif(n)
  aux <- function(z,u) (u-p2pnorm(z,m,s1,s2))^2
  apply(as.matrix(1:n),1,function(z)uniroot(function(x) p2pnorm(x,m,s1,s2)-u[z],lower=-100,upper=100)$root)
})

##############################################
# qs.2pnorm - Compute QS of 2-piece normal predictive density
# Inputs:
# - rlz, scalar, realizing value
# - m,s1,s2, distribution parameters (see p2pnorm above)
# Output:
# - QS (the greater the better)
##############################################
qs.2pnorm <- cmpfun(function(rlz,m,s1,s2){
  2*d2pnorm(rlz,m,s1,s2) - (s1+s2)/(sqrt(pi)*(s1+s2)^2)
})

##############################################
# crps.2pnorm - Compute CRPS of 2-piece normal predictive density
# Same as qs.2pnorm above, but for CRPS
##############################################
crps.2pnorm <- cmpfun(function(rlz,m,s1,s2){
  aux1 <- function(rlz,m,s1,s2){
    (4*(s1^2)/(s1+s2)) * (  ((rlz-m)/s1)*pnorm((rlz-m)/s1)+dnorm((rlz-m)/s1)      )    
  }
  aux2 <- function(m,s1,s2){
    (2/(sqrt(pi)*(s1+s2)^2)) * ( sqrt(2)*s2*(s2^2-s1^2)-(s1^3+s2^3) )    
  }
  sc <- ifelse(rlz <= m,aux1(rlz,m,s1,s2)-(rlz-m)+aux2(m,s1,s2),aux1(rlz,m,s2,s1)+(rlz-m)*( ((s1-s2)^2-4*s2^2) / ((s1+s2)^2) )+aux2(m,s2,s1))
  return(-sc)
})

##############################################
# main.2pnorm - compute score of 2-piece normal 
# Input:
# - rlz, scalar, realizing value
# - m,s1,s2, distribution parameters (see p2pnorm above)
# - sel, 1 = LS, 2 = QS, 3 = CRPS
# Output:
# - appropriate score
##############################################
main.2pnorm <- cmpfun(function(x,m,s1,s2,sel){
if (sel==1){ 
d2pnorm(x,m,s1,s2,lg=TRUE) 
} else if (sel == 2){
qs.2pnorm(x,m,s1,s2)
} else if (sel == 3){ 
crps.2pnorm(x,m,s1,s2)
} else if (sel == 4){ 
-(x-m)^2
}})

##############################################
# ls.mixn - Compute log score for a Gaussian mixture
# Same as `main' above, except that scoring rule is fixed to log score
##############################################
ls.mixn <- cmpfun(function(y,w,m,s){
return(log(sum(w*dnorm(y,m,s))))
})	

##############################################
# qs.mixn - Compute QS for Gaussian mixture density
# Same as main.mixn above, but fixed to QS
##############################################
qs.mixn <- cmpfun(function(y,w,m,s){
mat1 <- outer(m,m,function(x,y) (x-y)^2)
mat2 <- outer(s,s,function(x,y) (x^2+y^2))
aux2 <- exp(-0.5*mat1/mat2)*(1/sqrt(mat2))/sqrt(2*pi)
w <- matrix(w,length(w),1)
aux2 <- t(w) %*% aux2 %*% w
aux1 <- sum(w*dnorm(y,m,s))
return(as.numeric(2*aux1-aux2))
})

##############################################
# auxcrps - Auxiliary function for CRPS computations
##############################################
auxcrps <- cmpfun(function(mu,s){
return(2*s*dnorm(mu/s)+mu*(2*pnorm(mu/s)-1))
})

##############################################
# crps.mixn - Compute CRPS for Gaussian mixture density
# Same as main.mixn above, but fixed to CRPS
##############################################
crps.mixn <- cmpfun(function(y,w,m,s){
mat1 <- outer(m,m,function(x,y) x-y)
mat2 <- outer(s,s,function(x,y) sqrt(x^2+y^2))
w <- matrix(w,length(w),1)
aux2 <- t(w) %*% auxcrps(mat1,mat2) %*% w
return(as.numeric(-sum(w*auxcrps(y-m,s))+0.5*aux2))
})

##############################################
# serr.mixn - Compute negSE for a Gaussian mixture
# Same as `main' above, except that scoring rule is fixed to negSE
##############################################
serr.mixn <- cmpfun(function(y,w,m,s){
  return(-(y-sum(w*m))^2 )
})

##############################################
# main.mixn - Compute score for a Gaussian mixture density
# Inputs: 
# y, scalar (!) realization
# w, n x 1 vector of weights (nonnegative, sum to one)
# m, n x 1 vector of means
# s, n x 1 vector of standard deviations
# sel, scalar, which scoring rule to use (1 - log score, 2 - quadratic score, 3 - CRPS, 4 - Squared Error)
# Output:
# Scalar, value of scoring rule (the greater the better)
##############################################
main.mixn <- cmpfun(function(y,w,m,s,sel){
if (sel==1){
 return(ls.mixn(y,w,m,s))
} else if (sel==2){
return(qs.mixn(y,w,m,s))
} else if (sel==3){
return(crps.mixn(y,w,m,s))
} else if (sel==4){
return(serr.mixn(y,w,m,s))
}})

##############################################
# linms - Compute mean and standard deviation of linear prediction pool
# Inputs:
# w, (n x 1) vector of weights (positive, sum to one)
# m, (T x n) matrix of forecast means (rows - time periods, cols - individual methods)
# s, (T x n) matrix of forecast standard deviations
# Output:
# list with two elements for mean and sd of linear pool
# 1 - (T x 1) vector of means, 2 - (T x 1) vector of sds
##############################################
linms <- cmpfun(function(w,m,s){
T <- dim(m)[1]
n <- dim(m)[2]
w <- t(matrix(w,n,T))
mnew <- apply(m*w,1,sum)
snew <- sqrt(apply(w*(s^2),1,sum) + apply(w*((m-mnew)^2),1,sum))
return(list(m=mnew,s=snew))  
})

##############################################
# logms - Compute mean and standard deviation of logarithmic prediction pool
# Inputs/outputs analogous to ``linms'' above
# !!! Note:
# Output formula assumes that individual forecast densities are Gaussian !!!
##############################################
logms <- cmpfun(function(w,m,s){
T <- dim(m)[1]
n <- dim(m)[2]
w <- t(matrix(w,n,T))
a <- w/(s^2)
mnew <- apply(a*m,1,sum)/apply(a,1,sum)
snew <- sqrt(1/apply(a,1,sum))  
return(list(m=mnew,s=snew))
})

##############################################
# mlog - Multivariate inverse logistic transformation
# Input:  
# x, m-variate vector (real-valued)
# i-th element is log( prob(baseline)/prob(event i) )
# Output
# (m+1) variate vector of probabilities, baseline in last row (m+1)
##############################################
mlog <- cmpfun(function(x) c(1,exp(x))/(1+sum(exp(x))))

##############################################
# optweights.mixn - Estimate optimal weights in Gaussian pool, based on training sample
# Inputs:
# y, T x 1 vector of realizations (training sample)
# m, T x n matrix of means
# s, T x n matrix of standard deviations
# sel, scalar, which scoring rule to use (1 - log score (default), 2 - quadratic score, 3 - CRPS, 4 - Squared Error)
# type, either "linear" (default) or "log", type of desired pooling method
# om, string, optimization method (defaults to "BFGS", see 'optim' docu for other options)
# tol, scalar, convergence tolerance in optimization (defaults to 1e-7; only relevant if n > 2)
# Output:
# List of two elements: 1 - vector of optimal weights, 2 - average score in optimum
# Notes:
# - Function uses optim w/ inverse log transform of inputs to avoid need for constraints
# - Works well in most cases, but can give problems for optimal weights close to zero or one;
#   try to (carefully!) increase `tol' in this case 
##############################################
optweights.mixn <- cmpfun(function(y,m,s,sel=1,type="linear",om = "BFGS", tol=1e-7){
T <- dim(m)[1]
n <- dim(m)[2]
# Define objective function
if (type == "linear"){
of <- function(theta){
w <- mlog(theta)
aux <- apply(as.matrix(1:T),1,function(z) main.mixn(y=y[z],w=w,m=m[z,],s=s[z,],sel=sel))
return(-sum(aux))
}
} else {
of <- function(theta){
aux <- logms(mlog(theta),m,s)
mnew <- aux$m
snew <- aux$s
aux <- apply(as.matrix(1:T),1,function(z) main.mixn(y=y[z],w=1,m=mnew[z],s=snew[z],sel=sel))
return(-sum(aux))
}
}
# Optimize objective function
if (n > 2){
  wstar <- optim(rep(0,n-1),of,method="BFGS",control=list(reltol=tol))
  return(list(wstar=mlog(wstar$par),score=-wstar$value/T))
} else {
  wstar <- optimize(of,c(-25,25))
  return(list(wstar=mlog(wstar$minimum),score=-wstar$objective/T))
}
})