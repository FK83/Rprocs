##############################################
# optweights - Estimate optimal weights in Gaussian pool, based on training sample
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
optweights <- function(y,m,s,sel=1,type="linear",om = "BFGS", tol=1e-7){
T <- dim(m)[1]
n <- dim(m)[2]
# Define objective function
if (type == "linear"){
of <- function(theta){
w <- mlog(theta)
aux <- apply(as.matrix(1:T),1,function(z) main(sel=sel,y=y[z],w=w,m=m[z,],s=s[z,]))
return(-sum(aux))
}
} else {
of <- function(theta){
aux <- logms(mlog(theta),m,s)
mnew <- aux$m
snew <- aux$s
aux <- apply(as.matrix(1:T),1,function(z) main(sel=sel,y=y[z],w=1,m=mnew[z],s=snew[z]))
return(-sum(aux))
}
}
if (n > 2){
  wstar <- optim(rep(0,n-1),of,method="BFGS",control=list(reltol=tol))
  return(list(wstar=mlog(wstar$par),score=-wstar$value/T))
} else {
  wstar <- optimize(of,c(-25,25))
  return(list(wstar=mlog(wstar$minimum),score=-wstar$objective/T))
}
}

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
linms <- function(w,m,s){
T <- dim(m)[1]
n <- dim(m)[2]
w <- t(matrix(w,n,T))
mnew <- apply(m*w,1,sum)
snew <- sqrt(apply(w*(s^2),1,sum) + apply(w*((m-mnew)^2),1,sum))
return(list(m=mnew,s=snew))  
}

##############################################
# logms - Compute mean and standard deviation of logarithmic prediction pool
# Inputs/outputs analogous to ``linms'' above
# !!! Note:
# Output formula assumes that individual forecast densities are Gaussian !!!
##############################################
logms <- function(w,m,s){
T <- dim(m)[1]
n <- dim(m)[2]
w <- t(matrix(w,n,T))
a <- w/(s^2)
mnew <- apply(a*m,1,sum)/apply(a,1,sum)
snew <- sqrt(1/apply(a,1,sum))  
return(list(m=mnew,s=snew))
}

##############################################
# main - Compute score for a Gaussian mixture density
# Inputs: 
# y, scalar (!) realization
# w, n x 1 vector of weights (nonnegative, sum to one)
# m, n x 1 vector of means
# s, n x 1 vector of standard deviations
# sel, scalar, which scoring rule to use (1 - log score, 2 - quadratic score, 3 - CRPS, 4 - Squared Error)
# Output:
# Scalar, value of scoring rule (the greater the better)
##############################################
main <- function(sel,y,w,m,s){
if (sel==1) return(logsc(y,w,m,s))
if (sel==2) return(qs(y,w,m,s))
if (sel==3) return(crps(y,w,m,s))
if (sel==4) return(serr(y,w,m,s))
}

##############################################
# qs - Compute quadratic score for a Gaussian mixture
# Same as `main' above, except that scoring rule is fixed to quadratic score
##############################################
qs <- function(y,w,m,s){
n <- length(w)
aux1 <- 0
aux2 <- 0
for (i in 1:n){
	aux1 <- aux1 + 2*w[i]*dnorm(y,mean=m[i],sd=s[i])
	for (j in 1:n){
		a <- (m[i]-m[j])/s[j]
		b <- s[i]/s[j]
		aux2 <- aux2 + (w[i]*w[j]/s[j])*((2*pi*(1+b^2))^(-0.5))*exp(-(0.5*a^2)*(1/(1+b^2)))	
	}
}
return(aux1-aux2)
}

##############################################
# auxcrps - Auxiliary function for CRPS computations
##############################################
auxcrps <- function(mu,s){
return(2*s*dnorm(mu/s)+mu*(2*pnorm(mu/s)-1))
}

##############################################
# crps - Compute CRPS for a Gaussian mixture
# Same as `main' above, except that scoring rule is fixed to CRPS
##############################################
crps <- function(y,w,m,s){
n <- length(w)
aux1 <- 0
aux2 <- 0
for (i in 1:n){
	aux1 <- aux1 + w[i]*auxcrps(y-m[i],s[i])
	for (j in 1:n){
		aux2 <- aux2 + w[i]*w[j]*auxcrps(m[i]-m[j],sqrt(s[i]^2+s[j]^2))
	}
}
return(-aux1+0.5*aux2)
}


##############################################
# logsc - Compute log score for a Gaussian mixture
# Same as `main' above, except that scoring rule is fixed to log score
##############################################
logsc <- function(y,w,m,s){
n <- length(w)
aux1 <- 0
for (i in 1:n){
	aux1 <- aux1 + w[i]*dnorm(y,mean=m[i],sd=s[i])
}
return(log(aux1))
}	

##############################################
# serr - Compute negSE for a Gaussian mixture
# Same as `main' above, except that scoring rule is fixed to negSE
##############################################
serr <- function(y,w,m,s){
  return(-(y-sum(w*m))^2 )
}

##############################################
# mlog - Multivariate inverse logistic transformation
# Input:  
# x, m-variate vector (real-valued)
# i-th element is log( prob(baseline)/prob(event i) )
# Output
# (m+1) variate vector of probabilities, baseline in last row (m+1)
##############################################
mlog <- function(x) c(1,exp(x))/(1+sum(exp(x)))