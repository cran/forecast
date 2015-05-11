#Bootstrap functions
# Written by Fotios Petropoulos and Rob J Hyndman

stlboot <- function(x, n=1, include.original=TRUE, lambda=BoxCox.lambda(x, lower=0, upper=1),
  method=c("MBB","LPB")) {
  # The following is an implementation of the time series bootstrapping procedure,
  # as described by Bergmeir et al. (2014, Monash University Working Paper).
  # It allows bootstrapping of non-stationary data using an STL decomposition.
  # Inputs:
  #   x: a time series object or a vector
  #   n: the number of time series to be generated
  #   include.original: change to FALSE if the original should not be returned
  #   lambda: Box-Cox parameter
  #   method: either Moving-Block-Bootstrap or Linear-Process-Boostrap applied 
  #           to the originals.
  # Output:
  #   xs: a 2-dimensional array where each row represents a boostrap
  
  method <- match.arg(method)

  x <- as.ts(x)
  freq <- max(1, frequency(x))
  n <- max(n, 1)
  
  xs <- array(0, c(n, length(x)))
  if (include.original) {
    xs[1,] <- x # the first series is the original one
  }  
  
  if (n>1) {
    # Box-Cox transformation
    x.bc <- BoxCox(x, lambda)
    
    if (freq>1) {
      # STL decomposition
      x.stl <- stl(x.bc, s.window="periodic")$time.series
      seasonal <- x.stl[,1]
      trend <- x.stl[,2]
      remainder <- x.stl[,3]
    } else {      
      # TL decomposition
      x.tl <- tl(x.bc)$time.series
      seasonal <- rep(0, length(x))
      trend <- x.tl[,1]
      remainder <- x.tl[,2]
    }

    bstrap <- ifelse(method=="MBB", MBB, lpb)

    
    # Bootstrap some series
    for (i in (1 + include.original):n) {
      xs[i,] <- InvBoxCox(trend + seasonal + bstrap(remainder, (if(freq>1){freq} else {4})), lambda)
    }
  }
  
  return(xs)
}


# Trend estimation like STL without seasonality.
# Non-robust version
tl <- function(x, ...)
{
  x <- as.ts(x)
  tspx <- tsp(x)
  n <- length(x)
  tt <- 1:n
  fit <- supsmu(tt, x)
  out <- ts(cbind(trend=fit$y, remainder=x-fit$y))
  tsp(out) <- tsp(x)

  out <- structure(list(time.series=out),class="stl")
  return(out)
}



# Function to return some bootstrap samples of x using the Linear Process Bootstrap
# based on McMurry & Politis (2010).
lpb <- function(x, nsim=100)
{
  n <- length(x)
  meanx <- mean(x)
  y <- x - meanx
  gamma <- wacf(y, lag.max=n)$acf[,,1]
  s <- length(gamma)
  Gamma <- matrix(1, s, s)
  d <- row(Gamma) - col(Gamma)
  for(i in 1:(s-1))
    Gamma[d==i | d==(-i)] <- gamma[i+1]
  L <- t(chol(Gamma))
  W <- solve(L) %*% matrix(y,ncol=1)
  out <- ts(L %*% matrix(sample(W, n*nsim, replace=TRUE), nrow=n, ncol=nsim) + meanx)
  tsp(out) <- tsp(x)
  return(out)
}


MBB <- function(x, freq) {
  # The following is an implementation of the moving blocks bootstrap (MBB) algorithm,
  # as described by Bergmeir et al. (2014, Monash University Working Paper).
  # Inputs:
  #   x: a vector (the values of the time series)
  #   freq: the frequency of the data
  # Output:
  #   bx: a vector of the bootstrapped data
    
  bx <- array(0, (floor(length(x)/(2*freq))+2)*2*freq)
  for (i in 1:(floor(length(x)/(2*freq))+2)){
    c <- sample(1:(length(x)-2*freq+1),1)
    bx[((i-1)*2*freq+1):(i*2*freq)] <- x[c:(c+2*freq-1)]
  }
  start_from <- sample(0:(2*freq-1),1) + 1
  return(bx[start_from:(start_from+length(x)-1)])
}
