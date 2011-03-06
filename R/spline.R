###############################################
##### Forecasting Using Smoothing Splines #####
#####   Last updated by RJH. 29/5/02      #####
###############################################

################# FUNCTIONS ##################

## Set up Sigma of order (n x n)
make.Sigma <- function(n,n0=0)
{
    nn <- n + n0
    Sigma <- matrix(0,nrow=nn,ncol=nn)
    for (i in 1:nn)
        Sigma[i,i:nn] <- Sigma[i:nn,i] <- (i*i*(3*(i:nn)-i))/6
    return(Sigma / (n^3))
}

## Compute spline matrices
spline.matrices <- function(n,lambda,cc=1e2,n0=0)
{
    nn <- n+n0
    Sigma <- make.Sigma(n,n0)
    s <- cbind(rep(1,nn),(1:nn)/n)
    Omega <- cc * s %*% t(s) + Sigma/lambda + diag(nn)
    max.Omega <- max(Omega)
    inv.Omega <- solve(Omega/max.Omega,tol=1e-10)/max.Omega
    P <- chol(inv.Omega)
    return(list(s=s,Sigma=Sigma,Omega=Omega,inv.Omega=inv.Omega,P=P))
}

## Compute smoothing splines
## Return -loglikelihood
# Lambda multiplied by 1e4 to avoid numerical difficulties in optimization
spline.loglik <- function(lambda,y,cc=1e2)
{
    n <- length(y)
    mat <- spline.matrices(n,lambda/1e6,cc=cc)
    y.star <- mat$P %*% matrix(y)
    return(-log(det(mat$P)) + 0.5*n*log(sum(y.star^2)))
}

# Spline forecasts
splinef <- function(x, h=10, level=c(80,95), fan=FALSE)
{
    if(!is.ts(x))
        x <- ts(x)
    n <- length(x)
    freq <- frequency(x)
    xname <- deparse(substitute(x))

    # Scale x to avoid numerical problems.
    mean.x <- mean(x)
    stdev.x <- sqrt(var(x))
    newx <- (x-mean.x)/stdev.x

    # Find optimal lambda
    if(n > 100) # Use only last 100 observations to get lambda
        xx <- newx[(n-99):n]
    else
        xx <- newx
    lambda.est <- optimize(spline.loglik, interval=c(1e-6,1e7),y = xx)$minimum/1e6

    # Compute fitted values
    r <- 256 * smooth.spline(1:n,newx,spar=0)$lambda
    lss <- lambda.est*n^3 /(n-1)^3
    spar <- (log(lss/r) / log(256) + 1) /3
    splinefit <- smooth.spline(1:n,newx,spar=spar)
    sfits <- splinefit$y*stdev.x+mean.x

    # Compute matrices for optimal lambda
    mat <- spline.matrices(n,lambda.est)
    newmat <- spline.matrices(n,lambda.est,n0=h)

    # Get one-step predictors
    yfit <- e <- rep(NA,n)
    if(n > 1000)
        warning("Series too long to compute in-sample fits and residuals")
    else  # This is probably grossly inefficient but I can't think of a better way right now
    {
        for(i in 1:(n-1))
        {
            U <- mat$Omega[1:i,i+1]
            Oinv <- solve(mat$Omega[1:i,1:i]/1e6)/1e6
            yfit[i+1] <- t(U) %*% Oinv %*% newx[1:i]
            sd <- sqrt(mat$Omega[i+1,i+1] - t(U) %*% Oinv %*% U)
            e[i+1] <- (newx[i+1]-yfit[i+1])/sd
        }
    }
    # Compute sigma^2
    sigma2 <- mean(e^2,na.rm=TRUE)

    # Compute mean and var of forecasts
    U <- newmat$Omega[1:n,n+(1:h)]
    Omega0 <- newmat$Omega[n+(1:h),n+(1:h)]
    Yhat <- t(U) %*% mat$inv.Omega %*% newx
    sd <- sqrt(sigma2*diag(Omega0 - t(U) %*% mat$inv.Omega %*% U))

    # Compute prediction intervals.
    if(fan)
        level <- seq(51,99,by=3)
    else
    {
        if(min(level) > 0 & max(level) < 1)
            level <- 100*level
        else if(min(level) < 0 | max(level) > 99.99)
            stop("Confidence limit out of range")
    }
    nconf <- length(level)
    lower <- upper <- matrix(NA,nrow=h,ncol=nconf)
    for(i in 1:nconf)
    {
        conf.factor <- qnorm(0.5 + 0.005*level[i])
        upper[,i] <- Yhat + conf.factor*sd
        lower[,i] <- Yhat - conf.factor*sd
    }
    lower <- ts(lower,start=tsp(x)[2]+1/freq,f=freq)
    upper <- ts(upper,start=tsp(x)[2]+1/freq,f=freq)

    return(structure(list(method="Cubic Smoothing Spline",level=level,x=x,mean=ts(stdev.x*Yhat+mean.x,f=freq,start=tsp(x)[2]+1/freq),
            upper=ts(upper*stdev.x + mean.x,start=tsp(x)[2]+1/freq,f=freq),
            lower=ts(lower*stdev.x + mean.x,start=tsp(x)[2]+1/freq,f=freq),
            model=list(lambda=lambda.est*n^3,call=match.call()),
            fitted =ts(sfits,start=start(x),f=freq), residuals=ts(e,start=start(x),f=freq),
            onestepf = ts(yfit*stdev.x+mean.x,start=start(x),f=freq)),
            class=c("splineforecast","forecast")))
}

plot.splineforecast <- function(x,fitcol=2,...)
{
    plot.forecast(x,type="o",pch=19,...)
    lines(x$fitted,col=fitcol)
}
