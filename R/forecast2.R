# Mean forecast
meanf <- function(x,h=10,level=c(80,95),fan=FALSE)
{
    xname <- deparse(substitute(x))
    n <- length(x)
    if(!is.ts(x))
        x <- ts(x)
    start.x <- tsp(x)[1]
    f=ts(rep(mean(x),h),start=start.x,f=frequency(x))
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
    s <- sd(x)
    for(i in 1:nconf)
    {
        tfrac <- qt( 0.5 - level[i]/200, n-1)
        w <- -tfrac * s*sqrt(1+1/n)
        lower[,i] = f-w
        upper[,i] = f+w
    }
    freq=frequency(x)
    lower <- ts(lower,start=tsp(x)[2]+1/freq,f=freq)
    upper <- ts(upper,start=tsp(x)[2]+1/freq,f=freq)
    colnames(lower) <- colnames(upper) <- paste(level,"%",sep="")
    fits <- rep(NA,n-1)
    for(i in 1:(n-1))
        fits[i] <- mean(x[1:i])
    res <- x[2:n] - fits

    junk <- list(method="Mean",level=level,x=x,xname=xname,mean=ts(f,start=tsp(x)[2]+1/freq,f=freq),lower=lower,upper=upper,
        model=list(mu=f[1],mu.se=s/sqrt(length(x)),sd=s,residuals=ts(c(x-f[1]),start=start.x,f=freq)),
        fitted =ts(fits,start=start.x+1/freq,f=freq), residuals=ts(res,start=start.x+1/freq,f=freq))
    junk$model$call <- match.call()

    return(structure(junk,class="forecast"))
}

thetaf <- function(x,h=10,level=c(80,95),fan=FALSE)
{
    if(fan)
        level <- seq(51,99,by=3)
    else
    {
        if(min(level) > 0 & max(level) < 1)
            level <- 100*level
        else if(min(level) < 0 | max(level) > 99.99)
            stop("Confidence limit out of range")
    }
    fcast <- ses(x,h=h)
    tmp2 <- lsfit(0:(length(x)-1),x)$coef[2]/2
    alpha <- fcast$model$par["alpha"]
    n <- length(x)
    fcast$mean <- fcast$mean + tmp2*(0:(h-1) + (1-(1-alpha)^n)/alpha)
    fcast.se <- sqrt(fcast$model$sigma) * sqrt((0:(h-1))*alpha^2+1)
    nconf <- length(level)
    fcast$lower <- fcast$upper <- matrix(NA,nrow=h,ncol=nconf)
    for(i in 1:nconf)
    {
        zt <- -qnorm( 0.5 - level[i]/200)
        fcast$lower[,i] <- fcast$mean - zt*fcast.se
        fcast$upper[,i] <- fcast$mean + zt*fcast.se
    }
    fcast$level <- level
    fcast$method <- "Theta"
    fcast$model <- list(alpha=alpha,drift=tmp2,sigma=fcast$model$sigma)
    fcast$model$call <- match.call()
    return(fcast)
}

# Random walk
rwf <- function(x,h=10,drift=FALSE,level=c(80,95),fan=FALSE)
{
    xname <- deparse(substitute(x))
    n <- length(x)
    nn <- 1:h
    if(!is.ts(x))
        x <- ts(x)
    if(drift)
    {
        fit <- summary(lm(diff(x) ~ 1))
        b <- fit$coefficients[1,1]
        b.se <- fit$coefficients[1,2]
        s <- fit$sigma
        res <- residuals(fit)
        method <- "Random walk with drift"
    }
    else
    {
        b <- b.se <- 0
        s <- sd(diff(x))
        res <- diff(x)
        method <- "Random walk"
    }
    f <- x[n] + nn*b
    se <- sqrt((nn*s^2) + (nn*b.se)^2)

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
    z <- qnorm(.5 + level/200)
    lower <- upper <- matrix(NA,nrow=h,ncol=nconf)
    for(i in 1:nconf)
    {
        lower[,i] <- f - z[i]*se
        upper[,i] <- f + z[i]*se
    }
    freq=frequency(x)
    lower <- ts(lower,start=tsp(x)[2]+1/freq,f=freq)
    upper <- ts(upper,start=tsp(x)[2]+1/freq,f=freq)
    colnames(lower) <- colnames(upper) <- paste(level,"%",sep="")
    junk <- list(method=method,level=level,x=x,xname=xname,mean=ts(f,start=tsp(x)[2]+1/freq,f=freq),lower=lower,upper=upper,
        model=list(drift=b,drift.se=b.se,sd=s,residuals=ts(c(x-f[1]),start=tsp(x)[1],f=freq)),
        fitted = ts(x[-n],start=tsp(x)[1]+1/freq,f=freq),
        residuals = diff(x))
    junk$model$call <- match.call()

    return(structure(junk,class="forecast"))
}

BoxCox <- function(x,lambda)
{
    if(lambda==0)
        log(x)
    else
        (x^lambda - 1)/lambda
}

InvBoxCox <- function(x,lambda)
{
    if(lambda==0)
        exp(x)
    else
        (x*lambda + 1)^(1/lambda)
}

forecast.StructTS <- function(object,h=ifelse(object$call$type=="BSM", 2*object$xtsp[3], 10),level=c(80,95),fan=FALSE,...)
{
    xname <- deparse(substitute(x))
    x <- object$data
    pred <- predict(object,n.ahead=h)
    if(fan)
        level <- seq(51,99,by=3)
    else
    {
        if(min(level) > 0 & max(level) < 1)
            level <- 100*level
        else if(min(level) < 0 | max(level) > 99.99)
            stop("Confidence limit out of range")
    }
    nint <- length(level)
    lower <- matrix(NA,ncol=nint,nrow=length(pred$pred))
    upper <- lower
    for(i in 1:nint)
    {
        qq <- qnorm(0.5*(1+level[i]/100))
        lower[,i] <- pred$pred - qq*pred$se
        upper[,i] <- pred$pred + qq*pred$se
    }
    colnames(lower) = colnames(upper) = paste(level,"%",sep="")
    if(object$call$type=="BSM")
        method <- "Basic structural model"
    else if(object$call$type=="level")
        method <- "Local level structural model"
    else if(object$call$type=="trend")
        method <- "Local linear structural model"
    return(structure(list(method=method,model=object,level=level,mean=pred$pred,lower=lower,upper=upper,
        x=x,xname=xname,fitted=x-residuals(object),residuals=residuals(object)),
        class="forecast"))
}

forecast.HoltWinters <- function(object,h=ifelse(frequency(object$x)>1,2*frequency(object$x),10),level=c(80,95),fan=FALSE,...)
{
    xname <- deparse(substitute(x))
    x <- object$x
    pred <- predict(object,n.ahead=h,prediction.interval=TRUE,level=level[1]/100)
    pmean <- pred[,1]
    if(fan)
        level <- seq(51,99,by=3)
    else
    {
        if(min(level) > 0 & max(level) < 1)
            level <- 100*level
        else if(min(level) < 0 | max(level) > 99.99)
            stop("Confidence limit out of range")
    }
    nint <- length(level)
    upper <- lower <- matrix(NA,ncol=nint,nrow=length(pred[,1]))
    se <- (pred[,2]-pred[,3])/(2*qnorm(0.5*(1+level[1]/100)))
    for(i in 1:nint)
    {
        qq <- qnorm(0.5*(1+level[i]/100))
        lower[,i] <- pmean - qq*se
        upper[,i] <- pmean + qq*se
    }
    colnames(lower) = colnames(upper) = paste(level,"%",sep="")
    method <- "HoltWinters"
    return(structure(list(method=method, model=object, level=level,
        mean=pmean, lower=lower, upper=upper,
        x=x, xname=xname, fitted=object$fitted, residuals=residuals(object)),
        class="forecast"))
}

## CROSTON

croston <- function(x,h=10,alpha=0.1)
{
    if(sum(x<0) > 0)
        stop("Series should not contain negative values")
    out <- croston2(x,h,alpha)
    out$x <- x
    out$xname <- deparse(substitute(x))
    out$residuals <- x-out$fitted
    out$method <- "Croston's method"
    return(structure(out,class="forecast"))
}

croston2 <- function(x,h=10,alpha=0.1,nofits=FALSE)
{
    x <- as.ts(x)
    y <- x[x>0]
    tt <- diff(c(0,(1:length(x))[x>0]))
    if(length(y)==1 & length(tt)==1)
        return(rep(y/tt,h))
    else if(length(y)<=1 | length(tt)<=1)
        return(rep(NA,h))
    y.f <- forecast(HoltWinters(y,beta=FALSE,gamma=FALSE,alpha=alpha),h=h)
    p.f <- forecast(HoltWinters(tt,beta=FALSE,gamma=FALSE,alpha=alpha),h=h)
    freq <- frequency(x)
    tsp.x <- tsp(x)
    if (!is.null(tsp.x)) 
    {
        start.f <- tsp.x[2] + 1/tsp.x[3]
        freq.x <- tsp.x[3]
    }
    else 
    {
        start.f <- length(x) + 1
        freq.x <- 1
    }
    ratio <- ts(y.f$mean/p.f$mean,start=start.f, f = freq.x)
    if(nofits)
        return(ratio)
    else
    {
        n <- length(x)
        junk <- x*NA
        for(i in 1:(n-1))
            junk[i+1] <- croston2(x[1:i],h=1,alpha=alpha,nofits=TRUE)
        return(list(mean = ratio, fitted = junk, model=list(demand=y.f,period=p.f)))
    }
}
