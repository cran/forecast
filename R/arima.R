search.arima <- function(x, d=NA, D=NA, max.p=5, max.q=5,
    max.P=2, max.Q=2, max.order=5, stationary=FALSE, ic=c("aic","aicc","bic"),
    trace=FALSE,approximation=FALSE,xreg=NULL,offset=offset,allowdrift=TRUE)
{
    ic <- match.arg(ic)
    m <- frequency(x)
         
    oldwarn <- options()$warn
    options(warn=-1)
    on.exit(options(warn=oldwarn))

    if(allowdrift)
        maxK <- (d+D <= 1)
    else
        maxK <- ((d+D) == 0)

    # Choose model orders
    best.ic <- 1e20
    for(i in 0:max.p)
    {
        for(j in 0:max.q)
        {
            for(I in 0:max.P)
            {
                for(J in 0:max.Q)
                {
                    if(i+j+I+J <= max.order)
                    {
                        for(K in 0:maxK)
                        {
                            fit <- myarima(x,order=c(i,d,j),seasonal=c(I,D,J),constant=(K==1),trace=trace,ic=ic,approximation=approximation,offset=offset,xreg=xreg)
                            if(fit$ic < best.ic)
                            {
                                best.ic <- fit$ic
                                bestfit <- fit
                                constant <- (K==1)
                            }
                        }
                    }
                }
            }
        }
    }
 
    if(exists("bestfit"))
    {
        # Refit using ML if approximation used for IC
        if(approximation)
        {
            #constant <- length(bestfit$coef) - ncol(xreg) > sum(bestfit$arma[1:4])
            newbestfit <- myarima(x,order=bestfit$arma[c(1,6,2)],
                seasonal=bestfit$arma[c(3,7,4)],constant=constant,ic,trace=FALSE,approximation=FALSE,xreg=xreg)
            if(newbestfit$ic > 1e19)
            {
                options(warn=oldwarn)
                warning("Unable to fit final model using maximum likelihood. AIC value approximated")
            }
            else
                bestfit <- newbestfit            
        }
    }
    else
        stop("No ARIMA model able to be estimated")

    bestfit$x <- x
    bestfit$series <- deparse(substitute(x))
    bestfit$ic <- NULL
    bestfit$call <- match.call()
#    bestfit$xreg <- xreg

    if(trace)
        cat("\n\n")

    return(bestfit)
}

ndiffs <- function(x,alpha=0.05,test=c("kpss","adf"))
{
    test <- match.arg(test)
    require(tseries)
    x <- c(na.omit(c(x)))
    d <- 0
    oldwarn <- options(warn=-1)
    if(test=="kpss")
        dodiff <- kpss.test(x)$p.value < alpha
    else
        dodiff <- adf.test(x)$p.value > alpha
    if(is.na(dodiff))
    {
        options(warn=oldwarn$warn)
        return(d)
    }
    while(dodiff & d<2)
    {
        x <- diff(x)
        d <- d+1
        if(test=="kpss")
            dodiff <- kpss.test(x)$p.value < alpha
        else
            dodiff <- adf.test(x)$p.value > alpha
        if(is.na(dodiff))
            return(d-1)
    }
    options(warn=oldwarn$warn)
    return(d)
}

# Set up seasonal dummies using Fourier series
SeasDummy <- function(x)
{
    n <- length(x)
    m <- frequency(x)
    if(m==1)
        stop("Non-seasonal data")
    tt <- 1:n
    fmat <- matrix(NA,nrow=n,ncol=2*m)
    for(i in 1:m)
    {
        fmat[,2*i] <- sin(2*pi*i*tt/m)
        fmat[,2*(i-1)+1] <- cos(2*pi*i*tt/m)
    }
    return(fmat[,1:(m-1)])
}

# CANOVA-HANSEN TEST
# Largely based on uroot package code for CH.test()
SD.test <- function (wts, s=frequency(wts))
{
    if(s==1)
        stop("Not seasonal data")
    t0 <- start(wts)
    N <- length(wts)
    if(N <= s)
        stop("Insufficient data")
    frec <- rep(1, as.integer((s+1)/2))
    ltrunc <- round(s * (N/100)^0.25)
    R1 <- as.matrix(SeasDummy(wts))
    lmch <- lm(wts ~ R1)   # run the regression : y(i) = mu+f(i)'gamma(i)+e(i)
    Fhat <- Fhataux <- matrix(nrow = N, ncol = s-1)
    for (i in 1:(s-1))
        Fhataux[, i] <- R1[,i] * lmch$residuals
    for (i in 1:N)
    {
        for (n in 1:(s - 1))
            Fhat[i, n] <- sum(Fhataux[1:i, n])
    }
    wnw <- 1 - seq(1, ltrunc, 1)/(ltrunc + 1)
    Ne <- nrow(Fhataux)
    Omnw <- 0
    for (k in 1:ltrunc)
        Omnw <- Omnw + (t(Fhataux)[, (k + 1):Ne] %*% Fhataux[1:(Ne - k), ]) * wnw[k]
    Omfhat <- (crossprod(Fhataux) + Omnw + t(Omnw))/Ne
    sq <- seq(1, s-1, 2)
    frecob <- rep(0,s - 1)
    for (i in 1:length(frec))
    {
       if (frec[i] == 1 && i == as.integer(s/2))
           frecob[sq[i]] <- 1
       if (frec[i] == 1 && i < as.integer(s/2))
           frecob[sq[i]] <- frecob[sq[i] + 1] <- 1
    }
    a <- length(which(frecob == 1))
    A <- matrix(0, nrow = s - 1, ncol = a)
    j <- 1
    for (i in 1:(s - 1)) if (frecob[i] == 1)
    {
        A[i, j] <- 1
        ifelse(frecob[i] == 1, j <- j + 1, j <- j)
    }
    stL <- (1/N^2) * sum(diag(solve(t(A) %*% Omfhat %*% A, tol=1e-25) %*% t(A) %*% t(Fhat) %*% Fhat %*% A))
    return(stL)
}

# Number of seasonal differences
nsdiffs <- function(x, m=frequency(x))
{
    if(m<=1)
        stop("Non seasonal data")
    chstat <- SD.test(x, m)
    crit.values <- c(0.4617146,0.7479655,1.0007818,1.2375350,1.4625240,1.6920200,1.9043096,2.1169602,
        2.3268562,2.5406922,2.7391007)
    if(m <= 12)
        D <- as.numeric(chstat > crit.values[m-1])
    else if (m == 24)
        D <- as.numeric(chstat > 5.098624)
    else if (m ==52)
        D <- as.numeric(chstat > 10.341416)
    else if (m ==365)
        D <- as.numeric(chstat > 65.44445)
    else
        D <- as.numeric(chstat > 0.269 * m^(0.928))
    return(D=D)
}


forecast.Arima <- function (object, h = ifelse(object$arma[5] > 1, 2 * object$arma[5], 10),
    level = c(80, 95), fan=FALSE, xreg = NULL,...)
{
#    use.constant <- is.element("constant",names(object$coef))
    use.drift <- is.element("drift", names(object$coef))
    if (is.element("x", names(object)))
        x <- object$x
    else
        x <- object$x <- eval.parent(parse(text = object$series))
    usexreg <- (!is.null(xreg) | use.drift)# | use.constant)
#    if(use.constant)
#        xreg <- as.matrix(rep(1,h))
    if(!is.null(xreg))
    {
        xreg <- as.matrix(xreg)
        h <- nrow(xreg)
    }
    if (use.drift)
    {
        n <- length(x)
        if(!is.null(xreg))
            xreg <- cbind(xreg,(1:h)+n)
        else
            xreg <- as.matrix((1:h)+n)
    }
    if(usexreg)
    {
        if(!is.null(object$xreg))
            object$call$xreg <- object$xreg
        else # object from arima() rather than Arima()
        {
            xr <- object$call$xreg
            object$call$xreg <- if (!is.null(xr))         
                                    eval.parent(xr)
                                else NULL
        }
        pred <- predict(object, n.ahead = h, newxreg = xreg)
    }
    else
        pred <- predict(object, n.ahead = h)

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
    lower <- matrix(NA, ncol = nint, nrow = length(pred$pred))
    upper <- lower
    for (i in 1:nint)
    {
        qq <- qnorm(0.5 * (1 + level[i]/100))
        lower[, i] <- pred$pred - qq * pred$se
        upper[, i] <- pred$pred + qq * pred$se
    }
    colnames(lower) = colnames(upper) = paste(level, "%", sep = "")
    method <- arima.string(object)
    return(structure(list(method = method, model = object, level = level,
        mean = pred$pred, lower = lower, upper = upper, x = x,
        xname = deparse(substitute(x)), fitted = fitted(object), residuals = residuals(object)),
        class = "forecast"))
}


forecast.ar <- function(object,h=10,level=c(80,95),fan=FALSE,...)
{
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
    method <- paste("AR(",object$order,")",sep="")
    x <- as.ts(eval.parent(parse(text=object$series)))
    f=frequency(x)
    res <- ts(object$resid[-(1:object$order)],start=tsp(x)[1]+object$order/f,f=f)
    fits <- x-res

    return(structure(list(method=method,model=object,level=level,mean=pred$pred,lower=lower,upper=upper,
        x=x, xname = deparse(substitute(x)), fitted=fits,residuals=res)
        ,class="forecast"))
}

# Extract errors from ARIMA model (as distinct from residuals)
arima.errors <- function(z)
{
    if(!is.list(z))
        stop("z must be a list")
    if(is.element("x",names(z)))
        x <- z$x
    else
    {
        series.name <- z$series
        if(is.null(series.name))
            stop("missing component series in argument z\n")
        x <- eval.parent(parse(text = series.name))
    }
    if(!is.element("xreg",names(z)))
    {
        if(!is.element("xreg",names(z$coef)))
            return(x)
        else
            xreg <- eval.parent(z$coef$xreg)
    }
    else
        xreg <- z$xreg
    norder <- sum(z$arma[1:4])
    if(is.element("intercept",names(z$coef)))
        xreg <- cbind(rep(1,length(x)),xreg)
    return(ts(x - xreg %*% as.matrix(z$coef[(norder+1):length(z$coef)]),f=frequency(x),s=start(x)))
}

# Return one-step fits
fitted.Arima <- function(object,...)
{
    if(is.element("x",names(object)))
        x <- object$x
    else
        x <- eval.parent(parse(text=object$series))

    return(x - object$residuals)
}

# Calls arima from stats package and adds data to the returned object
# Also allows refitting to new data
# and drift terms to be included.
Arima <- function(x, order = c(0, 0, 0),
      seasonal = list(order = c(0, 0, 0), period = NA),
      xreg = NULL, include.mean = TRUE, include.drift = FALSE, transform.pars = TRUE,
      fixed = NULL, init = NULL, method = c("CSS-ML", "ML", "CSS"),
      n.cond, optim.control = list(), kappa = 1e6, model=NULL)
{
    if(!is.null(model))
        tmp <- arima2(x,model,xreg)
    else
    {
        if (!is.null(xreg))
        {
            nmxreg <- deparse(substitute(xreg))
            xreg <- as.matrix(xreg)
            if (is.null(colnames(xreg))) 
                colnames(xreg) <- if (ncol(xreg) == 1) nmxreg
                                  else paste(nmxreg, 1:ncol(xreg), sep = "")
        }
        if(include.drift)
        {
            drift <- 1:length(x)
            xreg <- cbind(xreg,drift=drift)
        }
        if(is.null(xreg))
            tmp <- stats:::arima(x=x,order=order,seasonal=seasonal,include.mean=include.mean,
                transform.pars=transform.pars,fixed=fixed,init=init,method=method,n.cond=n.cond,optim.control=optim.control,kappa=kappa)
        else
            tmp <- stats:::arima(x=x,order=order,seasonal=seasonal,xreg=xreg,include.mean=include.mean,
                   transform.pars=transform.pars,fixed=fixed,init=init,method=method,n.cond=n.cond,optim.control=optim.control,kappa=kappa)
    }
    tmp$x <- x
    tmp$series <- deparse(substitute(x))
    tmp$xreg <- xreg
    tmp$call <- match.call()
    return(tmp)
}

# Refits the model to new data x
arima2 <- function (x, model, xreg)
{
    use.drift <- is.element("drift",names(model$coef))
    use.intercept <- is.element("intercept",names(model$coef))
    use.xreg <- is.element("xreg",names(model$call))
    if(use.drift)
    {
        time.offset <- (tsp(x)[1] - tsp(model$x)[1])*tsp(x)[3]
        xreg <- cbind(xreg, time.offset+(1:length(x)))
        names(xreg)[ncol(xreg)] <- "drift"    
    }

    if(model$arma[5]>1 & sum(abs(model$arma[c(3,4,7)]))>0) # Seasonal model
    {
        if(use.xreg)
            refit <- Arima(x,order=model$arma[c(1,6,2)],seasonal=list(order=model$arma[c(3,7,4)],period=model$arma[5]),
                fixed=model$coef,include.mean=use.intercept,xreg=xreg)
        else
            refit <- Arima(x,order=model$arma[c(1,6,2)],seasonal=list(order=model$arma[c(3,7,4)],period=model$arma[5]),
                fixed=model$coef,include.mean=use.intercept)
    }
    else if(length(model$coef)>0) # Nonseasonal model with some parameters
    {
        if(use.xreg)
           refit <- Arima(x,order=model$arma[c(1,6,2)],fixed=model$coef,xreg=xreg,include.mean=use.intercept)
        else
            refit <- Arima(x,order=model$arma[c(1,6,2)],fixed=model$coef,include.mean=use.intercept)
    }
    else # No parameters
            refit <- Arima(x,order=model$arma[c(1,6,2)],include.mean=FALSE)

    refit$var.coef <- matrix(0,length(refit$coef),length(refit$coef))
    return(refit)
}

# Modified version of function in stats package

print.Arima <- function (x, digits = max(3, getOption("digits") - 3), se = TRUE,
    ...)
{
#    if (!is.element("x", names(x)))
#        x$x <- eval.parent(parse(text = x$series))

    cat("Series:",x$series,"\n")
    cat(arima.string(x),"\n")
    cat("\nCall:", deparse(x$call, width.cutoff = 75), "\n", sep = " ")
#    if(!is.null(x$xreg))
#    {
#        cat("\nRegression variables fitted:\n")
#        xreg <- as.matrix(x$xreg)
#        for(i in 1:3)
#            cat("  ",xreg[i,],"\n")
#        cat("   . . .\n")
#        for(i in 1:3)
#            cat("  ",xreg[nrow(xreg)-3+i,],"\n")
#    }
    if (length(x$coef) > 0) {
        cat("\nCoefficients:\n")
        coef <- round(x$coef, digits = digits)
        if (se && nrow(x$var.coef)) {
            ses <- rep(0, length(coef))
            ses[x$mask] <- round(sqrt(diag(x$var.coef)), digits = digits)
            coef <- matrix(coef, 1, dimnames = list(NULL, names(coef)))
            coef <- rbind(coef, s.e. = ses)
        }
        print.default(coef, print.gap = 2)
    }
    cm <- x$call$method
    if (is.null(cm) || cm != "CSS")
    {
        cat("\nsigma^2 estimated as ", format(x$sigma2, digits = digits),
            ":  log likelihood = ", format(round(x$loglik, 2)),"\n",sep="")
        npar <- length(x$coef) + 1
        nstar <- length(x$residuals) - x$arma[6] - x$arma[7]*x$arma[5]
        bic <- x$aic + npar*(log(nstar) - 2)
        aicc <- x$aic + 2*npar*(nstar/(nstar-npar-1) - 1)
        cat("AIC = ", format(round(x$aic, 2)), sep = "")
        cat("   AICc = ", format(round(aicc, 2)), sep = "")
        cat("   BIC = ", format(round(bic, 2)), "\n",sep = "")
    }
    else cat("\nsigma^2 estimated as ", format(x$sigma2, digits = digits),
        ":  part log likelihood = ", format(round(x$loglik, 2)),
        "\n", sep = "")
    invisible(x)
}

# Modified version of function in stats package

predict.Arima <- function(object, n.ahead = 1, newxreg = NULL, se.fit = TRUE, ...) 
{
    myNCOL <- function(x) if (is.null(x)) 
        0
    else NCOL(x)
    rsd <- object$residuals
    ## LINES ADDED
    if(!is.null(object$xreg))
        object$call$xreg <- object$xreg
    ## END ADDITION
    xr <- object$call$xreg
    xreg <- if (!is.null(xr)) 
        eval.parent(xr)
    else NULL
    ncxreg <- myNCOL(xreg)
    if (myNCOL(newxreg) != ncxreg) 
        stop("'xreg' and 'newxreg' have different numbers of columns: ", ncxreg, " != ", myNCOL(newxreg))
    class(xreg) <- NULL
    xtsp <- tsp(rsd)
    n <- length(rsd)
    arma <- object$arma
    coefs <- object$coef
    narma <- sum(arma[1:4])
    if (length(coefs) > narma) {
        if (names(coefs)[narma + 1] == "intercept") {
            xreg <- cbind(intercept = rep(1, n), xreg)
            newxreg <- cbind(intercept = rep(1, n.ahead), newxreg)
            ncxreg <- ncxreg + 1
        }
        xm <- if (narma == 0) 
            drop(as.matrix(newxreg) %*% coefs)
        else drop(as.matrix(newxreg) %*% coefs[-(1:narma)])
    }
    else xm <- 0
    if (arma[2] > 0) {
        ma <- coefs[arma[1] + 1:arma[2]]
        if (any(Mod(polyroot(c(1, ma))) < 1)) 
            warning("MA part of model is not invertible")
    }
    if (arma[4] > 0) {
        ma <- coefs[sum(arma[1:3]) + 1:arma[4]]
        if (any(Mod(polyroot(c(1, ma))) < 1)) 
            warning("seasonal MA part of model is not invertible")
    }
    z <- KalmanForecast(n.ahead, object$model)
    pred <- ts(z[[1]] + xm, start = xtsp[2] + deltat(rsd), frequency = xtsp[3])
    if (se.fit) {
        se <- ts(sqrt(z[[2]] * object$sigma2), start = xtsp[2] + 
            deltat(rsd), frequency = xtsp[3])
        return(list(pred = pred, se = se))
    }
    else return(pred)
}
