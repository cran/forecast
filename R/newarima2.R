auto.arima <- function(x, d=NA, D=NA, max.p=5, max.q=5,
    max.P=2, max.Q=2, max.order=5,
    start.p=2, start.q=2, start.P=1, start.Q=1,
    stationary=FALSE, ic=c("aic","aicc","bic"),
    stepwise=TRUE, trace=FALSE,
    approximation=length(x)>100 | frequency(x)>12,xreg=NULL,test=c("kpss","adf","pp"),allowdrift=TRUE)
{
    ic <- match.arg(ic)
    test <- match.arg(test)
    m <- frequency(x)

    # Choose order of differencing
    if(!is.null(xreg))
    {
        nmxreg <- deparse(substitute(xreg))
        xreg <- as.matrix(xreg)
        if(ncol(xreg)==1 & length(nmxreg) > 1)
            nmxreg <- "xreg"
        if (is.null(colnames(xreg))) 
            colnames(xreg) <- if (ncol(xreg) == 1) nmxreg
                              else paste(nmxreg, 1:ncol(xreg), sep = "")
        xx <- ts(residuals(lm(x ~ xreg)))
        tsp(xx) <- tsp(x)
    }
    else
        xx <- x
    if(stationary)
        d <- D <- 0
    if(m == 1)
        D <- max.P <- max.Q <- 0
    else if(is.na(D))
        D <- nsdiffs(xx)
    if(D > 0)
        dx <- diff(xx,D)
    else
        dx <- xx
    if(is.na(d))
        d <- ndiffs(dx,test=test)

    if(m > 1)
    {
        if(max.P > 0)
            max.p <- min(max.p, m-1)
        if(max.Q > 0)
            max.q <- min(max.q, m-1)
    }

    # Find constant offset for AIC calculation using simple AR(1) model
    if(approximation)
    {
        if(D==0)
            fit <- try(arima(x,order=c(1,d,0),xreg=xreg))
        else
            fit <- try(arima(x,order=c(1,d,0),seasonal=list(order=c(0,D,0),period=m,xreg=xreg)))
        if(class(fit) != "try-error")
            offset <- -2*fit$loglik - length(x)*log(fit$sigma2)
        else
        {
            warning("Unable to calculate AIC offset")
            offset <- 0
        }
    }
    else
        offset <- 0

    if(!stepwise)
    {
        bestfit <- search.arima(x,d,D,max.p,max.q,max.P,max.Q,max.order,stationary,ic,trace,approximation,xreg=xreg,offset=offset,allowdrift=allowdrift)
        bestfit$call <- match.call()
        return(bestfit)
    }

    # Starting model
    p <- start.p <- min(start.p,max.p)
    q <- start.q <- min(start.q,max.q)
    P <- start.P <- min(start.P,max.P)
    Q <- start.Q <- min(start.Q,max.Q)
    if(allowdrift)
        constant <- (d+D <= 1)
    else 
        constant <- ((d+D) == 0)
        
    results <- matrix(NA,nrow=100,ncol=8)
   
    oldwarn <- options()$warn
    options(warn=-1)
    on.exit(options(warn=oldwarn))

    bestfit <- myarima(x,order=c(p,d,q),seasonal=c(P,D,Q),constant=constant,ic,trace,approximation,offset=offset,xreg=xreg)
    results[1,] <- c(p,d,q,P,D,Q,constant,bestfit$ic)
    # Null model
    fit <- myarima(x,order=c(0,d,0),seasonal=c(0,D,0),constant=constant,ic,trace,approximation,offset=offset,xreg=xreg)
    results[2,] <- c(0,d,0,0,D,0,constant,fit$ic)
    if(fit$ic < bestfit$ic)
    {
        bestfit <- fit
        p <- q <- P <- Q <- 0
    }
    # Basic AR model
    if(max.p > 0 | max.P > 0)
    {
        fit <- myarima(x,order=c(max.p>0,d,0),seasonal=c((m>1)&(max.P>0),D,0),constant=constant,ic,trace,approximation,offset=offset,xreg=xreg)
        results[3,] <- c(1,d,0,m>1,D,0,constant,fit$ic)
        if(fit$ic < bestfit$ic)
        {
            bestfit <- fit
            p <- (max.p>0)
            P <- (m>1) & (max.P>0)
            q <- Q <- 0
        }
    }
    # Basic MA model
    if(max.q > 0 | max.Q > 0)
    {
        fit <- myarima(x,order=c(0,d,max.q>0),seasonal=c(0,D,(m>1)&(max.Q>0)),constant=constant,ic,trace,approximation,offset=offset,xreg=xreg)
        results[4,] <- c(0,d,1,0,D,m>1,constant,fit$ic)
        if(fit$ic < bestfit$ic)
        {
            bestfit <- fit
            p <- P <- 0
            Q <- (m>1) & (max.Q>0)
            q <- (max.q>0)
        }
    }
     k <- 4

    startk <- 0
    while(startk < k & k < 94)
    {
        startk <- k
        if(P > 0 & newmodel(p,d,q,P-1,D,Q,constant,results[1:k,]))
        {
            k <- k + 1
            fit <- myarima(x,order=c(p,d,q),seasonal=c(P-1,D,Q),constant=constant,ic,trace,approximation,offset=offset,xreg=xreg)
            results[k,] <- c(p,d,q,P-1,D,Q,constant,fit$ic)
            if(fit$ic < bestfit$ic)
            {
                bestfit <- fit
                P <- (P-1)
            }
        }
        if(P < max.P & newmodel(p,d,q,P+1,D,Q,constant,results[1:k,]))
        {
            k <- k + 1
            fit <- myarima(x,order=c(p,d,q),seasonal=c(P+1,D,Q),constant=constant,ic,trace,approximation,offset=offset,xreg=xreg)
            results[k,] <- c(p,d,q,P+1,D,Q,constant,fit$ic)
            if(fit$ic < bestfit$ic)
            {
                bestfit <- fit
                P <- (P+1)
            }
        }
        if(Q > 0 & newmodel(p,d,q,P,D,Q-1,constant,results[1:k,]))
        {
            k <- k + 1
            fit <- myarima(x,order=c(p,d,q),seasonal=c(P,D,Q-1),constant=constant,ic,trace,approximation,offset=offset,xreg=xreg)
            results[k,] <- c(p,d,q,P,D,Q-1,constant,fit$ic)
            if(fit$ic < bestfit$ic)
            {
                bestfit <- fit
                Q <- (Q-1)
            }
        }
        if(Q < max.Q & newmodel(p,d,q,P,D,Q+1,constant,results[1:k,]))
        {
            k <- k + 1
            fit <- myarima(x,order=c(p,d,q),seasonal=c(P,D,Q+1),constant=constant,ic,trace,approximation,offset=offset,xreg=xreg)
            results[k,] <- c(p,d,q,P,D,Q+1,constant,fit$ic)
            if(fit$ic < bestfit$ic)
            {
                bestfit <- fit
                Q <- (Q+1)
            }
        }
        if(Q > 0 & P > 0 & newmodel(p,d,q,P-1,D,Q-1,constant,results[1:k,]))
        {
            k <- k + 1
            fit <- myarima(x,order=c(p,d,q),seasonal=c(P-1,D,Q-1),constant=constant,ic,trace,approximation,offset=offset,xreg=xreg)
            results[k,] <- c(p,d,q,P-1,D,Q-1,constant,fit$ic)
            if(fit$ic < bestfit$ic)
            {
                bestfit <- fit
                Q <- (Q-1)
                P <- (P-1)
            }
        }
        if(Q < max.Q & P < max.P & newmodel(p,d,q,P+1,D,Q+1,constant,results[1:k,]))
        {
            k <- k + 1
            fit <- myarima(x,order=c(p,d,q),seasonal=c(P+1,D,Q+1),constant=constant,ic,trace,approximation,offset=offset,xreg=xreg)
            results[k,] <- c(p,d,q,P+1,D,Q+1,constant,fit$ic)
            if(fit$ic < bestfit$ic)
            {
                bestfit <- fit
                Q <- (Q+1)
                P <- (P+1)
            }
        }

        if(p > 0 & newmodel(p-1,d,q,P,D,Q,constant,results[1:k,]))
        {
            k <- k + 1
            fit <- myarima(x,order=c(p-1,d,q),seasonal=c(P,D,Q),constant=constant,ic,trace,approximation,offset=offset,xreg=xreg)
            results[k,] <- c(p-1,d,q,P,D,Q,constant,fit$ic)
            if(fit$ic < bestfit$ic)
            {
                bestfit <- fit
                p <- (p-1)
            }
        }
        if(p < max.p & newmodel(p+1,d,q,P,D,Q,constant,results[1:k,]))
        {
            k <- k + 1
            fit <- myarima(x,order=c(p+1,d,q),seasonal=c(P,D,Q),constant=constant,ic,trace,approximation,offset=offset,xreg=xreg)
            results[k,] <- c(p+1,d,q,P,D,Q,constant,fit$ic)
            if(fit$ic < bestfit$ic)
            {
                bestfit <- fit
                p <- (p+1)
            }
        }
        if(q > 0 & newmodel(p,d,q-1,P,D,Q,constant,results[1:k,]))
        {
            k <- k + 1
            fit <- myarima(x,order=c(p,d,q-1),seasonal=c(P,D,Q),constant=constant,ic,trace,approximation,offset=offset,xreg=xreg)
            results[k,] <- c(p,d,q-1,P,D,Q,constant,fit$ic)
            if(fit$ic < bestfit$ic)
            {
                bestfit <- fit
                q <- (q-1)
            }
        }
        if(q < max.q & newmodel(p,d,q+1,P,D,Q,constant,results[1:k,]))
        {
            k <- k + 1
            fit <- myarima(x,order=c(p,d,q+1),seasonal=c(P,D,Q),constant=constant,ic,trace,approximation,offset=offset,xreg=xreg)
            results[k,] <- c(p,d,q+1,P,D,Q,constant,fit$ic)
            if(fit$ic < bestfit$ic)
            {
                bestfit <- fit
                q <- (q+1)
            }
        }
        if(q > 0 & p > 0 & newmodel(p-1,d,q-1,P,D,Q,constant,results[1:k,]))
        {
            k <- k + 1
            fit <- myarima(x,order=c(p-1,d,q-1),seasonal=c(P,D,Q),constant=constant,ic,trace,approximation,offset=offset,xreg=xreg)
            results[k,] <- c(p-1,d,q-1,P,D,Q,constant,fit$ic)
            if(fit$ic < bestfit$ic)
            {
                bestfit <- fit
                q <- (q-1)
                p <- (p-1)
            }
        }
        if(q < max.q & p < max.p & newmodel(p+1,d,q+1,P,D,Q,constant,results[1:k,]))
        {
            k <- k + 1
            fit <- myarima(x,order=c(p+1,d,q+1),seasonal=c(P,D,Q),constant=constant,ic,trace,approximation,offset=offset,xreg=xreg)
            results[k,] <- c(p+1,d,q+1,P,D,Q,constant,fit$ic)
            if(fit$ic < bestfit$ic)
            {
                bestfit <- fit
                q <- (q+1)
                p <- (p+1)
            }
        }
        if(allowdrift | (d+D)==0)
        {
            if(newmodel(p,d,q,P,D,Q,!constant,results[1:k,]))
            {
                k <- k + 1
                fit <- myarima(x,order=c(p,d,q),seasonal=c(P,D,Q),constant=!constant,ic,trace,approximation,offset=offset,xreg=xreg)
                results[k,] <- c(p,d,q,P,D,Q,!constant,fit$ic)
                if(fit$ic < bestfit$ic)
                {
                    bestfit <- fit
                    constant <- !constant
                }
            }
        }
    }

    # Refit using ML if approximation used for IC
    if(approximation)
    {
        #constant <- length(bestfit$coef) > sum(bestfit$arma[1:4])
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

    if(bestfit$ic > 1e19)
    {
        cat("\n")
        stop("No suitable ARIMA model found")
    }
    
    # Return best fit
    bestfit$x <- x
    bestfit$series <- deparse(substitute(x))
    bestfit$ic <- NULL
    bestfit$call <- match.call()

    if(trace)
        cat("\n\n Best model:",arima.string(bestfit),"\n\n")

    return(bestfit)
}


# Calls arima from stats package and adds data to the returned object
# Also allows refitting to new data
# and drift terms to be included.
myarima <- function(x, order = c(0, 0, 0), seasonal = c(0, 0, 0), constant=TRUE, ic="aic", trace=FALSE,approximation=FALSE,offset=0,xreg=NULL)
{
    n <- length(x)
    m <- frequency(x)
    use.season <- (sum(seasonal)>0) & m>0
    diffs <- order[2]+seasonal[2]
    if(approximation)
        method <- "CSS"
    else
        method <- "CSS-ML"
    if(diffs==1 & constant)
    {
        xreg <- cbind(1:length(x),xreg)
        if(use.season)
            fit <- try(stats:::arima(x=x,order=order,seasonal=list(order=seasonal,period=m),xreg=xreg,method=method),silent=TRUE)
        else
            fit <- try(stats:::arima(x=x,order=order,xreg=xreg,method=method),silent=TRUE)
    }
    else
    {
        if(use.season)
            fit <- try(stats:::arima(x=x,order=order,seasonal=list(order=seasonal,period=m),include.mean=constant,method=method,xreg=xreg),silent=TRUE)
        else
            fit <- try(stats:::arima(x=x,order=order,include.mean=constant,method=method,xreg=xreg),silent=TRUE)
    }
    if(class(fit) != "try-error")
    {
        nstar <- n - order[2] - seasonal[2]*m
        if(diffs==1 & constant)
        {
            fitnames <- names(fit$coef)
            fitnames[length(fitnames)] <- "drift"
            names(fit$coef) <- fitnames
            fit$xreg <- xreg
        }
        npar <- length(fit$coef) + 1
        if(approximation)
            fit$aic <- offset + nstar * log(fit$sigma2) + 2 * npar
        if(!is.na(fit$aic))
        {
            fit$bic <- fit$aic + npar*(log(nstar) - 2)
            fit$aicc <- fit$aic + 2*npar*(nstar/(nstar-npar-1) - 1)
            fit$ic <- switch(ic,bic=fit$bic,aic=fit$aic,aicc=fit$aicc)
        }
        else
            fit$aic <- fit$bic <- fit$aicc <- fit$ic <- 1e20
        # Check for unit roots
        minroot <- 2
        if(order[1] + seasonal[1] > 0)
        {
            testvec <- fit$model$phi
            last.nonzero <- max(which(abs(testvec)>1e-8)) 
            if(last.nonzero > 0)
            {
                testvec <- testvec[1:last.nonzero]
                if(last.nonzero > 48)
                    warning("Unable to check for unit roots")
                else
                    minroot <- min(minroot,abs(polyroot(c(1,-testvec))))
            }
        }
        if(order[3] + seasonal[3] > 0)
        {
            testvec <- fit$model$theta
            last.nonzero <- max(which(abs(testvec)>1e-8)) 
            if(last.nonzero > 0)
            {
                testvec <- testvec[1:last.nonzero]
                if(last.nonzero > 48)
                    warning("Unable to check for unit roots")
                else
                    minroot <- min(minroot,abs(polyroot(c(1,testvec))))
            }
        }
        if(minroot < 1 + 1e-3)
            fit$ic <- 1e20 # Don't like this model
        if(trace)
            cat("\n",arima.string(fit),":",fit$ic)
        return(fit)
    }
    else
    {
        if(trace)
        {
            cat("\n ARIMA(",order[1],",",order[2],",",order[3],")",sep="")
            if(use.season)
                cat("(",seasonal[1],",",seasonal[2],",",seasonal[3],")[",m,"]",sep="")
            if(constant & (order[2]+seasonal[2] == 0))
                cat(" with non-zero mean")
            else if(constant & (order[2]+seasonal[2] == 1))
                cat(" with drift        ")
            else if(!constant & (order[2]+seasonal[2] == 0))
                cat(" with zero mean    ")
            else
                cat("         ")
            cat(" :",1e20,"*")
        }
        return(list(ic=1e20))
    }
}

newmodel <- function(p,d,q,P,D,Q,constant,results)
{
    n <- nrow(results)
    for(i in 1:n)
    {
        if(identical(c(p,d,q,P,D,Q,constant),results[i,1:7]))
            return(FALSE)
    }
    return(TRUE)
}

arima.string <- function(object)
{
    order <- object$arma[c(1,6,2,3,7,4,5)]
    result <- paste("ARIMA(",order[1],",",order[2],",",order[3],")",sep="")
    if(order[7]>1 & sum(order[4:6]) > 0)
        result <- paste(result,"(",order[4],",",order[5],",",order[6],")[",order[7],"]",sep="")
    if(is.element("constant",names(object$coef)) | is.element("intercept",names(object$coef)))
        result <- paste(result,"with non-zero mean")
    else if(is.element("drift",names(object$coef)))
        result <- paste(result,"with drift        ")
    else if(order[2]==0 & order[5]==0)
        result <- paste(result,"with zero mean    ")
    else
        result <- paste(result,"                  ")
    return(result)
}

summary.Arima <- function(object,...)
{
    print(object)
    cat("\nIn-sample error measures:\n")
    print(accuracy(object))
}
