## Measures of forecast accuracy
## Forecasts in f. This may be a numerical vector or the output from arima or ets or derivatives.
## Actual values in x
## test enables a subset of x and f to be tested.
forecasterrors <- function(f,x,test=1:length(x),lambda=NULL)
{
    data.x <- dx <- NULL
    if(is.list(f))
    {
        if(is.element("x",names(f)))
            data.x <- f$x
        if(is.element("mean",names(f)))
            f = f$mean
        else
            stop("Unknown list structure")
    }
    n <- length(x)
    if(is.null(lambda))
    {
        ff <- f
        xx <- x
        dx <- data.x
    }
    else
    {
        ff <- InvBoxCox(f,lambda)
        xx <- InvBoxCox(x,lambda)
        if(!is.null(data.x))
            dx <- InvBoxCox(data.x,lambda)
    }
    error <- (xx-ff[1:n])[test]
    me <- mean(error)
    mse <- mean(error^2)
    mae <- mean(abs(error))
    mape <- mean(100*abs(error/xx[test]))
    mpe <-  mean(100*error/xx[test])
    junk <- c(me,sqrt(mse),mae,mpe,mape)
    names(junk) <- c("ME","RMSE","MAE","MPE","MAPE")
    if(!is.null(dx))
    {
        scale <- mean(abs(diff(dx)),na.rm=TRUE)
        mase <- mean(abs(error/scale))
        junk <- c(junk,mase)
        names(junk)[6] <- "MASE"
    }
    if(n>1)
    {
        fpe <- (c(ff[2:n])/c(xx[1:(n-1)]) - 1)[test-1]
        ape <- (c(xx[2:n])/c(xx[1:(n-1)]) - 1)[test-1]
        theil <- sqrt(sum((fpe - ape)^2)/sum(ape^2))
        r1 <- acf(error,plot=FALSE,lag.max=2)$acf[2,1,1]
        nj <- length(junk)
        junk <- c(junk,r1,theil)
        names(junk)[nj+(1:2)] <- c("ACF1","Theil's U")
    }
    return(junk)
}


accuracy <- function(f,x,test=1:length(x),lambda=NULL)
{
    if(!missing(x))
        return(forecasterrors(f,x,test,lambda))
    if(class(f)=="Arima" & !is.element("x", names(f)))
        f$x <- eval(parse(text = f$series))
    if(!is.null(lambda))  # undo Box-Cox transformation
    {
        ff <- InvBoxCox(f$x,lambda)
        fits <- InvBoxCox(fitted(f),lambda)
    }
    else
    {
        ff <- f$x
        fits <- fitted(f)    # Don't use f$resid as this may contain multiplicative errors.
    }
    res <- ff-fits

    pe <- res/ff * 100 # Percentage error
    scale <- mean(abs(diff(ff)),na.rm=TRUE)
    out <- c(mean(res,na.rm=TRUE), sqrt(mean(res^2,na.rm=TRUE)), mean(abs(res),na.rm=TRUE), mean(pe,na.rm=TRUE), mean(abs(pe),na.rm=TRUE),
        mean(abs(res/scale),na.rm=TRUE))
    names(out) <- c("ME","RMSE","MAE","MPE","MAPE","MASE")
    return(out)
}
