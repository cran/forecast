## Generic forecast functions
## Part of forecast and demography packages

forecast <- function(object,...) UseMethod("forecast")

forecast.default <- function(object,...) forecast.ts(object,...)

forecast.ts <- function(object, h=ifelse(frequency(object)>1, 2*frequency(object), 10), level=c(80,95), fan=FALSE, ...)
{
    n <- length(object)
    if(n > 3)
    {
        if(frequency(object) < 13)
            forecast(ets(object,...),h=h,level=level,fan=fan)
        else
            stlf(object,h=h,level=level,fan=fan,...)
    }
    else 
        meanf(object,h=h,level=level,fan=fan,...)
}

as.data.frame.forecast <- function(x,...)
{
    nconf <- length(x$level)
    out <- matrix(x$mean, ncol=1)
    ists <- is.ts(x$mean)
    if(ists)
    {
        out <- ts(out)
        attributes(out)$tsp <- attributes(x$mean)$tsp
    }
    names <- c("Point Forecast")
    if (!is.null(x$lower) & !is.null(x$upper) & !is.null(x$level)) 
    {
        x$upper <- as.matrix(x$upper)
        x$lower <- as.matrix(x$lower)
        for (i in 1:nconf) 
        {
            out <- cbind(out, x$lower[, i], x$upper[, i])
            names <- c(names, paste("Lo", x$level[i]), paste("Hi", x$level[i]))
        }
    }
    colnames(out) <- names
    rownames(out) <- time(x$mean)
    # Rest of function borrowed from print.ts(), but with header() omitted
    if(!ists)
        return(as.data.frame(out))

    x <- as.ts(out)
    fr.x <- frequency(x)
    calendar <- any(fr.x == c(4, 12)) && length(start(x)) ==  2L
    Tsp <- tsp(x)
    if (is.null(Tsp)) 
    {
        warning("series is corrupt, with no 'tsp' attribute")
        print(unclass(x))
        return(invisible(x))
    }
    nn <- 1 + round((Tsp[2L] - Tsp[1L]) * Tsp[3L])
    if (NROW(x) != nn) 
    {
        warning(gettextf("series is corrupt: length %d with 'tsp' implying %d", NROW(x), nn), domain = NA, call. = FALSE)
        calendar <- FALSE
    }
    if (NCOL(x) == 1) 
    {
        if (calendar) 
        {
            if (fr.x > 1) 
            {
                dn2 <- if (fr.x == 12) 
                  month.abb
                else if (fr.x == 4) 
                  c("Qtr1", "Qtr2", "Qtr3", "Qtr4")
                else paste("p", 1L:fr.x, sep = "")
                if (NROW(x) <= fr.x && start(x)[1L] == end(x)[1L]) 
                {
                  dn1 <- start(x)[1L]
                  dn2 <- dn2[1 + (start(x)[2L] - 2 + seq_along(x))%%fr.x]
                  x <- matrix(format(x, ...), nrow = 1L, byrow = TRUE, 
                    dimnames = list(dn1, dn2))
                }
                else 
                {
                  start.pad <- start(x)[2L] - 1
                  end.pad <- fr.x - end(x)[2L]
                  dn1 <- start(x)[1L]:end(x)[1L]
                  x <- matrix(c(rep.int("", start.pad), format(x, ...), rep.int("", end.pad)), ncol = fr.x, 
                    byrow = TRUE, dimnames = list(dn1, dn2))
                }
            }
            else 
            {
                tx <- time(x)
                attributes(x) <- NULL
                names(x) <- tx
            }
        }
        else
            attr(x, "class") <- attr(x, "tsp") <- attr(x, "na.action") <- NULL
    }
    else 
    {
        if (calendar && fr.x > 1) 
        {
            tm <- time(x)
            t2 <- 1 + round(fr.x * ((tm + 0.001)%%1))
            p1 <- format(floor(zapsmall(tm)))
            rownames(x) <- if (fr.x == 12) 
                paste(month.abb[t2], p1, sep = " ")
            else paste(p1, if (fr.x == 4) 
                c("Q1", "Q2", "Q3", "Q4")[t2]
            else format(t2), sep = " ")
        }
        else
            rownames(x) <- format(time(x))
        attr(x, "class") <- attr(x, "tsp") <- attr(x, "na.action") <- NULL
    }
    return(as.data.frame(x))
}

print.forecast <- function(x ,...)
{
    print(as.data.frame(x))
}


summary.forecast <- function(object,...)
{
    cat(paste("\nForecast method:",object$method))
#    cat(paste("\n\nCall:\n",deparse(object$call)))
    cat(paste("\n\nModel Information:\n"))
    print(object$model)
    cat("\nIn-sample error measures:\n")
    print(accuracy(object))
    if(is.null(object$mean))
        cat("\n No forecasts\n")
    else
    {
        cat("\nForecasts:\n")
        print(object)
    }
}

plot.forecast <- function(x, include, plot.conf=TRUE, shaded=TRUE, shadebars=(length(x$mean)<5),
        shadecols=NULL, col=1, fcol=4, pi.col=1, pi.lty=2, ylim=NULL, main=NULL, ylab="",xlab="",...)
{
    if(is.element("x",names(x))) # Assume stored as x
        data=as.ts(x$x)
    else
    {
        data=NULL
        include=0
    }

    if(missing(include))
        include <- length(data)

    if(is.null(x$lower) | is.null(x$upper) | is.null(x$level))
        plot.conf=FALSE
        
    if(!shaded)
        shadebars <- FALSE

    # Extract components of predict if it exists
    # This occurs for the pegels functions.
    if(is.null(main))
        main=paste("Forecasts from ",x$method,sep="")

    freq <- frequency(data)
    strt <- start(data)
    n <- length(data)
    if(plot.conf)
    {
        upper <- as.matrix(x$upper)
        lower <- as.matrix(x$lower)
    }
    pred.mean <- x$mean
    # if(!is.null(lambda))  # undo Box-Cox transformation
    # {
        # pred.mean <- InvBoxCox(x$mean,lambda)
        # if(plot.conf)
        # {
            # lower <- InvBoxCox(lower,lambda)
            # upper <- InvBoxCox(upper,lambda)
        # }
        # xx <- InvBoxCox(data,lambda)
        # if(lambda<0 & plot.conf)
        # {
            # junk <- upper
            # upper <- lower
            # lower <- junk
        # }
    # }
    # else
        xx <- data
    # Remove final missing values
    nx <- max(which(!is.na(xx)))
    xxx <- xx[1:nx]
    include <- min(include,nx)
    if(is.null(ylim))
    {
        ylim <- range(c(xx[(n-include+1):n],pred.mean),na.rm=TRUE)
        if(plot.conf)
            ylim <- range(ylim,lower,upper,na.rm=TRUE)
    }
    npred <- length(pred.mean)
    plot(ts(c(xxx[(nx-include+1):nx], rep(NA, npred)), end=tsp(xx)[2] + (nx-n)/freq + npred/freq, f = freq),
        xlab=xlab,ylim=ylim,ylab=ylab,main=main,col=col,...)

    if(plot.conf)
    {
        xxx <- tsp(pred.mean)[1] - 1/freq + (1:npred)/freq            
        idx <- rev(order(x$level))
        nint <- length(x$level)
        for(i in 1:nint)
        {
			if(is.null(shadecols))
				shadecols <- heat.colors(length(x$level)+2)[switch(1+(length(x$level)>1),2,length(x$level):1)+1]
            if(shadebars)
			{
                for(j in 1:npred)
                {
                    polygon(xxx[j] + c(-0.5,0.5,0.5,-0.5)/freq, c(rep(lower[j,idx[i]],2),rep(upper[j,idx[i]],2)),
                        col = shadecols[i], border=FALSE)
                }
			}
            else if(shaded)
            {
                polygon(c(xxx,rev(xxx)), c(lower[,idx[i]],rev(upper[,idx[i]])),
                        col = shadecols[i], border=FALSE)
            }
            else if(npred == 1)
            {
                lines(xxx+c(-0.5,0.5)/freq,rep(lower[,idx[i]],2),col=pi.col,lty=pi.lty)
                lines(xxx+c(-0.5,0.5)/freq,rep(upper[,idx[i]],2),col=pi.col,lty=pi.lty)
            }
            else
            {
                lines(xxx,lower[,idx[i]],col=pi.col,lty=pi.lty)
                lines(xxx,upper[,idx[i]],col=pi.col,lty=pi.lty)
            }
        }
    }
    if(npred > 1 & !shadebars)
        lines(pred.mean, lty = 1,col=fcol)
    else
        points(pred.mean, col=fcol, pch=19)
    if(plot.conf)
        invisible(list(mean=pred.mean,lower=lower,upper=upper))
    else
        invisible(list(mean=pred.mean))
}

predict.default <- function(object, ...)
{
    forecast(object, ...)
}
