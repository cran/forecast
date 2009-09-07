### Functions to handle seasonality

monthdays <- function(x)
{
    if(!is.ts(x))
        stop("Not a time series")
    f <- frequency(x)
    if(f==12)
        days <- c(31,28,31,30,31,30,31,31,30,31,30,31)
    else if(f==4)
        days <- c(90,91,92,92)
    else
        stop("Not monthly or quarterly data")
    nyears <- round(length(x)/f+1)+1
    years <- (1:nyears) + (start(x)[1] - 1)
    leap.years <- ((years %% 4 == 0) & !(years %% 100 ==0 & years %% 400 != 0))[1:nyears]
    dummy <- t(matrix(rep(days,nyears),nrow=f))
    if(f==12)
        dummy[leap.years,2] <- 29
    else
        dummy[leap.years,1] <- 91
    xx <- c(t(dummy))[start(x)[2]-1+(1:length(x))]
    return(ts(xx,start=start(x),f=f))
}

sindexf <- function(object,h)
{
    if(class(object)=="stl")
    {
        ss <- object$time.series[,1]
        m <- frequency(ss)
        ss <- ss[length(ss)-(m:1)+1]
        tsp.x <- tsp(object$time.series)
    }
    else if(class(object)=="decomposed.ts")
    {
        ss <- object$figure
        m <- frequency(object$seasonal)
        n <- length(object$trend)
        ss <- rep(ss,n/m+1)[1:n]
        ss <- ss[n-(m:1)+1]
        tsp.x <- tsp(object$seasonal)
    }
    else
        stop("Object of unknown class")
    out <- ts(rep(ss,h/m+1)[1:h], f=m, start=tsp.x[2]+1/m)

    return(out)
}

seasadj <- function(object)
{
    if(class(object)=="stl")
        return(object$time.series[,2]+object$time.series[,3])
    else if(class(object)=="decomposed.ts")
            return(object$trend+object$random)
    else
        stop("Object of unknown class")
}

seasonaldummy <- function(x)
{
    if(!is.ts(x))
        stop("Not a time series")
    else
        fr.x <- frequency(x)
    if(fr.x==1)
        stop("Non-seasonal time series")
    dummy <- as.factor(cycle(x))
    dummy.mat <- matrix(0,ncol=frequency(x)-1,nrow=length(x))
    nrow <- 1:length(x)
    for(i in 1:(frequency(x)-1))
        dummy.mat[dummy==paste(i),i] =1
    colnames(dummy.mat) <- if (fr.x == 12)
                month.abb[1:11]
            else if(fr.x == 4)
                c("Q1", "Q2", "Q3")
            else paste("S",1:(fr.x-1),sep="")

    return(dummy.mat)
}

seasonaldummyf <- function(x, h)
{
    if(!is.ts(x))
        stop("Not a time series")
    f=frequency(x)
    return(seasonaldummy(ts(rep(0,h),start=tsp(x)[2]+1/f,freq=f)))
}
