simulate.ets <- function(object, nsim=length(object$x), seed=NULL, initstate=object$state[1,], bootstrap=TRUE,...)
{
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
        runif(1)
    if (is.null(seed)) 
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }

    if(bootstrap)
        e <- sample(object$residuals,nsim,replace=TRUE)
    else
        e <- rnorm(nsim,0,sqrt(object$sigma))
    if(object$components[1]=="M")
        e <- pmax(-1,e)
    tmp <- ts(.C("etssimulate",
            as.double(initstate),
            as.integer(object$m),
            as.integer(switch(object$components[1],"A"=1,"M"=2)),
            as.integer(switch(object$components[2],"N"=0,"A"=1,"M"=2)),
            as.integer(switch(object$components[3],"N"=0,"A"=1,"M"=2)),
            as.double(object$par["alpha"]),
            as.double(ifelse(object$components[2]=="N",0,object$par["beta"])),
            as.double(ifelse(object$components[3]=="N",0,object$par["gamma"])),
            as.double(ifelse(object$components[4]=="FALSE",1,object$par["phi"])),
            as.integer(nsim),
            as.double(numeric(nsim)),
            as.double(e),
        PACKAGE="forecast")[[11]],f=object$m,s=1)
    if(is.na(tmp[1]))
        stop("Problem with multiplicative damped trend")
    tmp
}

#
#simulate.Arima <- function(object, nsim=length(object$x), seed=NULL, ...)
#{
#    if (!exists(".Random.seed", envir = .GlobalEnv))
#        runif(1)
#    if (is.null(seed))
#        RNGstate <- .Random.seed
#    else
#    {
#        R.seed <- .Random.seed
#        set.seed(seed)
#        RNGstate <- structure(seed, kind = as.list(RNGkind()))
#        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
#    }
#    model <- list(ar=object$model$phi,ma=object$model$theta,sd=sqrt(object$sigma2))
#    tmp <- arima.sim(model,nsim,...)
#
##    attr(tmp, "seed") <- RNGstate
#    tmp
#}
