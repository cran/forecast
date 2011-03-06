simulate.ets <- function(object, nsim=length(object$x), seed=NULL, future=TRUE, bootstrap=FALSE,...)
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
	if(is.null(tsp(object$x)))
		object$x <- ts(object$x,f=1,s=1)
	
	if(future) 
		initstate <- object$state[length(object$x)+1,] 
		#initstate <- object$state[1,] ??
	else # choose a random starting point
		initstate <- object$state[sample(1:length(object$x),1),]
	
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
        PACKAGE="forecast")[[11]],f=object$m,s=tsp(object$x)[2]+1/tsp(object$x)[3])
    if(is.na(tmp[1]))
        stop("Problem with multiplicative damped trend")
    tmp
}


# Simulate ARIMA model starting with observed data x

myarima.sim <- function (model, n, x, e, ...) 
{
    data <- x
	model$x <- x
	if(is.null(tsp(data)))
		data <- ts(data,f=1,s=1)
    if (!is.list(model)) 
        stop("'model' must be list")
    p <- length(model$ar)
    if (p) 
    {
        minroots <- min(Mod(polyroot(c(1, -model$ar))))
        if (minroots <= 1) 
            stop("'ar' part of model is not stationary")
    }
    q <- length(model$ma)
    n.start <- length(x)
    d <- 0
    if (!is.null(ord <- model$order)) 
    {
        if (length(ord) != 3) 
            stop("'model$order' must be of length 3")
        if (p != ord[1]) 
            stop("inconsistent specification of 'ar' order")
        if (q != ord[3]) 
            stop("inconsistent specification of 'ma' order")
        d <- ord[2]
        if (d != round(d) || d < 0) 
            stop("number of differences must be a positive integer")
    }
    start.innov <- residuals(model)
    innov <- e
    x <- ts(c(start.innov, innov), start = 1 - n.start)
    if (length(model$ma)) 
        x <- filter(x, c(1, model$ma), sides = 1)
    if (length(model$ar)) 
        x <- filter(x, model$ar, method = "recursive")
    if(d == 0) # Adjust to ensure end matches approximately
    {
        # Last 20 diffs
        if(n.start >= 20)
            xdiff <- (model$x - x[1:n.start])[n.start-(19:0)]
        else
            xdiff <- model$x - x[1:n.start]
        # If all same sign, choose last
        if(all(sign(xdiff)==1) | all(sign(xdiff)==-1))
            xdiff <- xdiff[length(xdiff)]
        else # choose mean.
            xdiff <- mean(xdiff)
        x <- x + xdiff
    }
    if (n.start > 0) 
        x <- x[-(1:n.start)]
    if (d > 0)
        x <- diffinv(x, differences = d,xi=data[length(data)-(d:1)+1])[-(1:2)]
    x <- ts(x[1:n],f=frequency(data),s=tsp(data)[2]+1/tsp(data)[3])
    return(x)    
}

simulate.Arima <- function(object, nsim=length(object$x), seed=NULL, xreg=NULL, future=TRUE, bootstrap=FALSE, ...)
{
	if(sum(object$arma[c(3,4,7)])>0)
		stop("Seasonal ARIMA simulation is not yet implemented")
	if (!exists(".Random.seed", envir = .GlobalEnv))
        runif(1)
    if (is.null(seed))
        RNGstate <- .Random.seed
    else
    {
        R.seed <- .Random.seed
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    order <- object$arma[c(1, 6, 2)]
    if(order[1]>0)
        ar <- object$model$phi[1:order[1]]
    else
        ar <- NULL
    if(order[3]>0)
        ma <- object$model$theta[1:order[3]]
    else
        ma <- NULL
    if(object$arma[2] != length(ma))
        stop("MA length wrong")
	else if(object$arma[1] != length(ar))
		stop("AR length wrong")
    model <- list(order=object$arma[c(1, 6, 2)],ar=ar,ma=ma,sd=sqrt(object$sigma2),residuals=residuals(object))
    if (is.element("x", names(object))) 
        x <- object$x
    else 
        x <- object$x <- eval.parent(parse(text = object$series))
	if(is.null(tsp(x)))
		x <- ts(x,f=1,s=1)

    n <- length(x)
	d <- order[2]
	if(bootstrap)
		e <- sample(model$residuals,nsim+d,replace=TRUE)
	else
		e <- rnorm(nsim+d, 0, model$sd)
		
	use.drift <- is.element("drift", names(object$coef))
    usexreg <- (!is.null(xreg) | use.drift)
    if (!is.null(xreg)) 
	{
        xreg <- as.matrix(xreg)
        if(nrow(xreg) < nsim)
			stop("Not enough rows in xreg")
		else
			xreg <- xreg[1:nsim,]
    }
    if (use.drift) 
	{
		dft <- as.matrix(1:nsim) + n
        xreg <- cbind(xreg, dft)
    }
	narma <- sum(object$arma[1L:4L])
	if(length(object$coef) > narma)
	{
		if (names(object$coef)[narma + 1L] == "intercept") 
		{
			xreg <- cbind(intercept = rep(1, nsim), xreg)
			object$xreg <- cbind(intercept = rep(1, n), object$xreg)
		}
		if(!is.null(xreg))
		{
			xm <- if (narma == 0) 
					drop(as.matrix(xreg) %*% object$coef)
				else 
					drop(as.matrix(xreg) %*% object$coef[-(1L:narma)])
			oldxm <- if(narma == 0)
						drop(as.matrix(object$xreg) %*% object$coef)
					else 
						drop(as.matrix(object$xreg) %*% object$coef[-(1L:narma)])
		}
	}
	else 
		xm <- oldxm <- 0
		
	if(future)
		sim <- myarima.sim(model,nsim,x-oldxm,e=e) + xm
	else
		sim <- arima.sim(model,nsim,innov=e) + xm
	return(sim)
}

simulate.ar <- function(object, nsim=object$n.used, seed=NULL, future=TRUE, bootstrap=FALSE, ...)
{
    if (!exists(".Random.seed", envir = .GlobalEnv))
        runif(1)
    if (is.null(seed))
        RNGstate <- .Random.seed
    else
    {
        R.seed <- .Random.seed
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    model <- list(ar=object$ar,sd=sqrt(object$var.pred),residuals=object$resid)
    x.mean <- object$x.mean
    if (!is.element("x", names(object))) 
        object$x <- eval.parent(parse(text = object$series))
	if(is.null(tsp(object$x)))
		object$x <- ts(object$x,f=1,s=1)
    object$x <- eval.parent(parse(text = object$series)) - x.mean
	if(bootstrap)
		e <- sample(model$residuals,nsim,replace=TRUE)
	else
		e <- rnorm(nsim, 0, model$sd)
	if(future)
		return(myarima.sim(model,nsim,x=object$x,e=e) + x.mean)
	else
		return(arima.sim(model,nsim,innov=e) + x.mean)
}

simulate.fracdiff <- function(object, nsim=object$n, seed=NULL, future=TRUE, bootstrap=FALSE, ...)
{
	if (is.element("x", names(object))) 
        x <- object$x
    else 
		x <- object$x <- eval.parent(parse(text = as.character(object$call)[2]))
	if(is.null(tsp(x)))
		x <- ts(x,f=1,s=1)
    
    # Strip initial and final missing values
    xx <- na.ends(x)
    n <- length(xx)
    
    # Remove mean
    meanx <- mean(xx)
    xx <- xx - meanx

    y <- undo.na.ends(x,diffseries(xx, d = object$d))
    fit <- arima(y, order = c(length(object$ar), 0, length(object$ma)), 
        include.mean = FALSE, fixed = c(object$ar, -object$ma))
	# Simulate ARMA
	ysim <- simulate(fit,nsim,seed,future=future,bootstrap=bootstrap)
	# Undo differencing
    return(unfracdiff(xx,ysim,n,nsim,object$d))
	# bin.c <- (-1)^(0:(n + nsim)) * choose(object$d, (0:(n + nsim)))
    # b <- numeric(n)
    # xsim <- LHS <- numeric(nsim)
    # RHS <- cumsum(ysim)
    # bs <- cumsum(bin.c[1:nsim])
    # b <- bin.c[(1:n) + 1]
    # xsim[1] <- RHS[1] <- ysim[1] - sum(b * rev(xx))
    # for (k in 2:nsim) 
	# {
        # b <- b + bin.c[(1:n) + k]
        # RHS[k] <- RHS[k] - sum(b * rev(xx))
        # LHS[k] <- sum(rev(xsim[1:(k - 1)]) * bs[2:k])
        # xsim[k] <- RHS[k] - LHS[k]
    # }
	# tspx <- tsp(x)
	# return(ts(xsim,f=tspx[3],s=tspx[2]+1/tspx[3]))
}
