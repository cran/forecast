tslm <- function(formula,data,...)
{
    if(missing(data)) # Grab first variable
    {
        x <- get(as.character(formula)[2])
        data <- NULL
    }
    else
        x <- data[,1]
    
    if(!is.ts(x))
        stop("Not time series data")
    tspx <- tsp(x)
    # Add trend and seasonal to data frame
    trend <- time(x)
    season <- as.factor(cycle(x))
    if(is.null(data))
        data <- data.frame(trend,season)
    else
        data <- data.frame(data,trend,season)
    fit <- lm(formula,data=data,...)
    if(!is.null(fit$call$subset))
    {
        j <- eval(fit$call$subset)
        data <- data[j,]
        # Try to figure out times for subset. Assume they are contiguous.
        timesx <- timesx[j]
        tspx <- c(min(timesx),max(timesx),tspx[3])
    }
    fit$data <- ts(data)
    fit$residuals <- ts(fit$residuals)
    fit$fitted.values <- ts(fit$fitted.values)
    tsp(fit$data) <- tsp(fit$residuals) <- tsp(fit$fitted.values) <- tspx
    return(fit)
}

forecast.lm <- function(object, newdata, level=c(80,95), fan=FALSE, h=10, ...)
{
  if (fan) 
    level <- seq(51, 99, by = 3)
  else 
  {
    if (min(level) > 0 & max(level) < 1) 
      level <- 100 * level
    else if (min(level) < 0 | max(level) > 99.99) 
      stop("Confidence limit out of range")
  }
  if(!is.null(object$data))
    origdata <- object$data
  else
    origdata <- eval(object$call$data)
  if(is.element("ts",class(origdata)))
  {
    tspx <- tsp(origdata)
    timesx <- time(origdata)
  }
  else
    tspx <- NULL
  if(!is.null(object$call$subset))
  {
    j <- eval(object$call$subset)
    origdata <- origdata[j,]
    if(!is.null(tspx))
    {
        # Try to figure out times for subset. Assume they are contiguous.
        timesx <- timesx[j]
        tspx <- tsp(origdata) <- c(min(timesx),max(timesx),tspx[3])
    }
  }
  # Add trend and seasonal to data frame
  if(!missing(newdata))
    h <- nrow(newdata)
  if(!is.null(tspx))
  {
    x <- ts(1:h, start=tspx[2]+1/tspx[3], frequency=tspx[3])
    trend <- time(x)
    season <- as.factor(cycle(x))
    if(!missing(newdata))
        newdata <- data.frame(as.data.frame(newdata),trend,season)
    else
        newdata <- data.frame(trend,season)
  }  
  out <- list()
  nl <- length(level)
  for(i in 1:nl)
    out[[i]] <- predict(object, newdata=newdata, se.fit=TRUE, interval="prediction", level=level[i]/100, ...)
  fcast <- list(model=object,mean=out[[1]]$fit[,1],lower=out[[1]]$fit[,2],upper=out[[1]]$fit[,3],level=level,x=object$model[,1])
  fcast$method <- "Linear regression model"
  fcast$residuals <- residuals(object)
  fcast$fitted <- fitted(object)
  if(!is.null(tspx))
  {
     fcast$x <- ts(fcast$x)
     fcast$residuals <- ts(fcast$residuals)
     fcast$fitted <- ts(fcast$fitted)
     tsp(fcast$x) <- tsp(fcast$residuals) <- tsp(fcast$fitted) <- tspx
  }
  if(nl > 1)
  {
      for(i in 2:nl)
      {
          fcast$lower <- cbind(fcast$lower,out[[i]]$fit[,2])
          fcast$upper <- cbind(fcast$upper,out[[i]]$fit[,3])
      }
  }
  if(!is.null(tspx))
  {
    fcast$mean <- ts(fcast$mean, start=tspx[2]+1/tspx[3],frequency=tspx[3])
    fcast$upper <- ts(fcast$upper, start=tspx[2]+1/tspx[3],frequency=tspx[3])
    fcast$lower <- ts(fcast$lower, start=tspx[2]+1/tspx[3],frequency=tspx[3])
  }
  return(structure(fcast,class="forecast"))
}