forecast.lm <- function(object, newdata, level=c(80,95), fan=FALSE, ...)
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
    
  newdata <- as.data.frame(newdata)
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