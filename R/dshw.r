###Double seasonal Holt winters method with or without AR adjustment
#period1 is the short period
#period2 is the long period

dshw <- function(y, period1, period2, h=2*max(period1,period2), alpha=NULL, beta=NULL, gamma=NULL, omega=NULL, phi=NULL, lambda=NULL, armethod=TRUE)
{
  if(!armethod)
  {
    phi <- 0
  }
  if(period1 > period2)
  {
    tmp <- period2
    period2 <- period1
    period1 <- tmp
  }
  if(period1 < 1 | period1 == period2)
  stop("Inappropriate periods")

  if (!is.null(lambda))
  {
    origy <- y
    y <- BoxCox(y, lambda)
  }

  pars <- rep(NA,5)
  if(!is.null(alpha))
    pars[1] <- alpha
  if(!is.null(beta))
    pars[2] <- beta
  if(!is.null(gamma))
    pars[3] <- gamma
  if(!is.null(omega))
    pars[4] <- omega
  if(!is.null(phi))
    pars[5] <- phi

  # Estimate parameters
  if(sum(is.na(pars)) > 0)
  {
    pars <- parameters_optim_dshw(y,period1,period2,pars)
    alpha <- pars[1]
    beta <- pars[2]
    gamma <- pars[3]
    omega <- pars[4]
    phi <- pars[5]
  }

  ## Allocate space
  n <- length(y)
  yhat <- numeric(n)

  ## Starting values
  I <- sindex(y,period1)
  wstart <- sindex(y,period2)
  wstart <- wstart / rep(I,(length(wstart)/length(I)))
  w <- wstart
  x <- c(0,diff(y[1:period2]))
  t <- t.start <- mean(((y[1:period2]- y[(period2+1):(2*period2)])/period2 ) + x )/2
  s <- s.start <- (mean(y[1:(2*period2)])-(period2+0.5)*t)

  ## In-sample fit
  for(i in 1: n)
  {
    yhat[i] <- (s+t) * I[i]*w[i]
    snew <- alpha*(y[i]/(I[i]*w[i]))+(1-alpha)*(s+t)
    tnew <- beta*(snew-s)+(1-beta)*t
    I[i+period1] <- gamma*(y[i]/(snew*w[i])) + (1-gamma)*I[i]
    w[i+period2] <- omega*(y[i]/(snew*I[i])) + (1-omega)*w[i]
    s <- snew
    t <- tnew
  }

	# Forecasts
  fcast <- (s + (1:h)*t) * rep(I[n+(1:period1)],h/period1 + 1)[1:h] * rep(w[n+(1:period2)],h/period2 + 1)[1:h]
  
  # Calculate MSE and MAPE
  yhat <- ts(yhat)
  tsp(yhat) <- tsp(y)
  e <- y - yhat
  if(armethod)
	{
		yhat <- yhat + phi * c(0,e[-n])
		e <- y - yhat
	  fcast <- fcast + phi^(1:h)*e[n]
	}
	mse <- mean(e^2)
	mape <- mean(abs(e)/y)*100
	
  fcast <- ts(fcast,f=frequency(y),s=tsp(y)[2]+1/tsp(y)[3])

  if(!is.null(lambda))
  {
    y <- origy
    fcast <- InvBoxCox(fcast,lambda)
    yhat <- InvBoxCox(yhat,lambda)
  }

  return(structure(list(mean=fcast,method="DSHW",x=y,residuals=e,fitted=yhat,
              model=list(mape=mape,mse=mse,alpha=alpha,beta=beta, gamma=gamma,omega=omega,phi=phi,
                     l0=s.start,b0=t.start,s10=wstart,s20=I)),class="forecast"))

}

###------------------------------------------------------------------------------------------------
###Double seasonal holt winters parameter optimisation

parameters_optim_dshw <- function(y, period1, period2, pars)
{
  start <- c(0.1,0.01,0.001,0.001,0.0)[is.na(pars)]
  out <- optim(start, dshw.mse, y=y, period1=period1, period2=period2, pars=pars)
  pars[is.na(pars)] <- out$par
  return(pars)
}

dshw.mse <- function(par, y, period1, period2, pars)
{
  pars[is.na(pars)] <- par
  if(max(pars) > 0.99 | min(pars) < 0 | pars[5] > .9)
    return(1e10)
  else
    return(dshw(y, period1, period2, h=1, pars[1], pars[2], pars[3], pars[4], pars[5], armethod=(abs(pars[5]) >1e-7))$model$mse)
}

###------------------------------------------------------------------------------------------------
###Calculating seasonal indexes

sindex <- function(y,p)
{
  require(zoo)
  n <- length(y)
  n2 <- 2*p
  y2 <- y[1:n2]
  average <- numeric(n)
  simple_ma <- rollmean(y2, p)
  if (identical(p%%2,0))
  {
    centered_ma <- rollmean(simple_ma[1:(n2-p+1)],2)
    average[p/2 + 1:p] <- y2[p/2 + 1:p]/centered_ma[1:p]
    si <- average[c(p+(1:(p/2)),(1+p/2):p)]
  }
  else
  {
    average[(p-1)/2 + 1:p] <- y2[(p-1)/2 + 1:p]/simple_ma[1:p]
    si <- average[c(p+(1:((p-1)/2)),(1+(p-1)/2):p)]
  } 
  return(si)
}
