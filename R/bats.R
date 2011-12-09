# Author: srazbash
###############################################################################
#setwd("/Volumes/NO\ NAME/SlavaBATS/")
#setwd("E:/SlavaBATS")
#source("fitBATS.R", verbose=TRUE)
#source("forecastBATS.R", verbose=TRUE)
#source("makeMatrices.R", verbose=TRUE)
#source("checkAdmissibility.R", verbose=TRUE)
#source("makeParamVector.R", verbose=TRUE)
#source("adjustSeasonalSeeds.R", verbose=TRUE)
#source("getBATS.R", verbose=TRUE)

filterSpecifics <- function(y, box.cox, trend, damping, seasonal.periods, use.arma.errors, ...) {
	if((trend == FALSE) & (damping == TRUE)) {
		return(list(AIC=Inf))
	}
	#printCASE(box.cox, trend, damping, seasonal.periods, NULL, NULL, 0, 0)
	first.model <- fitSpecificBATS(y, use.box.cox=box.cox, use.beta=trend, use.damping=damping, seasonal.periods=seasonal.periods)
	if(!is.null(seasonal.periods)) {
		#printCASE(box.cox, trend, damping, NULL, NULL, NULL, 0, 0)
		non.seasonal.model <- fitSpecificBATS(y, use.box.cox=box.cox, use.beta=trend, use.damping=damping, seasonal.periods=NULL)
		if(first.model$AIC > non.seasonal.model$AIC) {
			seasonal.periods <- NULL
			first.model <- non.seasonal.model
		}
	}
	if(use.arma.errors) { 
		##Turn off warnings
		old.warning.level <- options()$warn
		options(warn=-1)
		arma <- auto.arima(as.numeric(first.model$errors), d=0, ...)
		###Re-enable warnings
		options(warn=old.warning.level)
		p <- arma$arma[1]
		q <- arma$arma[2]
		if((p != 0) | (q != 0)) { #Did auto.arima() find any AR() or MA() coefficients?
			if(p != 0) {
				ar.coefs <- numeric(p)
			} else {
				ar.coefs <- NULL
			}
			if(q != 0) {
				ma.coefs <- numeric(q)
			} else {
				ma.coefs <- NULL
			}
			starting.params <- first.model$parameters
			#printCASE(box.cox, trend, damping, seasonal.periods, ar.coefs, ma.coefs, p, q)
			second.model <- fitSpecificBATS(y, use.box.cox=box.cox, use.beta=trend, use.damping=damping, seasonal.periods=seasonal.periods, ar.coefs=ar.coefs, ma.coefs=ma.coefs)
			if(second.model$AIC < first.model$AIC) {
				return(second.model)
			} else {
				return(first.model)
			}
		} else { #Else auto.arima() did not find any AR() or MA()coefficients
			return(first.model)
		}
	} else {
		return(first.model)
	}
}


bats <- function(y, use.box.cox=NULL, use.trend=NULL, use.damped.trend=NULL, seasonal.periods=NULL, use.arma.errors=TRUE, ...) {
	if(any((y <= 0))) {
		stop("BATS requires positive data")
	}
  origy <- y
	if(any(class(y) == "msts")) {
		seasonal.periods <- attr(y,"msts")
	} else if(class(y) == "ts") {
		seasonal.periods <- frequency(y)
	}
  y <- as.numeric(y)
	best.aic <- NULL
	if(is.null(use.box.cox)) {
		use.box.cox <- c(FALSE, TRUE)
	} 
	if(is.null(use.trend)) {
		use.trend <- c(FALSE, TRUE)
	} else if(use.trend == FALSE) {
		use.damped.trend <- FALSE
	}
	if(is.null(use.damped.trend)) {
		use.damped.trend <- c(FALSE, TRUE)
	}
	for(box.cox in use.box.cox) {
		for(trend in use.trend) {
			for(damping in use.damped.trend) {
				current.model <- filterSpecifics(y, box.cox=box.cox, trend=trend, damping=damping, seasonal.periods=seasonal.periods, use.arma.errors=use.arma.errors, ...)
				if(!is.null(best.aic)) {
					if(current.model$AIC < best.aic) {
						best.aic <- current.model$AIC
						best.model <- current.model
					}
				} else {
					best.model <- current.model
					best.aic <- best.model$AIC
				}
			}
		}
	}
	best.model$call <- match.call()
	if(best.model$optim.return.code != 0) {
		warning("optim() did not converge.")
	}
  
  # Add ts attributes
  if(!any(class(origy) == "ts"))
  {
    if(is.null(seasonal.periods))
      origy <- ts(origy,s=1,f=1)
    else 
      origy <- msts(origy,seasonal.periods)
  }
  attributes(best.model$fitted.values) <- attributes(best.model$errors) <- attributes(origy)
  best.model$y <- origy
  
	return(best.model)
}

print.bats <- function(x,...) {
	cat("\n")
	cat(makeText(x))
	cat("\n")
#	cat("BATS( {")
#	if(!is.null(x$lambda)) {
#		cat(x$lambda)
#	} else {
#		cat("1")
#	}
#	cat("}, {")
#	if(!is.null(x$ar.coefficients)) {
#		cat(length(x$ar.coefficients))
#	} else {
#		cat("0")
#	}
#	cat(", ")
#	if(!is.null(x$ma.coefficients)) {
#		cat(length(x$ma.coefficients))
#	} else {
#		cat("0")
#	}
#	cat("}, {")
#	if(!is.null(x$damping.parameter)) {
#		cat(x$damping.parameter)
#	} else {
#		cat("0")
#	}
#	
#	if(!is.null(x$seasonal.periods)) {
#		cat("}, { ")
#		for(i in x$seasonal.periods) {
#			cat(i)
#			if(i != x$seasonal.periods[length(x$seasonal.periods)]) {
#				cat(", ")
#			} else {
#				cat("})")
#			}
#		}
#	} else {
#		cat("})\n\n")	
#	}
	cat("\nCall: ")
	print(x$call)
	cat("\nParameters:\n")
	cat("\nBox-Cox Parameter: ")
	cat(x$lambda)
	cat("\nAlpha: ")
	cat(x$alpha)
	cat("\nBeta: ")
	cat(x$beta)
	cat("\nDamping Parameter: ")
	cat(x$damping.parameter)
	cat("\nGamma Values: ")
	cat(x$gamma.values)
	cat("\nAR() Coefficients: ")
	cat(x$ar.coefficients)
	cat("\nMA() Coefficients: ")
	cat(x$ma.coefficients)
	cat("\n\n")
	cat("\nSeed States:\n")
	print(x$seed.states)
	
	cat("\nSigma: ")
	cat(sqrt(x$variance))
	
	cat("\n\nAIC: ")
	cat(x$AIC)
	cat("\n\n")
	
}

residuals.bats <- function(object, ...) {
	object$errors
}
