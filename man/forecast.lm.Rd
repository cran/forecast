\name{forecast.lm}
\alias{forecast.lm}
\title{Forecasting using lm objects}
\usage{\method{forecast}{lm}(object, newdata, level=c(80,95), fan=FALSE, ...)
}

\arguments{
\item{object}{An object of class "\code{lm}". Usually the result of a call to \code{\link[stats]{lm}}.}
\item{newdata}{An optional data frame or matrix in which to look for variables with which to predict.}
\item{level}{Confidence level for prediction intervals.}
\item{fan}{If TRUE, level is set to seq(50,99,by=1). This is suitable for fan plots.}
\item{...}{Other arguments passed to \code{\link[stats]{predict.lm}()}.}
}

\description{Returns forecasts obtained from linear models with explanatory variables specified by \code{newdata}.}

\details{\code{forecast.lm} produces predicted values, obtained by evaluating the regression function in \code{newdata}. It is largely a wrapper to \code{\link[stats]{predict.lm}} and returns the output in the form of a \code{forecast} object. If the original data are time series, the function will make the forecasts time series with the appropriate parameters.}

\value{An object of class "\code{forecast}".

The function \code{summary} is used to obtain and print a summary of the
results, while the function \code{plot} produces a plot of the forecasts and prediction intervals.

The generic accessor functions \code{fitted.values} and \code{residuals} extract useful features of
the value returned by \code{forecast.lm}.

An object of class \code{"forecast"} is a list containing at least the following elements:
\item{model}{A list containing information about the fitted model}
\item{method}{The name of the forecasting method as a character string}
\item{mean}{Point forecasts as a time series}
\item{lower}{Lower limits for prediction intervals}
\item{upper}{Upper limits for prediction intervals}
\item{level}{The confidence values associated with the prediction intervals}
\item{x}{The original time series (either \code{object} itself or the time series used to create the model stored as \code{object}).}
\item{residuals}{Residuals from the fitted model. That is x minus fitted values.}
\item{fitted}{Fitted values (one-step forecasts)}
}

\seealso{\code{\link[stats]{predict.lm}}, \code{\link[stats]{lm}}.}

\author{Rob J Hyndman}

\examples{
#To be completed.
}
\keyword{stats}
