\name{bats}
\alias{bats}

\title{BATS model (Exponential smoothing state space model with Box-Cox transformation, ARMA errors, Trend and Seasonal components)}
\usage{bats(y, use.box.cox=NULL, use.trend=NULL, use.damped.trend=NULL, 
  seasonal.periods=NULL, use.arma.errors=TRUE, ...)}

\arguments{
\item{y}{The time series to be forecast. Can be \code{numeric}, \code{msts} or \code{ts}. Only univariate time series is supported.}
\item{use.box.cox}{\code{TRUE/FALSE} indicates whether to use the Box-Cox transformation or not. If \code{NULL} then both are tried and the best fit is selected by AIC.}
\item{use.trend}{\code{TRUE/FALSE} indicates whether to include a trend or not. If \code{NULL} then both are tried and the best fit is selected by AIC.}
\item{use.damped.trend}{\code{TRUE/FALSE} indicates whether to include a damping parameter in the trend or not. If \code{NULL} then both are tried and the best fit is selected by AIC.}
\item{seasonal.periods}{If \code{y} is a numeric then seasonal periods can be specified with this parameter.}
\item{use.arma.errors}{\code{TRUE/FALSE} indicates whether to include ARMA errors or not. If \code{NULL} then both are tried and the best fit is selected by AIC.}
\item{...}{Additional parameters to be passed to \code{auto.arima} when choose an ARMA(p, q) model for the errors.}
}

\description{Fits a BATS model applied to \code{y}.}

\value{An object of class "\code{bats}".

The generic accessor functions \code{fitted.values} and \code{residuals} extract useful features of
the value returned by \code{ets} and associated functions.}

\references{De Livera, A.M., Hyndman, R.J., & Snyder, R. D. (2012), Forecasting time series with complex seasonal patterns using exponential smoothing, \emph{Journal of the American Statistical Association}, to appear. 
}

\author{Rob J Hyndman and Slava Razbash}

\examples{
fit <- bats(USAccDeaths)
plot(forecast(fit))
\dontrun{
taylor.fit <- bats(taylor)
plot(forecast(taylor.fit))
}
}

\keyword{ts}

