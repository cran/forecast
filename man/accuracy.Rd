\name{accuracy}
\alias{accuracy}
\title{Accuracy measures for forecast model}
\usage{accuracy(f, x, test=1:length(x),lambda=NULL)
}

\arguments{
\item{f}{An object of class \code{"forecast"}, or a numerical vector containing forecasts.}
\item{x}{An optional numerical vector containing actual values of the same length as object.}
\item{test}{Indicator of which elements of x and f to test.}
\item{lambda}{Box-Cox parameter. If present, function backtransforms data and forecasts using \code{\link{InvBoxCox}}.}
}

\description{Returns range of summary measures of the forecast accuracy. If \code{x} is provided, the function measures out-of-sample forecast accuracy
based on x-f. If \code{x} is not provided, the function produces in-sample accuracy measures of the one-step forecasts based on f["x"]-fitted(f).
All measures are defined and discussed in Hyndman and Koehler (2006).}

\value{Vector giving forecast accuracy measures.
}

\references{Hyndman, R.J. and Koehler, A.B. (2006) "Another look at measures of forecast accuracy". \emph{International Journal of Forecasting},
\bold{22}(4).}

\author{Rob J Hyndman}
\examples{
fit1 <- rwf(EuStockMarkets[1:200,1],h=100)
fit2 <- meanf(EuStockMarkets[1:200,1],h=100)
accuracy(fit1)
accuracy(fit2)
accuracy(fit1,EuStockMarkets[201:300,1])
accuracy(fit2,EuStockMarkets[201:300,1])
plot(fit1)
lines(EuStockMarkets[1:300,1])
}
\keyword{ts}
