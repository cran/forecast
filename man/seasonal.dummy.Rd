\name{seasonaldummy}
\alias{seasonaldummy}
\alias{seasonaldummyf}
\title{Seasonal dummy variables}
\usage{seasonaldummy(x)
seasonaldummyf(x,h)
}

\arguments{
\item{x}{Seasonal time series}
\item{h}{Number of periods ahead to forecast}
}
\description{Returns matrix of dummy variables suitable for use in \code{\link[stats]{arima}} or 
\code{\link[stats]{lm}}. The last season is omitted and used as the control.}

\value{Numerical matrix with number of rows equal to the \code{length(x)} and number of columns equal to \code{frequency(x)-1}.
}

\author{Rob J Hyndman}
\examples{plot(ldeaths)
plot(ldeaths)
month <- seasonaldummy(ldeaths)
deaths.lm  <- lm(ldeaths ~ month)
tsdisplay(residuals(deaths.lm))
ldeaths.fcast <- predict(deaths.lm,
   data.frame(month=I(seasonaldummyf(ldeaths,36))),
   se.fit=TRUE)
}
\keyword{ts}
