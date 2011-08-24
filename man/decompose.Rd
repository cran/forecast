\name{decompose}
\alias{decompose}
\alias{plot.decomposed.ts}
\title{
  Classical Seasonal Decomposition by Moving Averages
}
\description{
  Decompose a time series into seasonal, trend and irregular components using moving averages.  Deals with additive or multiplicative
  seasonal component.
}
\usage{
decompose(x, type = c("additive", "multiplicative"), filter = NULL)
}
\arguments{
  \item{x}{A time series.}
  \item{type}{The type of seasonal component. Can be abbreviated.}
  \item{filter}{A vector of filter coefficients in reverse time order (as for
    AR or MA coefficients), used for filtering out the seasonal
    component. If \code{NULL}, a moving average with symmetric window is
  performed.}
}
\details{
  The additive model used is:
  \deqn{Y_t = T_t + S_t + e_t}{Y[t] = T[t] + S[t] + e[t]}
  The multiplicative model used is:
  \deqn{Y_t = T_t\,S_t\, e_t}{Y[t] = T[t] * S[t] * e[t]}
  
  The function first determines the trend component using a moving
  average (if \code{filter} is \code{NULL}, a symmetric window with equal
  weights is used), and removes it from the time series.  Then, 
  the seasonal figure is computed by averaging, for each time unit,
  over all periods. The seasonal figure is then centered.
  Finally, the error component is determined by
  removing trend and seasonal figure (recycled as needed) from the original
  time series.
}
\note{
This function is a modification of the \code{\link[stats]{decompose}} function in the \code{stats} package.  In this version, the seasonal component is not incorrectly truncated, the trend is correctly calculated for odd frequencies, and the data is returned as part of the object. It is hoped that this function will replace the version in \code{stats} package, and then it will be removed from the \code{forecast} package.
}
\value{
  An object of class \code{"decomposed.ts"} with following components:
  \item{seasonal}{The seasonal component (i.e., the repeated seasonal figure)}
  \item{figure}{The estimated seasonal figure only}
  \item{trend}{The trend component}
  \item{random}{The remainder part}
  \item{type}{The value of \code{type}}
}
\references{
  M. Kendall and A. Stuart (1983)
  The Advanced Theory of Statistics, Vol.3,
  \emph{Griffin}, 410--414.
}
\author{
  David Meyer \email{David.Meyer@wu-wien.ac.at}. Revised by Rob J Hyndman \email{Rob.Hyndman@monash.edu}.
}
\seealso{\code{\link[stats]{decompose}}}

\examples{
m <- decompose(co2)
plot(m)
}
\keyword{ts}