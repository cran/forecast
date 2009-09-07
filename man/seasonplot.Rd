\name{seasonplot}
\alias{seasonplot}
\title{Seasonal plot}
\usage{seasonplot(x, s, season.labels = NULL, year.labels = FALSE,
                 year.labels.left = FALSE, type = "o", main, ylab = "",
                 xlab = NULL, col = 1, ...)
}
\arguments{
\item{x}{a numeric vector or time series.}
\item{s}{seasonal frequency of x}
\item{season.labels}{Labels for each season in the "year"}
\item{year.labels}{Logical flag indicating whether labels for each year of data should be plotted on the right.}
\item{year.labels.left}{Logical flag indicating whether labels for each year of data should be plotted on the left.}
\item{type}{plot type (as for \code{\link[graphics]{plot}})}
\item{main}{Main title.}
\item{ylab}{Y-axis label}
\item{xlab}{X-axis label}
\item{col}{Colour}
\item{\dots}{additional arguments to \code{\link[graphics]{plot}}.} }
\description{Plots a seasonal plot as described in Makridakis, Wheelwright and Hyndman (1998, chapter 2).
}

\value{  None.
}
\references{  Makridakis, Wheelwright and Hyndman (1998) \emph{Forecasting:
methods and applications}, Wiley: New York. \url{http://www.robhyndman.info/forecasting/}
}
\author{Rob J Hyndman}
\seealso{\code{\link[stats]{monthplot}}}
\examples{seasonplot(AirPassengers)
}
\keyword{ts}
