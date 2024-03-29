% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/checkresiduals.R
\name{checkresiduals}
\alias{checkresiduals}
\title{Check that residuals from a time series model look like white noise}
\usage{
checkresiduals(object, lag, test, plot = TRUE, ...)
}
\arguments{
\item{object}{Either a time series model, a forecast object, or a time
series (assumed to be residuals).}

\item{lag}{Number of lags to use in the Ljung-Box or Breusch-Godfrey test.
If missing, it is set to \code{min(10,n/5)} for non-seasonal data, and
\code{min(2m, n/5)} for seasonal data, where \code{n} is the length of the series,
and \code{m} is the seasonal period of the data. It is further constrained to be
at least \code{df+3} where \code{df} is the degrees of freedom of the model. This
ensures there are at least 3 degrees of freedom used in the chi-squared test.}

\item{test}{Test to use for serial correlation. By default, if \code{object}
is of class \code{lm}, then \code{test="BG"}. Otherwise, \code{test="LB"}.
Setting \code{test=FALSE} will prevent the test results being printed.}

\item{plot}{Logical. If \code{TRUE}, will produce the plot.}

\item{...}{Other arguments are passed to \code{\link{ggtsdisplay}}.}
}
\value{
None
}
\description{
If \code{plot=TRUE}, produces a time plot of the residuals, the
corresponding ACF, and a histogram. If \code{test} is not \code{FALSE},
the output from either a Ljung-Box test or Breusch-Godfrey test is printed.
}
\examples{

fit <- ets(WWWusage)
checkresiduals(fit)

}
\seealso{
\code{\link{ggtsdisplay}}, \code{\link[stats]{Box.test}},
\code{\link[lmtest]{bgtest}}
}
\author{
Rob J Hyndman
}
