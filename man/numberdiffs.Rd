\name{ndiffs}
\alias{ndiffs}
\title{Number of differences}
\usage{ndiffs(x, alpha=0.05, test=c("kpss","adf", "pp"))
}
\arguments{
\item{x}{A univariate time series}
\item{alpha}{Level of the test}
\item{test}{Type of unit root test to use}
}

\description{Uses a unit root test to determine the number of differences required for time series \code{x}. If \code{test="kpss"},
the KPSS test is used with the null hypothesis that \code{x} has a stationary root against a unit-root alternative. Then the test returns 
the least number of differences required to pass the test at the level \code{alpha}. If \code{test="adf"}, the Augmented Dickey-Fuller test is used 
and if \code{test="pp"} the Phillips-Perron test is used. In both of these cases, the null hypothesis is that \code{x} has a unit root against a 
stationary root alternative. Then the test returns the least number of differences required to fail the test at the level \code{alpha}.}

\value{An integer.}

\author{Rob J Hyndman}

\examples{ndiffs(WWWusage)}

\keyword{ts}
