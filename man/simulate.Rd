\name{simulate.ets}
\alias{simulate.ets}
\title{Simulation from an ETS model}
\usage{
\method{simulate}{ets}(object, nsim=length(object$x), seed=NULL, initstate=object$state[1,], bootstrap=TRUE, ...)
}

\arguments{
\item{object}{An object of class "\code{ets}". Usually the result of a call to \code{\link{ets}}.}
\item{nsim}{Number of periods for the simulated series}
\item{seed}{Either NULL or an integer that will be used in a call to \code{\link{set.seed}} before simulating the time seriers.
    The default, NULL will not change the random generator state.}
\item{initstate}{State at time 0.}
\item{bootstrap}{If TRUE, simulation uses resampled errors rather than normally distributed errors.}
\item{...}{Other arguments.}
}

\description{Returns a time series based on the model object \code{object}.}

\value{An object of class "\code{ts}".}

\seealso{\code{\link{ets}}}

\author{Rob J Hyndman}

\examples{fit <- ets(USAccDeaths)
plot(simulate(fit, 100))
}
\keyword{ts}
