% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggplot.R
\docType{data}
\name{StatForecast}
\alias{StatForecast}
\alias{GeomForecast}
\alias{geom_forecast}
\title{Forecast plot}
\format{
An object of class \code{StatForecast} (inherits from \code{Stat}, \code{ggproto}, \code{gg}) of length 3.

An object of class \code{GeomForecast} (inherits from \code{Geom}, \code{ggproto}, \code{gg}) of length 7.
}
\usage{
StatForecast

GeomForecast

geom_forecast(
  mapping = NULL,
  data = NULL,
  stat = "forecast",
  position = "identity",
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE,
  PI = TRUE,
  showgap = TRUE,
  series = NULL,
  ...
)
}
\arguments{
\item{mapping}{Set of aesthetic mappings created by \code{\link[ggplot2:aes]{aes()}}. If specified and
\code{inherit.aes = TRUE} (the default), it is combined with the default mapping
at the top level of the plot. You must supply \code{mapping} if there is no plot
mapping.}

\item{data}{The data to be displayed in this layer. There are three options:

If \code{NULL}, the default, the data is inherited from the plot data as
specified in the call to \code{\link[ggplot2]{ggplot}}.

A \code{data.frame}, or other object, will override the plot data. All
objects will be fortified to produce a data frame. See \code{\link[ggplot2]{fortify}}
for which variables will be created.

A \code{function} will be called with a single argument, the plot data. The
return value must be a \code{data.frame}, and will be used as the layer
data.}

\item{stat}{The stat object to use calculate the data.}

\item{position}{Position adjustment, either as a string, or the result of a
call to a position adjustment function.}

\item{na.rm}{If \code{FALSE} (the default), removes missing values with a
warning.  If \code{TRUE} silently removes missing values.}

\item{show.legend}{logical. Should this layer be included in the legends?
\code{NA}, the default, includes if any aesthetics are mapped. \code{FALSE}
never includes, and \code{TRUE} always includes.}

\item{inherit.aes}{If \code{FALSE}, overrides the default aesthetics, rather
than combining with them. This is most useful for helper functions that
define both data and aesthetics and shouldn't inherit behaviour from the
default plot specification, e.g. \code{\link[ggplot2]{borders}}.}

\item{PI}{If \code{FALSE}, confidence intervals will not be plotted, giving
only the forecast line.}

\item{showgap}{If \code{showgap=FALSE}, the gap between the historical
observations and the forecasts is removed.}

\item{series}{Matches an unidentified forecast layer with a coloured object
on the plot.}

\item{...}{Additional arguments for \code{\link{forecast.ts}}, other
arguments are passed on to \code{\link[ggplot2]{layer}}. These are often aesthetics,
used to set an aesthetic to a fixed value, like \code{color = "red"} or
\code{alpha = .5}. They may also be parameters to the paired geom/stat.}
}
\value{
A layer for a ggplot graph.
}
\description{
Generates forecasts from \code{forecast.ts} and adds them to the plot.
Forecasts can be modified via sending forecast specific arguments above.
}
\details{
Multivariate forecasting is supported by having each time series on a
different group.

You can also pass \code{geom_forecast} a \code{forecast} object to add it to
the plot.

The aesthetics required for the forecasting to work includes forecast
observations on the y axis, and the \code{time} of the observations on the x
axis. Refer to the examples below. To automatically set up aesthetics, use
\code{autoplot}.
}
\examples{

\dontrun{
library(ggplot2)
autoplot(USAccDeaths) + geom_forecast()

lungDeaths <- cbind(mdeaths, fdeaths)
autoplot(lungDeaths) + geom_forecast()

# Using fortify.ts
p <- ggplot(aes(x=x, y=y), data=USAccDeaths)
p <- p + geom_line()
p + geom_forecast()

# Without fortify.ts
data <- data.frame(USAccDeaths=as.numeric(USAccDeaths), time=as.numeric(time(USAccDeaths)))
p <- ggplot(aes(x=time, y=USAccDeaths), data=data)
p <- p + geom_line()
p + geom_forecast()

p + geom_forecast(h=60)
p <- ggplot(aes(x=time, y=USAccDeaths), data=data)
p + geom_forecast(level=c(70,98))
p + geom_forecast(level=c(70,98),colour="lightblue")

#Add forecasts to multivariate series with colour groups
lungDeaths <- cbind(mdeaths, fdeaths)
autoplot(lungDeaths) + geom_forecast(forecast(mdeaths), series="mdeaths")
}

}
\seealso{
\code{\link[generics]{forecast}}, \code{\link[ggplot2]{ggproto}}
}
\author{
Mitchell O'Hara-Wild
}
\keyword{datasets}
