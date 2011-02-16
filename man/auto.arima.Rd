\name{auto.arima}
\alias{best.arima}
\alias{auto.arima}
\title{Fit best ARIMA model to univariate time series}
\usage{
auto.arima(x, d = NA, D = NA, max.p = 5, max.q = 5,
           max.P = 2, max.Q = 2, max.order = 5,
           start.p=2, start.q=2, start.P=1, start.Q=1,
           stationary = FALSE, ic = c("aic","aicc", "bic"), 
           stepwise=TRUE, trace=FALSE, 
           approximation=length(x)>100 | frequency(x)>12, xreg=NULL,
           test=c("kpss","adf","pp"), allowdrift=TRUE)
}

\arguments{
\item{x}{a univariate time series}
\item{d}{Order of first-differencing. If missing, will choose a value based on KPSS test.}
\item{D}{Order of seasonal-differencing. If missing, will choose a value based on CH test.}
\item{max.p}{Maximum value of p}
\item{max.q}{Maximum value of q}
\item{max.P}{Maximum value of P}
\item{max.Q}{Maximum value of Q}
\item{max.order}{Maximum value of p+q+P+Q if model selection is not stepwise.}
\item{start.p}{Starting value of p in stepwise procedure.}
\item{start.q}{Starting value of q in stepwise procedure.}
\item{start.P}{Starting value of P in stepwise procedure.}
\item{start.Q}{Starting value of Q in stepwise procedure.}
\item{stationary}{If \code{TRUE}, restricts search to stationary models.}
\item{ic}{Information criterion to be used in model selection.}
\item{stepwise}{If \code{TRUE}, will do stepwise selection (faster). Otherwise, 
                it searches over all models. Non-stepwise selection can be very slow,
                especially for seasonal models.}
\item{trace}{If \code{TRUE}, the list of ARIMA models considered will be reported.}
\item{approximation}{If \code{TRUE}, estimation is via conditional sums of squares and
    the information criteria used for model selection are approximated. The final model is 
    still computed using maximum likelihood estimation. Approximation should be used
    for long time series or a high seasonal period to avoid excessive computation
    times.}
\item{xreg}{Optionally, a vector or matrix of external regressors, which must have the same number of rows as x.}
\item{test}{Type of unit root test to use. See \code{\link{ndiffs}} for details.}
\item{allowdrift}{If \code{TRUE}, models with drift terms are considered.}
}
\description{Returns best ARIMA model according to either AIC, AICc or BIC value.
The function conducts a search over
possible model within the order constraints provided.}

\details{Non-stepwise selection can be slow, especially for seasonal data. Non-seasonal differences chosen using the KPSS test. 
Seasonal differences chosen using a variation on the Canova-Hansen test. Stepwise algorithm outlined in Hyndman and Khandakar (2008).
}
\value{Same as for \code{\link{arima}}}

\references{Hyndman, R.J. and Khandakar, Y. (2008) "Automatic time series forecasting: The forecast
package for R", \emph{Journal of Statistical Software}, \bold{26}(3).}


\seealso{\code{\link{Arima}}}

\author{Rob J Hyndman}
\examples{fit <- auto.arima(WWWusage)
plot(forecast(fit,h=20))
}
\keyword{ts}
