## This is a virtual copy of the decompose() function in the stats package due to David Meyer.
## The stats package should be updated to include this version in R2.14.0, and then I'll delete
## this fix from the forecast package.

decompose <- function (x, type = c("additive", "multiplicative"), filter = NULL) 
{
    type <- match.arg(type)
    l <- length(x)
    f <- frequency(x)
    if (f <= 1 || length(na.omit(x)) < 2 * f) 
        stop("time series has no or less than 2 periods")
    if (is.null(filter)) 
        filter <- if (!f%%2) 
            c(0.5, rep(1, f - 1), 0.5)/f
        else rep(1, f)/f
    trend <- filter(x, filter)
    season <- if (type == "additive") 
        x - trend
    else x/trend
    periods <- l%/%f
    index <- seq(1L, l, by=f) - 1
    figure <- numeric(f)
    for (i in 1L:f) figure[i] <- mean(season[index + i], na.rm=TRUE)
    figure <- if (type == "additive") 
        figure - mean(figure)
    else figure/mean(figure)
    seasonal <- ts(rep(figure, periods+1)[seq_len(l)], start = start(x), frequency = f)
    structure(list(x=x, seasonal = seasonal, trend = trend, random = if (type == 
        "additive") x - seasonal - trend else x/seasonal/trend, 
        figure = figure, type = type), class = "decomposed.ts")
}
