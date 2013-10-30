\name{array.yourcast}
\alias{array.yourcast}

\title{Array generation tool for YourCast}

\description{Creates array from YourCast output for each
geographical unit.}

\usage{
array.yourcast(x, unlog=FALSE)
}
\arguments{

\item{x}{A \code{\link{yourcast}} output object.}

\item{unlog}{Logical. Should the dependent variable be unlogged? Default: \code{FALSE}.}
}  %%%end of arguments

\value{Creates array(s) from \code{\link{yourcast}} output. One array is created per geographic area. If there is one geographic area, an array is returned; else, a named list of arrays is returned. Each array is of size \eqn{T} (number of times) by \eqn{A} (number of ages) by 3, where the last dimension captures the type of data: \code{y} (observed values), \code{yhat} (predicted values), and \code{comb} (observed values in-sample and predicted values out-of-sample).}

\author{Konstantin Kashin \email{kkashin@fas.harvard.edu}.}

\references{\url{http://gking.harvard.edu/yourcast}}

\seealso{\code{\link{yourcast}} function and documentation
(\code{help(yourcast)})}

\examples{
# Run Lee-Carter model for Figure 2.6 in Demographic Forecasting
data(chp.2.6.1)
ff.allc <- log(allc2/popu2) ~  time	
ylc.allc <- yourcast(formula=ff.allc, dataobj=chp.2.6.1, model="LC",
                       elim.collinear=FALSE,
                       sample.frame=c(1950,2000,2001,2060))
                       
yc.array <- array.yourcast(ylc.allc)
dimnames(yc.array)

# predicted mortality rates (observed in-sample and predicted out-of-sample)
yc.array[,,"comb"]
} %%% end of examples