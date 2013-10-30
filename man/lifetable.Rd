\name{lifetable}
\alias{lifetable}

\title{Life table calculation for YourCast}

\description{Creates a list of observed and predicted period life tables across from YourCast output for each geographical unit. This function should only be used when forecasting all-cause mortality rates.}

\usage{
lifetable(obj,ax=0.5,a0=0.07,nx=NULL,dv.log=NA)
}

\arguments{

\item{obj}{A \code{\link{yourcast}} output object containing all-cause mortality rates.}

\item{ax}{A scalar, numeric vector of length \eqn{A}, or a numeric \eqn{T} by \eqn{A} matrix specifying the average years lived within a given age interval \eqn{[x,x+n)} by individuals who die within that age interval. If a scalar, \code{ax} is assumed to be constant across all age intervals (except for the first age interval; see \code{a0}) and across time. If a vector, \code{ax} is assumed to be constant across time. See \sQuote{Details}. Default: 0.5.}

\item{a0}{A scalar indicating the average years lived within the first age interval by individuals who die within that age interval. Ignored if \code{ax} is a vector or a matrix. Default: 0.07.}

\item{nx}{A scalar or vector of length \eqn{A} defining the age interval widths. If a scalar, the age intervals have identical widths. If \code{NULL}, the age interval widths are calculated using the ages in the \code{\link{yourcast}} object. The last age interval is assumed to be open-ended so its length is set to \code{NA}. Default: \code{NULL}.}

\item{dv.log}{Logical. Is the dependent variable logged? If \code{NA}, the function parses the formula in the \code{\link{yourcast}} object to determine if the dependent variable was logged. Default: \code{NA}.}

}  %%%end of arguments

\details{The input must be a \code{\link{yourcast}} object containing forecasts of all-cause mortality rates. The death rates do not have to be for one-year age intervals, albeit \code{ax} and \code{a0} should be adjusted accordingly based on the age interval size. The rates may either be logged or unlogged. If \code{dv.log} is set to \code{NULL}, the function parses the formula automatically to determine whether the rates are logged. Setting \code{dv.log} to \code{TRUE} or \code{FALSE} manually sets whether the rates are logged or not.

Specifying \code{ax} is necessary to convert mortality rates into conditional probabilities of death. Suppose that \eqn{m_x} is the mortality rate for the \eqn{n}-year age interval \eqn{[x,x+n)}. \eqn{a_x} is then defined as the average years lived within the age interval \eqn{[x,x+n)} by individuals who die at that age. The default value of \eqn{a_x} is 0.5, which means that individuals on average live half the year within the age that they die (for single-year age intervals). The default value of \eqn{a_0} is set to 0.07, following the reported values within the Human Mortality Database for adult females in the United States for the \eqn{[0,1)} age interval in the last decade.

Alternatively, one can specify \code{ax} as a vector of length \eqn{A}, where \eqn{A} is the number of age intervals. In this case, \eqn{a_x} is assumed to be constant across time for each age interval. One can also pass a \eqn{T} by \eqn{A} matrix, in which case \eqn{a_x} is allowed to vary over time.

The last age interval is open-ended: \eqn{[x,\infty)}. The conditional probability of death is set to 1 and \eqn{a_x} is automatically set to \eqn{1/m_x}.
} %%% end of details

\value{If a single geographic area, a list with names \code{obs} and \code{pred} is returned, where the elements are life table arrays calculated using observed and predicted mortality rates, respectively. If multiple geographic areas, a named list of geographic areas is returned, where each element contains a list of observed and predicted life tables.

Each life table array has dimensions \eqn{T} by \eqn{A} by \eqn{10}, where the third dimension contains the following variables:
\item{x}{Age interval.}
\item{nx}{Length of age interval.}
\item{ax}{Average years lived within age interval by individuals who die during that age interval.}
\item{qx}{Conditional probability of death for the age interval.}
\item{px}{Conditional probability of survival for the age interval.}
\item{lx}{Number of survivors at the start of the age interval. Set to 100,000 for the first age interval.}
\item{dx}{Number of individuals dying within the age interval.}
\item{Lx}{Person-years lived by individuals in the age interval.}
\item{Tx}{Person-years remaining at start of age interval.}
\item{ex}{Remaining life expectancy at start of age interval.}
	}

\author{Konstantin Kashin \email{kkashin@fas.harvard.edu}}

\references{
Preston, S.H., Patrick Heuveline, Michel Guillot. (2000) \emph{Demography: Measuring and Modeling Population Processes.} Blackwell Publishing.

Wilmoth, J.R. et al. \dQuote{Methods Protocol for the Human Mortality Database} \url{http://www.mortality.org/}.

\url{http://gking.harvard.edu/yourcast}
}

\seealso{\code{\link{yourcast}} function and documentation
(\code{help(yourcast)})}

\examples{
# Run Lee-Carter model for Figure 2.6 in Demographic Forecasting
data(chp.2.6.1)
ff.allc <- log(allc2/popu2) ~  time	
ylc.allc <- yourcast(formula=ff.allc, dataobj=chp.2.6.1, model="LC",
                       elim.collinear=FALSE,
                       sample.frame=c(1950,2000,2001,2060))

# calculate lifetable function assuming ax=2.5 and a0=1.907
# a0 calculated for both males and females from 
# Human Mortality Database for New Zealand for 2000
lt <- lifetable(ylc.allc,ax=2.5,a0=1.907)

# observed lifetable for csid 5150 (New Zealand) for 2000
lt[["5150"]]$obs["2000",,]

# predicted lifetable for csid 5150 (New Zealand) for 2015
lt[["5150"]]$pred["2015",,]
}
