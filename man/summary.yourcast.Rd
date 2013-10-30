\name{summary.yourcast}
\alias{summary.yourcast}

\title{Summarize yourcast output}

\description{yourcast class summary function}

\usage{
\method{summary}{yourcast}(object, ...)
}

\arguments{
  \item{object}{\code{\link{yourcast}} output object to be summarized}
  \item{...}{Arguments to be passed to or from other methods.}
}

\value{A list of summary objects of class \sQuote{summary.yourcast}:
  
\item{sample.frame}{Four element vector containing, in order,
        the start and end time periods to be used for the observed
        data and the start and end time periods to be forecast}

      \item{params}{Vector. Smoothing parameters used in model.}

      \item{model}{String. Estimation technique used.}

      \item{formula}{Language. Formula used in analysis.}

      \item{numcs}{Integer. Number of cross sections used in analysis.}

      \item{cntry.codes}{Vector. Unique geographical codes of cross
	sections used in analysis.}

      \item{cntry.names}{Vector. If \code{G.names} dataframe submitted
	to \code{\link{yourcast}}, the names of all unique geographical
	codes of cross sections used in analysis.}

      \item{coef}{List. List of matrices with estimated coefficents for
	each cross section. For Lee-Carter models, only one set of
	coefficents calculated. Standard errors are not calculated for
	these coefficients.}
}

\references{\url{http://gking.harvard.edu/yourcast}}

\author{Jon Bischof \email{jbischof@fas.harvard.edu}}

\seealso{\code{\link{yourcast}},\code{\link{print.summary.yourcast}} }

