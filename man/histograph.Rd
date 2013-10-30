\name{histograph}
\alias{histograph}

\title{Histograms for model ebayes}

\description{Draws histograms for priors calculated by model ebayes.}

\usage{
histograph(d1.a, d1.t, dt.da, SD, depvar=" ",
           model="ebayes", graphics.file=NA)
  }

  
\arguments{
  \item{d1.a}{Numeric vector. First derivative respect to age.}
  
  \item{d1.t}{Numeric vector. First derivative respect to time.}
  
  \item{dt.da}{Numeric vector. Second derivative respect to age
    and time.}
  
  \item{SD}{Numeric vector. Standard deviation.}
  
  \item{depvar}{String with the name of the dependent
    variable. Default: \code{" "}}
  
  \item{model}{String with the name of the model. Default: \code{"ebayes"}.}
  
  \item{graphics.file}{String or \code{NA}. If string the name of a file to be 
    appended to a directory path where the graphics will be saved, if \code{NA}
    it is displayed in the screen and not saved.  Default: \code{NA}.}
 }

\value{

  Histograms of the vectors \code{d1.a}, \code{d1.t}, \code{dt.da}, and \code{SD};
  see demos \code{chp.11.7}, \code{chp.11.8}, or \code{chp.11.9}.
}

\references{ \url{http://gking.harvard.edu/yourcast}  }
\author{Federico Girosi \email{girosi@rand.org}}

\keyword{file}% at least one, from doc/KEYWORDS

