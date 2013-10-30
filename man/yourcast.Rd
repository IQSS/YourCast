\name{yourcast}
\alias{yourcast}

\title{Time-series cross-sectional Forecasting}
\description{
  Runs a set of regression models to forecast time-series cross-sectional data by either
  considering independent regressions in each cross-sectional unit or by
  using a variety of techniques to smooth across units.}

\usage{
yourcast(formula=NULL, dataobj=NULL,sample.frame=c(1950,2000,2001,2030), 
                     standardize=TRUE, elim.collinear=FALSE,
                     tol=0.9999, solve.tol = 1.e-10,svdtol=10^(-10),
                     userfile=NULL, savetmp = T, model.frame=FALSE,
                     debug = F,  rerun= "yourcast.savetmp", 
          ### specific to models
                     model="OLS",zero.mean=FALSE,
          #### smooth over ages           
                     Ha.sigma = 0.3,
                     Ha.sigma.sd= 0.1, Ha.deriv=c(0,0,1),
                     Ha.age.weight=0, Ha.time.weight=0,
          #### smooth over time
                     Ht.sigma= 0.3,
                     Ht.sigma.sd=0.1,  Ht.deriv=c(0,0,1),
                     Ht.age.weight=0, Ht.time.weight=0,
          #### smooth over age-time
                     Hat.sigma=0.2,
                     Hat.sigma.sd=0.1,Hat.a.deriv=c(0,1),Hat.t.deriv=c(0,1), 
                     Hat.age.weight=0,Hat.time.weight=0,
          #### smooth over cntry-time
                     Hct.sigma=0.3, Hct.sigma.sd =0.1,
                     Hct.t.deriv=1, Hct.time.weight = 0,
                     LI.sigma.mean=0.2,LI.sigma.sd = 0.1, nsample= 500,
                     low.pow=T, verbose=TRUE)
}

\arguments{
  \item{formula}{ A standard \R formula of the form \eqn{y \sim x_1 +
      x_2}, except that an explanatory variable is included for a
    particular cross-section only if it is both listed in the formula
    and available in that cross-section's data set (see
    \code{dataobj}). Explanatory variables in the formula but not available for
    a cross-section (or in a cross-sectional dataset but not in the
    formula) are excluded. (For mortality forecasting, the specification
    looks like log(deaths/population) \eqn{ \sim x_1 + x_2}, with deaths
    and population stored as separate variables in each dataframe.) (May
    be set to \code{NULL} if \code{savetmp} was set to \code{TRUE} on
    the last run, in which case the value of formula will come from the
    saved file.)}
  
  \item{dataobj}{A object of class \sQuote{yourcast} or equivalent. See
    \code{help(yourprep)} for more details.

    The \code{dataobj} may be supplied in one of four ways. Most
    commonly, the argument will specify (1) an
    object (in working memory) or (2) a string with the name of a file
    in the working directory. However, if (3) \code{dataobj} is a string
    referring to a directory on disk, then
    each element of the list above should be stored in a file in that
    directory, with element \sQuote{data} consisting of a subdirectory
    containing separate ASCII data files. (If this option is chosen, a
    complete data object, called \file{dataobj.Rdata}, will be stored in the
    directory named, and it will be loaded automatically if
    \code{yourcast} is run again with this chosen option.)
    (4) The last option is for dataobj to be set to \code{NULL}, after
    which the function will look for a \file{yourcast.savetmp} file in the
    working directory from a previous run of the function where the
    argument \code{savetmp} was set to \code{TRUE}.

    The function \code{\link{yourprep}} is available to help construct
    the \code{dataobj} in the proper format from individual cross section
    files in the working directory or the workspace. This function also
    performs a number of diagnostics to
    ensure that the data is entered properly and can be read by
    \code{yourcast}. See \code{help(yourprep)} for more information}
  
  \item{sample.frame}{Vector. A four element vector containing, in order, the start
    and end time periods to be used for the observed data and the start
    and end time periods to be forecast. Years identified here that are
    not available for a cross-section are ignored. Default:
    \code{c(1950,2000,2001,2030)}.}
  
  \item{standardize}{Boolean. Should the covariates in each
    cross-sectional unit be standardized (to zero mean and standard
    deviation of  1)? Standardization is performed for both the in-
    and out-of-sample periods. Default: \code{TRUE}.}
  
  \item{elim.collinear}{Boolean. Whether collinearity among covariates
    should be tested and those that are collinear shoul be eliminated.
    Default: \code{FALSE}.}
  
  \item{tol}{Double scalar. Tolerance to find collinearities among
    covariates. Default: \code{0.9999}.}
  
  \item{solve.tol}{A real number smaller than one that is used in the
    argument of the R-function \code{solve} to invert matrices (see
    description for \code{tol}). Default: \eqn{1^{-10}}.}
  
  \item{svdtol}{A scalar; the tolerance used in inverting a matrix by
    SVD. Default: \eqn{10^{-10}}.}
  
  \item{userfile}{A string with the name of a file that contains your
    values for some or all of \code{yourcast}'s arguments. This file
    contains \R code that changes default values of arguments. E.g.,
    the file might contain:
\preformatted{
    index.code <- 30
    data <- "WHOmortalityData"
  }
  If an option is specified in \code{userfile}, it takes precidence over
  command line options, so it is normally best to specify each option in
  either the \code{userfile} or the command line but not both. Default: \code{NULL}}

\item{savetmp}{If \code{TRUE}, \code{yourcast} saves a file in the default directory
  (called \file{yourcast.savetmp}) with preliminary calculations. If the value
  of \code{formula} or \code{dataobj} is missing when \code{yourcast} is
  called, \code{yourcast}
  will get their values from this file, if it exists. This saves a
  minute or so of computing time for large data sets and is useful for
  multiple runs on the same data with different formulas specified or
  different prior values. If \code{FALSE}, no file is saved. (The structure of
  \file{yourcast.savetmp} is for the convenience of \code{yourcast} and is not
  intended to be read by the user or saved for more than one run.)
  Default: \code{TRUE}.}

  \item{model.frame}{If \code{TRUE}, include entire input dataobj in the output
    object. Default: \code{FALSE}.}
  
  \item{debug}{Boolean.  It puts the
    environment that contains parameters and arguments of the
    simulation in the user workspace. Default \code{FALSE}.}
  
  \item{rerun}{String. The name of the file that is saved in the default
    directory with preliminary calculations; see
    \code{savetmp}. Default: \code{yourcast.savetmp}} 
  
  \item{model}{A string indicating the forecasting method, including:
    Bayes maximum a posteriori (\code{map}), Bayes with Gibbs sampling
    (\code{bayes}), Ordinary Least Squares (\code{ols}), Poisson
    (\code{poisson}), and Lee-Carter (\code{LC}). Default: \code{ols}.
    (We usually recommend \code{map}.)

    \code{\link{yourcast}} also includes a procedure to help users set the sigma
    parameters below automatically for the case of model=\code{map}, and
    smoothing over age, time, or age and time, but for only one
    country. You may do this by running a preprocessing instance of
    \code{\link{yourcast}} first by setting this parameter to \code{ebayes} and using
    either the data to be analyzed or a larger data set which is likely
    to have similar or related parameter values. When \code{ebayes} is chosen,
    the \code{\link{yourcast}} output object will contain only the parameter values to
    feed into the next run of \code{\link{yourcast}}.}
  
  \item{zero.mean}{A boolean or named vector with a value of \eqn{ \bar\mu}
    for each age group. If \code{TRUE}, the prior has zero mean. If \code{FALSE}, the
    prior has nonzero mean centered around the observed mean age profile
    (i.e., the average of \eqn{Y} over time and levels of the geographic
    index for each age group). Default: \code{FALSE}.}
  
  \item{Ha.sigma}{This can be set in one of three ways: (1) a scalar
    which sets \eqn{ \sigma_a}, the prior standard deviation of \eqn{ E(Y)},
    indicating how much to smooth \eqn{ E(Y)} over age groups (which may
    vary over geographic areas and time periods, and with the standard
    deviations averaged over age groups). A larger standard deviation
    represents more prior uncertainty, which allows the data to play a
    greater role. (2) \code{NA} to not smooth in this way. (3) To have \code{\link{yourcast}}
    search for a good value based on a target value of the derivative of
    \eqn{ E(Y)} with respect to age, set to a vector of elements containing
    the start and end of a range in sigma in which to look (such as 0.05
    and 1.5), the number of values to look at within this range (such as
    5), and the target value of the derivative of \eqn{ E(Y)} with respect
    to age (such as 0.05). The vector may also include a fifth element,
    which is the target value of the total standard deviation of \eqn{E(Y)}
    over all dimensions of the prior (such as 0.1). (You may choose to
    run \code{\link{yourcast}} with model=\code{ebayes} on a related data set to find an
    approximate target value of the derivative and standard deviation
    automatically.) Default: \code{0.30}.}
  
  \item{Ha.sigma.sd}{A scalar; the standard deviation of parameter
    Ha.sigma (for Gibbs sampling only). Default: \code{0.1}.}
  
  \item{Ha.deriv}{A numeric vector, each element of which is \eqn{n},the
    degree of a (discrete) derivative of the
    smoothness functional with respect to the age group. Element \eqn{k} of
    this vector refers to the \eqn{(k-1)}th derivative, where 0 excludes
    the derviative, 1 includes it, and values in between include the
    derivative but weight it down proportionally. The first element of
    the vector corresponds to the weight on the derivative with respect
    to age of order 0 (the identity operator), the second to the weight
    on the derivative of order 1 (the 1st derivative), etc. For example,
    c(0, 1, 1) corresponds to a mixed functional that penalizes the
    first and second derivatives equally. The higher the order of
    derivative, the more local smoothness over age groups; and lowest
    specified derivative controls the form of prior
    indifference. Default: \code{c(0, 0, 1)}, which usually works well.}

  \item{Ha.age.weight}{A scalar or a numeric vector with weights that
    determine how much smoothing occurs for different age groups. If set
    to 0 or NA, age groups are weighted equally; if set to a nonzero
    scalar, the weight for age group \eqn{a} is set proportional to
    \eqn{a^Ha.age.weight};
    if a vector of length A, the \eqn{ a}th element is the
    weight of age group \eqn{ a}. Default: \code{0}.}
  
  \item{Ha.time.weight}{A scalar or a numeric vector with weights that
    determine how much smoothing occurs for different time periods when
    smoothing over age groups. If \code{0} or \code{NA}, time periods are weighted
    equally; if set to a nonzero scalar value, the weight for time
    period \eqn{ t} in smoothing age groups is proportional to
    \eqn{t^Ha.time.weight}; if the argument is a vector of length T, the 
    \eqn{t}th element is the weight of time period \eqn{t}. Default: \code{0}.}
  
  \item{Ht.sigma}{This can be set in one of three ways: (1) a scalar
    which sets \eqn{ \sigma_t}, the prior standard deviation of \eqn{E(Y)},
    indicating how much to smooth \eqn{E(Y)} over time periods (which may
    vary over geographic areas and age groups, and with the standard
    deviations averaged over time periods). A larger standard deviation
    represents more prior uncertainty, which allows the data to play a
    greater role. (2) NA to not smooth in this way. (3) To have \code{\link{yourcast}}
    search for a good value based on a target value of the derivative of
    \eqn{ E(Y)} with respect to time, set to a vector of elements containing
    the start and end of a range in sigma in which to look (such as 0.05
    and 1.5), the number of values to look at within this range (such as
    5), and the target value of the derivative of \eqn{ E(Y)} with respect
    to time (such as 0.05). The vector may also include a fifth element,
    which is the target value of the total standard deviation of \eqn{ E(Y)}
    over all dimensions of the prior (such as 0.1). (You may choose to
    run \code{\link{yourcast}} with model=\code{ebayes} on a related data set to find an
    approximate target value of the derivative and standard deviation
    automatically.) Default: \code{0.30}.}
  
  \item{Ht.sigma.sd}{A scalar; the standard deviation of parameter
    \code{Ht.sigma} (for Gibbs sampling only). Default: \code{0.1}.}
  
  \item{Ht.deriv}{A numeric vector, each element of which is
    \eqn{n}, the degree of a (discrete) derivative of the
    smoothness functional with respect to time. Element \eqn{k} of this
    vector refers to the \eqn{(k-1)}th derivative, where 0 excludes the
    derviative, 1 includes it, and values in between include the
    derivative but weight it down proportionally. The first element of
    the vector corresponds to the weight on the derivative with respect
    to time of order 0 (the identity operator), the second to the weight
    on the derivative of order 1 (the 1st derivative), etc. For example,
    \code{c(0, 1, 1)} corresponds to a mixed functional that penalizes the
    first and second derivatives equally. The higher the order of
    derivative, the more local smoothness over time; and lowest
    specified derivative controls the form of prior
    indifference. Default: \code{c(0, 0, 1)}, which usually works well.}

  \item{Ht.age.weight}{A scalar or a numeric vector with weights that
    determine how much smoothing occurs for different age groups when
    smoothing over time. If set to \code{0} or \code{NA}, age groups are weighted
    equally in smoothing over time; if set to a nonzero scalar, the
    weight for age group \eqn{a} is set proportional to \eqn{a^Ht.age.weight};
    if a vector of length A, the \eqn{ath} element is the weight of age
    group \eqn{a}. Default: 0.}
  
  \item{Ht.time.weight}{A scalar or a numeric vector with weights that
    determine how much smoothing occurs for different time periods when
    smoothing over time. If \code{0} or \code{NA}, time periods are weighted equally;
    if set to a nonzero scalar value, the weight for time period \eqn{t} in
    smoothing time periods is proportional to \eqn{t^Ht.time.weight}; if
    the argument is a vector of length T, the \eqn{t}th element is the
    weight of time period \eqn{t}. Default: 0.}
  
  \item{Hat.sigma}{This can be set in one of three ways: (1) a
    scalar which sets \eqn{ \sigma_{at}}, the prior standard deviation
    of \eqn{E(Y)}, indicating how much to smooth the time trend in \eqn{E(Y)} over
    age groups. A larger standard deviation represents more prior
    uncertainty, which allows the data to play a greater role. (2) NA to
    not smooth in this way. (3) To have \code{\link{yourcast}} search for a good value
    based on a target value of the derivative of \eqn{E(Y)} with respect to
    age and time, set to a vector of elements containing the start and
    end of a range in sigma in which to look (such as 0.05 and 1.5), the
    number of values to look at within this range (such as 5), and the
    target value of the derivative of \eqn{ E(Y)} with respect to age and
    time (such as 0.05). The vector may also include a fifth element,
    which is the target value of the total standard deviation of \eqn{E(Y)}
    over all dimensions of the prior (such as 0.1). (You may choose to
    run \code{\link{yourcast}} with model=\code{ebayes} on a related data set to find an
    approximate target value of the derivative and standard deviation
    automatically.) Default: \code{0.2}. }
  
  \item{Hat.sigma.sd}{A scalar; the standard deviation of parameter
    \code{Hat.sigma} (for Gibbs sampling only). Default: \code{0.1}.}
  
  \item{Hat.a.deriv}{A numeric vector, each element of which is \eqn{n}, the degree of a (discrete) derivative of the
    smoothness functional of time trends with respect to age
    groups. Element \eqn{k}  of this vector refers to the \eqn{(k-1)}th
    derivative of the time trend v with respect to age, where 0 excludes
    the derviative, 1 includes it, and values in between include the
    derivative but weight it down proportionally. The first element of
    the vector corresponds to the weight on the derivative of the time
    trend with respect to age of order 0 (the identity operator), the
    second to the weight on the derivative of order 1 (the 1st
    derivative), etc. For example, \code{c(0, 1, 1)} corresponds to a mixed
    functional that penalizes the first and second derivatives
    equally. The higher the order of derivative, the more local
    smoothness over time; and lowest specified derivative controls the
    form of prior indifference. Default: \code{c(0, 0, 1)}, which usually works
    well.}
  
  \item{Hat.t.deriv}{A numeric vector, each element of which is \eqn{n}, the degree of a (discrete) derivative of the
    smoothness functional of age derivative with respect to
    time. Element \eqn{k} of this vector refers to the \eqn{(k-1)}th
    derivative of the age derivative with respect to time, where 0
    excludes the derviative, 1 includes it, and values in between
    include the derivative but weight it down proportionally. The first
    element of the vector corresponds to the weight on the age
    derivative with respect to time of order 0 (the identity operator),
    the second to the weight on the derivative of order 1 (the 1st
    derivative), etc. For example, \code{c(0, 1, 1)} corresponds to a mixed
    functional that penalizes the first and second derivatives
    equally. The higher the order of derivative, the more local
    smoothness over time; and lowest specified derivative controls the
    form of prior indifference. Default: \code{c(0, 0, 1)}, which usually works
    well.}
  
  \item{Hat.age.weight}{A scalar or a numeric vector with weights that
    determines how much smoothing occurs for different age groups when
    smoothing over age and time. If set to \code{0} or \code{NA}, age groups are
    weighted equally in smoothing over time; if set to a nonzero scalar,
    the weight for age group \eqn{a} is set proportional to
    \eqn{a^Ht.age.weight}; if a vector of length A, the \eqn{a}th element is the
    weight of age group \eqn{a}. Default: \code{0}.}
  
  \item{Hat.time.weight}{A scalar or a numeric vector with weights that
    determine how much smoothing occurs for different time periods when
    smoothing over age and time. If \code{0} or \code{NA}, time periods are weighted
    equally; if set to a nonzero scalar value, the weight for time
    period \eqn{t} in smoothing time periods is proportional to 
    \eqn{t^Ht.time.weight}; if the argument is a vector of length T, the 
    \eqn{t}th element is the weight of time period \eqn{t}. Default: \code{0}.}
  
  \item{Hct.sigma}{A scalar which sets \eqn{\sigma_t}, the prior standard
    deviation of \eqn{E(Y)}, which indicates how to smooth  \eqn{E(Y)} over
    geographic areas, or NA to not smooth in this way. The parameter 
    \eqn{\sigma_ct} is the expected prior standard deviation of \eqn{E(Y)} for a
    geographic area (varying over time periods and age groups, and with
    the standard deviations averaged over geographic areas). (A larger
    standard deviation represents more prior uncertainty, which allows
    the data to play a greater role.) Default: \code{0.3}.}
  
  \item{Hct.sigma.sd}{A scalar; the standard deviation of parameter
    Ht.sigma (for Gibbs sampling only). Default: \code{0.1}.}
  
  \item{Hct.t.deriv}{A numeric vector; controls whether smoothing the
    level or the time trend of \eqn{E(Y)} over geographic areas (both
    cannot presently be done simultaneously). To smooth the level of 
    \eqn{E(Y)} over geographic areas, set to 1, the identity. To smooth the
    time trend, set this (as in \code{Hat.t.deriv}) to the weight of the
    partial derivative taken with respect to time in the standard
    smoothness functional for the prior. The use of the first or higher
    order partial derivatives are supported. Default: \code{1}.}
  
  \item{Hct.time.weight}{A scalar or a numeric vector with weights
    that determine how much smoothing occurs for different time periods
    when smoothing over geographic areas. If \code{0} or \code{NA}, time periods are
    weighted equally; if set to a nonzero scalar value, the weight for
    time period \eqn{t} in smoothing over areas is proportional to 
    \eqn{t^Hct.time.weight}; if the argument is a vector of length T, the 
    \eqn{t}th element is the weight of time period \eqn{t}. Default: \code{0}. }
  
  
  \item{LI.sigma.mean}{A scalar; used in the likelihood and in the
    calculation of the priors in conjunction with \code{Ha.sigma.sd},
    \code{Hat.sigma.sd}, \code{Ht.sigma.sd}, and \code{Hct.sigma.sd}.
    Default: \code{0.2}.}
  
  \item{LI.sigma.sd}{A scalar; the standard deviation of
    \code{LI.sigma.mean} used in the calculation of the priors. Default: \code{0.1}.}
  
  \item{nsample}{A scalar; represents the number of iterations in the
    Gibbs algorithm \code{bayes}. Default: \code{500}.}
  
  \item{low.pow}{Boolean. Whether to include lower-power of explanatory
    variables in the simulation as derived from \code{formula}. For
    example \eqn{ y \sim x^4}, if \code{low.pow} = \code{TRUE}, then
    \eqn{x, x^2, x^3, x^4} will be included. Default: \code{TRUE}.}

  \item{verbose}{Boolean. Suppress verbose output. Default: \code{FALSE}}
}

\value{
  Returns a list of class \sQuote{yourcast} containing the following components: 

  \item{call}{The full call, including all command line options when
    yourcast was called.}
  
  \item{userfile}{The full userfile if it was specified.}
  
  \item{yhat}{A list with the same cross-sectional elements as the input
    data, but with two columns: \sQuote{y} for the observed dependent
    variable and \sQuote{yhat} for the predicted values. These include both
    in-sample and out-of-sample values, as distinguished by the values of
    \code{sample.frame}.}
  
  \item{coeff}{A list with the same cross-sectional elements as the
    input data, elements of which are the estimated coefficients if
    calculated by the chosen model.}
  
    \item{sigma}{A list with the same cross-sectional elements as the
    input data, elements of which are the estimated standard error of
    the estimate of the regression (the standard deviation of the
    dependent variable given the explanatory variables).}

  \item{aux}{List. A list of summary information about the yourcast analysis
    used by \code{\link{plot.yourcast}}}

  \item{params}{Vector. Smoothing parameters used in model.}
  
}
\references{\url{http://gking.harvard.edu/yourcast} }

\author{Federico Girosi \email{girosi@rand.org}; 
  Elena Villalon \email{evillalon@iq.harvard.edu};
  Gary King \email{king@harvard.edu}}

\keyword{file}% at least one, from doc/KEYWORDS
\keyword{datasets}% __ONLY ONE__ keyword per line
