\name{yourprep}
\alias{yourprep}

\title{Data object creation wizard for YourCast}

\description{Builds the data object for \code{\link{yourcast}}
     function from files in working directory or other specified
     directory and checks for errors}

   \usage{yourprep(dpath=getwd(),tag="csid",index.code="ggggaa",
                     datalist=NULL,G.names=NULL,A.names=NULL,
                     T.names=NULL,proximity=NULL,year.var=FALSE,
                     sample.frame=NULL,summary=FALSE,verbose=FALSE,

                     #lagging utility
                     lag=NULL,formula=NULL,vars.nolag=NULL)}


\arguments{

\item{dpath}{String. Name of the directory where data files are
          stored.  If \code{NULL} then defaults to working directory.
          Default: \code{NULL}}

\item{tag}{String. Group of characters placed before CSID code in
          filenames to indicate which files in \code{dpath} function
          should load. The \code{tag} can also be used to differentiate
	  between different groups to be considered in separate
	  analysis; for example, \sQuote{m} for male deaths and
	  \sQuote{f} for female deaths. Default: \code{"csid"}}

\item{index.code}{String indicating how the CSID index variable is
          coded in the input data. Between 0 and 4 of the following two
	  characters are used in this order: \code{g} for the geographic
	  index (such
          as country) and \code{a} for a grouped continuous variable
          like
          an age group. For example, \code{ggggaa} would have the
          function interpret
          \sQuote{245045} by using \sQuote{2450} as the country code and
          \sQuote{45} as the age group. Default: \code{"ggggaa"}}



\item{datalist}{A list of cross section dataframes already loaded into
          the workspace to be added to the \code{dataobj}. Names of
          list elements should be the numerical CSID code for each
          cross section, and dataframes should be formated identically
          to files loaded from an external directory (see Details)}

\item{A.names, G.names, T.names}{String. Filename of optional
          two-column data files that list all valid numerical codes
          (in the first column) and corresponding alphanumeric names
          (optionally in the second column) for the indices
          corresponding to geographic areas in \code{G.names}, age
          groups in \code{A.names}, and time periods in
          \code{T.names}. Function will search \code{dpath} for file
          with specified name; please include column labels. The
          optional alphanumeric identifiers are most commonly only
          used for geographic areas since numerical values for age
          groups and time periods are usually meaingful on their
          own. However, if other grouped continuous variable used in
          place of ages, for example, specifying these labels will be
          important for output to be meaningful. NOTE: Auxiliary files
	  will loaded automatically by \code{yourprep()} if they are
	  saved in the \code{dpath} and labeled with the \code{tag}
	  specified by
	  the user. See \sQuote{Details} section for more
	  infromation. Default: \code{NULL}}


\item{proximity}{Data file with codes to construct the symmetric
          matrix (geographic region by geographic region) of proximity
          scores for geographic smoothing used by the \sQuote{map} and \sQuote{bayes}
          methods. The larger the relative score, the more proximate
          that pair of countries is in the prior; a zero element means
          the two geographic areas are unrelated (the diagonal is
          ignored).  Each row of the \code{proximity} file has three
          columns, consisting of geographic codes for two countries
          and a score indicating the proximity or similarity of the
          two geographic regions; please include column labels. For
          convenience, geographic regions that are unrelated (and
          would have zero entries in the symmetric matrix) may be
          omitted from \code{proximity}. In addition, \code{proximity}
          may include rows corresponding to geographic regions not
          included in the present analysis. Default: \code{NULL}}

\item{year.var}{Boolean. Should be \code{TRUE} if \code{year} coded as
          separate variable rather than as rowname for cross section
          data files.  Function will look for \code{year} variable to
          use as rownames and then drop it from the dataframe. Change
          will only be made to dataframe if it does not already have
          rownames or if exisiting rownames are merely a
          \sQuote{1...N} index of row numbers, so it is possible to
          apply correction even if some cross sections do not have a
          \code{year} variable and already have the correct
          rownames. Default: \code{FALSE}}

\item{sample.frame}{Optional four element vector containing, in order,
          the start and end time periods to be used for the observed
          data and the start and end time periods to be forecast. All
          cross sections do not have to begin at starting date, but
          must contain all years after the first observed
          value. Variables to be forecasted should be coded as
          \code{NA} in the out-of-sample period. Note that this makes
          it easy to reserve a range of values of the dependent
          variable for out-of-sample forecasting evaluation; our
          \code{summary} and \code{plot} functions in
          \code{\link{yourcast}} will make these comparisons
          automatically if the out-of-sample data are
          included. \code{yourprep()} uses this information only to
          verify that cross sections are correctly
          constructed, but it should also be included if one wants to
	  use the lag utility. Default: \code{NULL}}

\item{summary}{Boolean. If \code{TRUE}, means for available
          observations on each variable are displayed for the cross
          sections read by \code{yourprep()}. Default: \code{FALSE}}

\item{verbose}{Boolean. If \code{TRUE}, function prints name of each
          cross section or auxiliary file as it is read into the
          \code{dataobj}. Default: \code{FALSE}}

\item{lag}{Number of years covariate data needs to be lagged from
  current position is cross section files. See \sQuote{Details} for more
  information. Default: \code{NULL}}

\item{formula}{Formula. The formula that one will use in the subsequent run of
  \code{yourcast()}. This helps the lagging utility distinguish between
  the response variable (which will not be shifted between cross
  sections) and the covariates of interest that should be lagged and
  included in the final cross sections of the dataobj. If the covariate
  \sQuote{index} is included in the formula, the lagging utility will
  include a variable in the cross sections that starts from 1 and counts
  the number of time periods since the start of the cross section.
  If a lag is requested, the formula argument must be specified.
  Default: \code{NULL}}

\item{vars.nolag}{Vector of strings. Vector of variables to be included
  in the dataobj but not lagged. These variables do not need to be
  included in the formula, and if found there will not ignored when the
  other covariates are lagged.}

      }

\details{Creates \code{dataobj} input for \code{\link{yourcast}} from
     files in working directory or other specified directory. Checks
     that all cross sections in \code{data} list titled properly and
     if all years up to last predicted year included in the dataframes
     (if \code{sample.frame} argument specified). Please note, however,
     that all cross sections from the same geographic area must have the
     same observation and prediction years in the dataframe (even if
     \code{NA}) for the graphing software \code{\link{plot.yourcast}} to
     work.

     The cross section files must be named according to the CSID
     identifiers for country code and age group, preceeded by the
     specified tag (default: \code{"csid"}) so that \code{yourprep()} can
     identify the file from other files in the dpath. For example, for
     the USA (country code 2450) time series of 45 year old
     individuals, the file name should be \file{csid245045.txt} if the
     tag is left as the default. Files must have an extension so that
     the program can recognize how the data is coded. Currently, fixed
     width text files (\file{*.txt}), comma-separated values
     (\file{*.csv}), and Stata v.5-10 (\file{*.dta}) files are
     supported, and multiple file types may be used in the same run of
     the program. \file{*.Rdata} objects can be included with the
     \code{datalist} option after they are loaded to a list in the
     workspace. \code{yourprep()} includes diagnostics to ensure that
     objects are properly named and not included accidentally, but
     users should examine the specified \code{dpath} before running
     \code{yourprep()} to minimize errors.

     Each cross section file should be labeled columns of time-series
     data for the dependent variable(s) (e.g., disease, pop) and the
     covariates that will be used in the forecast. The rownames for
     the dataframe should be the observation year (if the year is
     coded as a separate variable, set \code{year.var=TRUE}). The
     files must contain the full time series that will be specified in
     the \code{sample.frame} argument in \code{\link{yourcast}} after
     the first observed year. For instance, if
     \code{sample.frame=c(1950,2000,2001,2030)}, then files would have
     observations that start between 1950 and 2000 and include all
     other years (even if the entries are \code{NA}) up to the last
     year of prediction, i.e., 2030.

     Optional auxiliary files such as \code{G.names} should be named
     according to the filename specified in the respective
     arguments. If specified, these files must have extensions and be
     coded in one of the three supported file types. However, these
     files will be automatically loaded by \code{yourprep()} if they are
     saved in the \code{dpath} and labeled with the tag specified by the
     user. The default names for these files must be used (e.g.,
     \sQuote{G.names} and \sQuote{proximity}). For example,
     if the \code{tag} is left as the default and there is a file in the
     \code{dpath} labeled \file{csid.G.names.txt}, \code{yourprep()}
     will load this
     automatically and save the input as the \code{G.names} element of
     the \sQuote{dataobj} list. \code{yourprep()} arguments such as
     \code{G.names} take precedence over \file{TAG.*} files in
     the\code{dpath}.

     \code{yourprep()} also includes a lagging utility (activated once
     one specifies a lag length with the \sQuote{lag} argument). This
     utility is useful for when the data in each cross section is, for
     example, the response and covariates for 50 year olds in each year
     but the desired content for each cross section is the response for
     50 year olds and the covariates for 25 year olds 25 years prior to
     each year (implying a lag of 25 years). In order to have
     \code{yourprep()} perform this lagging automatically, include cross
     sections for each age group with data starting the same number of
     years before the first observation year as the requested lag
     period. Thus if \code{lag=25} and the first observation year is
     1950, then the cross sections should all start at 1925. Age groups
     younger than the length of the lag will not retain covariate data
     (except perhaps an \sQuote{index} variable) in the output
     object. The covariates lagged are the predictor variables specified
     in the formula argument.

     If data for a cohort 25 years (in this case) younger is not
     available for some cohort over age 25, \code{yourprep()} will look
     for the closest cohort available and issue a warning message. 
   }



\value{
\item{dataobj}{A list with several components:

\describe{ 

\item{data}{A list with the cross-sectional data matrices as elements.}

\item{proximity}{A three-column matrix of proximity scores for geographic smoothing used by the \sQuote{map}
and \sQuote{bayes} methods. For each row, the first two columns indicate the country pair. The third column indicates the proximity score. The larger each score, the more
proximate that pair of countries is in the prior; a zero element means
the two geographic areas are unrelated (the diagonal is ignored).}

\item{G.names, A.names, T.names}{Optional two-column dataframes that
list all valid numerical codes (in the first column, labeled codes)
and corresponding alphanumeric names (optionally in the second
column, labeled name) for the indices corresponding to the geographic
areas in \code{G.names}, age groups in \code{A.names}, and time
periods in \code{T.names}.}

\item{index.code}{A string indicating how the index variable is coded
in the input data.}}}}

\author{Jon Bischof \email{jbischof@fas.harvard.edu}}

\references{\url{http://gking.harvard.edu/yourcast}}

\seealso{\code{\link{yourcast}} function and documentation
(\code{help(yourcast)})}

\examples{
\dontrun{
# Working directory automatically set to directory with cross
# section and auxiliary files to begin. Files for this example
# in 'data' folder of YourCast library.

#Old working directory to be restored later
oldwd <- getwd()
# Now setting wd to 'data' folder in YourCast library
setwd(system.file("data",package="YourCast"))

# Simple run of the function, using option that turns year variable
# into label in each cs. Use sample.frame argument for all diagnostics
# to work
 
dta <- yourprep(G.names="cntry.codes.txt", proximity="proximity.txt",
year.var=TRUE,verbose=TRUE,sample.frame=c(1950,2000,2001,2030))


# With summary output (means of variables in each cross section) 

%\dontrun{
dta <- yourprep(G.names="cntry.codes.txt", proximity ="proximity.txt",
year.var=TRUE,summary=TRUE)
%}

# Function can also add datafiles already loaded into R as objects in
# the workspace with "datalist" option if put into a list and properly
# labeled. All diagnostics still performed 
# 'csid204545', etc., are dataframes in workspace

# Labels changed to nonsense ones so as not to confuse with other files

data(csid204545)
data(csid204550)
data(csid204555)

datalist <- list("123456"=csid204545,"234567"=csid204550,
"345678"=csid204555) 

# Verbose option turned on and datalist argument added 

dta <- yourprep(G.names="cntry.codes.txt", proximity="proximity.txt",
year.var=TRUE,verbose=TRUE,datalist=datalist)

# Setting working directory back
setwd(oldwd)
rm(oldwd)
}
}
