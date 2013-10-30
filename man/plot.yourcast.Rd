\name{plot.yourcast}
\alias{plot.yourcast}

\title{Plot generation tool for YourCast}

\description{Creates graphics from YourCast output for each
geographical unit and prints to the device window or a .pdf file in
the specified directory.}

\usage{
\method{plot}{yourcast}(x, plots=c("age","time"), title=NULL, subtitle=NULL, 
     age.opts=list(), time.opts=list(), threedim.opts=list(),totcount.opts=list(), 
     ex0.opts=list(),print="device", print.args=list(), dv.log=NA, ...)
}

\arguments{

\item{x}{A \code{\link{yourcast}} output object.}

\item{plots}{Vector of strings specifying the plot types to be generated. The available plots are age profile (\code{"age"}), time profile (\code{"time"}), 3D age/time profile (\code{"threedim"}), total counts by year (\code{"totcounts"}), and life expectancy at birth (\code{"ex0"}). The order in which the plots are specified is the order in which they are generated in the device. A maximum of 3 plot types is allowed. Default: \code{c("age","time")}.}

\item{title}{String. Main title for the plots. Concatenated with country name if \code{dataobj$G.names} exists. Default: \code{NULL}.}

\item{subtitle}{String. Subtitle for the plots. Default: \code{NULL}.}

\item{age.opts}{A list of options for the age profile plot. See \sQuote{Details}.}

\item{time.opts}{A list of options for the time profile plot. See \sQuote{Details}.}

\item{threedim.opts}{A list of options for the 3D plot. See \sQuote{Details}.}

\item{totcount.opts}{A list of options for the total counts plot. See \sQuote{Details}.}

\item{ex0.opts}{A list of options for the life expectancy at birth plot. See \sQuote{Details}.}

\item{print}{String specifying whether graphical output should be
    displayed sequentially on a device window (\code{"device"}) or saved
    directly to a \file{.pdf} file in the \code{dpath} (\code{"pdf"}). 
    Default: \code{"device"}.}

\item{print.args}{A list of options for the device. See \sQuote{Details}.}

\item{dv.log}{Logical. Is the dependent variable logged? If \code{NA}, the function parses the formula in the \code{\link{yourcast}} object to determine if the dependent variable is logged. Default: \code{NA}.}

\item{...}{Further arguments passed to or from other methods.}

}  %%%end of arguments

\details{Plots age profiles, time profiles, 3D age/time profiles, total counts (by time), and life expectancy at birth from \code{\link{yourcast}} output. Total count (by time) plots should only be used if forecasting counts (eg. population or death counts) and life expectancy at birth should only be used if forecasting mortality rates - both are common quantities of interest calculated from these forecasts and are thus included in the function. See \code{\link{lifetable}} for details of life expectancy at birth calculations from mortality rates. All plots are created using the \code{\link{ggplot2}} package, except for the 3D age/time profile, which is created using the \code{wireframe} function from the \code{lattice} library.
	
The function supports multiple plots in each call (up to 3) and arranges them horizontally in the order specified in \code{plots} argument. The device window is automatically resized to accomodate different number of plots, albeit device width and height may be specified manually through the \code{print.args} option (see below).

Plots are titled with the \code{title} and the \code{G.names} dataframe if it was supplied to \code{\link{yourcast}} in the \code{dataobj}. For example if \code{title="Respiratory Infections"} and the geographic identifier for that region is matched with \code{"Belize"}, the plot will be titled \dQuote{Respiratory Infections, Belize}. One or both labels will be utilized by the function if available. A subtitle may also be specified through \code{subtitle}.

It is important to note that \code{plot.yourcast} will only work if
all cross sections within the same geographic unit are of the same
dimensions. If, for example, a cross section for one age group has
fewer yearly observations than another from the same group, these
missing years must be filled in with \code{NA}, even if they occur in
the beginning of the sample period. This does not hold across
geographic units, however.

Finally, note that \code{plot.yourcast} opens a new device window for each new
plot. This is done so that the size of the device can be controlled to
keep the side-by-side plots from appearing distorted when launched.

\strong{Options for \sQuote{age}, \sQuote{time}, \sQuote{totcount}, and \sQuote{ex0} plots}

Options for the age profile, time profile, total count, and life expectancy at birth plots may be specified using the \code{age.opts}, \code{time.opts}, \code{totcount.opts}, and \code{ex0.opts} arguments, respectively, which are lists with any of following components:

\describe{
\item{\code{xlab}}{String specifying the x-axis label. Default is \code{"Age"} for \sQuote{age} plots and \code{"Time"} for \sQuote{time}, \sQuote{totcount}, and \sQuote{ex0} plots.}

\item{\code{ylab}}{String specifying the y-axis label. Default is \code{"Data and
    Forecasts"} for \sQuote{age} plots and \code{"Forecasts"} for \sQuote{time} plots.}
    
\item{\code{insamp.obs}}{Logical. Should the observed data be plotted? Default for \sQuote{age}, \sQuote{totcounts}, \sQuote{ex0} is \code{TRUE}; default for  \sQuote{time} is \code{FALSE}.}

\item{\code{insamp.predict}}{Logical. Should the predicted data be plotted for the in-sample period? Default: \code{TRUE}.}
}

Additionally, the following options may be passed to \code{age.opts} and \code{time.opts} lists:
\describe{
\item{\code{age.select}}{A numeric vector listing the ages to be plotted. If a scalar is supplied, it is used as the step size in a sequence of ages between the minimum and the maximum age. If \code{NULL}, all ages are plotted. Default: \code{NULL}.}

\item{\code{time.select}}{A numeric vector listing the times to be plotted. If a scalar is supplied, it is used as the step size in a sequence of times between the minimum and the maximum time. If \code{NULL}, all times are plotted. Default: \code{NULL}.}

\item{\code{unlog}}{Logical. Should the dependent variable be unlogged for \sQuote{age} and \sQuote{time} plots? Ignored if \code{dv.log} is set to \code{FALSE}. Default: \code{FALSE}.}
} % end of describe envir

Option \code{age.select} may be included in the \code{totcount.opts} list, but it functions differently than for the \sQuote{age} and \sQuote{time} plots. Instead of designating the ages to plot, \code{age.select} should be a numeric vector listing the ages to be summed over in the calculation of total counts over time. 

For \sQuote{ex0} plots, options \code{ax}, \code{a0}, and \code{nx} that control the calculation of life tables from mortality rates in \code{\link{yourcast}} objects may be included in the \code{ex0.opts} list. By default, \code{ax} is set to 0.5, \code{a0} is set to 0.07, and \code{nx} is set to \code{NULL}. See \code{\link{lifetable}} for additional information.

\strong{Options for \sQuote{threedim} plot}

Options for the 3D age/time profile may be specified through \code{threedim.opts}. The argument must be a list with any of the following components:
\describe{
\item{\code{xlab}}{String specifying the x-axis label. This corresponds to the time dimension. Default: \code{"Time"}.}

\item{\code{ylab}}{String specifying the y-axis label. This corresponds to the age dimension. Default: \code{"Age"}.}

\item{\code{zlab}}{String specifying the z-axis label. This corresponds to the dependent variable. Default: \code{"Forecasts"}}

\item{\code{insamp.pred}}{Logical. If \code{TRUE}, the forecasted data is plotted for the in-sample period. If \code{FALSE}, the observed data is plotted instead. Default: \code{TRUE}.}

\item{\code{unlog}}{Should the dependent variable be unlogged? Ignored if \code{dv.log} is set to \code{FALSE}. Default: \code{FALSE}.}

\item{\code{args.wireframe}}{A list of arguments (must be labeled) to be
  passed to \code{wireframe}. Note that some arguments to \code{wireframe}
  such as \code{zlab} and \code{screen} should be made as arguments to
  \code{threedim.opts}; they will be overwritten if found in
  \code{args.wireframe}.}

\item{\code{screen}}{List with three elements \sQuote{x}, \sQuote{y},
    and \sQuote{z} that rotate the viewing angle for the 3D plot. Default: \code{list(z=-40, x=-60, y=0)}.}
} % end of describe


\strong{Print options}

The function prints sequentially to the device or saves \file{.pdf} files
with the requested plot for each geographic unit in the sample. If
requested, \file{.pdf} files will be saved in a specified directory or
the working directory. Use \code{print.args} to set print options:
\describe{
\item{\code{height}}{Height of device. Default is 5 inches.}

\item{\code{width}}{Width of device. Default width adjusts automatically according to number of plots.}

\item{\code{filename}}{A string or vector of strings denoting the filename(s) of \file{.pdf}s created for each geographical area if \code{print="pdf"}. If a string, the \code{csid} code is automatically appended if multiple geographies so as to ensure unique filenames. If the vector is of the length equal to the number of geographies, the \code{csid} is not appended, albeit the order of the labels should be the same as the order of the areas in \code{dataobj$data}. Default: \code{NULL}, in which case plots are titled as \code{plotgggg.pdf}, where \code{gggg} denotes the \code{csid}.}

\item{\code{dpath}}{A string giving the directory where \file{.pdf} files will be output if \code{print="pdf"}. Default: \code{getwd()}.}
} % end of describe
} %%% end of details

\value{Device windows with requested plots or \file{.pdf} files saved
in the \code{dpath}.}

\author{Konstantin Kashin \email{kkashin@fas.harvard.edu}, based on an earlier version by Jon Bischof.}

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
# default
plot(ylc.allc,title="All causes of death")

\dontrun{
# output plots to pdf
plot(ylc.allc,title="All causes of death", print="pdf")

# output plots to pdf with names "acd"+csid
plot(ylc.allc,title="All causes of death", print="pdf",
print.args=list(filename="acd"))

# order matters
plot(ylc.allc,plots=c("time","age"),title="All causes of death")
plot(ylc.allc,plots=c("age","time"),title="All causes of death")}

# remove in-sampled predicted data from age profile plot
plot(ylc.allc,title="All causes of death",age.opts=list(insamp.predict=FALSE))

# plot every 10th age instead of every 5th age for age profile plot
plot(ylc.allc,title="All causes of death",
age.opts=list(insamp.predict=FALSE,age.select=10))

\dontrun{
# plot every 5th year instead of every year for time series plot
plot(ylc.allc,title="All causes of death",
age.opts=list(insamp.predict=FALSE),time.opts=list(time.select=5))

# only plot data from 2000+
plot(ylc.allc,title="All causes of death",
age.opts=list(insamp.predict=FALSE, time.select=seq(2001,2060,5)),
time.opts=list(time.select=seq(2001,2060,5)))

# plot unlogged mortality rates
plot(ylc.allc,title="All causes of death",
age.opts=list(insamp.predict=FALSE, unlog=TRUE),time.opts=list(unlog=TRUE))

# plot `threedim' plot in addition to age and time profiles
plot(ylc.allc,plots=c("threedim","age","time"),
title="All causes of death",age.opts=list(insamp.predict=FALSE))

# plot life expectancy at birth in addition to age and time profiles
# since these are 5-year age groups, assume ax=2.5 and a0=1.907
plot(ylc.allc,plots=c("ex0","age","time"),title="All causes of death",
age.opts=list(insamp.predict=FALSE), ex0.opts=list(a0=1.907,ax=2.5))}

} %%% end of examples