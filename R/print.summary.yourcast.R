print.summary.yourcast <- function(x, ...) {

if(!is.null(x$sample.frame)) {
cat("Model:",x$model,"\n")
cat("Number of cross sections:",x$numcs,"\n")
cat("Formula:",deparse(x$formula),"\n")

cat("\n")
cat("Observed period: ",x$sample.frame[1],"-",
    x$sample.frame[2],"\n",sep="")
cat("Forecast period: ",x$sample.frame[3],"-",
    x$sample.frame[4],"\n",sep="")

cat("\n")
cat("Smoothing parameters:\n")
print(round(x$params,4))

cat("\n")
cat("Geo units included:\n")
print(x$cntry.codes)

cat("\n")
cat("See 'help(plot.yourcast)' for instructions on how to plot observed and predicted 'y' values\n")
}

else {cat("This is a special output object from the 'ebayes' model that provides empirical estimates of the the smoothing parameters to be used in the 'map' model. There were no predicted values produced in this analysis. Please see the YourCast documentation for more details (especially the sections on smoothing).\n")
cat("\n")
cat("Estimated smoothing parameters:\n")
print(unlist(x$params))}
}
