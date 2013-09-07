print.yourcast <- function(x, ...) {

if(!is.null(x$call)) {
sample.frame <- x$aux$sample.frame
  
cat("Model:",x$call$model,"\n")
cat("Number of cross sections:",length(x$yhat),"\n")
cat("Formula:",deparse(x$call$formula),"\n")
cat("\n")
cat("Observed period: ",sample.frame[1],"-",sample.frame[2],"\n",
    sep="")
cat("Forecast period: ",sample.frame[3],"-",sample.frame[4],"\n",
    sep="")
cat("\n")
cat("Smoothing parameters:\n")
names(x$params) <- c("Ha.sigma","Ht.sigma","Hat.sigma")
print(round(x$params,4))

cat("\n")
cat("See 'help(plot.yourcast)' for instructions on how to plot observed and predicted 'y' values\n") }

else {cat("This is a special output object from the 'ebayes' model that provides empirical estimates of the the smoothing parameters to be used in the 'map' model. There were no predicted values produced in this analysis. Please see the YourCast documentation for more details (especially the sections on smoothing).\n")
cat("\n")
cat("Estimated smoothing parameters:\n")
print(unlist(x$summary))}
}
