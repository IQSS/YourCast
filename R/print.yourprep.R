print.yourprep <- function(x, ...) {

# Print number of cross sections
cat("Total number of cross sections:",length(x$data),"\n")
cat("\n")

cat("CSID coding rule:",x$index.code,"\n")
cat("\n")
  
cat("This is a yourprep() output object to be used in the 'dataobj'
arugment of yourcast() to supply the function with the original data
and labeling information. To examine the content of this output
object, use names().\n")

}
