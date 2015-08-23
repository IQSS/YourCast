## ----setup, echo=FALSE---------------------------------------------------
# global chunk options
knitr::opts_chunk$set(cache = FALSE, collapse=TRUE, fig.path = "figures/", dev = "pdf")

## ----install_cran, eval=FALSE--------------------------------------------
#  install.packages("YourCast")

## ----install_git, echo=TRUE, eval=FALSE----------------------------------
#  devtools::install_github("IQSS/YourCast")

## ----load, echo=TRUE, eval=FALSE-----------------------------------------
#  library(YourCast)

## ----loaddata, eval=FALSE------------------------------------------------
#  ydata <- yourprep(tag="cancerMales")

## ----forecast, eval=FALSE------------------------------------------------
#  ylc <- yourcast(formula=log(rspi2/popu2) ~ time, dataobj=dta, model="LC", verbose=FALSE)

## ----seedata1, eval=FALSE------------------------------------------------
#  data(package="YourCast")

## ----seedata2, eval=TRUE-------------------------------------------------
system.file("data", package = "YourCast")
#data(package="YourCast")$results[,"Item"]

