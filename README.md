# YourCast

[![Travis-CI Build Status](https://travis-ci.org/IQSS/YourCast.png?branch=master)](https://travis-ci.org/IQSS/YourCast)

[YourCast][] is an R package that makes forecasts by running sets of linear regressions together in a variety of sophisticated ways. YourCast avoids the bias that results when stacking datasets from separate cross-sections and assuming constant parameters, and the inefficiency that results from running independent regressions in each cross-section.

Read the YourCast [vignette][] for additional details.

## How to install

### Installation requirements
`YourCast` requires [R][] version 3.2.1 or higher.

### Installation from CRAN
```R
install.packages("YourCast")
```

[vignette]: http://cran.r-project.org/web/packages/YourCast/vignettes/YourCast.pdf
[YourCast]: http://gking.harvard.edu/yourcast
[R]: http://cran.r-project.org
