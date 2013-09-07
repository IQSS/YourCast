

chapter11.4 <- function()
{
    message("Figure 11.4 in 'Demographic Forecasting'...")
    message("Lung cancer for Peru...")

    cat("\n")
    cat("Loading data...\n")
    data(chp.11.4)
    cat("> data(chp.11.4)\n")
    
    cat("\n")
    message("Formula for male lung disease...")
    ff <- log(lung2/popu2) ~ time + log(time -1876)
    cat("> ff <- log(lung2/popu2) ~ time + log(time -1876)\n")
    message("Running yourcast with OLS model...")
    user.prompt()
    
    cat("\n")
    cat("> yols <- yourcast(formula=ff, dataobj=chp.11.4)\n")
    yols <- yourcast(formula=ff, dataobj=chp.11.4)
    
    message("Generating the graphics for OLS...")
    user.prompt()
    
    cat("\n")
    cat("> plot(yols, title=\"Lung Cancer\",age.opts=list(insamp.predict=FALSE))\n")
    plot(yols, title="Lung Cancer",age.opts=list(insamp.predict=FALSE))
   
  }

chapter11.4()
