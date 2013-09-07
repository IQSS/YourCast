

chapter11.11 <- function()
{
  message("Figure 11.11 in 'Demographic Forecasting'...")
  message("Transportation accidents for Chile and Argentina...")

  cat("\n")
  cat("Loading data...\n")
  data(chp.11.11)
  cat("> data(chp.11.11)\n")
   
  message("Formula for male transportation accidents...")
  ff <- log(trns2/popu2) ~ log(gdp) + time^2
  cat("> ff <- log(trns2/popu2) ~ log(gdp) + time^2\n")
  
  message("Running yourcast with OLS model...")
  user.prompt()
  cat("> yols <- yourcast(formula=ff, dataobj=chp.11.11, model=\"OLS\",
                   elim.collinear=FALSE)\n")
  yols <- yourcast(formula=ff, dataobj=chp.11.11, model="OLS",
                   elim.collinear=FALSE)
  
  message("Generating the graphics for OLS...")
  user.prompt()
  cat("> plot(yols, title=\"Transportation Accidents\")\n")
  plot(yols, title="Transportation Accidents")
 
  }

chapter11.11()
