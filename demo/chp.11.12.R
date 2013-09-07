


chapter11.12 <- function()
{
  message("Figure 11.12 in 'Demographic Forecasting'...")
  message("Transportation accidents for Argentina...")

  cat("\n")
  cat("Loading data...\n")
  data(chp.11.12)
  cat("> data(chp.11.12)\n")
  
  message("Formula for male transportation accidents...")
  cat("> ff <- log((trns2+ 0.5)/popu2) ~ log(gdp) + time + (time-1900)^2\n")
  ff <- log((trns2+ 0.5)/popu2) ~ log(gdp) + time + (time-1900)^2
  
  message("Running yourcast with MAP model...")
  user.prompt()

  cat("\n")
  cat("Setting means for priors...\n")
  cat("> zmean <- c(-8.334608,-7.848482,-7.940896,-8.022037,-8.062267,
             -8.077525,-8.041380,-8.028983,-7.990504, -7.919344,
             -7.874932, -7.708156, -7.557175, -7.330945)\n")
  zmean <- c(-8.334608,-7.848482,-7.940896,-8.022037,-8.062267,-8.077525,-8.041380,-8.028983,
             -7.990504, -7.919344, -7.874932, -7.708156, -7.557175, -7.330945)
  names(zmean) <- 3:16*5

  cat("\n")
  cat("> ymap <- yourcast(formula=ff, dataobj=chp.11.12, model=\"map\",
                     elim.collinear=TRUE,zero.mean=zmean,  
                     Ha.sigma=1.5, Ht.sigma=0.94, Hat.sigma=0.34,
                     Hct.sigma=NA, Ha.deriv=c(0,1,0),
                     Hat.a.deriv=c(0,1),
                     Hat.t.deriv=c(0,0,1),
                     Ht.deriv=c(0,0,1), low.pow=F)\n")
  ymap <- yourcast(formula=ff, dataobj=chp.11.12, model="map",
                     elim.collinear=TRUE,zero.mean=zmean,  
                     Ha.sigma=1.5, Ht.sigma=0.94, Hat.sigma=0.34,
                     Hct.sigma=NA, Ha.deriv=c(0,1,0),
                     Hat.a.deriv=c(0,1),
                     Hat.t.deriv=c(0,0,1),
                     Ht.deriv=c(0,0,1), low.pow=F)

  
  message("Generating the graphics for map...")
  user.prompt()
  cat("> plot(ymap, title=\"Transportation Accidents\")\n")
  plot(ymap, title="Transportation Accidents")

  }

chapter11.12()
