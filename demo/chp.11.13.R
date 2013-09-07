

chapter11.13 <- function()
{
  message("Figure 11.13 in 'Demographic Forecasting'...")
  message("Transportation accidents for Argentina, Chile, Canada, Colombia,
       Costa Rica, Cuba, and USA...")

  cat("\n")
  cat("Loading data for Argentina, Chile, Canada, Colombia,
       Costa Rica, Cuba, and USA...\n")
  data(chp.11.13)
  cat("> data(chp.11.13)\n")

  cat("\n")
  message("Formula for male transportation accidents...")
  cat("> ff <- log(trns2/popu2) ~ log(gdp) + (time-1900)^2  + time\n")
  ff <- log(trns2/popu2) ~ log(gdp) + (time-1900)^2  + time

  
  message("Running BAYES model, this may take over 10 minutes; final
counter to reach 500...")
  user.prompt()

  cat("Setting means for priors...\n")
  cat("> zmean <- c(-8.334608,-7.848482,-7.940896,-8.022037,-8.062267,-8.077525,
             -8.041380,-8.028983, -7.990504, -7.919344, -7.874932,
             -7.708156, -7.557175, -7.330945)\n")
  zmean <- c(-8.334608,-7.848482,-7.940896,-8.022037,-8.062267,-8.077525,
             -8.041380,-8.028983, -7.990504, -7.919344, -7.874932,
             -7.708156, -7.557175, -7.330945)
  names(zmean) <- 3:16*5

  cat("\n")
  cat("> ybayes <- yourcast(formula=ff, dataobj=chp.11.13, model=\"bayes\",
                     elim.collinear=T,zero.mean=zmean, low.pow=F,
                     nsample=500, Ha.sigma=0.3,Ha.sigma.sd=0,
                     Ht.sigma=1, Ht.sigma.sd=0,Hat.sigma=0.01,
                     Hat.sigma.sd=0,Hct.sigma=0.01,Hct.sigma.sd=0,
                     Hct.t.deriv=c(0,0,1), Ha.deriv=c(0,1),
                     Hat.a.deriv=c(0,1),Ht.deriv=c(0,0,1),
                     Hat.t.deriv=c(0,0,1))\n")
  ybayes <- yourcast(formula=ff, dataobj=chp.11.13, model="bayes",
                     elim.collinear=T,zero.mean=zmean, low.pow=F,
                     nsample=500, Ha.sigma=0.3,Ha.sigma.sd=0,
                     Ht.sigma=1, Ht.sigma.sd=0,Hat.sigma=0.01,
                     Hat.sigma.sd=0,Hct.sigma=0.01,Hct.sigma.sd=0,
                     Hct.t.deriv=c(0,0,1), Ha.deriv=c(0,1),
                     Hat.a.deriv=c(0,1),Ht.deriv=c(0,0,1),
                     Hat.t.deriv=c(0,0,1))
  
  message("Generating the graphics for BAYES...")
  user.prompt()
  plot(ybayes, title="Transportation Accidents")
  
  }

chapter11.13()
