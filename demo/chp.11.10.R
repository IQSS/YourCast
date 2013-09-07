
chapter11.10 <- function()
{
  message("Figure 11.10 in 'Demographic Forecasting'...")
  message("Breast cancer, Chile, Cuba, Belgium, and the Netherlands...")

  cat("\n")
  cat("Loading data...\n")
  data(chp.11.10)
  cat("> data(chp.11.10)\n")
  
  message("Formula for female breast cancer...")
  ff <- log(brst3/popu3) ~ log(hc) + log(gdp) + log(tobacco3) + log(fat) + time
  cat("ff <- log(brst3/popu3) ~ log(hc) + log(gdp) + log(tobacco3) + log(fat) + time\n")

  cat("\n")
  cat("Setting means for priors...\n")
  cat("> z.mean <- c(-11.391495,-10.147588, -9.305023, -8.692574, -8.232481,
              -7.953798, -7.798399, -7.678475, -7.577912, -7.433581,
              -7.293615, -6.926301)\n")
  z.mean <- c(-11.391495,-10.147588, -9.305023, -8.692574, -8.232481,
              -7.953798, -7.798399, -7.678475, -7.577912, -7.433581,
              -7.293615, -6.926301)
  names(z.mean) <- 5:16*5
  
  message("Running yourcast with MAP model...")
  user.prompt()
  cat("> ymap <- yourcast(formula=ff, dataobj=chp.11.10, model=\"MAP\",
                   elim.collinear=FALSE,
                   Ha.sigma=1.22,Ht.sigma=0.94,Hat.sigma=0.38,
                   zero.mean=z.mean,  Ha.deriv=c(0,0,1),
                   Hat.a.deriv=c(0,1), Hat.t.deriv=c(0,0,1),
                   Ht.deriv=c(0,0,1))\n")
  ymap <- yourcast(formula=ff, dataobj=chp.11.10, model="MAP",
                   elim.collinear=FALSE,
                   Ha.sigma=1.22,Ht.sigma=0.94,Hat.sigma=0.38,
                   zero.mean=z.mean,  Ha.deriv=c(0,0,1),
                   Hat.a.deriv=c(0,1), Hat.t.deriv=c(0,0,1),
                   Ht.deriv=c(0,0,1))
                   
  message("Generating the graphics for MAP...")
  user.prompt()

  cat("> plot(ymap, title=\"Breast Cancer\",plots=\"age\")\n")
  plot(ymap, title="Breast Cancer",plots="age")
  
  }

chapter11.10()
