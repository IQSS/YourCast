
chapter11.9 <- function()
{
   message("Figure 11.9 in 'Demographic Forecasting'...")
   message("Breast cancer, Croatia...")
   
   cat("\n")
   message("First loading dataobj with 43 related countries (including
country of interest) to estimate smoothing parameters. Only cross
sections with more than 20 yearly observations used...")
   data(chp.11.9.1)
   cat("> data(chp.11.9.1)\n")

   cat("\n")
   message("Formula for female breast cancer...")
   ff <- log(brst3/popu3) ~ log(hc) + log(gdp) + log(tobacco3) + log(fat) + time
   cat("> ff <- log(brst3/popu3) ~ log(hc) + log(gdp) + log(tobacco3) + log(fat) + time\n")
  
   message("Preliminary run of yourcast with model 'ebayes'...")
   user.prompt()

   cat("\n")
   cat("Setting means for priors...\n")
   cat("> z.mean <- c(-11.391495, -10.147588, -9.305023, -8.692574,
        -8.232481, -7.953798, -7.798399,-7.678475, -7.577912,
        -7.433581, -7.293615, -6.926301)\n")
   z.mean <- c(-11.391495,-10.147588, -9.305023, -8.692574, -8.232481,
               -7.953798, -7.798399,-7.678475, -7.577912, -7.433581,
               -7.293615, -6.926301)
   names(z.mean) <- 5:16*5

   cat("\n")
   cat("> yebayes <- yourcast(formula=ff, dataobj=chp.11.9.1,
                      elim.collinear= FALSE,
                      model=\"ebayes\",Ha.deriv=c(0,0,1), Ha.sigma=0.2,
                      Hat.sigma=NA, Ht.sigma=NA,zero.mean=z.mean)\n")
   yebayes <- yourcast(formula=ff, dataobj=chp.11.9.1, elim.collinear= FALSE,
                      model="ebayes",Ha.deriv=c(0,0,1), Ha.sigma=0.2,
                      Hat.sigma=NA, Ht.sigma=NA,zero.mean=z.mean)
  
  auto <- unlist(yebayes$summary)
  distvec   <- yebayes$summary.vec
  d1.a.vec <- distvec$d1.a
  d1.t.vec <- distvec$d1.t
  dtda.vec <- distvec$dtda
  SD.vec <- distvec$SD
   
  cat("\n")
  message("Plotting histogram of estimated priors...")
  user.prompt()
  depvar <- YourCast:::leftside.formula(ff)$numerator
      histograph(d1.a.vec, d1.t.vec,dtda.vec,SD.vec, depvar,
                 model="ebayes", graphics.file=NA)
    
  d1.a <- auto["d1.a"]
  d1.t <- auto["d1.t"]
  dtda <- auto["dtda"]
  SD   <- auto["SD"]


  cat("\n")
  cat("Now plugging in estimated smoothing parameters to 'map' model...\n")
  message("Breast cancer, Croatia...")
  
  cat("\n")
  message("Loading Croatia data only; only can run one country at a
time when using 'ebayes'...")
  data(chp.11.9.2)
  cat("> data(chp.11.9.2)\n")
   
  message("Formula for female breast cancer...")
  cat("> ff <- log(brst3/popu3) ~ log(hc) + log(gdp) + log(tobacco3) + log(fat) + time\n")
  
  
  message("First running yourcast with model OLS for comparison...")
  user.prompt()

  cat("yols <- yourcast(formula=ff, dataobj=chp.11.9.2, model=\"OLS\",
       elim.collinear=FALSE)\n")
  yols <- yourcast(formula=ff, dataobj=chp.11.9.2, model="OLS",
                   elim.collinear=FALSE)

  cat("\n") 
  message("Generating the graphics for OLS...")
  user.prompt()
  cat("> plot(yols, title=\"Breast cancer\")\n")
  plot(yols, title="Breast cancer")

  cat("\n")
  message("Now running yourcast with 'map' model; this may take several minutes")
  user.prompt()
   
  cat("\n")
  cat("Setting means for priors...\n")
  cat("> z.mean <- c(-11.391495, -10.147588, -9.305023, -8.692574,
        -8.232481, -7.953798, -7.798399,-7.678475, -7.577912,
        -7.433581, -7.293615, -6.926301)\n")
  z.mean <- c(-11.391495,-10.147588, -9.305023, -8.692574, -8.232481,
               -7.953798, -7.798399,-7.678475, -7.577912, -7.433581,
               -7.293615, -6.926301)
  names(z.mean) <- 5:16*5
  sims <- c(sims=100)

  cat("ymap <- yourcast(model=\"map\", Ha.sigma=c(0.1,1.5,6,d1.a, SD=NULL),
                   Ht.sigma=c(0.1,1.5,6,d1.t, sims),
                   Hat.sigma=c(0.1,1.5,6,dtda),
                   zero.mean=z.mean,  Ha.deriv=c(0,0,1),
                   Hat.a.deriv=c(0,1), Hat.t.deriv=c(0,0,1),
                   Ht.deriv=c(0,0,1))\n")
  ymap <- yourcast(model="map", Ha.sigma=c(0.1,1.5,6,d1.a, SD=NULL),
                   Ht.sigma=c(0.1,1.5,6,d1.t, sims),
                   Hat.sigma=c(0.1,1.5,6,dtda),
                   zero.mean=z.mean,  Ha.deriv=c(0,0,1),
                   Hat.a.deriv=c(0,1), Hat.t.deriv=c(0,0,1),
                   Ht.deriv=c(0,0,1))

  message("Generating the graphics for MAP...")
  user.prompt()
  cat("plot(ymap,title=\"Breast cancer\")\n")
  plot(ymap,title="Breast cancer")
 
  }

chapter11.9() 
