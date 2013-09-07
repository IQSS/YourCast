

chapter11.7 <- function()
  
{
   message("Figure 11.7 in 'Demographic Forecasting'...")
   message("Lung cancer, Peru...")
   
   cat("\n")
   message("First loading dataobj with 51 related countries (including
country of interest) to estimate smoothing parameters. Only cross
sections with more than 20 yearly observations used...")
   data(chp.11.7.1)
   cat("> data(chp.11.7.1)\n")

   cat("\n")
   message("Formula for male lung disease...")
   ff <- log(lung2 /popu2) ~ time + log(time - 1876) + log(gdp) + log(tobacco2)
   cat("ff <- log(lung2 /popu2) ~ time + log(time - 1876) + log(gdp) + log(tobacco2)\n")
   
   message("Preliminary run of yourcast with model 'ebayes'...")
   user.prompt()

   cat("\n")
   cat("Setting means for priors...\n")
   cat("zmean <- c(-11.123375, -10.188764,  -9.239278,  -8.419195,  -7.699627,
              -7.141269,  -6.702073,-6.400106,  -6.223774,  -6.153787,  -6.231790)\n")
   zmean <- c(-11.123375, -10.188764,  -9.239278,  -8.419195,  -7.699627,
              -7.141269,  -6.702073,-6.400106,  -6.223774,  -6.153787,  -6.231790)
   names(zmean) <- 6:16*5

   cat("\n")
   cat("yebayes <- yourcast(formula=ff, dataobj=chp.11.7.1,
                       elim.collinear= FALSE,
                       model=\"ebayes\", zero.mean=zmean,
                       Ha.deriv=c(0,0,1), Ha.sigma=0.3,
                       Hat.sigma=NA, Ht.sigma=NA)\n")
   yebayes <- yourcast(formula=ff, dataobj=chp.11.7.1,
                       elim.collinear= FALSE,
                       model="ebayes", zero.mean=zmean,
                       Ha.deriv=c(0,0,1), Ha.sigma=0.3,
                       Hat.sigma=NA, Ht.sigma=NA)
   
   auto <- unlist(yebayes$summary)
   distvec   <- yebayes$summary.vec
   d1.a.vec  <- distvec$d1.a
   d1.t.vec  <- distvec$d1.t
   dtda.vec  <- distvec$dtda
   SD.vec <- distvec$SD

   cat("\n")
   message("Plotting histogram of estimated priors...")
   user.prompt()
   depvar <- YourCast:::leftside.formula(ff)$numerator
   fs <- histograph(d1.a.vec, d1.t.vec,dtda.vec,SD.vec, depvar,
                 model="ebayes", graphics.file=NA)
  
   d1.a <- auto["d1.a"]
   d1.t <- auto["d1.t"]
   dtda <- auto["dtda"]
   SD   <- auto["SD"]
   sims <- c(sims=50)

   cat("\n")
   cat("Now plugging in estimated smoothing parameters to 'map' model...\n")
   message("Lung cancer, Peru...")
   user.prompt()

   cat("\n")
   message("Loading Peru data only...")
   data(chp.11.7.2)
   cat("> data(chp.11.7.2)\n")

   cat("\n")
   message("Formula for male lung disease...")
   ff <- log(lung2/popu2) ~ time + log(time - 1876)
   cat("ff <- log(lung2/popu2) ~ time + log(time - 1876)\n")
   
   message("Running yourcast with 'map' model; this may take several minutes")
   user.prompt()
   cat("ymap <- yourcast(formula=ff, dataobj=chp.11.7.2, model=\"map\",
                    elim.collinear=FALSE, 
                    Ha.sigma=c(0.1,1.5,5,d1.a, SD=NULL),
                    Ht.sigma=c(0.1,1.5,5,d1.t, sims),
                    Hat.sigma=c(0.05,1.5,5,dtda),
                    zero.mean=zmean,
                    Ha.deriv=c(0,0,1),
                    Hat.a.deriv=c(0,1), Hat.t.deriv=c(0,0,1),
                    Ht.deriv=c(0,0,1))\n")
   ymap <- yourcast(formula=ff, dataobj=chp.11.7.2, model="map",
                    elim.collinear=FALSE, 
                    Ha.sigma=c(0.1,1.5,5,d1.a, SD=NULL),
                    Ht.sigma=c(0.1,1.5,5,d1.t, sims),
                    Hat.sigma=c(0.05,1.5,5,dtda),
                    zero.mean=zmean,
                    Ha.deriv=c(0,0,1),
                    Hat.a.deriv=c(0,1), Hat.t.deriv=c(0,0,1),
                    Ht.deriv=c(0,0,1))
  
   
   cat("\n")
   message("Generating the graphics for MAP...")
   user.prompt()
   
   cat("> plot(ymap, title=\"Lung Cancer\")\n")
   plot(ymap, title="Lung Cancer")
       
   
 }

chapter11.7()
