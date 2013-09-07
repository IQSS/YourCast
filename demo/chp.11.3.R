

chapter11.3 <- function()
{
    message("Figure 11.3 in 'Demographic Forecasting'...")
    message("Respiratory infections, Bulgaria...")

    cat("\n")
    cat("Loading data...\n")
    data(chp.11.3)
    cat("> data(chp.11.3)\n")
   
    cat("\n")
    message("Formula for male respiratory infections...")
    ff <- log(rspi2/popu2) ~ time 
    cat("> ff <- log(rspi2/popu2) ~ time\n")
    message("Running yourcast with model 'map'...")
    user.prompt()

    cat("Setting means for priors...\n")
    cat("> zmean <- c(-7.474950, -10.391050, -10.745170, -10.511022,
        -10.450573, -10.321841, -10.066730, -9.721626,  -9.362865,
        -8.995520,  -8.607914,  -8.233437,  -7.752187,  -7.240793,
        -6.626354,  -6.019082,  -4.938154)\n")
    zmean <- c(-7.474950, -10.391050, -10.745170, -10.511022, -10.450573, -10.321841, -10.066730, 
               -9.721626,  -9.362865,  -8.995520,  -8.607914,  -8.233437,  -7.752187,  -7.240793, 
               -6.626354,  -6.019082,  -4.938154)
    names(zmean) <- 0:16*5

    cat("\n")
    cat("ymap <- yourcast(formula=ff, dataobj=chp.11.3, model=\"map\",
                     Ha.sigma=0.2, Ht.sigma=NA,Hat.sigma=0.05,
                     Hat.a.deriv=c(0,1), Hat.t.deriv=c(0,1),
                     zero.mean=zmean)\n")
    ymap <- yourcast(formula=ff, dataobj=chp.11.3, model="map",
                     Ha.sigma=0.2, Ht.sigma=NA,Hat.sigma=0.05,
                     Hat.a.deriv=c(0,1), Hat.t.deriv=c(0,1),
                     zero.mean=zmean)
    
    cat("\n")
    message("Generating the graphics for MAP...")
    user.prompt()
    
    cat("> plot(ymap, title=\"Respiratory Infections\",age.opts=list(insamp.predict=FALSE))\n")
    plot(ymap,title="Respiratory Infections",age.opts=list(insamp.predict=FALSE))

  }

chapter11.3()
