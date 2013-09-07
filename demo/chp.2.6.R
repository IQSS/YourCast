chapter2.6 <- function()
{
    message("Figure 2.6 in 'Demographic Forecasting'...")
    message("Part 1: All causes of death, New Zealand and Hungary...")
    
    cat("\n")
    cat("Loading data...\n")
    data(chp.2.6.1)
    cat("> data(chp.2.6.1)\n")

    message("Formula for all causes of death in male population...")
    cat("> ff <- log(allc2/popu2) ~  time\n")
    ff.allc <- log(allc2/popu2) ~  time
    
    message("Running yourcast with model LC...")
    cat("> ylc.allc <- yourcast(formula=ff.allc, dataobj=chp.2.6.1, model=\"LC\",
                       elim.collinear=FALSE,
                       sample.frame=c(1950,2000,2001,2060))\n")
    ylc.allc <- yourcast(formula=ff.allc, dataobj=chp.2.6.1, model="LC",
                       elim.collinear=FALSE,
                       sample.frame=c(1950,2000,2001,2060))

    message("Generating the graphics for LC...")
    user.prompt()

    cat("> plot(ylc.allc,title==\"All causes of death\",age.opts=list(insamp.predict=FALSE))\n")
    plot(ylc.allc,title="All causes of death",age.opts=list(insamp.predict=FALSE))
   

    message("Part 2: Transportation accidents for males, Portugal")
    user.prompt()
    
    cat("\n")
    cat("Loading data...\n")
    data(chp.2.6.2)
    cat("> data(chp.2.6.2)\n")
    
    message("Formula for transportation accidents in male population...")
    cat("> ff.trns <- log(trns2/popu2) ~ time\n")
    ff.trns <- log(trns2/popu2) ~ time
   
    message("Running yourcast with model LC...")
    user.prompt()
    cat("> ylc.trns <- yourcast(formula=ff.trns, dataobj=chp.2.6.2, model=\"LC\",
                      sample.frame=c(1950,2000,2001,2060))\n")
    ylc.trns <- yourcast(formula=ff.trns, dataobj=chp.2.6.2, model="LC",
                      sample.frame=c(1950,2000,2001,2060))
    
    message("Generating the graphics for LC...")
    user.prompt()
    cat("> plot(ylc.trns,title=\"Transportation accidents\",age.opts=list(insamp.predict=FALSE))\n")
    plot(ylc.trns,title="Transportation accidents", age.opts=list(insamp.predict=FALSE))
 
  }

chapter2.6()
