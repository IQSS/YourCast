

chapter2.7 <- function()
  {
    message("Figure 2.7 in 'Demographic Forecasting'...")
    message("Part 1: Male suicide deaths, USA...")
    
    cat("\n")
    cat("Loading data...\n")
    data(chp.2.7.1)
    cat("> data(chp.2.7.1)\n")

    message("Formula for male suicide deaths...")
    cat("> ff.suic <- log(suic2/popu2) ~ time\n")
    ff.suic <- log(suic2/popu2) ~ time
    
    message("Running yourcast with model LC...")
    user.prompt()
    cat("> ylc.suic <- yourcast(formula=ff.suic, dataobj=chp.2.7.1, model=\"LC\",
                      sample.frame=c(1950,2000,2001,2060))\n")
    ylc.suic <- yourcast(formula=ff.suic, dataobj=chp.2.7.1, model="LC",
                      sample.frame=c(1950,2000,2001,2060))

    message("Generating the graphics for LC...")
    user.prompt()

    cat("> plot(ylc.suic,title=\"Male suicide\",age.opts=list(insamp.predict=FALSE))\n")
    plot(ylc.suic,title="Male suicide",age.opts=list(insamp.predict=FALSE))

      
    message("Part 2: Female digestive disease, Hungary...\n")
    user.prompt()
    
    cat("\n")
    cat("Loading data...\n")
    data(chp.2.7.2)
    cat("> data(chp.2.7.2)\n")
    
    message("Formula for female digestive disease...")
    cat("> ff.dgst <- log(dgst3/popu3) ~ time\n")
    ff.dgst <- log(dgst3/popu3) ~ time
    
    message("Running yourcast with model LC...\n")
    user.prompt()
    cat("> ylc.dgst <- yourcast(formula=ff.dgst, dataobj=chp.2.7.2, model=\"LC\",
                      sample.frame=c(1950,2000,2001,2060))\n")
    ylc.dgst <- yourcast(formula=ff.dgst, dataobj=chp.2.7.2, model="LC",
                      sample.frame=c(1950,2000,2001,2060))

    message("Generating the graphics for LC...")
    user.prompt()

    cat("> plot(ylc.dgst,title=\"Female digestive disease\",age.opts=list(insamp.predict=FALSE))\n")
    plot(ylc.dgst,title="Female digestive disease",age.opts=list(insamp.predict=FALSE))    
    
    message("Part 3: Female cervix cancer, United Kingdom")
    user.prompt()
    
    cat("\n")
    cat("Loading data...\n")
    data(chp.2.7.3)
    cat("> data(chp.2.7.3)\n")
         
    message("Formula for female cervix cancer...")
    cat("> ff.cerv <- log(cerv3/popu3) ~ time\n")
    ff.cerv <- log(cerv3/popu3) ~ time
   
    message("Running yourcast with model LC...")
    user.prompt()
    cat("> ylc.cerv <- yourcast(formula=ff.cerv, dataobj=chp.2.7.3, model=\"LC\",
                      sample.frame=c(1950,2000,2001,2060))\n")
    ylc.cerv <- yourcast(formula=ff.cerv, dataobj=chp.2.7.3, model="LC",
                         sample.frame=c(1950,2000,2001,2060))

    message("Generating the graphics for LC...")
    user.prompt()

    cat("> plot(ylc.cerv,title=\"Female cervix cancer\",age.opts=list(insamp.predict=FALSE))\n")
    plot(ylc.cerv,title="Female cervix cancer",age.opts=list(insamp.predict=FALSE))

  }


chapter2.7()
