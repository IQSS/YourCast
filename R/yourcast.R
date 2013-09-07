
##
## DESCRIPTION:  The driver function, which the user calls with his choice
##               of parameters or globals that initializes the simulation, and the
##               model for the predictions. We save two environmnets, the
##               env.base created by yourcast and env.who created by
##               the initialization file WHO(). The globals or static variables
##               are shared by all the functions. 
##               It reads arguments of yourcast and pass them to WHO() so that the names
##               are changed internally. It checks dataobj, which is a list   
##               of three elements; the data list whose elements are identified by csid
##               or cntry+code combination and each element is a matrix of rows corresponding
##               to the year series and columns the covariates or dependemnt variable data. 
##               The index.code that explicity indicates the digits for cntry, age and year
##               say ggggaatttt has 4 digits for cntry, two for age and 4 for year.
##               Optionally a matrix with the code-cntry correspondance one for every cntry.
##               For our programs we use dataobj$G.names = geo.names(),
##               defined in build.from.formula. Parameter sample.frame is a 4 element
##               vector that contains the years of data
##               and prediction, with second element being the yrest dividing the periods 
##               for insample (depvar data use to forecast) and outsample (forecastin years). 
##               Updates env.who with the input globals, if any and then runs the simulation.
##               The formula to apply to the covariates and dependent variable.
##               Optionally, you may just give a number for sample.frame indicating the
##               yrest dividing insample and outsample since the time series is part of
##               the list element data of dataobj.
##               The userfile where the user may stored the parameters instead of explicitly
##               changing them in the function call.
##               The order of preference is first the arguments stored in userfile and if no
##               file is specified then the arguments of yourcast.  If both dataobj and formula
##               are missing then we will check if the file yourcast.savetmp exists in the current
##               working directory and upload it with the resulst of previous simulation. Those
##               results are preprocessing lists with the depvar insample and outsample and the
##               covariates.  The user may want to use previous results with different models or
##               different parameters related to models only.
##               In addition, dataobj may refer to a file and the data, index.code and, optionally,
##               G.names will be uploaded from the file. Also, the file dataobj may contain just
##               the preprocessing of any previous simulation and then,
##               yourcast will just run the model part. 
##               It chooses the model to forecast with model from either argument of yourcast
##               or the userfile. 
##               Produces the foreocast output and plots it with plot.forecast.
##              
## INPUT: formula, i.e. log(allc2/popu2 + 0.5) ~ gdp^2 + tobacco2 + time
##        (important DO NOT FORGET  ADD 0.5 if using the log)
##        dataobj, a list with 4 elements: "data" "index.code" "proximity"  "G.names"
##        dataobj$data is a list of csid units with a time serises for covaraites and diseases
##        dataobj$index.code is a string "ggggaatttt", 4 g's for geographic index
##        2 a's for age group, and 4 t's for time
##        dataobj$proximity the adjacency file
##        dataobj$G.names the file giving names for country codess
##        sample.frame, a vector c(time=1950, 2000, 2001, 2020); year first obv (for depvar),
##        year year last insample, year first outsample, year last forecast
##
## FORMAT: yourcast(userfile="myfile.R",formula=ff, dataobj=myobject,model="poisson",.... ) 
##
## VALUE:  The forecast object and the graphics
##
## WRITTEN BY: Elena Villalon 
##             evillalon@iq.harvard.edu
##             IQSS, Harvard University
##
## LAST MODIFIED: 10/14/2005
## 
## 
## ************************************************************************
## ************************************************************************
### loading the libraries first

yourcast <- function(formula=NULL, dataobj=NULL,
                     sample.frame=c(1950, 2000, 2001, 2030), 
                     standardize=TRUE, elim.collinear=FALSE,
                     tol=0.9999, solve.tol = 1.e-10,svdtol=10^(-10),
                     userfile=NULL, savetmp = T,model.frame=FALSE, 
                     debug = F,  rerun= "yourcast.savetmp", 
### specific to models
                     model="OLS",zero.mean=FALSE,
###                     ols.sigma.param=list(sigma.bar=1,use.deaths=FALSE,average=TRUE,model.based=FALSE), 
          #### smooth over ages           
                     Ha.sigma = 0.3,
                     Ha.sigma.sd= 0.1, Ha.deriv=c(0,0,1),
                     Ha.age.weight=0, Ha.time.weight=0,
          #### smooth over time
                     Ht.sigma= 0.3,
                     Ht.sigma.sd=0.1,  Ht.deriv=c(0,0,1),
                     Ht.age.weight=0, Ht.time.weight=0,
          #### smooth over age-time
                     Hat.sigma=0.2,
                     Hat.sigma.sd=0.1,Hat.a.deriv=c(0,1),Hat.t.deriv=c(0,1), 
                     Hat.age.weight=0,Hat.time.weight=0,
          #### smooth over cntry-time
                     Hct.sigma=0.3, Hct.sigma.sd =0.1,
                     Hct.t.deriv=1, Hct.time.weight = 0,
                     LI.sigma.mean=0.2,LI.sigma.sd = 0.1, nsample= 500,
                     low.pow=T, verbose=TRUE)
{
 #USER INPUT ERROR CHECKING
 # sma is not really needed (Ferdi) 
 #Load lpSolve package and print error if it is not available
  #if (!require(sma, quietly = TRUE))  { 
  #  stop("The sma package is required.\nInstall it from CRAN.")
  #}
  guirun <- NULL ###provision to run the simulation with a GUI
### done with data input check model  
   
### some initialization and clean up
  try(rm(env.base, env.who, inherits=T), silent=T);
### I'm defining base environment
  env.base <- environment();
 

  if(is.character(model))
    model <- trim.blanks(toupper(model))
  if(length(dataobj) > 0 )
    dataobj0 <- dataobj
  if(is.character(dataobj))
    dataobj0 <- trim.blanks(dataobj0) ### saving the initial input
  if(length(dataobj) >0)
    checkdataobj(dataobj$data, sample.frame)
  rerun <- trim.blanks(rerun)
  if(length(rerun) > 0 && !is.character(rerun))
    stop(rerun)
  
  if(length(userfile) > 0)
    messout(paste("Inputs in userfile ", userfile, "overwrite yourcast arguments"), verbose)
          
### previous run of yourcast may be stored in file
###  print(formals(yourcast)$Ha.sigma)
###  print(match.call(yourcast))
###  print(str(match.call(yourcast)))
 
  if((length(dataobj$data) <= 0 || length(formula) <= 0) && 
     !file.exists(rerun) && length(userfile) <= 0)
    stop("Insufficient data")

   
### for debug = T then make it available as global
  if(debug){
    messout("Setting to mode debugging with the global env.base", verbose)
    env.base <<- env.base
  }
  
##parsing for insample and outsample periods year
  yrest <- NULL
  if(length(sample.frame) > 0)
    yrest <- find.yrest(sample.frame)
    
  year.digits <- nchar(yrest)

### get the number of digist for each of cntry,yourcast.savetmpages and years
### from the time cross sectional identifiers are stored in
### dataobj$index.code, parse it and assign them to env.base
  
  if(length(dataobj$index.code) > 0){
    ix.code <- parse.index.code(icode=dataobj$index.code, yrest)
    assign.number.digits(ix.code, env.base)    
  }

  Hct.c.deriv <-  NULL
  if(length(dataobj$proximity) >0)
    Hct.c.deriv <- dataobj$proximity
  
### check if the user has provided other inputs with reflection
### give the name of the driver and call for the simulation
  
  driver <- match.call()

### print(driver)
  driver <- as.character(driver)
### name of the calling function 
  mortality.driver <- driver[1]
### print(mortality.driver)
### give me its arguments
 
  args  <- names(formals(mortality.driver))
  names(args) <- args
  
###  print(args) 
###
### is the call coming from the global environment or from a GUI
### if from command line; creates env.who and stored all input data there
  
  if(length(guirun) <= 0  || !is.environment(guirun)){
    env.who <- WHO.yourcast(args, env.base)
    sims <- get("N", envir=env.who)
    
  }else 
  env.who <- guirun

### find out the arguments of yourcast to be added to output

  callst <- call.yourcast(formlst=formals(mortality.driver),
                          callmatch= match.call(), model.frame, rerun)

### exclude dataobj (way too big) if model.frame =F,i.e. include just its name
   callst <- call.dataobj1(callst, model.frame, env.who)
 
### exclude dataobj (way too big)  

### tobe added to output
  userfilelst <- NA
  if(length(userfile) > 0){
    userfilelst <- args.userfile(userfile,env.base)
    messout("Retreiving userfile parameters for...",verbose)
    messout(paste(names(userfilelst),sep=" "),verbose)
    
     
  }
### finds out the inputs for dataobj, if any.  Process them
### It migth be either a file or an object in memory; assign to env.base

  if(length(dataobj) > 0 ) {

    dataobj <- parse.dataobj(dataobj,env.base)
  
    
        
  }
###last check the userfile and if it is not null, then its inputs
### parameters will overwrite the arguments of yourcast
  if(any(!is.na(userfilelst)) )
    args1 <- try(get("userfile", env.who), silent=T)
  else
    args1 <- NULL

### update args with the user supplied file, if any
  if (class(args1) != "try-error" && !is.null(args1)){
    messout("Updating global with the userfile parameters", verbose)
    
    update.args(args, userfile, env.base);
    
    env.who <- WHO.yourcast(args, env.base, env.who)
    
  }
 
### rownames instead of tagging with years change to
### time cross-sectional id's for consistency with old code
 if(length(dataobj$data) > 0) {
     dataobj <- sort.dataobj(dataobj)
     ix <- 1:length(dataobj$data)
     names(ix) <- names(dataobj$data)
     data <- lapply(ix, FUN="elimallna.csid",dataobj$data, verbose)
     data <- conversion.rownames(dataobj$data, year.digits)

     dataobj$data <- data
     aux <- dataobj[-grep("data", names(dataobj))]
     aux$sample.frame <- sample.frame
     assign("aux", aux, envir=env.who)

   }
###sanity checks
  if(is.character(userfile) && !file.exists(userfile))
    if((length(unlist(dataobj)) <= 0 || length(formula) <=0)
       && !file.exists(rerun))
    stop("Insufficient data")
 
### if no dataobject and formula use previous run
### which might be stored in the file yourcast.savetmp
  objread <- FALSE
  if((length(unlist(dataobj)) <= 0 || length(formula) <=0) && !is.character(userfile))
     {
       dataobj <- parse.dataobj(rerun,env.base)
       objread <- TRUE
       load(rerun)
       if(is.na(callst$formula)){
         callst$formula <- formula
       }
     }

### check to make sure that formula is either NULL or formula class
  if(class(formula) %in% c("formula")==FALSE){
    stop("Please make sure that formula is either a 'formula' object (or set to 'NULL' if you are using a formula from a saved file).")
  }


 
### done with data input check model  
  if(is.character(model))
    model <- trim.blanks(toupper(model))  
 ### for consistency with past code; taking the log with model poisson
  if(length(formula) > 0)
    {
      log.poiss <- check.depvar(formula)    
      assign("log.poiss", log.poiss, envir=env.who)
    
    }
   
### parameters args are stored correctly in env.who;
### bring them locally for easy of computation.
### This will only run if we have dataobj$data and formula
###################################################################
### RUN THE PREPROCESSING, which output are the lists variables
### Estimate of whoinsampy, whoutsampy, whocov, whoinsampx, whoutsampx
  prepross <- NULL ##preprocessing 

  if(length(dataobj) > 0 && !is.character(dataobj)) ##preprocessing
    {
    
      if(length(ind <-grep("lag", formula))> 0 )
        {
          lst <-  lagdataobj(formula, dataobj)
          datamat <- lst$datamat
          formula <- lst$ff
          dataobj$data <- datamat
 ###       print(formula)
        }
 
      prepross <- input.to.model(datamat=dataobj$data,ff=formula,
                                 all.pow=low.pow, sample.frame,
                                 index.code=dataobj$index.code,
                                 Gnames= dataobj$G.names,  
                                 standard=standardize, elim.collinear,
                                 tol=tol, solve.tol=solve.tol,
                                 log.poiss=log.poiss,verb=verbose)
 
      if(class(Hct.c.deriv)!= "try-error" && length(Hct.c.deriv) > 0)
        hc <- Hct.c.deriv
      else 
        hc <- dataobj$proximity
      
      assign.to.env(prepross, model, prox = hc, ewho=env.who);
### if savetmp =T then saves the prepross list return by input.to.model
### into a file in the working directory whose name is given by parameter rerun
      
      if(savetmp)
        save.yourcast.file(rerun, env.base)
            
      filexists <- try(file.exists(dataobj0),silent=T)
      if(class(filexists) != "try-error")
        if(filexists && !identical(dataobj0, rerun))
          dataobj <- parse.dataobj(rerun, env.base)
         
    }### done with preprocessing
##########################################################
    
### RUN MODEL REGRESSION ANALYSIS 
 
     
### sanity checks 
  if(model=="LS") model <- "OLS";
  if(length(formula) > 0) 
    depvar <- leftside.formula(formula)$numerator
  
  if((model=="MAP" || model=="EBAYES") && !is.na(Hct.sigma) ){
    messout("Setting Hct.sigma = NA...",verbose)
    Hct.sigma <- NA
    assign("who.Hct.sigma", Hct.sigma, envir=env.who)
  }
    
    
  if((model=="BAYES" || model=="MAP" || model=="EBAYES") && !is.na(Hct.sigma))
    who.Hct.c.deriv <- proximity.fill(dataobj,get("who.Hct.c.deriv", envir=env.who))
    
   
  if(model=="MAP" || model=="EBAYES"){
 
    autoset <- c(Ha.sigma["d1.a"], Ht.sigma["d1.t"], Hat.sigma["dtda"])
 
    autoset <- find.all.args()
    autoset <- unlist(autoset)
    assign("autoset", autoset, envir=env.who)
    summary.measures <- NULL
    if(length(autoset) > 0){
      summary.measures <- names(autoset)
      assign("summary.measures", summary.measures, envir=env.who)
      messout("Printing autoset....",verbose, obj=unlist(autoset))
     
    }
    
  }
   
  ### which model are you running ??
 
  indbayes <- grep(model, c("MAP", "BAYES", "EBAYES"))
### consistency check for zero.mean
### and reads a file with vector value for zero.mean
  if(length(indbayes) > 0)
    zero.mean <- zero.mean.age.profile(indbayes, env.who)
  
  smooth1 <- get.param.smooth(env.who)### c(Ha.sigma, Ht.sigma, Hat.sigma)
  if(length(Ha.sigma) >= 3 || length(Ht.sigma) >= 3 || length(Hat.sigma) >= 3)
    smooth <-  NULL
  else
    smooth <- smooth1
 
 
  dfirst  <- get("who.digit.first", envir=env.who)
  cdigits <- get("who.cntry.digits", envir=env.who)
  whoinsampy  <- get("whoinsampy", envir=env.who)
 
  count.cntry <- unique.default(unlist(sapply(names(whoinsampy), substr, dfirst, cdigits)))
  G <- rbind(NULL, smooth1)

  final.sigmas <- NULL
  if(length(indbayes) > 0 && length(smooth) < 3){
      Ha.sigma.vec <- try(get("Ha.sigma.vec", envir=env.who), silent=T)
      if(class(Ha.sigma.vec) == "try-error")
        Ha.sigma.vec <- try(get("Ha.sigma.vec", envir=env.base), silent=T)
      Ht.sigma.vec <- try(get("Ht.sigma.vec", envir=env.who), silent=T)
      if(class(Ht.sigma.vec) == "try-error")
        Ht.sigma.vec <- try(get("Ht.sigma.vec", envir=env.base), silent=T)
      Hat.sigma.vec <- try(get("Hat.sigma.vec", envir=env.who), silent=T)
      if(class(Hat.sigma.vec) == "try-error")
        Hat.sigma.vec <- try(get("Hat.sigma.vec", envir=env.base), silent=T)
      
      G <- build.prior.param(smooth,count.cntry, Ha.sigma.vec, Ht.sigma.vec, Hat.sigma.vec, sims,
                             stats=c(d1.a=d1.a.stat,d1.t=d1.t.stat,dtda=dt.da.stat),
                             autoset, model, env.who)
      final.sigmas <- G[1, ]
      messout("\n",verbose, obj=final.sigmas)
      
    }
  modeltopass <- model
   
  if(identical(model,"EBAYES"))
    modeltopass <- "MAP"
  
  ix <- grep(model, c("OLS", "LS", "POISSON", "LC", "BAYES", "EBAYES"))
  bool <- (length(Ha.sigma) >= 1 || length(Ht.sigma) >= 1 || length(Hat.sigma) >= 1)
 

  if(length(count.cntry) >  1 || length(ix) > 0 || bool){
   
    m <- try(switch(modeltopass, OLS=ols(env.base),LC=lc(env.base),POISSON=glm.poisson(env.base),
                       MAP=cxc(ebase=env.base),
                       BAYES=gibbs.sampler(), model="", NULL), silent=F);
    params <- NULL
  }else{
    ### runs one for one cntry and model MAP
   
    lst <- onecntry.bayes.empirical(G)
    m <- lst$m
    params <- lst$params
  
  }

  if(class(m) =="try-error")
    messout("Data is inconsistent.",verbose)

  if (is.null(m))
    stop(paste("Model ",model," is not available"));
       
### m contains obj= yhatin, yhatout, insampy, outsampy but listed with csid
### conversion.cntry.mat uses list.by.cntry to convert into country list matrices
### so it is a more useful format for the graphics output.
  
### EBAYES specific; calling Federico code
  if(identical(model, "EBAYES")){
 
    formula <- get("formula", envir=env.who)
   
    mebayes <- model.ebayes(autoset,formula,m,summary.measures,
                            depvar,graphics.file=NA,env.who)
   
    names(mebayes$summary) <- NULL ### to avoid problems with outputs
    messout("\n",verbose,obj=unlist(mebayes$summary))
    class(mebayes) <- "yourcast"
    return(mebayes)
  }
  
  model <- model.string()
  
  outy <- yhat.mat(ebase=env.base)
  coeff <- m$coeff
  std <- m$std

  outputlist <- list(yhat=outy,call=callst,userfile=userfilelst,
                     coeff=coeff, sigma=std, aux=aux)
   if(length(params) <= 0)
     params <- c(Ha.sigma=Ha.sigma, Ht.sigma=Ht.sigma, Hat.sigma=Hat.sigma)
  
   outputlist <- c(outputlist, list(params=params))
   
  if(get("save.output", envir=env.who)  && !is.na(model)){
    messout("Saving output file (option save.output = TRUE", verbose)
    depvar <- leftside.formula(formula)
    build.file.output(outputlist,model, depvar);
  }
     
### some formatting for yhat
  
  outy <- format.yhat(outy,year.digits)
  outputlist$yhat  <- outy
  if(length(final.sigmas) > 0)
    outputlist$final.sigmas <- final.sigmas
  class(outputlist) <- "yourcast"
     
  return(invisible(outputlist)) 
}

### DESCRIPTION: Takes the list data and ydigits a number
###              The elements of the list data are identified with csid
###              or cross country and age group combo,i.e. 245045 for USA
###              2450 and age = 45. The elements are matrices  
###              of numbers, whose rows names are the years of sample.frame
###              or observation years, i. e. 1985. It combined the the csid
###              tags with the year numbers, i.e. 2450451985, and tags rownames
###              of elements o data with the cross time sectinal idetifiers
###
### AUTHOR Elena villalon
###        IQSS Harvard Univ
###        evillalon@iq.harvard.edu
##########################################################

conversion.rownames <- function(data, ydigits){
  ind <- 1:length(data)
  names(ind) <- names(data)
  data <- lapply(ind, function(n)
                 {
                   mat <- as.matrix(data[[n]])
                   rname <- rownames(mat)
                   csid <- names(ind)[n]
                   nc <- unique.default(nchar(rname))
                   if(length(nc) > 1)
                     stop("Rows of data need to be named consistently")
                   
                   if (nc <= ydigits){
                     rname <- paste(csid, rname, sep="")
                     rownames(mat) <- rname
                   }
                   return(mat)})
}

### DESCRIPTION: helper function to yourcast; builds the call list output 
### It finds out if there is a dataobj in the argument inputs of yourcast
### If model.frame is True then it adds the dataobj to the call otherwise it adds
### just its name.
### INPUT: callst, a list with the arguments inputs to yourcast
####       model.frame a boolean indicating if the dataobj should be added to callst
###        env.who an environment
### OUTPUT: the callst with the element dataobj either added to it or the name of the object
### Elena Villalon
### evillalon@iq.harvard.edu
######################################################################################
call.dataobj1 <- function(callst, model.frame, env.who)
  {
    ebase <- try(get("env.base", envir=parent.frame()),silent=T)
    ix <- grep("dataobj", names(callst))
 ###   verbose <- get("verbose", envir=parent.frame())
    indc <- 1:length(callst)
    names(indc) <- names(callst) 
    if(length(ix) <= 0)
      return(callst)
    
    dbj <- callst[ix]
      
    if(!model.frame){
      callst[[ix]] <- dbj$dataobj
      dobj <- callst[[ix]]
    }else{
      callst[[ix]] <- dbj$dataobj$val
      dobj <- dbj$dataobj$dataname
    }
    if(length(dobj) > 0){  
      assign("dobj", dobj, envir=env.who)
      assign("dobj", dobj, envir=ebase)
    }
    return(callst)  
    }
### DESCRIPTION It takes the inputs of arguments of yourcast and a finds
###             the contents, assigning them to the environmnet that
###             yourcast is using to store the globals.
###
### AUTHOR Elena Villalon
###        IQSS Harvard Univ
###        evillalon@iq.harvard.edu
###################################################

update.args <- function(input.lst,userfile= NULL, ebase=NULL){
  
  ebase <- get("env.base", envir=parent.frame())
  ewho <- get("env.who", envir=ebase)
  sample.frame <- get("sample.frame", envir=ebase)
  whoyrest <- get("yrest", envir=ebase)
  env.base <- ebase;
  env.who <- ewho
  if(length(userfile) <= 0)
    {
      userfile <- get("userfile", envir=ewho)
      filename <- userfile
    }else 
    filename <- read.file.input(userfile, ebase);
  
###  print(filename) 
  if(length(filename) <= 0)
    return(list())

  source(file=filename, local=T)
     
  n.input <-  length(input.lst)
  for(i in 1:n.input )
    {
      val <- try(eval(as.symbol(input.lst[i])), silent=T)
      inlsti <- trim.blanks(input.lst[i])
 ###     print(inlsti)
      if(class(val) != "try-error" && inlsti=="sample.frame")
        { 
          yrest <- find.yrest(sample.frame)
          year.digits <- as.character(nchar(yrest))
       
          assign("yrest", as.numeric(yrest), envir=ebase)
          assign("whoyrest", as.numeric(whoyrest), envir=ewho)
          assign("year.digits", as.numeric(year.digits), envir=ebase)
          assign("who.year.digits", as.numeric(year.digits), envir=ewho)
          next; 
        }
      if(class(val) != "try-error" && inlsti =="dataobj")
        {
 ###         print(dataobj)
          dataobj <- parse.dataobj(val, ebase)
          assign("dataobj", dataobj, envir=ebase)
         
          if(length(dataobj$index.code) > 0 )
            {
              ix.code <- parse.index.code(icode=dataobj$index.code)
              assign.number.digits(ix.code, env.base)
            }
          if(length(dataobj$proximity) >0)
            {
              Hct.c.deriv <- dataobj$proximity
              assign("Hct.c.deriv", Hct.c.deriv, envir=ebase)
              assign("who.Hct.c.deriv", Hct.c.deriv, envir=ewho)
            }
        
          next;
        }
      if (class(val)!=  "try-error" && inlsti=="model")
        {
          assign("model", val, envir=ebase) 
          assign("whomodel", val, envir=ewho)
          next; 
        }
      if(class(val) != "try-error" && inlsti=="formula")
        {
          assign("formula", val, envir=ebase)
          assign("formula", val, envir=ewho) 
          next;
        }
      if(class(val) != "try-error" && inlsti=="savetmp"){
          assign(inlsti, val, envir=ebase)
       
          next; 
        }
      if (class(val)!=  "try-error" && !is.null(val))
        {
         
          ind <- grep(inlsti, ls(envir=ewho))
          ch <- ls(envir=ewho)[ind]
          assign(inlsti, val, envir=ebase)
        }
    }
  
 
}
 args.userfile <- function(userfile, ebase){  
  
   userfile <- read.file.input(userfile, ebase);
  
###  print(filename) 
  if(length(userfile) <= 0)
    return(list())

  ev <- environment()
   
  source(file=userfile, local=T)
  param <- ls(envir=ev)
 
  lst <- lapply(param, function(ch, ev) {
    val <- try(get(ch, envir=ev), silent=T)
    if(class(val)!= "try-error")
      return(val)
    else
        return(NULL)
  }, ev)
   names(lst) <- param
   ixenv <- grep("ev",names(lst))
   ixuser <- grep("userfile",names(lst))
   ix <- c(ixenv, ixuser)
   lst <- lst[-ix]
   ix1 <- grep("ebase", names(lst))
   ix2 <- grep("env.base", names(lst))
   ix <- c(ix1, ix2)
   if(length(ix) > 0)
     lst <- lst[-ix]
   return(lst)}
    
## DESCRIPTION A vector ix.code with 4 elements that are named with the
##             digits parameters names.  Using those names we
##             extract their values and assign them to the environmnets
##             It works in conjuction with parse.index.code
##
## AUTHOR: Elena Villalon
##         IQSS, Harvard Univ
##         evillalon@iq.harvard.edu
##############################################
   assign.number.digits <- function(ix.code, ebase)
    {
      ebase <- get("env.base", envir=parent.frame())
      ewho <- try(get("env.who", envir=ebase), silent=T)
      env.base <- ebase
      cdigits <- ix.code["cntry.digits"]
      try(assign("cntry.digits", cdigits, envir=ebase), silent=T)
      adigits <- ix.code["age.digits"]
      try(assign("age.digits", adigits, envir=ebase), silent=T)
      ydigits <- ix.code["year.digits"]
      try(assign("year.digits", ydigits, envir=ebase), silent=T)
      dfirst <- ix.code["digit.first"]
      try(assign("digit.first", dfirst, envir=ebase), silent=T)
      
      if(class(ewho) != "try-error"){
        assign("who.cntry.digits", as.numeric(cdigits), envir=ewho)
        assign("who.age.digits", as.numeric(adigits), envir=ewho)
        assign("who.year.digits",as.numeric(ydigits), envir=ewho)
        assign("who.digit.first", as.numeric(dfirst), envir=ewho)
      }
        
        
      return(list())
    }
### DESCRIPTION if dataobj is either an object or a character string
###             it calls two different functions
###
### AUTHOR Elena Villalon
###        IQSS Harvard Univ
###        evillalon@iq.harvard.edu
###
################################################################

parse.dataobj <- function(dataobj,ebase)
  {
    rerun <- get("rerun", envir=ebase)
    verbose <- get("verbose", envir=ebase)
    if(!is.character(dataobj)){
   
      dataobj <- assign.dataobj(dataobj, ebase)
      return(dataobj)
     }         
    if(is.character(dataobj) && (identical(dataobj,rerun) || !file.info(dataobj)$isdir)){
   
      dataobj <- check.for.files(dataobj, ebase)
      return(dataobj)
    }
    
    if(!identical(dataobj,rerun) && file.exists(dataobj) && file.info(dataobj)$isdir){
      messout(paste("Reading dataobj from directory ", dataobj,"...",sep=""),verbose)
      dcats <- dataobj
      dataobj <- yourprep(dpath= dcats,verbose=verbose)
      return(dataobj)
    }
    stop("Not a valid dataobj")
  }
### DESCRIPTION It takes different elements of dataobj and defined them
##              or set them to null if they are not defind. Note that
##              both data and index.code need to have some values
## AUTHOR Elena Villalon
##        evillalon@iq.harvard.edu
##
#####################################################################
assign.dataobj <- function(dataobj, ebase)
  {
    
    ewho <- try(get("env.who", envir=ebase), silent=T)
    data <- try(dataobj$data, silent=T)
    index.code <- try(dataobj$index.code, silent=T)
    G.names <- try(dataobj$G.names, silent=T)
    A.names <-  try(dataobj$A.names, silent=T)
    T.names <- try(dataobj$T.names, silent=T)
    Hct.c.deriv <- try(dataobj$proximity, silent=T)
    
    if(class(data)=="try-error")
      data <- NULL
    if(class(index.code)=="try-error")
      index.code <- NULL
    if(class(G.names)=="try-error")
      G.names <- NULL
    if(class(A.names)=="try-error")
      A.names <- NULL
    if(class(T.names)=="try-error")
      T.names <- NULL
    if(class(Hct.c.deriv)=="try-error")
      Hct.c.deriv <- NULL
###    print(dim(Hct.c.deriv))
    dataobj$data <- data
    dataobj$index.code <- index.code
    dataobj$G.names <- G.names
    dataobj$T.names <- T.names
    dataobj$A.names <- A.names
    dataobj$proximity <- Hct.c.deriv
    
    return(dataobj)
} 
            
### DESCRIPTION: It finds if userfile is in the working directory
###              or somewhere else
### AUTHOR Elena Villalon
################################################

read.file.input <- function(userfile=NULL, ebase=NULL,
                            datapath=getwd()){
  
  if(length(userfile) <= 0)
    userfile <- get("userfile", envir=ebase)
  
  chu0 <- grep("/", userfile) ### for unix

  chw <- NULL
  if(length(grep(":",userfile)) > 0)
    
  chw <- grep("([\\]+)", userfile) ### for windows
  chu <- c(chu0, chw)
  
  look <- (length(chu) <= 0 )
  if(length(userfile) > 0 && look && length(chu0) > 0){
    dir <- datapath
    filename <- paste(dir,"/",userfile, sep="");
    
  }else if(length(userfile) > 0 && look && length(chw) > 0){
    dir <- datapath
    filename <- paste(dir,":\\",userfile, sep="");
 }else if(length(userfile) > 0) 
      filename <- userfile 
  else
    filename <- NULL
  return(filename) 
}

### DESCRIPTION Helper function to yourcast that is used after
###             the preprocessing part of yourcast, input.to.model
###             to put in the environmnets the list elements return
###             with input.to.model and assign to the argument lst.
###             mod stands for model string, and ewho is the environmnet.
### OUTPUT a bunch of assigmnets of variables and parameters needed to
###        run the different regression models.
###
### AUTHOR Elena Villalon
###        IQSS, Harvard Univ
###        evillalon@iq.harvard.edu
#########################################################
assign.to.env <- function(lst, mod, prox, ewho){
 
  whocov <- lst$whocov
  assign("whocov", whocov, envir=ewho)
  
  whoinsampx <- lst$whoinsampx
  assign("whoinsampx", whoinsampx, envir=ewho)
                
  whoutsampx <- lst$whoutsampx
  assign("whoutsampx", whoutsampx, envir=ewho)
  
  whoinsampy <- lst$whoinsampy                 
  whoutsampy <- lst$whoutsampy
  
  assign("whoinsampy", whoinsampy, envir=ewho)
  assign("whoutsampy", whoutsampy, envir=ewho)
  whopopul <- lst$whopopul
  assign("whopopul", whopopul, envir=ewho)
                
  whopopulos <- lst$whopopulos
  assign("whopopulos", whopopulos, envir=ewho)
                
  cov.lst <- lst$cov.lst
  assign("cov.lst", cov.lst, envir=ewho)
                
  age.vec <- lst$age.vec
  assign("age.vec", age.vec, envir=ewho)
                
  cntry.vec <- lst$cntry.vec
  assign("cntry.vec", cntry.vec, envir=ewho)
  
  geo.index <- lst$geo.index
  assign("geo.index", geo.index, envir=ewho)

  c.vec <- as.character(cntry.vec)
  cntry.names.lst <- sapply(c.vec, function(ch,geo.index) {geo.index[ch]}, geo.index)
  names(cntry.names.lst) <- c.vec
  assign("cntry.names.lst", cntry.names.lst, envir=ewho)
  
  log.poiss <- lst$log.poiss
 
  assign("log.poiss", log.poiss, envir=ewho)
  
  formula <- lst$formula
  assign("formula", formula, envir=ewho)
 
  assign("who.Hct.c.deriv", prox, envir=ewho)
  who.Hct.c.deriv <- get("who.Hct.c.deriv", envir=ewho)  
   
}
### DESCRIPTION:helper function to yourcast to update values of smoothness parameters

get.param.smooth <- function(env.who)
  {
    smooth <- c(Ha.sigma=get("who.Ha.sigma", envir=env.who),
                Ht.sigma=get("who.Ht.sigma", envir=env.who),
                Hat.sigma=get("who.Hat.sigma", envir=env.who))
    return(smooth)
  }
### DESCRIPTION Given a formula for the simulation such as
###             log(depvar) ~ tobacco.txt + time + gdp.txt^3 or
###             depvar ~ tobacco.txt + time + gdp.txt^3 or
###             It splits it into the r-h-s and l-h-s
###             Finds out if the left-hand-side is the log of depvar
###             and assign values T or F to log.poiss accordingly
###
### AUTHOR Elena Villalon
###        IQSS Harvard Univ
###        evillalon@iq.harvard.edu
##################################################################
check.depvar <- function(formula)
  {
    fc  <- as.character(formula)
    ls  <- fc[2]
    vec <- strsplit(ls, NULL)[[1]]
    
    log.poiss <- F
    ix <- unlist(sapply(c("l", "o", "g"), match, vec))
    if (all(ix == 1:3))
      log.poiss <- T
    
    ch <- sub("log\\(([a-z A-Z / \\. *]+)+\\)","\\1", ls) ### from "log(gdp)" to "gdp"
 ###   ch <-  sub("([a-z A-Z \\.]+)\\(([a-z A-Z / * \\+ \\. 0-9]+)+\\)","\\2",ls)
 
    log.poiss0 <- ifelse(nchar(ls) > nchar(ch), T, F)
    if (log.poiss0 && log.poiss)
      return(log.poiss)
    else{
       v <- FALSE
      return(v)}

  }
    
### DESCRIPTION It takes a filename and an environment
###             and load the filename locally. It checks what is
###             in the file.  The vector input.lst contains the lists
###             of the preprocessing part of yourcast.  If the information is in
###             obj or, unless, enough information for the basic lists needed
###             by the programs we just pass that to the environments
###             and return obj (the input).  If we do not have all information but we
###             provide for whoinsampy, whoutsampy, whoinsampx, whoutsampx
###             whopopul, whopopulos, then we reconstruct with them the rest
###             of thew lists needed by yourcast.
###             If noneof those lists are provided we do some error checking
###             to make sure that dataobj is in the file and that the names of
###             the list elements correspond to "data" and "index.code" as a minimum.
###             After  all error checking is completed then we assign dataobj
###             and index.code to the environments and returns the object dataobj.
###
### INPUT       the obj which is a filename with inputs list and the environmnet
###
### OUTPUT      If obj file contains some of the preprocessing stuff that
###             yourcast needs, returns the file name itself obj after assigning
###             all lists in obj to the environments.
###             If obj contains the dataobj list input argument of yourcast, then\
###             it parses it and assign "data", and "index.code" to environments,
###             and return the object list dataobj.
###
### AUTHOR      Elena Villalon
###             IQSS, Harvard Univ
###             evillalon@iq.harvard.edu
##############################################################################


check.for.files <- function(obj, ebase)
  {
    ewho <- get("env.who", envir=ebase)
    verbose <- get("verbose", envir=ebase)
    ev <- environment()
    load(obj, envir=ev)
###    print(ls(envir=ev))
   
    input.lst <- c("whoinsampx", "whoutsampx", "whoinsampy",
                   "whoutsampy", "whopopul", "whopopulos","whocov",
                   "cov.lst","age.vec", "cntry.vec", "geo.index",
                   "log.poiss", "formula", "who.Hct.c.deriv")
 
   who.Hct.c.deriv <- get("who.Hct.c.deriv", envir=ewho)
    ln <- length(input.lst)
    for(n in 1:ln)
      {     
        ch <- input.lst[n]
        val <- try(eval(as.symbol(ch)), silent=T)
        
        if(class(val)!="try-error"  && n <=6 )
          assign(ch, val, envir=ewho)
        else if (n <= 6)
          break;
        
        if(class(val)=="try-error"  && n ==7){
          whocov <- build.whocov(ewho)
          assign("whocov", whocov, envir=ewho)
         
          next;
        }else if (n==7){
          assign(ch, val, envir=ewho)
          next;
       }
        
        if (class(val)=="try-error"  && n ==8){
          whoinsampx <- get("whoinsampx", envir=ewho)
          cov.lst <- list.covariates(whoinsampx)
          assign("cov.lst", cov.lst, envir=ewho)
          next;
        }else if(n==8){
          assign(ch, val, envir=ewho)
          next;
        }
        
        if(class(val)=="try-error" && n == 9){
          age.vec <- build.age.cntry.vec(ebase)$age.vec
          assign("age.vec", age.vec, envir=ewho)
          next;
        }else if(n==9){
          assign(ch, val, envir=ewho)
          next; 
        }
        
        if (class(val)=="try-error" && n == 10){
          cntry.vec <- build.age.cntry.vec(ebase)$cntry.vec
          assign("cntry.vec", cntry.vec, envir=ewho)
          next;
        }else if(n==10){
          assign(ch, val, envir=ewho)
          cntry.vec <- val
          next;
       }
     
        if (class(val)=="try-error" && n == 11){
     
          val <- NULL
          assign(ch, val, envir=ewho)
          next;
       }else if(n==11 ){
     
          c.vec <- as.character(cntry.vec)
          
          assign(ch, val, envir=ewho)
          geo.index <- val
          if(length(geo.index) > 0){
            cntry.names.lst <- sapply(c.vec,
                                      function(ch,geo.index) {geo.index[ch]}, geo.index)
          
          }else
            cntry.names.lst <- rep(" ", length(c.vec))
          
          names(cntry.names.lst) <- NULL
          names(cntry.names.lst) <- as.character(c.vec)
          assign("cntry.names.lst", cntry.names.lst, envir=ewho)
          next; 
        }
        if(class(val) != "try-error")       
          assign(ch, val, envir=ewho)
        else
          assign(ch, NULL, envir=ewho)
       
      }
   
    if(n >= ln){
      messout(paste("Reading parameters from file... ", obj, sep=""),verbose)
###                "\nIgnoring yourcast arguments that are not model specific and ", 
###                "are not in the userfile", "\n"))
###      print(ls(envir=ev))
###      print(dim(who.Hct.c.deriv))
     
      digits <- age.cntry.digits(age.vec, cntry.vec, ebase)
     
      assign("age.digits", as.numeric(digits$age.digits), envir=ebase)
      assign("cntry.digits", as.numeric(digits$cntry.digits), envir=ebase)
      assign("digit.first", as.numeric(digits$digit.first), envir=ebase)
      assign("who.age.digits", as.numeric(digits$age.digits), envir=ewho)
      assign("who.cntry.digits", as.numeric(digits$cntry.digits), envir=ewho)
      assign("who.digit.first", as.numeric(digits$digit.first), envir=ewho)
      assign("who.Hct.c.deriv", who.Hct.c.deriv, envir=ewho)
      assign("Hct.c.deriv", who.Hct.c.deriv, envir=ebase)
      return(obj)
    }
    
    if(class(try(dataobj$proximity, silent=T))!= "try-error"){
        assign("who.Hct.c.deriv", dataobj$proximity, envir=ewho)
        assign("who.Hct.c.deriv", dataobj$proximity, envir=ebase)
      }
       
    
    dataobj <- try(eval(as.symbol("dataobj")), silent=T)
    
    if(class(dataobj) == "try-error")
      stop("Insufficient data")
### checks if data and index.code are in dataobj list         
    nmobj <- names(dataobj)
    if(length(nmobj) < 2)
      stop("Insufficient data")
    
    vec <- c("data", "index.code")
    logic <- all(is.element(vec, nmobj))
    if(!logic)
       stop("Insufficient data")
   
    dataobj <- assign.dataobj(dataobj, ebase)
    index.code <- dataobj$index.code
###parsing the index.code for the relevant digits
    ix.code <- parse.index.code(icode=dataobj$index.code, yrest=NULL)
    assign.number.digits(ix.code, ebase)    
            
    return(dataobj)     
  }

### DESCRIPTION finds the number of digist for age, country and years
###             using vectors age.vec and age.vec, and assign them to environmnets
###             It also creates a list with them to return with function call
### Author Elena Villalon
###        IQSS, Harvard Univ
###        evillalon@iq.havard.edu
##################################################################
age.cntry.digits <- function(age.vec, cntry.vec, ebase)
  {
    ewho <- get("env.who", envir=ebase)
    verbose <- get("verbose", envir=ebase)
    age.digits <- max(as.numeric(unlist(sapply(as.character(age.vec), nchar))))
    age.digits <- as.numeric(age.digits)
    assign("who.age.digits", as.numeric(age.digits), envir=ewho)
    cntry.digits <- max(as.numeric(unlist(sapply(as.character(cntry.vec), nchar))))
    cntry.digits <- as.numeric(cntry.digits)
    assign("who.cntry.digits", as.numeric(cntry.digits), envir=ewho)
    nm <- names(get("whoinsampx", envir=ewho))[1]
   
    digit.first <- nchar(nm) - age.digits - cntry.digits
    assign("who.digit.first", digit.first, envir=ewho)
    if( digit.first < 0)
      digit.first <- 0
    assign("who.age.digits", age.digits, envir=ewho)
    lst <- list(age.digits=age.digits, cntry.digits=cntry.digits, digit.first=digit.first)
    return(lst)
  }
### DESCRIPTION from sample.frame takes the year that divides
###             the insample and outsample perios called yrest and return it
### AUTHOR Elena Villalon
###        IQSS, Harvard Univ
###        evillalon@iq.harvard.edu
##########################################################

 find.yrest <- function(sample.frame)
    {
      if(length(sample.frame) == 1)
        yrest <- sample.frame
      else if(length(sample.frame)==4)
        yrest <- ifelse(try(sample.frame[2], silent=T) != "try-error", sample.frame[2], NULL)
      else
        stop("Wrong input data for sample.frame")
      return(yrest)
    }

### DESCRIPTION takes as argument the name rerun of the temp file to save
###             preprocessing lists and parameters of yourcast in the file rerun
###             Elements toi be saved are the insample and outsample periods
###             for dependent variable and covariates
###
### AUTHOR Elena Villalon
###        IQSS, Harvard Univ
###        evillalon@iq.harvard.edu
##################################################################
save.yourcast.file <- function(rerun, ebase = env.base){
  ebase <- get("env.base", envir=parent.frame());
  verbose <- get("verbose", envir=ebase)
  env.base <- ebase
 
### get the environment where data is located
  ewho <- get("env.who", envir=ebase)
 
### create the blank file to store data for directory: whooutpath
 
  filename  <- paste(getwd(),"/", rerun,sep="")
 
  messout(paste("Saving preprocessing in the file...",filename,sep=""),verbose)
  filename <- rerun
  if (file.exists(filename)){
    file.remove(filename)}
  file.create(filename)
### name for present enviroment:
  esave <- environment()
    
### what to store in filename:   
    what <- c("dobj",  "whocov","whoinsampx", "whoutsampx", "whoinsampy","whoutsampy", 
              "whopopul","whopopulos", "cov.lst", "age.vec", "cntry.vec","geo.index", 
              "cntry.names.lst", "log.poiss", "formula", "who.Hct.c.deriv", "aux") 

   
### write those values in present environment esave

    for (i in 1:length(what)){
     
        
      wi <- try(get(what[i], envir=ewho, inherits=T), silent=T)
      if(class(wi) == "try-error")
         wi <- try(get(what[i], envir=ebase, inherits=T), silent=T)
     
        assign(what[i],wi, envir=esave)
     
         
    }
  
   
 
### better to find what to exclude
### name of environments that start with "e":    
    ind.ex <- match(ls(envir=esave,pattern="^e"), ls(envir=esave))
    ind.ex <- c(ind.ex,  match("filename", ls(envir=esave)))
    ind.ex <- c(ind.ex, match("what", ls(envir=esave)))
    ind.ex <- c(ind.ex, match("ch", ls(envir=esave)))
    ind.ex <- c(ind.ex, match("^esave$", ls(envir=esave)))
    ind.ex <- c(ind.ex, match("^ebase", ls(envir=esave)))
    ind.ex <- c(ind.ex, match("^ev", ls(envir=esave)))
    ind.ex <- c(ind.ex, match("^ewho", ls(envir=esave)))
  
### indeces that start with "ind":
    ind.ex <- c(ind.ex, match(ls(envir=esave, pattern="^ind"), ls(envir=esave)))
    ind.ex <- na.omit(ind.ex)
### now save everything in envir= esave but variables in [indx.ex]:
###   save(list=ls(envir=esave)[-ind.ex], file=filename,compress=T)
### Safer to save what you need
      
     save("whocov","whoinsampx", "whoutsampx", "whoinsampy","whoutsampy", 
          "whopopul","whopopulos", "cov.lst", "age.vec", "who.Hct.c.deriv", 
          "cntry.vec","geo.index", "cntry.names.lst", "log.poiss","formula","dobj","aux",      
          file=filename,compress=T)
  
###load("yourcast.savetmp")
###> dobj
###[1] "datam"
###> dobj <- eval(as.symbol("datam"))
###> names(dobj)
###[1] "data"       "index.code" "proximity"  "G.names"   
###> dobj$index.code
###[1] "ggggaatttt"
  
  messout(paste("The size of ", filename, " is = ", file.info(filename)$size, sep=""),verbose)
}

#### DESCRIPTION Takes whoinsampx and whoutsampx and joins
###              them together to form whocov that contains
###              the insample and outsample periods for the covariates
###
### AUTHOR Elena Villalon
###        IQSS Harvard Univ
###        evillalon@iq.harvard.edu
###################################################################
build.whocov <- function(ewho=NULL, whoinsampx=NULL, whoutsampx=NULL){
  if(length(ewho)>0){ 
    whoinsampx <- get("whoinsampx", envir=ewho)
    whoutsampx <- get("whoutsampx", envir=ewho)
  }
  ind <- 1:length(whoinsampx)
  names(ind) <- names(whoinsampx)
  whocov <- lapply(ind, function(n){
    inmat <- whoinsampx[[n]]
    outmat <- whoutsampx[[n]]
    mat <- rbind(inmat, outmat)
    return(mat)})
  return(whocov)
}

#### DESCRIPTION Takes number of digits for cntry, age, and the insampy list
###              and reconstruct the vectors of age and cntry for sample data 
###
### AUTHOR Elena Villalon
###        IQSS Harvard Univ
###        evillalon@iq.harvard.edu
###################################################################
build.age.cntry.vec <- function(ebase){
  cdigits <- get("cntry.digits", envir=ebase)
  adigits <- get("age.digits", envir=ebase)
  fdigit <- get("digit.first", envir=ebase)
  st <- cdigits + 1 + adigits
  whoinsampx <- get("whoinsampx", envir=get("ewho", envir=ebase))
  tag <- names(whoinsampx)
  age.vec  <- unique.default(unlist(sapply(tag, substring, fdigit + cdigits+1, st)))
  cntry.vec <- unique.default(unlist(sapply(tag,substring, fdigit + 1,cdigits)))
  lst <- list(age.vec = age.vec, cntry.vec = cntry.vec)
  return(lst)
}

### DESCRIPTION Takes two arguments from a fuction call, such as
### yourcast(par1=1, par2=2,...,parn=n). The first argument formlst
### are the formals of the function call: formals(yourcast); and
### the second argument callmath is = match.call(yourcast).
### The formals arguments gives a list of the defaults parameters
### of the function call.  The match.call returns, in language, the
### actual function call for the run with the new values, if any,
### assign t the parameters of the function, yourcast(model="ols").
### We find the new values assign to the arguments of the function call
### and substitute the new values in the formals list.  Those that are
### redefined remains with their defaults values. Return the modified
### formals list with new arguments for those parameters included in
### the function call for the run.
###
### AUTHOR Elena Villalon
###        IQSS, Harvard Univ
###        evillalon@iq.harvard.edu
####################################################################


 call.yourcast <- function(formlst, callmatch, svtmp, rerun)
    {
      ebase <- get("env.base", envir=parent.frame())
      formu <- get("formula", envir=ebase)
      dbj <- get("dataobj", envir=ebase)
      hct <- get("Hct.c.deriv", envir=ebase)
      ev <- environment()
   
      what <- c("dobj", "formula", "who.Hct.c.deriv") 
      formlst <- lapply(formlst, function(mat) {
        if(length(mat) <= 0)
          return(NA)
        else
          return(mat)
      })
                 
  
     
      cm  <- as.character(callmatch)
      nmf <- names(formlst)
      nmc <- names(callmatch)
      ind <- 1:length(formlst)
     
      names(ind) <- names(formlst)
     
      modform <- lapply(ind, function(n){
      
        val <- formlst[[n]]
       
        nm  <- nmf[n]
        ix  <- grep(nm, nmc)
      
       if(!identical(nm, "dataobj")){
          
         if(length(ix) <= 0)
           return(val)
         else if (identical(nm, "formula")){
           return(formu)
         }else
         return(as.character(callmatch[ix]))
       }
     
     
        if(length(ix) <= 0)
          return(NA)
      
        if(!svtmp )
          return(as.character(callmatch[ix]))
        else{
          val <- eval(as.symbol(callmatch[[ix]]))
          lst <- list(val=val, dataname=as.character(callmatch[ix]))
          return(lst)
        }
      })
                        
      
      modform <- lapply(modform,function(ch){
        if(!is.character(ch))
          return(ch);
       
        ix <- grep("[0-9]+", ch)
        ix1 <- grep("[\\( , a-z A-Z \\)]+", ch)
        if(length(ix) > 0 && length(ix1) <=0){
        
           return(as.numeric(ch))
        }
### converting the string "c(1, 2, 3)" to regular vector
        if(length(ix) >0){
          ##getting rid of c and parenthesis
          noandcom <- gsub("[c \\(  \\)]+", "", ch)
          chno <- strsplit(noandcom,",")[[1]]
          if(all(is.numeric(chno)))
            return(as.numeric(chno))
          else
            return(chno)
        }
        return(ch);
      } )
    }
               
find.dataobject <- function(envir= NULL,
                            loc="http://gking.harvard.edu/yourcast/dataobj.dat"){
  ev <- envir
  if(length(envir) <= 0)
    ev <- parent.frame()
  
   load(url(loc), envir=ev)
}

### DESCRIPTION it takes the list yhat and rename the rownames with the
### year digits ony, say we have "--2450451989" then convert it to "1989"
### INPUT: list yhat with elements csid and ecah is a matrix of two columns
###        whose names are the csid cross time series identifiers
###        Convert it to time series, y substracting the cross-sectional id
### OUTPUT: the same list with the new rownames that only have years
### Elena Villalon
### evillalon@iq.harvard.edu    

    format.yhat <- function(y, ydg)
  {
    
    y <- lapply(y, function(mat, ydg) {
      vec <- rownames(mat)
      
   ##take away leading and trailing "-"   
      vec <- sapply(vec, function(mm) {
        mm <- sub("(^[-]*)", "", mm)
        mm <- sub("([-]*$)", "", mm)
        return(mm)})
   ##obtain the years digits   
      vec <- sapply(vec, function(x){
        nc <- nchar(x)
        ly <- nc - ydg + 1
        sb <- x 
        if( ly > 0)
          sb <- substring(x, ly)
        return(sb)})
   ##rename the rownames   
      vec <- unlist(vec)
      rownames(mat) <- vec
      return(mat)}, ydg)
     return(y)
  }
### DESCRIPTION 
###        It checks if any of the columns of mat contains only na's
###        If it is the case then removes the entry
###
### INPUT: n integer; dat list with csid matrices
###        mat = dat[[n]] 
###
### OUTPUT either matrix mat or NULL depending on if mat has any columns all NA's
###
### AUTHOR Elena Villalon
###        evillalon@iq.harvard.edu
###        March 20th, 2006
###
##################################################################################

elimallna.csid <- function(n, dat,verbose)
  {
    mat <- dat[[n]]
    nm <- names(dat)[n]
    nc0 <- ncol(mat)
    ix <- as.list(1:ncol(mat))
    names(ix) <- colnames(mat)
    ix <- sapply(ix, function(n){
      cc <- na.omit(mat[,n])
      
      if(length(cc) >0)
        return(n)
      else
        return(NULL)
      
    })
    ix <- unlist(ix)
    if(length(ix) < nc0){
         messout(paste("Elimination of unit ", nm, "because of NA's",sep=""),verbose) 
      return(list())
    }else
     return(mat)
    
  }

### DESCRIPTION: helper function to yourcast;
### it finds out if proximity matrix is in the data set
### and assigned it to the environmnet

proximity.fill <- function(dataobj, who.Hct.c.deriv)
  {
    ebase <- get("env.base", envir=parent.frame())
    ewho <- get("env.who", envir=ebase)
   
    
    if((length(who.Hct.c.deriv) <= 0 || class(who.Hct.c.deriv) == "try-error")
       && length(dataobj) >0){
      who.Hct.c.deriv <- dataobj$proximity
           
      assign("who.Hct.c.deriv", who.Hct.c.deriv, envir=ewho)
      assign("Hct.c.deriv", who.Hct.c.deriv, envir=ebase)
    }
         
    if(length(who.Hct.c.deriv) <= 0 ){       
      stop("You need to provide for proximity matrix")}
    return(who.Hct.c.deriv)
  }

### DESCRIPTION helper function to yourcast to make sure we have
###            all the arguments to run bayesian type models
###            Creates auset with the target for ebayes
### OUTPUT:autoset, vector with values for d1.a,d1.t, dtda,SD
### INPUT: environments to extract parameters
### USES: ebayesparam.vec 
### AUTHOR Elena Villalon
###        evillalon@iq.harvard.edu
###        May 22, 2006
###########################################################
find.all.args <- function(autoset=NULL)
  {
    env.base <- get("env.base", envir=parent.frame())
    env.who <- get("env.who", envir=env.base)
   
    args.ebayes <- c("who.Ha.sigma","who.Ht.sigma","who.Hat.sigma",
                     "graphics.file", "output.file")
    whomodel <- get("model", envir=env.who)
    
    
    for(ch in args.ebayes){
      nmch <- ch
     
      ch <- try(get(ch, envir=env.who), silent=T)
  
      if(class(ch)=="try-error")
        ch <- try(get(ch, envir=env.base), silent=T)
      if(class(ch)=="try-error")
        stop(message("missing variable ", ch))
    
      if(length(ch) >= 3){
        autoset <- ebayesparam.vec(ch,nmch, autoset, env.base)
      }
    }
   
   
    nmauto <- names(autoset)
    if(length(grep("SD", nmauto)) <= 0 && length(nmauto) >0){
      assign("SD", NA, envir=env.who)
      assign("SD", NA, envir=env.base)
    }
      
   return(autoset)     
  }
### DESCRIPTION helper function to find.all.args to build the target values
###             if length(who.Hx.sigma) >= 3;  
###             d1.a, d1.t, dtda, SD and autoset=c(d1.a, d1.t, dtda, SD); 
###             and Ha.sigma.vec, Ht.sigma.vec, Hat.sigma.vec
###             Assigns NULL values to Ha.sigma, Ht.sigma, Hat.sigma
##              
### OUTPUT:autoset, vector with values for d1.a,d1.t, dtda,SD
###        assignmnets to env.who, for Ha.(Ht.)(Hat.)sigma.vec,
###        and set Ha.(Ht.)(Hat.)sigma=NULL asiigning to env
###
### INPUT: environments to extract parameters
###        who.Hx.sigma with Hx=Ha, Ht, Hat; a vector numeric values
###        who.Ha.sigma=c(0.1, 1.5, 6, 0.4, NA) means
###        Ha.sigma.vec= seq(0.1, 1.5, lenght=6); d1.a=0.4, SD=NA
###        
### AUTHOR Elena Villalon
###        evillalon@iq.harvard.edu
###        May 22, 2006
###########################################################
ebayesparam.vec <- function(who.Hx.sigma,nmhx, autoset,env.base)
  {
    env.base <- get("env.base",envir=parent.frame())
    env.who <- get("env.who", envir=env.base)

    if(length(who.Hx.sigma)<= 2){
      stop(paste("Provide either a scalar or a vector of  3 or more elements for ",nmhx))
      
      return(autoset)
    }
    
 ###   cat("Print diagnosis...\n")
 ###   print(who.Hx.sigma)
 ###   print(nmhx)
 ###   print(autoset)
 ###   cat("diagnostic...\n")
 ###   nmHxsigma <- names(who.Hx.sigma)
 ###   ln   <- length(who.Hx.sigma)
 ###   lnnm <- length(nmHxsigma)
     
    if(identical("who.Ha.sigma", nmhx)){
      Hx.vec <- "Ha.sigma.vec"
      d1.x <- "d1.a"
      d1x <- d1.x
      Hx <- "Ha.sigma"
      
        
    }
    if(identical("who.Ht.sigma", nmhx)){
      Hx.vec <- "Ht.sigma.vec"
      d1.x <- "d1.t"
      d1x <- d1.x
      Hx <- "Ht.sigma"
         
    }
    if(identical("who.Hat.sigma", nmhx)){
      Hx.vec <- "Hat.sigma.vec"
      d1.x <- "dtda"
      d1x <- d1.x
      Hx <- "Hat.sigma"
      
    }
  ###  print(Hx.vec)
  ###  print(d1.x)
  ###  print(d1x)
  ###  print(Hx)
    
    tonull <- who.Hx.sigma
    hx <- unlist(who.Hx.sigma[1:3])
  
    Hx.sigma.vec <- seq(hx[1], hx[2], length=hx[3])
    assign(Hx.vec,Hx.sigma.vec, envir=env.who)
    assign(Hx.vec,Hx.sigma.vec, envir=env.base)
 
    if(length(who.Hx.sigma) <= 3){
      assign(nmhx, NULL, envir=env.who)
      assign(Hx, NULL, envir=env.base)
      return(autoset)
    }
    
    hxvec <- names(who.Hx.sigma)    

    if(length(hxvec) <= 0) ### length(who.Hx.sigma) >= 4
      {
        
        d1.x <- who.Hx.sigma[4]
        nmauto <- names(autoset)
        autoset <- c(autoset, d1.x)
        names(autoset) <- c(nmauto,d1x)
        if(length(who.Hx.sigma) == 4){
          assign(nmhx, NULL, envir=env.who)
          assign(Hx, NULL, envir=env.base)
          return(autoset)
        }

        SD <- who.Hx.sigma[5]
        nmauto <- names(autoset)
        autoset <- c(autoset, SD=SD)
        names(autoset) <- c(nmauto, "SD")
        if(length(who.Hx.sigma) ==5){
          assign(nmhx, NULL, envir=env.who)
          assign(Hx, NULL, envir=env.base)
          return(autoset)
        }

        sims <- who.Hx.sigma[6]
        assign("sims", sims, envir=env.base)
        assign("N", sims, envir=env.who)
        assign(nmhx, NULL, envir=env.who)
        assign(Hx, NULL, envir=env.base)
        return(autoset)
           
      }

      
    if(length(hxvec) > 0){
   
      hxixsg <- grep(Hx, hxvec)
      if(length(hxixsg) > 0){
        Hx.sigma.vec <- who.Hx.sigma[hxixsg]
        assign(Hx.vec,Hx.sigma.vec, envir=env.who)
        assign(Hx.vec,Hx.sigma.vec, envir=env.base)
      }
      
      if(length(grep(d1.x, hxvec) ) > 0){
        d1.x <- who.Hx.sigma[d1.x]
        nmauto <- names(autoset)
        autoset <- c(autoset, d1.x)
        names(autoset) <- c(nmauto,d1x) 
      }
       if(length(grep("SD", hxvec)) > 0){
        SD <- who.Hx.sigma["SD"]
        nmauto <- names(autoset)
        autoset <- c(autoset, SD=SD)
        names(autoset) <- c(nmauto, "SD")
      }
    
      if(length(grep("sims", hxvec)) > 0){
        
        sims <- who.Hx.sigma["sims"]
        assign("sims", sims, envir=env.base)
        assign("N", sims, envir=env.who)
       
      }
    }
    
    assign(nmhx, NULL, envir=env.who)
    assign(Hx, NULL, envir=env.base)
  
    return(autoset)
  }    
### DESCRIPTION:For bayesian and map models, it finds if zero.mean
###             is a vector of the same length as age.vec.  Also, it 
###             checks if the names of zero.mean are according to age.groups
######################################################################
zero.mean.age.profile <- function(indbayes, env.who){
    if(length(indbayes) <= 0)
      return(invisible(list()))
    
    who.zero.mean <- get("who.zero.mean", envir=env.who) 
    age.vec <- sort(get("age.vec", envir=env.who))
    if(!is.numeric(who.zero.mean))
      return(who.zero.mean)
    ln <- abs(length(who.zero.mean) - length(age.vec))
    if(ln > 1.e-10)
      stop("Number of elements in zero.mean are different from age groups")
    ages <- sort(as.numeric(names(who.zero.mean)))
    ages <- as.numeric(ages)
    age.vec <- as.numeric(age.vec)
    if(!all(ages == age.vec) && length(ages) > 0)
      stop("Provide zero.mean values for all age groups")
   
      names(who.zero.mean) <- age.vec
      assign("who.zero.mean", who.zero.mean, envir=env.who)
    
    return(who.zero.mean)
  
}

messout <- function(str=NULL, verb=TRUE, obj=NULL){
  if(verb && length(str))
  message(str);
  if(length(obj) >0 && verb)
    print(obj)
}
 
 checkdataobj <- function(dat=NULL, smpfrm=c(1950,2000,2001,2030)){
    ln <- length(smpfrm)
    noy <- ln +1
    smpvec <- seq(smpfrm[1],smpfrm[ln], by=1)
    for(n in 1:length(dat)){
      mat <- dat[[n]]
      nm <- rownames(mat)
      if(length(nm) <= 0){
        message("Can't check the data portion of dataobj")
        return(NULL)
      }
        
      nmnum <- as.numeric(nm)
      minentry <- min(nmnum)
      mn  <- paste("^",minentry,"$",sep="")
      tmp <- sapply(smpvec,function(x) {paste("^",x,"$",sep="")})
      tmp <- unlist(tmp)
      ind <- match(mn, tmp)
      if(is.na(ind)){
        message("Your data object is missing some years for csid ", names(dat)[n])
        return(NULL)
      }
      smp <- seq(smpvec[ind], smpvec[ln],by=1)
      if(all(smp %in% nmnum)) next;
     
        message("Your data object is missing some years for csid ", names(dat)[n])
    }
         
      
  }
