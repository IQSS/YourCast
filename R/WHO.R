##
## USED GLOBALS: envir=env.base, our base environment, which for working purposes we
##               designated it to be the .GlobalEnv.  It is now the environment for function main()
##               that wraps up all programs and functions.
##               Also, envir=env.who (an environment), which we create to stores all
##               static variables to be shared by all functions and programs.  
##               We assign env.who to the env.base as, assign("env.who", envir=env.base)
##               if you want to see what is in env.who:  >  ls(envir=get("env.who", envir=env.base)); 
##
## DESCRIPTION:  The programs static variables, are given default values with
##               the function WHO(); variables are stored in the environment env.who,
##               which belongs to env.base. Any variables belonging to env.who,
##               may be retreived from env.base as:
##               get("whodisease",envir=get("env.who", envir=env.base)); 
##               and you may change its value as
#               assign("whodisease",newdisease, envir=get("env.who", envir=env.base))
##
## FORMAT: WHO(ebase=env.base) 
##
## VALUE: all variables used in the programs with their default values in WHO()
##
## WRITTEN BY: Elena Villalon & Federico Girosi
##             evillalon@latte.harvard.edu,fgirosi@latte.harvard.edu
##             CBRSS, Harvard University
##
## LAST MODIFIED: 10/27/2003
## 
## 
## ************************************************************************
## ************************************************************************

### env.base is the environmnet for function main, which may be globalized 
### assign("env.base", env.base, envir= get(".GlobalEnv", envir=.GlobalEnv))
##ebase = get("env.base", envir=.GlobalEnv)

WHO.yourcast <- function(input.lst, ebase=NULL,env.who=NULL){
### we are not attaching ewho to env.base but returning it
### because we need to make it available to stand alone functions.
###  ebase <- get("env.base", envir=parent.frame())
### find the environment to store variables
### the parent.frame() is env.base (driver environment),
### but if driver does not exist then default is .GlobalEnv

  ebase <- get("env.base", envir=parent.frame())

  chk.env <- T ### checking if you are assigning to input env.who or to the local env.  
  if(length(env.who) <= 0){
   env.who <- environment()  ###env.who is the same as local environment 
   chk.env <- F
 }
  
  env.who <- args.of.yourcast(input.lst, ebase,ind=1, ewho=env.who)
  

### Enter all variables as vectors; such as whousercntrylist <- c(2450, 4280)
### cause of death we want to predict

whoyrest <- try(get("yrest",envir=ebase), silent=T)
delta.tol <- .Machine$double.eps

########################### DEFINE STRUCTURES OF cstsid ###############################
###

 who.cntry.digits <- try(get("cntry.digits", envir=ebase), silent=T)
 who.year.digits  <- try(get("year.digits", envir=ebase), silent=T)
 who.digit.first  <- try(get("digit.first", envir=ebase), silent=T)

  ## the covariates specified in who.cov.select will be deleted from the
### age groups specified in who.age.select
who.age.digits <- try(get("age.digits", envir=ebase), silent=T)


whomodel <- get("model",envir=ebase)

### if 1 the covariates (as stored in whocov) are standardized to zero mean and std 1
whostandardize <- get("standardize", envir=ebase)

### if 1 during the preprocessing we search for collinearities and delete redundant
### covariates if necessary
whoelim.collinear <- get("elim.collinear", envir=ebase) 


### BELOW WE HAVE QUANTITIES RELATED TO THE PRIORS
 
### If FALSE the prior  does not have zero mean
### and it is centered around who.mean.age.profile. Having a non-zero mean is achieved by
### subtracting the mean from the dependet variable and adding it back after the forecast
### has been done. If TRUE no action is taken.
who.zero.mean <- get("zero.mean", envir=ebase)

who.zero.mean <- check.zero.mean(who.zero.mean, whomodel)

    
############################ CHOICES FOR OLS ##################################
who.ols.sigma.param <- try(get("ols.sigma.param", envir=ebase), silent=T)
if(class(who.ols.sigma.param)== "try-error")
  who.ols.sigma.param <- list(sigma.bar=1,use.deaths=FALSE,average=TRUE,model.based=FALSE)

############################ PRIOR OVER AGE GROUPS  ###############################

###
who.Ha.deriv <- get("Ha.deriv",envir=ebase);

### If 0 or NA then all age groups are weighted equally. If scalar, then the weight
### of age group a is proportional to a^who.Ha.age.weight. If vector of length A
### then it is taken as the weight vector.
who.Ha.age.weight <- get("Ha.age.weight",envir=ebase);
who.Ha.time.weight <- get("Ha.time.weight",envir=ebase);

### the average standard deviation of the prior. If NA the prior is not used
### (it is like having an infinite standard deviation 
who.Ha.sigma <- try(get("Ha.sigma", envir=ebase),silent = T)
 
if(class(who.Ha.sigma)=="try-error")
  who.Ha.sigma <- NA


who.Ha.sigma.sd <- get("Ha.sigma.sd", envir=ebase)


############################ PRIOR OVER AGE AND TIME GROUPS  ###############################

### the "age part" of the age/time prior
who.Hat.a.deriv <- get("Hat.a.deriv", envir=ebase);

### the "time part" of the age/time prior
who.Hat.t.deriv <- get("Hat.t.deriv",envir=ebase);

### If 0 or NA then all age groups are weighted equally. If scalar, then the weight
### of age group a is proportional to a^who.Ha.age.weight. If vector of length A
### then it is taken as the weight vector.
who.Hat.age.weight <- get("Hat.age.weight", envir=ebase);

who.Hat.time.weight <- get("Hat.time.weight", envir=ebase);

who.Hat.sigma <- try(get("Hat.sigma", envir=ebase), silent=T)

if(class(who.Hat.sigma) == "try-error")
  who.Hat.sigma <- NA
  
who.Hat.sigma.sd <- get("Hat.sigma.sd", envir=ebase)

############################ PRIOR OVER  TIME GROUPS  ###############################

### 
who.Ht.deriv <- get("Ht.deriv", envir=ebase)

who.Ht.age.weight <- get("Ht.age.weight", envir=ebase)

who.Ht.time.weight <- get("Ht.time.weight",envir=ebase)

who.Ht.sigma <- try(get("Ht.sigma", envir=ebase), silent=T)

if(class(who.Ht.sigma) == "try-error")
  who.Ht.sigma <- NA
   
who.Ht.sigma.sd <- get("Ht.sigma.sd", envir=ebase)
                                       
############################ PRIOR OVER CNTRY GROUPS  ###############################
who.Hct.t.deriv <- get("Hct.t.deriv", envir=ebase); ### smooths trend over cntry's
### who.Hct.t.deriv <- c(1)  ### smooths mortality
### we do not offer the option of setting who.Hct.time.weight != NA
### because it gets messy when different countries have different
### length of time series.
### Uniform weighting is the only option we have now
### Not implemented in the code                  
##who.Hct.cntry.weight <- get("Hct.cntry.weight", envir=ebase)
##
who.Hct.c.deriv <- try(get("Hct.c.deriv", envir=ebase), silent=T)
if(class(who.Hct.c.deriv) == "try-error") who.Hct.c.deriv <- NA
who.Hct.sigma   <-  get("Hct.sigma", envir=ebase)
if(class(who.Hct.sigma) == "try-error")
  who.Hct.sigma <- NA

who.Hct.sigma.sd <- get("Hct.sigma.sd", envir= ebase)
who.Hct.time.weight <- get("Hct.time.weight", envir=ebase)
LI.sigma.mean <- try(get("LI.sigma.mean", envir=ebase), silent=T)
if(class(LI.sigma.mean)=="try-error")
  LI.sigma.mean <- 0.2
LI.sigma.sd <- try(get("LI.sigma.sd", envir=ebase), silent=T)
if(class(LI.sigma.sd)=="try-error")
  LI.sigma.sd <- 0.1

N <- try(get("sims", envir=ebase), silent=T)
if(class(N)=="try-error")
  N <- 120
  
n.row <- try(get("n.row", envir=ebase), silent=T)
if(class(n.row)=="try-error")
  n.row <- 1000
filename <- try(get("filename", envir=ebase), silent=T)
if(class(filename) =="try-error")
  filename <- "smalltable.tex"
digits <- try(get("digits", envir=ebase), silent=T)
if(class(digits)=="try-error")
  digits <- 3

maxprint <- 1
graphics.file <- NA
output.file <- FALSE
### directory where mortality data are
datapath <-  getwd()
data.path <- datapath
whodatapath <- datapath

### directory where covariates are
whocovpath <- datapath
 
### directory where priors are  
whopriorpath <- datapath
  
### directory where output files will be saved
whooutpath <- datapath
### initially was a global to save the output of the model into a file but now is only
#
  ## available hgere and need to be set up in this function.  
save.output <- FALSE
 
#### for poisson model only: no log of mortality or depvar
log.poiss <- FALSE
##option to invert a matrix in cxc.age.cntry
svdonly <- 0
if(chk.env){ ### env.who is not the local environmnet
  ev <- environment()
 
  lstfile <- ls(envir=ev)
  for(nm in lstfile){
    val <- get(nm, envir=ev)
    if(!is.environment(val))
    assign(nm, val, envir=env.who)
  }
}
  
return(env.who)
} ## end WHO.yourcast

### DESCRIPTION takes the input list of arguments of yourcast and the nevironmnets
###             assign the values to envs if they have any but for dataobj and formula
###
###             Elena Villalon
###             IQSS Harvard Univ
###             evillalon@iq.harvard.edu
###
######################################################################
args.of.yourcast <- function(input.lst, ebase=NULL, ind=1,ewho=NULL)
  {
    if(length(ewho) <= 0)
      ewho <- environment()
    n.input <-  length(input.lst)
    vec <- ind:n.input
   
    for(i in vec){
   
      if(trim.blanks(input.lst[[i]]) == "dataobj")
        next;
      if(trim.blanks(input.lst[[i]]) == "formula")
        next;
       
      val <- get(input.lst[[i]], envir=ebase, inherits=T)
     
      inlst <- input.lst[[i]]
      if(length(val) > 0 )
        assign(inlst, val, envir=ewho)
       
    }
   
    return(ewho)
  }


WHO.userfile <- function(ebase)
  {
### Enter all variables as vectors; such as whousercntrylist <- c(2450, 4280)
### cause of death we want to predict
ewho <- get("env.who", envir=ebase)

whoyrest <- try(get("yrest",envir=ebase), silent=T)
assign("whoyrest", whoyrest, envir=ewho)
delta.tol <- try(get("delta.tol", envir=ebase), silent=T)
assign("delta.tol", delta.tol, envir=ewho)

########################### DEFINE STRUCTURES OF cstsid ###############################
###

 who.cntry.digits <- try(get("cntry.digits", envir=ebase), silent=T)
 who.year.digits  <- try(get("year.digits", envir=ebase), silent=T)
 who.digit.first  <- try(get("digit.first", envir=ebase), silent=T)
 assign("who.cntry.digits", who.cntry.digits, envir=ewho)
 assign("who.year.digits", who.year.digits, envir=ewho)
 assign("who.digit.first", who.digit.first, envir=ewho)

### the covariates specified in who.cov.select will be deleted from the
### age groups specified in who.age.select
who.age.digits <- try(get("age.digits", envir=ebase), silent=T)
assign("who.age.digits", who.age.digits, envir=ewho)

whomodel <- get("model",envir=ebase)
assign("whomodel", whomodel, envir=ewho)
### if 1 the covariates (as stored in whocov) are standardized to zero mean and std 1
whostandardize <- get("standardize", envir=ebase)
assign("whostandardize", whostandardize, envir=ewho)

### if 1 during the preprocessing we search for collinearities and delete redundant
### covariates if necessary
whoelim.collinear <- get("elim.collinear", envir=ebase) 


### BELOW WE HAVE QUANTITIES RELATED TO THE PRIORS
 
### If FALSE the prior  does not have zero mean
### and it is centered around who.mean.age.profile. Having a non-zero mean is achieved by
### subtracting the mean from the dependet variable and adding it back after the forecast
### has been done. If TRUE no action is taken.
who.zero.mean <- get("zero.mean", envir=ebase)

who.zero.mean <- check.zero.mean(who.zero.mean, whomodel)
############################ CHOICES FOR OLS ##################################
who.ols.sigma.param <- try(get("ols.sigma.param", envir=ebase), silent=T)

############################ PRIOR OVER AGE GROUPS  ###############################

###
who.Ha.deriv <- get("Ha.deriv",envir=ebase);

### If 0 or NA then all age groups are weighted equally. If scalar, then the weight
### of age group a is proportional to a^who.Ha.age.weight. If vector of length A
### then it is taken as the weight vector.
who.Ha.age.weight <- get("Ha.age.weight",envir=ebase);
who.Ha.time.weight <- get("Ha.time.weight",envir=ebase);

### the average standard deviation of the prior. If NA the prior is not used
### (it is like having an infinite standard deviation 
who.Ha.sigma <- get("Ha.sigma", envir=ebase)
who.Ha.sigma.sd <- get("Ha.sigma.sd", envir=ebase)

############################ PRIOR OVER AGE AND TIME GROUPS  ###############################

### the "age part" of the age/time prior
who.Hat.a.deriv <- get("Hat.a.deriv", envir=ebase);

### the "time part" of the age/time prior
who.Hat.t.deriv <- get("Hat.t.deriv",envir=ebase);

### If 0 or NA then all age groups are weighted equally. If scalar, then the weight
### of age group a is proportional to a^who.Ha.age.weight. If vector of length A
### then it is taken as the weight vector.
who.Hat.age.weight <- get("Hat.age.weight", envir=ebase);

who.Hat.time.weight <- get("Hat.time.weight", envir=ebase);

who.Hat.sigma <- get("Hat.sigma", envir=ebase)
who.Hat.sigma.sd <- get("Hat.sigma.sd", envir=ebase)

############################ PRIOR OVER  TIME GROUPS  ###############################

### 
who.Ht.deriv <- get("Ht.deriv", envir=ebase)

who.Ht.age.weight <- get("Ht.age.weight", envir=ebase)

who.Ht.time.weight <- get("Ht.time.weight",envir=ebase)

who.Ht.sigma <- get("Ht.sigma", envir=ebase)
who.Ht.sigma.sd <- get("Ht.sigma.sd", envir=ebase)
                                       
############################ PRIOR OVER CNTRY GROUPS  ###############################
who.Hct.t.deriv <- get("Hct.t.deriv", envir=ebase); ### smooths trend over cntry's
### who.Hct.t.deriv <- c(1)  ### smooths mortality
### we do not offer the option of setting who.Hct.time.weight != NA
### because it gets messy when different countries have different
### length of time series.
### Uniform weighting is the only option we have now
### Not implemented in the code                  
##who.Hct.cntry.weight <- get("Hct.cntry.weight", envir=ebase)
##
who.Hct.c.deriv <- try(get("Hct.c.deriv", envir=ebase), silent=T)

who.Hct.sigma   <-  get("Hct.sigma", envir=ebase)
who.Hct.sigma.sd <- get("Hct.sigma.sd", envir= ebase)
who.Hct.time.weight <- get("Hct.time.weight", envir=ebase)

### directory where mortality data are
datapath <-  getwd()
 
}
### DESCRIPTION: Checks if zero.mean (zmean) is logical or
###              character string.  If it is a string, it
###              checks if the string is the name of a file.
###              It loads the file and evaluate its contents
### INPUT: zmean, a character string, logical, or numeric vector
###        model, a string with model name
### OUPUT: either stop the simulation or return the zero.mean
###         as logical or numeric vector
### AUTHOR:  Elena Villalon
###          evillalon@iq.harvard.edu
###          IQSS, Harvard University
###          05/11/2006
#####################################################
check.zero.mean <- function(zmean, model)
  {
  
    model <- toupper(model)
    ind <- grep(model, c("MAP","EBAYES", "BAYES"))
    
    if(length(ind) <= 0 ||is.logical(zmean))
      return(zmean)
    if(is.character(zmean)){
      if(!file.exists(zmean))
        stop("Provide for a valid input for zero.mean")
      else
        zmean <- eval(as.symbol(load(zmean)))
       
    }
    return(zmean)
  }
