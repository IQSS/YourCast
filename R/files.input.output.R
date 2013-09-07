

##
## USED GLOBALS: envir=env.base, our base environment,is now the environment for yourcast()
##               that wraps up all programs and functions, and also drives the simulation.
##               Also, envir=env.who (an environment), which we create to store all
##               static variables to be shared by functions and programs.  
##               We assign env.who to the env.base, assign("env.who", envir=env.base)
##               if you want to see what is in env.who:  >  ls(envir=get("env.who", envir=env.base));
##
## USED: make.output.filename()
##
## DESCRIPTION:  The output from different models, i.e. "LC", "OLS", "POISSON", and "CXC"
##               is stored in the list "lst.output", which is also in envr= env.who  
##               The length of the list and elements depend on the model used for the predictions.
##               We retreive with get the lst.output and put it in the environment inside
##               esave of build.file.output; and after that save lst.output in the filename. 
##
## FORMAT: build.file.output(ebase=env.base);   build.file.output()  
##
## VALUE:  store in filename the outputs from different models, i.e. yhatin, yhatout
##         for the predicted values of mortality in the insample and outsample periods
##         Also, insampy and outsampy for the actual data in the insample and outsample.
##         Other paramters whyich are modeled specific 
##
## WRITTEN BY: Elena Villalon 
##             evillalon@latte.harvard.edu,
##             CBRSS, Harvard University
##
## LAST MODIFIED: 11/12/2003
## 
## 
## ************************************************************************
## ************************************************************************
build.file.output <- function(lst= NULL, whomodel=NULL, depv=NULL){

  
  if( is.na(whomodel) ||  
     ( toupper(whomodel) != "LC" && toupper(whomodel) =="OLS" && 
      toupper(whomodel) == "POISSON" && toupper(whomodel) =="MAP" && toupper(whomodel) =="BAYES" ) )
      return(list())
    
      coeff <- lst$coeff
      yhat  <- lst$yhat  
      sigma <- lst$sigma
      num <- depv$numerator
      den <- depv$denominator
      formula <- depv$formula
### create the blank file to store data for directory: whooutpath
  num <- strsplit(num, "\\.")[[1]][1]
  den <- strsplit(den, "\\.")[[1]][1]
  filename  <-  make.output.filename(whomodel, num, den);
  filename  <- paste(getwd(),filename,sep="")
  if (file.exists(filename)){
    file.remove(filename)}
    file.create(filename)
### name for present enviroment:
    esave <- environment()
### 
### what to store in filename:
###
   n.lst <- length(names(lst))
   for(i in 1:n.lst)
     assign(names(lst)[i], lst[[i]],envir=esave)

   
     save("whomodel", "coeff", "yhat", "call","sigma","formula",  
        file=filename, compress=T);
  
### this also save the environmnet esave with all in it 
###  save(esave,file=paste(whooutpath,"evenv",sep=""), compress=T, envir=esave); 
}


 make.output.filename <- function(whomodel,num,den){
  return(paste(whomodel,"_",num,"_", den,".dat",sep=""));
}




## ************************************************************************
## ************************************************************************
##
## FUNCTION NAME:      conversion.cntry.mat
##
## 
## IMPORTED functions: list.by.cntry
##
## USED GLOBALS: environments (env.base, env.who or aliases ebase, ewho)
##               the object outputs from models (LC, OLS, POISSON, CXC)
##
## DESCRIPTION: The object outputs contain yhatin, yhatout,insampy,outsampy
##              They are listed by csid or cntry+age combinations,
##              so that they are 1 col matrices with rows=no of years.  
##              We use list.by.cntry to obtain the listing for every
##              country, which are TxA matrices with rows = no of years, and
##              cols= age groups.  This is the output which is needed
##              for the graphics.  
##
## FORMAT:  filename  <- conversion.cntry.mat(obj=lst.output)
##
## VALUE: the object output with yhatin, yhatout, insampy, outsampy
##        listed by cntry, i.e. list elements are TxA matrices for evry cntry
##
## WRITTEN BY: Elena Villalon
##             evillalon@latte.harvard.edu, 
##             CBRSS, Harvard University
##
## Date: 12/09/2003
## 
## ************************************************************************
## ************************************************************************

conversion.cntry.mat <- function(obj=NULL, ebase=env.base){
### I do not think that we need to check for obj and environmnets
### since obj will be passed when we call the function, but
    ebase <- get("env.base", envir=parent.frame())
    env.base <- ebase
    ewho <- get("env.who", envir=ebase)
    whomodel <- get("whomodel", envir=ewho)
    who.digit.first <- get("who.digit.first", envir=ewho)
    who.cntry.digits <- get("who.cntry.digits", envir=ewho)
    who.year.digits <- get("who.year.digits", envir=ewho)
    who.age.digits <- get("who.age.digits", envir=ewho)
    whopopul <- get("whopopul", envir=ewho)
    whopopulos <- get("whopopulos", envir=ewho)
    
  
    if(length(obj) <= 0)
      obj <- get("lst.output", envir=ewho)
### obj should come when the function is called.
  insampy  <- obj$insampy
  outsampy <- obj$outsampy
   
     obj$insampy  <- list.by.cntry(insampy, who.digit.first, who.cntry.digits,
                                   who.age.digits, who.year.digits)
   
###  obj$outsampy <- list.by.cntry(outsampy)
   
if(!is.na(whomodel)){
    yhatin   <- obj$yhatin
   
    yhatout  <- obj$yhatout
   
### list.by.cntry needs the env.who
    obj$yhatin   <- list.by.cntry(yhatin,  who.digit.first, who.cntry.digits,
                                   who.age.digits, who.year.digits)
   
    obj$yhatout  <-  list.by.cntry(yhatout, who.digit.first, who.cntry.digits,
                                   who.age.digits, who.year.digits)
   
  }else{
    obj$whopopul <- list.by.cntry(whopopul,  who.digit.first, who.cntry.digits,
                                   who.age.digits, who.year.digits)
   
    obj$whopopulos <- list.by.cntry(whopopulos, who.digit.first, who.cntry.digits,
                                   who.age.digits, who.year.digits)
   
  }
  assign("lst.output",obj,envir=ewho)
  return(obj)}

### DESCRIPTION Takes the insample and outsample for depvar and estimation
###             joins insample and outsample by row and then the depvar and
###             forecast by columns for every csid (cntry+ age id)

yhat.mat <- function(obj=NULL, ebase=env.base){
### I do not think that we need to check for obj and environmnets
### since obj will be passed when we call the function, but
    ebase <- get("env.base", envir=parent.frame())
    env.base <- ebase
    ewho <- get("env.who", envir=ebase)
       
    if(length(obj) <= 0)
      obj <- get("lst.output", envir=ewho)
### obj should come when the function is called.
    insampy  <- obj$insampy
    outsampy <- obj$outsampy
    yhatin <- obj$yhatin
    yhatout <- obj$yhatout

    
    y <- bind.list(x=insampy,y=outsampy,bycol=FALSE,namex=TRUE)
    nrowy <- lapply(y, nrow)
       
    yhat <- bind.list(x=yhatin,y=yhatout,bycol=FALSE,namex=TRUE)
    nrowyhat <- lapply(yhat, nrow)
   
    ind <- as.list(1:length(nrowy))
    ind <- lapply(ind, function(n){
      if( (nrowy[[n]] - nrowyhat[[n]]) != 0)
        stop('data years and estimation do not agree')
    })
    yyhat <- bind.list(x=y,y=yhat,bycol=T,namex=TRUE,colname=c("y","yhat"))
    
    return(yyhat)
    
  }

 leftside.formula <- function(ff)
      {
              
        fc  <- as.character(ff)
        rs <- fc[3]
        ls  <- fc[2]
        lsf <- as.formula(paste("1 ~", ls))
        vec <- all.vars(lsf)
        if(length(vec) > 1){
          num <- trim.blanks(vec[1])
          den <- trim.blanks(vec[2])
        }else{
         num <- trim.blanks(vec)
         den <- NULL
       }
       
        lst <- list(numerator=num, denominator=den, formula=ff)
        return(lst)  
      }
