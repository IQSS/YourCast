## ************************************************************************
## ************************************************************************
##
## FILE NAME:         cast.funcs
##
## PLACE:     Preprocessing module, funcs functions.
## index of functions:  disease.number,subregvec,
##                      drop.ages, count.observations,
##                      extend.datamat, cntryid, long.causes 
##
## ************************************************************************
## ************************************************************************
##


arguments <- function(driver="yourcast"){
  
  args <- names((formals(driver)))
              }
  
## ************************************************************************
## ************************************************************************
##
## FUNCTION NAME:      digitpull
##
## PLACE:     Preprocessing module, the old Gauss digitpull.g
##            function rewritten in R
##
## IMPORTED functions: none
##
## DESCRIPTION: digitpull() pulls connected digits from a numeric variable
##
## FORMAT:  x <- digitpull(v, start, stop)
##
## INPUT:   v           vector of numeric data, all of the same length,
##                      and with at least stop digits
##          
##          startdig    the first digit of v to return
##
##          stopdig     the last digit of v to return. 
##                      (stopdig must be >= startdig)
##
## OUTPUT:  x           startdig to stopdig digits of v
##
## REMARKS:    As currently written, digitpull() only works on the 
##             integer portion of v
##
## WRITTEN BY: Federico Girosi & Anita Gulyasne Goldpergel
##             fgirosi@latte.harvard.edu, anitag@latte.harvard.edu
##             CBRSS, Harvard University
## 
## Last modified: 11/14/2002
## 
## ************************************************************************
## ************************************************************************

digitpull <- function(v, startdig, stopdig) {
  v <- as.numeric(v)
  v <- matrix(v, ncol = 1)
  v <- trunc(v)

### error checking
  if (stopdig < startdig) {
    stop(message="stopdig < startdig in proc preproc.digitpull")
  }

  if (min(abs(v)) < 10^(stopdig-1)) {
    stop(message="narrowest element in v not as wide as stopdig in proc digitpull()")
  }

  if (startdig < 1) {
    stop(message="startdig < 1 in proc digitpull()")
  }
  
  if (ncol(v) > 1) {
    stop(message="v not column vector in proc digitpull()")
  }
### end error checking  

  n <- trunc(log10(v[1]))+1
  a1 <- trunc(v/10^(n-stopdig))
  a2 <- trunc(v/10^(n-startdig+1))
  ret <- a1 - a2*10^(stopdig - startdig +1)
  return(matrix(ret))
}

## ************************************************************************
## ************************************************************************
##
## FUNCTION NAME:      modify.mean.age.profile
##
## 
## IMPORTED functions: list.by.cntry, list.by.csid
##
## USED GLOBALS: none
## DESCRIPTION: given a cross-sectional time series in the usual list
##              format, it applies a given function to each age profile
##              and returns the modified cross-sectional time-series
##              
## FORMAT: new.y <- modify.mean.age.profile(y, func)
##
## INPUT: y: a cross-sectional time series in list format
##
##        func: a function of one argument, which takes an age profile and returns
##              a vector of the same length
##
## OUTPUT:  new.y: a cross-sectional time series obtained from y by replacing each age profile
##                 x by the vector func(x).
##     
##
## WRITTEN BY: Federico Girosi
##             fgirosi@latte.harvard.edu
##             CBRSS, Harvard University
## 
## 
## ************************************************************************
## ************************************************************************

 
 modify.age.profiles <- function(y,func,param){
  ebase <- get("env.base", envir=parent.frame())
  env.base <- ebase
  ewho <- get("env.who", envir=ebase)
  df <- get("who.digit.first", envir=ewho)
  dc <- get("who.cntry.digits", envir=ewho)
  da <- get("who.age.digits", envir=ewho)
  dy <- get("who.year.digits", envir=ewho)
  agev <- get("age.vec", envir=ewho)
  n.age <- length(agev)
  insampy <- list.by.cntry(y,who.digit.first=df,who.cntry.digits=dc, who.age.digits=da, who.year.digits=dy);
  
   insampy <-  lapply(insampy,FUN=function(x,func,param,agev){
     n.age <- length(agev)
     mat <- t(apply(x,1,FUN=func,param))
     dm <- dim(mat)
    
     if(length(dm) <= 0)
       mat <- as.matrix(mat)
     
     if(n.age <= 1 && dm[1]==1){##one age group
       mat <- t(mat)
       colnames(mat) <- agev}
 
     return(mat)},func,param,agev);
 
   new.whoinsampy <- list.by.csid(insampy);

   return(invisible(new.whoinsampy))
 }





## ######################################################################################
##
## FUNCTION NAME:      bind.list
##
## IMPORTED functions: none
##
## DESCRIPTION: given two lists of equal length, whose elements are matrices,
##              it returns one list of same length, whose elements are the concatenation
##              of the corresponding elements in the two lists
##              
##
## FORMAT: lst <- bind.list(x,y,bycol=bycol,namex=namex)
##
## INPUT: x,y:  lists of length n. x[[i]] and y[[i]] must have either the same
##              number of rows or the same number of columns
##        bycol: logical. If TRUE concatenation is performed along columns
##               otherwise along rows
##        namex: logical. If TRUE the resulting list inherits the names of x
##               otherwise it inherits the names of y
##
## OUTPUT: lst: a list of length n. if bycol=FALSE then lst[[i]] = rbind(x[[i]],y[[i]]),
##              otherwise lst[[i]] = cbind(x[[i]],y[[i]]).
##              if namex=TRUE then names(lst) = names(x), otherwise names(lst)=names(y)
##
##  Federico Girosi
##  CBRSS
##  Harvard University
##  34 Kirkland St. 
##  Cambridge, MA 02143
##  fgirosi@latte.harvard.edu
##  (C) Copyright 2003   


bind.list <-  function(x,y,bycol=FALSE,namex=TRUE,colname=NULL) {
 
  n <- length(x);
  z <- as.list(1:n);
  
  if(length(x) <= 0)
    return(y)
  if(length(y) <= 0)
    return(x)
 
  z <- lapply(z,FUN=function(n,x=x,y=y,bycol=bycol,colname=colname){
                             
 ###                             print(ncol(x[[n]]))
 ###                             print(ncol(y[[n]]))
 ###                             print(nrow(x[[n]]))
 ###                             print(nrow(y[[n]]))
                             
                              cbool <- (ncol(x[[n]]) > 0 && ncol(y[[n]]) > 0)
                              cbool <- cbool && (ncol(x[[n]]) != ncol(y[[n]]))
                              
                              if (!bycol && cbool)                                
                                stop("Different number of cols no binding possible")
                              
                              rbool <- (nrow(x[[n]]) > 0 && nrow(y[[n]]) > 0)
                              rbool <- rbool && (nrow(x[[n]]) != nrow(y[[n]]))
                              if (bycol && rbool)
                                stop("Different number of rows no binding possible")
                                   
                              if (!bycol)
                               return(rbind(x[[n]],y[[n]]))
                              if(length(colname) > 0 )
                               {
                                 mat <- cbind(x[[n]],y[[n]])
                                 colnames(mat) <- colname
                                 return(mat)
                               }else
                              return(cbind(x[[n]],y[[n]]))},
              x,y,bycol,colname);
  if (namex)
    names(z) <- names(x)
  else
    names(z) <- names(y)
  return(z);
}


## ************************************************************************
## ************************************************************************
##
## FUNCTION NAME:      make.average.age.profile
##
## 
## IMPORTED functions: list.by.cntry
##
## USED GLOBALS: none
##
## DESCRIPTION: given a cross-sectional time series in list format,
##              it returns the average  age profile (where the average
##              is taken over countries and years)
##
## FORMAT:  mean.age.profile <- make.average.age.profile(y)
##
## INPUT:   y: (list) a cross-sectional time series over C countries
##                      and A age groups, in list format
##
##
## OUTPUT: mean.age.profile: A x 1 vector. The average age profile, with elements named
##                           according to age groups
##
## WRITTEN BY: Federico Girosi
##             fgirosi@latte.harvard.edu
##             CBRSS, Harvard University
## 
## 
## ************************************************************************
## ************************************************************************

make.average.age.profile <- function(y,bycntry=FALSE, ebase=env.base){
  ebase <- get("env.base", envir=parent.frame())
  env.base <- ebase
  ewho <- get("env.who", envir=ebase)
  df <- get("who.digit.first", envir=ewho)
  dc <- get("who.cntry.digits", envir=ewho)
  da <- get("who.age.digits", envir=ewho)
  dy <- get("who.year.digits", envir=ewho)

  insampy <- list.by.cntry(y,who.digit.first=df,who.cntry.digits=dc, who.age.digits=da, who.year.digits=dy);
  age.profile.list <- lapply(insampy,FUN=function(x){return(colMeans(as.data.frame(x),na.rm=TRUE))});
  result <- age.profile.list;
 
  if (!bycntry){
    mean.age.profile <- 0*age.profile.list[[1]];
    for (i in 1:length(age.profile.list)){
      mean.age.profile <- mean.age.profile + age.profile.list[[i]];
    }
    result <- mean.age.profile/length(age.profile.list);
  }
  return(result);
}


## ************************************************************************
## ************************************************************************
##
## FUNCTION NAME:      list.by.cntry
##
## 
## IMPORTED functions: none
##
## USED GLOBALS: none
##
## DESCRIPTION: given a cross-sectional time series over C countries
## and A age groups, in list format, it returns a list of C matrices,
## one for each country. The matrices have dimension T x A, where T is
## the length of the time series and A is the number of age groups.
## Each column of a matrix is the time series for the corresponding age group
##
## FORMAT:  insampy <- list.by.cntry(whoinsampy)
##
## INPUT:   whoinsampy: (list) a cross-sectional time series over C countries
##                      and A age groups, in list format.
##
##
## OUTPUT: insampy: list of C matrices. The names of the elements of
##                  the list are the country codes. Each matrix is T x A. The names of
##                  the rows are the years of the time series, and the names of the
##                  columns are the names of the age groups.
##
## WRITTEN BY: Federico Girosi
##             fgirosi@latte.harvard.edu
##             CBRSS, Harvard University
## 
## 
## ************************************************************************
## ************************************************************************

list.by.cntry <- function(x, ebase=1,
                          who.digit.first=0, who.cntry.digits=4,
                          who.age.digits=2, who.year.digits=4)
{
  if (!identical(.GlobalEnv, parent.frame()) && length(ebase) > 0 ){
    
    ebase <- get("env.base", envir=parent.frame());
    env.base <- ebase;
    ewho <- get("env.who", envir=ebase)
    cntry.vec <- get("cntry.vec", envir=ewho)
    age.vec <- get("age.vec", envir=ewho)
      
    who.digit.first  <- get("who.digit.first", envir=ewho) 
    who.cntry.digits <- get("who.cntry.digits", envir=ewho)
    who.age.digits   <- get("who.age.digits",envir=ewho)
    who.year.digits  <- get("who.year.digits", envir=ewho)
  }
### Structure of dataset cstsid:

  digit.cntry.begin <- who.digit.first + 1
  digit.cntry.end   <- who.digit.first + who.cntry.digits
  digit.age.begin   <- digit.cntry.end + 1
  digit.age.end     <- digit.cntry.end + who.age.digits
  digit.year.begin  <- digit.age.end   + 1
  digit.year.end    <- digit.age.end   + who.year.digits
  
  cs.vec <- as.numeric(names(x));
 
  cs.cntry.vec <- digitpull(cs.vec,digit.cntry.begin,digit.cntry.end);  
  cntry.vec <- unique.default(cs.cntry.vec);
       
  age.vec <- unique.default(digitpull(cs.vec,digit.age.begin,digit.age.end));
  
  n.cntry <- length(cntry.vec);
  n.age <- length(age.vec);
  newlist <-  list();
  
  for (i in 1:n.cntry){
    cntry <- cntry.vec[i];
### get all the elements of the list x which correspond to country cntry 
    c.list <- x[cs.cntry.vec == cntry];
###    print(names(c.list))
    colnam <- digitpull(as.numeric(names(c.list)),digit.age.begin,digit.age.end)
    
   ### rownam <- digitpull(as.numeric(rownames(c.list[[1]])),digit.year.begin, digit.year.end);
    rownam <- as.numeric(rownames(c.list[[1]]))
    c.mat <- as.matrix(as.data.frame(c.list));
    colnames(c.mat) <- colnam;
    rownames(c.mat) <- rownam;
    newlist <-  c(newlist, list(c.mat));
  }
  names(newlist) <- cntry.vec;
  return(newlist);
}


## ************************************************************************
## ************************************************************************
##
## FUNCTION NAME:      list.by.csid
##
## 
## IMPORTED functions: none
##
## USED GLOBALS: none
##
## DESCRIPTION: it inverts the effect of list.by.cntry.
##
## FORMAT:  whoinsampy <- list.by.csid(insampy)
##
## INPUT:   insampy: (list) the output of list.by.cntry
##
##
## OUTPUT: whoinsampy: (list) ther cross-sectional time series which,
##                      when given as input to list.by.cntry, generates insampy
## 
## WRITTEN BY: Elena Villalon
##             evillalon@latte.harvard.edu
##             CBRSS, Harvard University
## 
## ***********************************************************************
## ************************************************************************

list.by.csid <- function(insampy.c,col.name="dth", ebase=1){
    if (identical(.GlobalEnv, parent.frame()) || length(ebase) <= 0){
### if you are working in the library no need to source
      message("Setting the number of age digits equal to 2") 
      who.age.digits   <- 2
    }else{
      ebase <- get("env.base", envir=parent.frame());
      env.base <- ebase;
      ewho <- get("env.who", envir=ebase)
  ### Structure of dataset cstsid: 
      who.age.digits   <- get("who.age.digits",envir=ewho)
    }
 
    cntry.vec <- as.numeric(names(insampy.c));
    age.vec <- as.numeric(colnames(insampy.c[[1]]));
    n.age <- length(age.vec);
 
    csid.vec  <- sort(kronecker(cntry.vec *10^(who.age.digits), age.vec, FUN="+"))
    n.csid <- length(csid.vec)
    indx <- as.list(1:n.csid)
    names(indx) <- as.character(csid.vec)
 
  newlist <- function(x){
  
    resto <- x %% n.age
    indx <-  x %/% n.age
 
    if(resto > 0) indx <- indx + 1
    naming <- as.character(csid.vec[x])
    nm.vec <-  rownames(insampy.c[[indx]]);
    nm.vec <-  paste(naming, nm.vec, sep="") 
    if(resto == 0){
      x <- insampy.c[[indx]][,n.age]
      x <- as.matrix(x);
      rownames(x) <- nm.vec;
      colnames(x) <- col.name;
    } else {
      x <- insampy.c[[indx]][,resto];
      x <- as.matrix(x);
      rownames(x) <- nm.vec;
      colnames(x) <- col.name;
    }
    return(x)
  }
  
  whoinsampy <- lapply(indx,FUN="newlist")
  return(invisible(whoinsampy))
}


## ************************************************************************
## ************************************************************************
##
## FUNCTION NAME:      split.matrix
##
## 
## IMPORTED functions: none
##
## USED GLOBALS: none
##
## DESCRIPTION: it splits a matrix in blocks along the rows.
##              
## FORMAT: Y <- split.matrix(X,blocks)
##
## INPUT: X: a matrix, vector, data frame or 2-dimensional array
##
##        blocks: a vector of integers, which sum up to the number of rows of X.
##                it specify the size of the blocks
##
## OUTPUT: Y a list of length equal to the length of blocks, with the blocks as elements
##         the blocks are stored as matrices.
##
## EXAMPLE: if X is a matrix with 10 rows and blocks = c(2,3,5) then Y has the 3 elements;
##           the first elements has the first 2 rows of X, the second has rows 3 to 5, and
##           the 3rd element has rows 6 to 10.
##
## WRITTEN BY: Federico Girosi
##             fgirosi@latte.harvard.edu
##             CBRSS, Harvard University
## 
## 
## ************************************************************************
## ************************************************************************


split.matrix <- function(x,blocks){
  x <- as.data.frame(x);
  if (nrow(x) != sum(blocks)) stop("Block structure is wrong in partition.vector");
  f <- factor(rep(1:length(blocks),blocks));
  res <- split(x,f);
  res <- lapply(res,FUN=as.matrix);
  return(res);  
}


## ************************************************************************
## ************************************************************************
##
## FUNCTION NAME:  make.forecast
##
## 
## IMPORTED functions: none
##
## USED GLOBALS: none
##
## DESCRIPTION: given a list of coefficients and covariates it returns
##              the corresponding list of predicted values 
##              
## FORMAT: yhat <- make.forecast(coeff,X)
##
## INPUT: coeff: list of regression coefficients, one element for each cross-section;
##
##            X: list  of covariate matrices, one element for each cross-section;
##
##
## OUTPUT: yhat: list of predicted values: yhat[[n]] =X[[n]] %*% coeff[[n]].
##               It inherits the names of X.
##
##
## WRITTEN BY: Federico Girosi
##             fgirosi@latte.harvard.edu
##             CBRSS, Harvard University
## 
## 
## ************************************************************************
## ************************************************************************

make.forecast <- function(coeff,X, verbose=T){
   verb <- try(get("verbose", envir=get("env.base", envir=parent.frame())),silent=T)
    if(class(verb)!="try-error")
      verbose <- verb
  indx <- seq(1,length(X));
  names(indx) <-  names(X);
  yhat <- lapply(indx,FUN=function(n,coeff){
    beta <- coeff[[n]];
    Z <- X[[n]];
    mult <- try(Z%*%beta, silent=T)
    if(class(mult)=="try-error"){
      messout(names(X[n]), verbose)
      messout(dim(beta), verbose)
      messout(dim(Z),verbose)
      return(0)
    }else
      return(mult)}, coeff);  
##    return(Z%*%beta)},coeff);
  return(yhat);
}

## ************************************************************************
## ************************************************************************
##
## FUNCTION NAME:  death.from.logmortality
##
## 
## IMPORTED functions: none
##
## USED GLOBALS: none
##
## DESCRIPTION: to compute the cross-sectioanl time series of deaths
##              given the cross-sectional time series of log mortality and population
##   
##              
## FORMAT:      d <- death.from.logmortality(logm,popu)
##
## INPUT:    logm: (list) log-mortality cross-sectional time series
##
##           popu: (list) population cross-sectional time series
##  
##
## OUTPUT:      d: (list) cross-sectional time series of deaths
##               
## WRITTEN BY: Federico Girosi
##             fgirosi@latte.harvard.edu
##             CBRSS, Harvard University
##
## Last modified: 08/7/2003
## 
## ************************************************************************
## ************************************************************************


death.from.logmortality <- function(logm,popu){
  ind <- as.list(1:length(logm));
  names(ind) <-  names(logm);
  func <- function(n,logm,popu){
    x <- logm[[n]];
    y <- popu[[n]];
    dth <- y*exp(x)
    dth[dth <= 0.5] <- 0.5
    return(dth)
  }
  d <-  lapply(ind,FUN="func",logm,popu);
  return(invisible(d));
}




##
## DESCRIPTION: to compute the cross-sectioanl time series of deaths
##              given the cross-sectional time series of log mortality and population
##   
##              
## FORMAT:      s <- model.string()
##
## INPUT:    none
##
## OUTPUT:      s: (string) a string which uniquely identifies the model being used,
##                 consisting of whomodel and some other strings specifying parameters value.
##                 Right now is equal to whomodel, except when whomodel is "OLS".
##                 In this case, if we use an heteroskedastic OLS, with observations
##                 weighted by number of deaths (a la Wilmoth) the string is "OLSH"
##                
## WRITTEN BY: Federico Girosi
##             fgirosi@latte.harvard.edu
##             CBRSS, Harvard University
##
## Last modified: 08/7/2003
## 
## ************************************************************************
## ************************************************************************


model.string <- function(ebase=env.base){
ebase <- get("env.base", envir=parent.frame());
  env.base <- ebase;
  ewho <- get("env.who", envir=ebase)
  whomodel <- get("whomodel", envir= ewho)
  who.ols.sigma.param <- get("who.ols.sigma.param", envir=ewho)
  s <-  whomodel;
  if (whomodel == "OLS"){
    if (who.ols.sigma.param$use.deaths == TRUE && who.ols.sigma.param$average == FALSE){
      s <- "OLSH";
    }
  }
  return(s);
}


## ************************************************************************
## ************************************************************************
##
## FUNCTION NAME:  strpad
##
## 
## IMPORTED functions: none
##
## USED GLOBALS: none
##
## DESCRIPTION: to convert a list of strings of different lenghts into
##              a list of strings of the same length, given by the maximum
##              string length in the list. This is achieved by padding the shorter
##              string with blanks, It is useful in printing statements.##   
##              
## FORMAT:      s <- strpad(x)
##
## INPUT:       x: (list) a list of strings
##
## OUTPUT:      s: (list) the same list as x, but with  right-padding with blanks, so
##                 that all the elements have the same length.
##                
## WRITTEN BY: Federico Girosi
##             fgirosi@latte.harvard.edu
##             CBRSS, Harvard University
##
## Last modified: 08/7/2003
## 
## ************************************************************************
## ************************************************************************


strpad <- function(x){
  l <- max(nchar(x));
  z <-  lapply(x,FUN=function(x,l){a <- paste(x,paste(rep(" ",l-nchar(x)),sep="",collapse=""),sep="");return(a)},l)
}




## ************************************************************************
## ************************************************************************
##
## FUNCTION NAME:  list.coeff.by.cntry 
##
## 
## IMPORTED functions: none
##
## USED GLOBALS: none
##
## DESCRIPTION: 
##              
## FORMAT:   s <- list.coeff.by.cntry(coeff,cntry.vec,age.vec)
##
## INPUT:    coeff: (list) a cross-sectional list of regression coefficients,
##                  indexed by country and age.
##
##       cntry.vec: (vector) vector of country codes represented in the list coeff
##
##         age.vec: (vector) vector of age groups represented in the list coeff
##
##         if either of cntry.vec or age.vec are not defined, get them from the
##         ewho environmnet: ewho= get("env.who", envir= env.base); ebase= env.base
##
## OUTPUT:       s: (list) a cross-sectional list indexed by countries.
##                  For each county the regression coefficients corresponding
##                  to different age groups are concatenated one after the
##                  the other, according to the order of the age groups in age.vec
##
##
##                
## WRITTEN BY: Federico Girosi
##             fgirosi@latte.harvard.edu
##             CBRSS, Harvard University
##
## Last modified: 08/7/2003
## 
## ************************************************************************
## ************************************************************************


list.coeff.by.cntry <- function(coeff,cntry.vec=NULL,age.vec=NULL, ebase = env.base){
ebase <- get("env.base", envir=parent.frame());
  env.base <- ebase;
  ewho <- get("env.who", envir=ebase)
  who.age.digits <- get("who.age.digits", envir=ewho)
  if(length(cntry.vec) <= 0)
    cntry.vec <- get("cntry.vec", envir=ewho)
  if(length(age.vec) <= 0)
    age.vec <- get("age.vec", envir=ewho)
  
  age.char <- formatC(age.vec,width= who.age.digits,format="d", flag="0")

### cov.list is the vector of covariates in the n-th cross-section;
  names(cntry.vec) <- cntry.vec;
  b <- lapply(cntry.vec,FUN=function(x,coeff,age.char){
    ind <- paste(x,age.char,sep="");
    return(unlist(coeff[ind]))},coeff,age.char);
  return(b);  
}



list.by.cntry.long <- function(x, ebase=env.base){
ebase <- get("env.base", envir=parent.frame());
  env.base <- ebase;
### In terms of the global, the structure of dataset cstsid:
  ewho <- get("env.who", envir=ebase)
  who.digit.first  <- get("who.digit.first", envir=ewho) 
  who.cntry.digits <- get("who.cntry.digits", envir=ewho)
  who.age.digits   <- get("who.age.digits",envir=ewho)
  who.year.digits  <- get("who.year.digits", envir=ewho) 
### Structure of dataset cstsid:
  digit.cntry.begin <- who.digit.first + 1 
  digit.cntry.end   <- who.digit.first + who.cntry.digits
  digit.age.begin   <- digit.cntry.end + 1
  digit.age.end     <- digit.cntry.end + who.age.digits
  digit.year.begin  <- digit.age.end   + 1 
  digit.year.end    <- digit.age.end   + who.year.digits 
  
  cs.vec <- as.numeric(names(x));
  cs.cntry.vec <- digitpull(cs.vec,digit.cntry.begin,digit.cntry.end);  
  cntry.vec <- unique.default(cs.cntry.vec);
  n.cntry <- length(cntry.vec);
  nam <- lapply(x,function(x){return(rownames(x))})
  ind <-  as.list(cntry.vec);
  names(ind) <- cntry.vec;

  func <- function(n,x,cs.cntry.vec,nam){
    cntry <- n;
    c.list <- x[cs.cntry.vec == cntry];
    nam.list <- nam[cs.cntry.vec == cntry];
    y.long <- unlist(c.list);
    y.long.names <- unlist(nam.list)
    names(y.long) <- y.long.names;
    y.long <- as.array(y.long);
    return(y.long)
  }


  l <- lapply(ind,FUN=func,x,cs.cntry.vec,nam);
 return(invisible(l))
}

######################################################################

unlist.by.cntry.long <- function(y, ebase= env.base){
ebase <- get("env.base", envir=parent.frame());
  env.base <- ebase;
### In terms of the global, the structure of dataset cstsid:
  ewho <- get("env.who", envir=ebase)
  who.digit.first  <- get("who.digit.first", envir=ewho) 
  who.cntry.digits <- get("who.cntry.digits", envir=ewho)
  who.age.digits   <- get("who.age.digits",envir=ewho)
  who.year.digits  <- get("who.year.digits", envir=ewho) 
### Structure of dataset cstsid: 
  digit.cntry.begin <- who.digit.first + 1 
  digit.cntry.end   <- who.digit.first + who.cntry.digits
  digit.age.begin   <- digit.cntry.end + 1
  digit.age.end     <- digit.cntry.end + who.age.digits
  digit.year.begin  <- digit.age.end   + 1 
  digit.year.end    <- digit.age.end   + who.year.digits 
  
  new.list <- NULL;
  for (i in 1:length(y)){
    x <- y[[i]];
    agevec <- digitpull(as.numeric(rownames(x)),digit.age.begin,digit.age.end);
    f <- factor(agevec);
    z <- split(as.data.frame(x),f);
    names(z) <- paste(names(y)[i],conv.char(as.numeric(names(z))),sep="");
    z <-  lapply(z,function(y){as.matrix(y)});
    new.list <- c(new.list,z);
  }
  return(invisible(new.list))
}



## ************************************************************************
## ************************************************************************
##
## PROGRAM NAME:      list.covariates
##
## PLACE:     Preprocessing module, covariates for all age groups
##
## IMPORTED functions: none
##
## DESCRIPTION: It obtains the names of the covariates contributing to list elemnets
##              for every csid  or cntry and age combination.  
##
## FORMAT:      list.covariates(lst)
##
## INPUT:      The covariates whocov, whoinsampx or whoutsampx,from make.mortality.data;
##             They are list whose elements are csid or country-age ccombinations. 
##
## OUTPUT:     Another list with same number of elements as the original list (or input list)
##             but with only one row with the names of contributing covariates.
##
## WRITTEN BY: Elena Villalon 
##             evillalon@latte.harvard.edu
##             CBRSS, Harvard University
## 
## Last modified: 12/11/2003
##
## ************************************************************************
## *************************

  list.covariates <- function(list){
    if( ncol(list[[1]])<=1 )
        cov.lst <- lapply(list, function(x){whocovlist <- rownames(x)})
    else
      cov.lst <- lapply(list, function(x){whocovlist <- colnames(x)})}

## ************************************************************************
### DESCRIPTION: function that builds covariates from files cov.txt;
###              cov="fat","hc", "gdp", "tfr". The files have one single entry (or row)
###              for each cntry, year combination and do not depend on age groups.
###              Expands the entries to the number of age.groups specify
###              in the function, by repeating as many times as age groups
###              any entry belonging to a single cntry+year combination.
###
### INPUT: age.groups, numeric vector of ages;
###        select,a numeric code representing the cntry we want to select or NA for all of them;
###        covname, a character string with the covariate name to expand (either of fat, hc, gdo, tfr);
###        covstring, string with location in the directory path to find files cov.txt;
###        save.FULL, a bolean either T or F, which will save some files.
###        ncol is the cols of original reduced file (2 0r 3 for gender independent or dependent);
###        cov.path, string with directory name to save the file
###
### OUTPUT: expanded matrix of three cols,first col is covariate values and 
### second cstsid = cntry+age+year, third is gender or NA (if gender independent)
### We include the choice of saving the matrices of expanded covariates 
### into files with the input parameter save.FULL = TRUE, which will generate
### the corresponding FULL_cov.txt file.
###
### USES:select.gender, save.file.WHO
###
## WRITTEN BY: Elena Villalon 
##             evillalon@iq.harvard.edu, 
##             IQS, Harvard University
## 
## Last modified: 04/27/2006
###


    expand.cov <- function(covname, covstring, ncol,type=NA,
                          select= NA, save.FULL=T,
                          cov.path=NULL, lstind=lstind,age.groups= (0:16)*5){
      
      todiscard <- lstind$digit.first
      cdigits <- lstind$cntry.digits
      cdigits <- cdigits + todiscard
      adigits <- lstind$age.digits
      ydigits <- lstind$year.digits
      
         if(length(cov.path) <= 0)
           cov.path <- getwd()
           cov <- NULL
           flag <- T
        if(length(covstring) > 0){
          cov <- scan(file=covstring,
                      ##what=c(numeric(0),integer(0),integer(0)), (You do not need to specify what;)
                      quiet=T,multi.line=T)
          cov <- matrix(cov,ncol=ncol,byrow=T)    
        }else{
         
           cov <- as.matrix(covname)
           covname <- colnames(cov)[1]
         }
      
    
### for covariates are gender independent, if given then gender=NA
         if(ncol >= 3 && is.character(covname))
           colnames(cov) <- c(covname, "cstsid","gender")
         else if(ncol >= 2 && is.character(covname))
           colnames(cov) <- c(covname, "cstsid")
        
### if only US
        if (length(na.omit(select)) > 0 && select=="US"){
          indx.US <- grep("^2450",as.character(cov[,2]))
          cov <- cov[indx.US, ]}
### else all the countries are included.         

        if(identical(all.equal(type,"age.independent"),T)){
        
           mm <-  array(0,dim=0)
           covfac <- split.data.frame(cov, cov[,3])
           mmlst <- lapply(covfac,
                           FUN= select.gender,
                           covname=covname,ncol=ncol, lstind=lstind, age.g=age.groups)
           for(mat in mmlst)
           mm <- rbind(mm, mat)
         }else{
          
           mm <- select.gender(cov,covname,ncol,lstind, age.g=age.groups)
          
         }
        if(save.FULL == T )
          save.file.WHO(mm, covname,ncol=3,type,cov.path)

### asssuming list of covariates is c("fat","gdp","hc", "tfr")
### creates output and text files for preprocessing routines
### Now, to be consistent with previous code eliminate col= 3 with gender=NA

        return(invisible(mm)) }

   select.gender <- function(cov,covname,ncol, lstind, age.g)
{
  todiscard <- lstind$digit.first
  cdigits <- lstind$cntry.digits
  cdigits   <- cdigits + todiscard
  adigits <- lstind$age.digits
  ydigits <- lstind$year.digits
  age.groups <- age.g

  wd  <- adigits
  
  age.n <- length(age.groups)
  time.series <- unique.default(as.numeric(cov[,2])%%10^ydigits)
  time.n    <- length(time.series) 
  cntry.vec <-  unique.default(as.numeric(cov[,2])%/%10^(ydigits+adigits))
  cntry.n   <- length(cntry.vec)
  cov.ex1 <- rep(cov[,1],age.n)
  cov.ex2 <- rep(cov[,2],age.n)
  if ( ncol >= 3)
    cov.ex3 <- rep(cov[,3],age.n)
  else if( ncol >= 2){
    cov.ex3 <- NA*cov.ex1
    
  }else
  stop("Insufficient cols")
         
  cov.expand <- cbind(cov.ex1,cov.ex2,cov.ex3)
       
  colnames(cov.expand) <- c(covname, "cstsid","gender")
       
### list by years of time series
  cov.lst <- split.data.frame(cov.expand,as.numeric(cov.expand[,2])%%10^ydigits)
### Given that all year groups have same covariates.
### Repeating the values of covs for age ="00" modify
### the time series to include all age groups  
       
  cov.lst <- lapply(cov.lst, function(x, age.groups){
    ord <- order(as.numeric(x[,2]))
    x     <- x[ord, ]
    yr    <- as.numeric(x[1,2])%%10^ydigits
    yr <- sapply(yr, FUN="complete.flags", ydigits)
    years.string <- formatC(age.groups, width =wd, format="d", flag="0")
    ctry.age <- kronecker(unique.default(as.numeric(x[,2])%/%10^(ydigits+adigits)),
                          years.string, FUN="paste", sep="")
    ctry.age <- sapply(ctry.age, FUN="complete.flags", (cdigits+adigits))
    ctry.age.yr <- as.numeric(paste(ctry.age,yr,sep=""))
    ctry.age.yr <- sort(ctry.age.yr)
    if(length(ctry.age.yr) != length(x[,2])){
      stop("Building the wrong matrix in expand.cov for make.covariates")
    }else
    x[,2] <- ctry.age.yr
    return(x)}, age.groups)
### end cov.lst
      
        nmc <- names(cov.lst)
        nmc <- sapply(nmc, FUN="complete.flags", adigits)
     
       names(nmc) <- NULL
       names(cov.lst) <- nmc
        cov.frame <- lapply(cov.lst, as.data.frame)
        st <- unlist(cov.frame, recursive=F)
        nmcov <- paste(names(cov.frame),".",covname,sep="")
        nmcst <- paste(names(cov.frame),".","cstsid",sep="")
        nmgen <- paste(names(cov.frame),".","gender",sep="")
### build a matrix with 3 cols= covariate, cstsid and gender
### with age groups=5*(0:16); and years=1920:2000.
### Covariates are repetition of one given age, say "45" or "00"
### for every single year in time series.
       
        mm <- matrix(,nrow=length(cov.lst)* age.n * cntry.n, ncol= 3)  
        colnames(mm) <- c(covname, "cstsid","gender")
        for(n in 1:3){  ##for2 build cols of datamat
          mmc  <- dim(0)
          coln <- sapply(1:length(st), function(i,st){
                  x <- n + (i-1) * 3
                  if(n <=1 && is.element(names(st)[x],nmcov)==T)
                    mmc <- c(mmc,as.vector(st[x]))
                  else if(n <=2 && is.element(names(st)[x],nmcst)==T)
                    mmc <- c(mmc,as.vector(st[x]))
                  else if(is.element(names(st)[x],nmgen)==T)
                    mmc <- c(mmc,as.vector(st[x]))
                  return(mmc)}, st)
          coln <- unlist(coln)
                    
          if(length(coln) < length(cov.lst)* age.n * cntry.n)
            warning("Wrong building jumbo expand.cov in make.covariates")
          else
            mm[,n] <- coln
        } ##end for2 build cols of datamat
### order matrix rows according to values of second col=cstsid
        indx <- order(mm[,2])
        mm   <- mm[indx, ]
        colnames(mm) <-  c(covname, "cstsid","gender")
    
  return(mm)}

 
### expand.cov creates for any age.groups-independent covariates, 
### the matrix of two cols (for strata independent) or three cols for
### gender dependent covariates, other cols are cov values and the cstsid.  
### If we want to save that matrix for further development
### the function expand.cov has the parameter save.WHO = (T or F). 
### If we set the parameter save.WHO= TRUE, expand.cov will call
### function save.file.WHO, which creates the file WHO_cov.txt
### cov= "fat", "gdp", "hc", "tfr"; in the data directory
### specified with whocovpath.
### save.file.WHO(mm, covname) takes two argument the matrix
### of two (three) cols with the covariate and cstsid, and the covariate name. 
### Depending on the covariate it opens the appropiate file name
### and save the matrix = mm in whocovpath directory.
### It also return mm as a list, but has no use in this file.
##
## WRITTEN BY: Elena Villalon 
##             evillalon@latte.harvard.edu, 
##             CBRSS, Harvard University
## 
## Last modified: 10/3/2003


save.file.WHO <- function(mm, covname, ncol, type,whocovpath, verbose=T){
   ebase <- try(get("env.base", envir=parent.frame()), silent=T)
   if(class(ebase)!="try-error")
     verbose <- get("verbose", envir=ebase)
 

###colnames(mat) <-  c(covname, "cstsid","gender")
### asssuming list of covariates is c("fat","gdp","hc", "tfr")
### creates output and text files for preprocessing routines
  Covs.lst <- list();
  cv.symb = as.symbol(covname)
    Covs.lst <- c(Covs.lst, cv.symb=list(mm))
    if(length(grep("txt", covname)) > 0)
      cov.path <- paste(whocovpath,"/FULL.",covname,sep="")
    else
       cov.path <- paste(whocovpath,"/FULL.",covname,".txt",sep="")
    
    file.remove(cov.path)
    file.create(cov.path)
    mm[,2] <- as.character(mm[,2])
### want to create some formatting:
### mm[,1] <- formatC(mm[,1],digits=5,width=11,flag= " ", format="f")
### any of the FULL. files have three cols, the last is strata
### 
    if(ncol(mm) <= 2)
      mm[,3] <- NA* mm[,1]
    write.table(mm, file=cov.path,quote=F, sep="\t",eol="\n",row.names=F, col.names=F, na="NA")
       messout(paste("Size FULL.",covname, " is = ",file.info(cov.path)$size,sep=""),verbose) 
       
   return(invisible(Covs.lst)) }
## ************************************************************************
 find.age.groups <- function(alldths, lstind){
     
     todiscard <- lstind$digit.first
    
     cdigits <- lstind$cntry.digits
     cdigits <- cdigits + todiscard
     adigits <- lstind$age.digits
     ydigits <- lstind$year.digits
     csid <- names(alldths)
     age.groups <- sapply(csid,substr, cdigits+1, cdigits+adigits)
     age.groups <- sapply(age.groups, trim.blanks)
     age.groups <- unique.default(age.groups)
     age.groups <- as.numeric(age.groups)
     return(age.groups)
     
   }

complete.flags <- function(m, digits){
  tag <- m
  if(nchar(m) < digits){
    npad <- digits - nchar(m)
    tag <- paste("",rep(0,npad),collapse="", sep="")
    tag <- paste(tag,m,sep="")
  }
  return(tag)
}
### DESCRIPTION: Orders the dataobj input to yourcast by country code
###              Two elements the data and proximity list elements
###              Needed to forecast with Bayes model
###
###              evillalon@iq.harvard.edu
###              02/22/2007
###
sort.dataobj <- function(dataobj){
  dat  <- dataobj$data
  prox <- dataobj$proximity
###print(names(dat))
### fixing data first
  ord  <- order(names(dat))
  dataobj$data  <- dat[ord]
  if(length(dataobj$proximity) <= 0)
    return(dataobj)
### ordering proximity
  ord  <- order(prox[,1])
  prox <- prox[ord, ]
  rwnm <- apply(prox,1, function(rw){
    nm <- paste(rw[1], rw[2], sep="")
    return(nm)})
  
  rwnm <- unlist(rwnm)
  names(rwnm) <- NULL
  if(length(rownames(prox)) >0 && any(rownames(prox) != rwnm))
    message("Bad rownaming of proximity matrix")
 
  ord <- order(rwnm)
  rownames(prox) <- rwnm
  prox <- prox[ord, ]
  dataobj$proximity <- prox
  return(dataobj)
}
              

## ************************************************************************
## ************************************************************************
##
## FUNCTION NAME:      cntry.names
##
## PLACE:     Preprocessing module, the old Gauss cntry_names.g (and the used
##            token2.g) functions rewritten in R.
## 
## IMPORTED functions: none
##
## USED GLOBALS: whodatapath from env.who
##
## DESCRIPTION: It returns a vector with all the names of the 
##              countries in the database, indexed by their
##              country codes.
##
## FORMAT:  c.names <- cntry.names()
##          c.names <- cntry.names(file)
##
## INPUT:   file: (optional) string, name of the file with path.
##                The named file is containing pairs of
##                country code and country name. Each line contains 
##                two entries, separated by at least one blank space: 
##                an integer number (the WHO code) and the name 
##                of the corresponding country. If filename is missing
##                then default <whodatapath>/cntry.odes.txt is used.
##                The environment ebase=env.base, which contains all progams and
##                the environmnet= env.who with all globals of the WHO code  
##
## OUTPUT: names = Kx1 string array. K is the maximum country code.
##                 names[j] contains the name of the country whose
##                 country code is j or "" otherwise.
##                 For ex. names[1025] = "Benin", names[1]=""
##
## WRITTEN BY: Federico Girosi & Anita Gulyasne Goldpergel & Elena Villalon 
##             fgirosi@latte.harvard.edu, anitag@latte.harvard.edu
##             evillalon@latte.harvard.edu 
##             CBRSS, Harvard University
## 
## Last modified: 12/12/2003
## 
## ************************************************************************
## ************************************************************************

cntry.names <- function(filename=NULL, flag= F, ebase=parent.frame()) {
 
  if(length(filename) <= 0 ){  
    ebase <- get("env.base", envir= parent.frame())
    env.base <- ebase
    ewho <- get("env.who", envir=ebase)
    whodatapath <- get("whodatapath", envir=ewho)
    codes.names <- get("codes.names", envir=ewho)

    filename <-  paste(whodatapath,codes.names,sep="")
  }
### It uses the NEW file from the data directory!!!
### (Which is separated with one tabulator character)

### quote is neccessary because of Cote d'Ivoire for ex.
  tmp <- scan(file=filename,what=c(integer(0),character(0)),quiet=T,multi.line=T,
              sep="\t",quote="\"")
  tmp <- matrix(tmp,ncol=2,byrow=T)
  codes <- as.integer(tmp[,1])
  maxcode <- max(codes)
  names <- array("",dim=maxcode)
  names[codes] <- tmp[,2]
  return(names)
  }
