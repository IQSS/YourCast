### Module for Poisson regression model with glm: driver glm.poisson()
### Contains following functions
### glm.poisson() Poisson linear regression model build upon
### glm.colinear.poiss(xw, yw,yoff, nmx, delta.tol=0.005, tol=0.999), which
### checks for colinearities using covdel

## ************************************************************************
## ************************************************************************
##
## FUNCTION NAME: glm.poisson  
##
## INPUT:     whoinsampx, whoutsampx matrix of covariates for given csid 
##            whoinsampy, whoutsampy (1 col) matrices = ln(mortality) 
##            mortality = dth/population, and dth are observations
##            delta.tol decrement to supress colinearities among cols of whoinsampx;
##
## DESCRIPTION: same as ols but uses glm to obtain the estimates for dth's 
##
## IMPORTED functions:
##                      covdel
##                      glm.colinear.poiss
##
## OUTPUT:    list of matrix elements for each csid; yhat = ln(dth/popu) 
##            with predicted (log of mortalities) yhatin,for years of insample
##            and yhatout for years of outsample
##            The fitting coefficients coeff and standard deviation. 
##           
##
##  WRITTEN BY: Federico Girosi & Elena Villalon 
##              fgirosi@latte.harvard.edu
##              CBRSS, Harvard University
## 
## Last modified: 05/18/2003
## 
## ************************************************************************
## ************************************************************************

glm.poisson <- function(ebase=env.base){
  ebase <- get("env.base", envir=parent.frame())
  env.base <- ebase
  ewho <- get("env.who", envir=ebase)
  verbose <- get("verbose", envir=ebase)
  messout("Running Poisson with glm-modified", verbose)
  ttmp <- proc.time()

  whoinsampy <- get("whoinsampy", envir=ewho)
  whoutsampy <- get("whoutsampy", envir=ewho)
  whopopul   <- get("whopopul", envir=ewho)
  whopopulos <- get("whopopulos", envir=ewho)
  whoinsampx <- get("whoinsampx", envir=ewho)
  whoutsampx <- get("whoutsampx", envir=ewho)
  whocov   <- get("whocov", envir=ewho)
  whoyrest <- get("sample.frame", envir=ewho)[2]
  age.vec  <- get("age.vec", envir=ewho)
  delta.tol <- get("delta.tol", envir=ewho)
  tol <- get("tol", envir=ewho)
  cntry.names.lst <- try(get("cntry.names.lst", envir=ewho), silent=T) 
  log.p <- try(get("log.poiss", envir=ewho), silent=T)
  ff <- try(get("formula", envir=ewho), silent=T)
  ls <- as.character(ff)[2]
  ixdiv <- grep("/", ls)
  linx <- (length(ixdiv) > 0)
  
  ly <- length(whoinsampy)
  yhatin  <- whoinsampy;
  yhatout <- whoutsampy;
  coeff <- vector(mode="list",length= ly);
  std <- vector(mode="list",length=ly);  
  names(coeff) <- names(whoinsampy);
  names(std) <- names(whoinsampy);
  
### if you not wish to eliminate colinears set delta.tol <- 0
### experience shows that delta.tol =0.005 is a good choice to get rid
### of all NA's coefficients from glm delta.tol <- 0.005
### Set it equal to delta.tol of eliminating correlation, which is also 0.005

### for values of dth that might be zero, dth <- zero.tol
  zero.tol <- 1.e-4 
  for (i in 1:length(whoinsampy)){
### assuming whotransform = 1 so whoinsampy = ln(mortality)
    if(log.p)
      y <- exp(whoinsampy[[i]])  ### y is mortality= dths/popul
    else
      y <- whoinsampy[[i]]
    
    ypop  <-  whopopul[[i]]  ### population insample
    yoss  <- whopopulos[[i]] ### population outsample
    
### testing:
    
    y <- round(y * ypop - 0.5) 
    ylen <-  length(y)
    y[y <= zero.tol] <- 0
    yna <- na.omit(y)
    n.obs <- nrow(yna);
    yoff <- log(ypop)
### not to get many warnings from glm round to integers    
    ind <- na.omit(seq(along=y)[y > 0.5])
    y[ind] <- round(y[ind])
### check if too many 0 in data
    y0 <- yna[yna <= zero.tol]
    ly0 <- length(y0)
### obtain useful quantities for testing
    ymean   <- colMeans(as.data.frame(ypop), na.rm=T)
    ypopbar <- rep(ymean, length(y))
### for testing
### yoff <- matrix(1, nrow=length(ypop))
    
    x <- whoinsampx[[i]];
    n.cov <- length(colnames(x))
    if(length(ypop) != nrow(x)) messout("Wrong ypop and x must be of same length", verbose)
### for testing
    xout <- whoutsampx[[i]]
    yw <- y;
    xw <-  x;
    nmx <- names(whoinsampx[i])
    obj <- try(glm.colinear.poiss(xw, yw, yoff, nmx, delta.tol= 0, tol=tol, verb=verbose), silent=T)
    if(class(obj)=="try-error"){
      csid <- names(whoinsampy)[i]
      mess <- paste("poisson error for...", csid)
      stop(mess)
    }
### yx is only needed if working with glm.fit
### (note that glm calls glm.fit; so NA's are removed)
### yx  <- na.omit(cbind(y,x))
### yna <- yx[,1]
### xna <- yx[,-1]
### if using glm.fit with no colinear checking 
### obj <- glm.fit(xna,yna,family=poisson(), offset=yoff)
### offset= log(p) when p != 1
### Poisson(p exp(Z \beta)) =Poisson(exp( Z \beta + log(p))
    ri <- obj$covind
    if (length(ri) > 0 ) {
        xout <- xout[,-ri]
        x <- x[,-ri]
        whocov[[i]] <- whocov[[i]][,-ri]
        
      }
        beta <- obj$coeff
### if too little data set beta = NA
### define degree of freedoms: dg.f
    
    dg.f <- n.obs - (n.cov + 1)
    if (( dg.f > 0 && (n.cov/dg.f > (1 - ly0/dg.f))) || dg.f <= 0)
      {
        messout(paste("Warning: for csid: ", names(whoinsampy[i]), sep=""),verbose)
        messout(paste("n.dth zeros= ",ly0, "; n.obs= ", n.obs,"; n.cov= ", n.cov, sep=""), verbose)
        beta <- matrix(,nrow=n.cov)
      }
    
    coeff[[i]] <- beta
    
    dthatin  <-  exp(x %*% beta) ### dth estimation for insamp
    dthatout <-  exp(xout %*% beta)  ### dth estimation for outsample

    if(log.p){      
      yhatin[[i]]  <- x %*% beta   ### log(mortality) insample
      yhatout[[i]] <- xout %*% beta  ### log(mortality) outsample
      std[[i]] <- sqrt(crossprod(na.omit(yhatin[[i]] - whoinsampy[[i]]))/(n.obs - 1));
      
    } else{
      yhatin[[i]]  <- dthatin   ### mortality insample
      yhatout[[i]] <- dthatout  ### mortality outsample
      std[[i]] <- sqrt(crossprod(na.omit(yhatin[[i]] - whoinsampy[[i]]))/(n.obs - 1));
        
    }
     
### take care of new globals after deleting colinears
    whoinsampx[[i]] <- x
    whoutsampx[[i]] <- xout
    }
###   print("Time spent with ols-glm: "); print(proc.time() - ttmp)
  rm(ttmp)
   assign("whocov", whocov, envir= ewho)
   assign("whoinsampx", whoinsampx, envir= ewho)
   assign("whoutsampx", whoutsampx, envir=ewho)
   model <- model.string()  
   lst <- list(yrest=whoyrest,
              model=model, age.vec=age.vec, cntry.lst=cntry.names.lst,
              coeff=coeff,yhatin=yhatin,yhatout=yhatout,std=std,
               insampy=whoinsampy,outsampy=whoutsampy)
   assign("lst.output", lst, envir=ewho)
   return(lst)
        
}
    

## ************************************************************************
## ************************************************************************
##
## FUNCTION NAME: glm.colinear.poiss 
##
## INPUT:     whoinsampx matrix of covariates for given csid 
##            with rows = years of observations in the insample matrix
##            and no. of cols= the covariates;
##            whoinsampy matrix of dth's observations for each csid
##            delta.tol tolerance decrement and tol factor for linear dependencies.
##
## DESCRIPTION: glm finds coeff with NA's (because of linear dependencies);  
##              i)  set is.na(coeff) = 0 if delta.tol <= 0; or, alternatively,
##              ii) if delta.tol > 0, call covdel which eliminates colinearities
##                  with values of tol that are decreased by delta.tol,
##                  until colinears (or NA's coefficients) are removed
##
## OUTPUT:    matrix whoinsampx modified,if necessary, eliminating colinearities 
##            the coefficients obtained with glm; and 
##            tol after decreasing in steps delta.tol
##           
##
##  WRITTEN BY: Elena Villalon
##              evillalon@latte.harvard.edu
##              CBRSS, Harvard University
##
## Last modified: 05/18/2003
## 
## ************************************************************************
## ************************************************************************
 glm.colinear.poiss <- function(xw, yw,yoff, nmx, delta.tol=0.0, tol=0.999, verb=T){
    n   <- floor( tol/delta.tol)
    ri  <- vector(,length=0)
    xw <- as.matrix(xw)
    x <- xw
    verbose <- verb
### options(show.error.messages=F)
    res <- try(glm(yw~x-1, family=poisson(), na.action=na.omit,offset=yoff),
               silent=T)
    
    if(length(res) <= 1 && class(res)=="try-error")
      messout("glm failure...",verbose)
      
    coeff <- res$coeff
    isna  <- seq(along= coeff)[is.na(coeff)]
    isnames <- names(coeff[isna])
      
  test <-  1
  if (length(isna) > 0 && delta.tol <= 0){  ### if1
    coeff[isna] <- 0
  }else  if (length(isna) > 0 && delta.tol >0){  ### if1
    messout(paste("Eliminating colinearities for glm-gaussian and csid= ", nmx, sep=""),verbose)
       while(n > 0){  ### while
           tol <- tol - delta.tol
           ret <- covdel(xw, tol)
           x   <- as.matrix(ret$covx)
           ri   <- ret$covind
           rn   <- ret$covnam
           rnx <- paste("x", rn, sep="")
           if (!any(is.element(rnx,names(coeff))))
              stop(paste("Wrong covariates deleted for csid= ", nmx))
           
           if((tol - delta.tol) > 0  &&  (length(rn) < length(isna))){  ### if2
                   n <- n - 1
### may eliminate more than (isna) with glm, if delta.tol is large
           }else if ((tol - delta.tol) > 0 && (length(rn) >= length(isna))) { ### if2
              res <- glm(yw~x-1, family=poisson(), na.action=na.omit,offset=yoff)
              coeff <- res$coeff
              n   <-  0
           }else {n <-  0
                  test <- F}  ### if2
               
            } ### while(n >0)
        } ### if1
        
      
    if (test==F)  
         messout(paste("Correlation not eliminated for glm-gauss csid= ", nmx, sep=""),verbose)
   
    nmc <- names(coeff)
### to remove x from glm on the names of covariates
    nmc <- substring(nmc, 2)
    coeff <- matrix(coeff)
    rownames(coeff) <- nmc 
    options(show.error.messages=T)
    return(list(covind= ri, coeff=coeff, tol= tol))
  }
      
 
   
  
    
