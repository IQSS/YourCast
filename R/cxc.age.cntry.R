  ##################################################################
##
## FUNCTION NAME:  cxc() & cxc.model() 
##
## INPUT: environments with globals and preprocessing matrices
##        user choices for the priors
##  
##
## IMPORTED:  globals and environmnets
##
## DESCRIPTION: Computes the beta for any of the priors age, time, age-time 
##              and for the cntry prior.  Estimates mortality rates
##              for insample and outsample according to model predictions)
##               
##
## OUTPUT : a list with the prediction of the models & and the environmnet
##          env.cxc with all matrices and paramters needed for gibbs model. 
##
## WRITTEN BY: Federico Girosi    
##             girosi@rand.org;  
##             CBRSS, Harvard University
##             (modified by EV to include cross-cntry smoothing)
## 
## Last modified: 10/01/2003
## 
## ***************************************                      
    

cxc <- function(Hct.only=F,  ebase=env.base){

 
  ebase <- get("env.base", envir=parent.frame())
  env.base <- ebase
  ewho <- get("env.who", envir=ebase)
  whoinsampx <- get("whoinsampx", envir=ewho)
  whoutsampx <- get("whoutsampx", envir=ewho)
  who.zero.mean <- get("who.zero.mean", envir=ewho)
  whoinsampy <- get("whoinsampy", envir=ewho)
  whoutsampy <- get("whoutsampy", envir=ewho)
  verbose <- get("verbose", envir=ebase)
  insampy  <- whoinsampy
  outsampy <- whoutsampy
  model <- model.string(); 

### the following are needed to output it yourcast

  whoyrest <- get("whoyrest", envir=ewho)
  age.vec <- get("age.vec", envir=ewho)
  cntry.names.lst <- get("cntry.names.lst", envir=ewho)
  whomodel <- get("whomodel", envir=ewho)
  whomodel <- toupper(trim.blanks(whomodel))
   
  ix <- NULL
##  if(identical(whomodel, "MAP") || identical(whomodel,"EBAYES"))
    messout("Running MAP model", verbose);   
### this is the main call, which computes the coefficients and other quantities of interest
    ecxc <- cxc.model(Hct.only,env.base);

    coeff <-  get("coeff", envir=ecxc)
      S.list <- get("S.list",envir=ecxc)
  
    ix <- check.coeff.na(coeff)
   
    ix <- unlist(ix)
  
  
   if(length(ix) <= 0){
    
     return(lst <- list(yrest=whoyrest,model=model,age.vec=age.vec,
                        cntry.lst=cntry.names.lst, 
                        coeff=list(),yhatin=list(),yhatout=list(),
                        insampy=insampy,outsampy=outsampy,ecxc=ecxc,std=NULL))
   }
  if(length(ix) > 0){
    coeff <- coeff[ix]
    whoinsampx <- whoinsampx[ix]
    whoutsampx <- whoutsampx[ix]
}
### now we make the forecasts of the dependent variable
### coeff <- res$coeff;
  yhatin <- make.forecast(coeff,whoinsampx);
  yhatout <- make.forecast(coeff,whoutsampx);  
  
### if the prior was not zero mean the dependent variable is not what
### we wish to forecast (it has zero mean): we have to add the mean age profile back to it

 if (!is.logical(who.zero.mean)) {
    who.mean.age.profile <- who.zero.mean;
    func <- function(x,param){x+param};
    yhatin <- modify.age.profiles(yhatin,func,who.mean.age.profile)
    yhatout <- modify.age.profiles(yhatout,func,who.mean.age.profile)  
   
  }else{
    if(!who.zero.mean){
      who.mean.age.profile <- make.average.age.profile(whoinsampy)
      func <- function(x,param){x+param};
      yhatin <- modify.age.profiles(yhatin,func,who.mean.age.profile)
      yhatout <- modify.age.profiles(yhatout,func,who.mean.age.profile)  
         }
  }
        
   ### here we compute the standard deviations of the predicted values
  if(identical(whomodel, "MAP"))
    std <- std.MAP(S.list)
  else
    std <- NA
  
  model <- model.string(); 
  lst <- list(yrest=whoyrest,model=model,age.vec=age.vec,
              cntry.lst=cntry.names.lst, 
              coeff=coeff,yhatin=yhatin,yhatout=yhatout,
              insampy=insampy,outsampy=outsampy, ecxc=ecxc, std=std);
### make the assignment to ewho to build output
  assign("lst.output", lst, envir=ewho)
###  lstout <<- lst
  return(lst);
}



######################################################################
######################################################################
######################################################################

cxc.model <- function(Hct.only=F, ebase=env.base){
  ebase <- get("env.base", envir=parent.frame())
  env.base <- ebase
  ewho <- get("env.who", envir=ebase)
  age.vec <- get("age.vec", envir=ewho)
  n.age <- length(age.vec)
  
  svdonly <- get("svdonly", envir=ewho)
  who.Ha.time.weight <- get("who.Ha.time.weight", envir=ewho)
  who.Ha.deriv <- get("who.Ha.deriv", envir=ewho)
  who.Ha.age.weight <- get("who.Ha.age.weight", envir=ewho)
  who.Hat.t.deriv <- get("who.Hat.t.deriv",envir=ewho)
  who.Hat.time.weight <- get("who.Hat.time.weight", envir=ewho)
  who.Hat.a.deriv <- get("who.Hat.a.deriv", envir=ewho)
  who.Hat.age.weight <- get("who.Hat.age.weight", envir=ewho)
  who.Ht.deriv <- get("who.Ht.deriv", envir=ewho)
  who.Ht.time.weight <- get("who.Ht.time.weight", envir=ewho)
  who.Ht.age.weight <- get("who.Ht.age.weight", envir=ewho)
  who.Hct.t.deriv <- get("who.Hct.t.deriv", envir=ewho)
  who.Hct.time.weight <- get("who.Hct.time.weight", envir=ewho)
  whoinsampy <- get("whoinsampy", envir=ewho)
  who.zero.mean <- get("who.zero.mean",envir=ewho)
  whoinsampx <- get("whoinsampx", envir=ewho)
 verbose <- get("verbose", envir=ebase)
  cntry.vec <- get("cntry.vec", envir=ewho)
  whocov <- get("whocov", envir=ewho)
  who.age.digits <- get("who.age.digits", envir=ewho)
  who.cntry.digits <- get("who.cntry.digits", envir=ewho)
  whomodel <- get("whomodel", envir=ewho)
  if(identical(toupper(whomodel),"EBAYES"))
    whomodel <- "MAP"
  
  who.Ha.sigma <- get("who.Ha.sigma", envir=ewho)
  who.Hat.sigma <- get("who.Hat.sigma", envir=ewho)
  who.Ht.sigma <- get("who.Ht.sigma", envir=ewho)
  who.Hct.sigma <- get("who.Hct.sigma", envir=ewho)
  param <- list(Ha.sigma=who.Ha.sigma,Hat.sigma=who.Hat.sigma,
                Ht.sigma=who.Ht.sigma, Hct.sigma=who.Hct.sigma);
  Ha.theta <- 0
  Ha.sigma <- param$Ha.sigma;
  if (is.na(Ha.sigma)) {Ha.theta <-  0;
  } else { Ha.theta <- 1/Ha.sigma^2}

  Hat.theta <- 0
  Hat.sigma <- param$Hat.sigma;
  if (is.na(Hat.sigma)) {Hat.theta <-  0;
  } else { Hat.theta <- 1/Hat.sigma^2}
  
  Ht.theta <- 0
  Ht.sigma <- param$Ht.sigma;
  if (is.na(Ht.sigma)) {Ht.theta <-  0;
  } else { Ht.theta <- 1/Ht.sigma^2}
  
  age.prior <- (!is.na(who.Ha.sigma) || !is.na(who.Ht.sigma) || !is.na(who.Hat.sigma))

  
  Hct.sigma <- param$Hct.sigma;
  if (is.na(Hct.sigma)) {Hct.theta <-  0;
  } else { Hct.theta <- 1/Hct.sigma^2}
### cxc does not need to calculate cross cntry smoothing even
### if we have Hct.sigma != NA and we calculate it in the Gibbs
   if( age.prior == F && Hct.theta != 0)
### indicates that cntry prior is the only prior
    Hct.only <- T

### Here we compute basic quantities used by the prior
if(Ha.theta != 0){    
  messout("Preparing for smoothing over age groups",verbose);
  time.prior.param <- list(time.der=1,time.weight=who.Ha.time.weight)
 
  W.age <- derivative.prior(n.age,who.Ha.deriv,who.Ha.age.weight);
 
  who.C.age <- build.C.age.time(W.age,time.prior.param, env.base);
 
  
 }
###  assign("W.age", W.age, envir=ewho)
###  assign("who.C.age", who.C.age, envir=ewho)
 if(Hat.theta != 0){ 
  messout("Preparing for smoothing of time trend over age groups", verbose);  
  time.prior.param <- list(time.der=who.Hat.t.deriv,time.weight=who.Hat.time.weight)
  W.age.time <- derivative.prior(n.age,who.Hat.a.deriv,who.Hat.age.weight)
  who.C.age.time <- build.C.age.time(W.age.time,time.prior.param, env.base);}

###  assign("W.age.time", W.age.time, envir=ewho)
###  assign("who.C.age.time", who.C.age.time, envir=ewho)
  if(Ht.theta != 0){
  messout("Preparing for smoothing over time", verbose);  
  time.prior.param <- list(time.der=who.Ht.deriv,time.weight=who.Ht.time.weight)
  W.time <- derivative.prior(n.age,1,who.Ht.age.weight);
  who.C.time <- build.C.age.time(W.time,time.prior.param, env.base);}

###  assign("W.time", W.time, envir=ewho)
###  assign("who.C.time", who.C.time, envir=ewho)
 
if(length(cntry.vec) && Hct.theta !=0 ){
  messout("Preparing for smoothing of time trend over countries", verbose);
  time.prior.param <- list(time.der=who.Hct.t.deriv,time.weight=who.Hct.time.weight)    
### omega.cntry <- build.adjacency(cntry.vec)
### W.cntry <- laplacian.from.adjacency(cntry.vec)
### W.cntry goes into calculations for who.C.cntry;
### but not needed for build.C.cntry.time,which calls laplacian.from.adjacency
### inside the function anyway.
 
  res <- build.C.cntry.time(cntry.vec,time.prior.param, env.base);
 
  if (length(res) > 0){
    who.C.cntry <- res$who.C.cntry
    W.cntry <- res$W.cntry}
  else{
    who.Hct.sigma <- NA
    Hct.theta <- 0}
 
}
 
###
### First we subtract the mean age profile if requested
  insampy <- whoinsampy;
  
  if (!is.logical(who.zero.mean)) {
    who.mean.age.profile <- who.zero.mean
    func <- function(x,param){x-param};
    insampy <- modify.age.profiles(whoinsampy,func,who.mean.age.profile)
###   print(who.zero.mean)
  }else{
    if(!who.zero.mean){
    
      who.mean.age.profile <- make.average.age.profile(whoinsampy,ebase=env.base)
      func <- function(x,param){x-param};
      insampy <- modify.age.profiles(whoinsampy,func,who.mean.age.profile)
    
    }
  }
  
  
### sigma is a list like whoinsampy

  ols.result <- ols(env.base);

  sigma <- estimate.standard.deviations(ols.result, ebase=env.base)
              
  
  sigma.ols <- sigma
 
  
  indx <- seq(1,length(insampy));
    names(indx) <- names(insampy);
### assign("ols.result", ols.result, envir=ewho)
### assign("sigma.ols", sigma, envir=ewho)
### likelihood.weights is a list like insampy, with the weights
### of each individual observation. Missing observations have 0 weight
  likelihood.weights <- lapply(indx,FUN=function(n,insampy,sigma){y <- insampy[[n]];
                                                    s <- 1/sigma[[n]];
                                                    s[is.na(y)] <- 0;
                                                    return(s)},insampy,sigma)
  likelihood.weights <- lapply(likelihood.weights, as.vector)
### Y is a list like insampy. It contains the dependent variable weighted with
### likelihood.weights

  Y <- lapply(indx,FUN=function(n,likelihood.weights,insampy){y <- insampy[[n]];
                                                      w <- likelihood.weights[[n]];
                                                      y[is.na(y)] <- 1;
                                                      y <- y*w;
                                                      return(y)},likelihood.weights,insampy);

### X is a list like whoinsampx.It contains insampx weighted with
### likelihood.weights

  X <- lapply(indx,FUN=function(n,likelihood.weights){x <- whoinsampx[[n]];
                                   w <- likelihood.weights[[n]];
                                   x <- w*x;
                                   return(x);
                                 },likelihood.weights);

### csid lists; each element= cntry+age  
  beta <- as.list(rep(NA,length(insampy)));
  names(beta) <- names(insampy);
  beta.hat.csid.lst  <- ols.result$coeff;
###  ols.beta  <<- ols.result$coeff
  
### cntry lists; each element=cntry
  clist <- vector(mode="list",length=length(cntry.vec));
  names(clist) <- as.character(cntry.vec);
  Ha.wc.list <- clist;
  Hat.wc.list <- clist;
  Ht.wc.list <- clist;
  D.list <- clist;
  S.list <- clist;
  beta.hat.list <- clist;
  D.age.lst <- clist;
  D.age.time.lst <- clist;
  D.time.lst <- clist;

  
#######################################################################################

#################### HERE WE START TO LOOP OVER COUNTRIES #############################

#######################################################################################
  
  
  age.char <- formatC(age.vec, width=who.age.digits, format="d", flag="0")
  age.prior <- (!is.na(who.Ha.sigma) || !is.na(who.Ht.sigma) || !is.na(who.Hat.sigma))
  beta.dim.all <- list()
  t.list.all <- list()
  if (age.prior ){## do not loop if all age, time, age-time priors are
    
  for (i in 1:length(cntry.vec)){
   
    cntry <-  cntry.vec[i];
   
    cntry.str <- as.character(cntry);
###    cat(cntry, "\n")
    
### beta.dim is a list like age.vec. Each element is the number of covariates
### in the correspoding age group
    
    nam <- paste(cntry,age.char,sep="");
    beta.dim  <- sapply(nam, function(n){nc <- ncol(whocov[[n]])});
    beta.dim.all[[i]] <- beta.dim 
    cov.name.list <-  lapply(nam, function(n){nc <- colnames(whocov[[n]])});
    names(cov.name.list) <- nam
   
### t.list is a list like age.vec. Each element is the length of the time series
### of the covariates in whocov in the correspoding age group
    t.list  <- sapply(nam, function(n){nc <- nrow(whocov[[n]])});
    t.list.all[[i]] <- t.list 
### We are assuming that the time series in whocov have, within each country, the same length
### for different age groups. If this is not the case the average length is taken.
### This affects the computation of wc.
    
### All of the sudden this did not work!!!  
###   if(!all.equal.numeric(t.list, rep(t.list[1],length(t.list))))
###      message("Lenght of covariates differ in different age groups");
  
### LIKELIHOOD and DATA MATRICES
    Q <- cntry.Xprime.X(cntry,X, env.base);
   
    v <- cntry.Xprime.Y(cntry,X,Y,env.base);
   
    ZZ <- cntry.Xprime.X(cntry,whocov,env.base);
  
### AGE PRIOR
###
  
    if( Ha.theta != 0){
     
      D.age <- make.age.time.prior.matrix(cntry,age.vec,W.age,who.C.age,env.base);
### this is the covariance of the improper prior (the pseudoinverse of D)
      D.age.pinv <- sample.improper.normal(D.age,1,1)$covar; 
      Ha.wc <- sum(diag(D.age.pinv%*%ZZ))/(length(age.vec)*mean(t.list));
      Ha.wc.list[[cntry.str]] <- Ha.wc;
      D.age.lst[[cntry.str]] <- D.age
      cola <-  colnames(D.age)
     
      age.block <- lapply(age.char,function(ch, cola){
            ch  <- paste("^", ch, sep="")
            ind <- grep(ch, cola)
            return(ind)},cola)
      names(age.block) <- age.char
    }
    
### AGE TIME PRIOR
###
  
    if (Hat.theta != 0){
      D.age.time <- make.age.time.prior.matrix(cntry,age.vec,W.age.time,who.C.age.time, env.base);
      D.age.time.pinv <- sample.improper.normal(D.age.time,1,1)$covar;
### this is the covariance of the improper prior
      Hat.wc <- sum(diag(D.age.time.pinv%*%ZZ))/(length(age.vec)*mean(t.list));
      Hat.wc.list[[cntry.str]] <- Hat.wc;
      D.age.time.lst[[cntry.str]] <- D.age.time
      colat <-  colnames(D.age.time)
      age.time.block <- lapply(age.char,function(ch, colat){
            ch  <- paste("^", ch, sep="")
            ind <- grep(ch, colat)
            return(ind)},colat)
      names(age.time.block) <- age.char
    }
    
### TIME PRIOR
###
  
    if (Ht.theta != 0){
      D.time <- make.age.time.prior.matrix(cntry,age.vec,W.time,who.C.time, env.base);
      D.time.pinv <- sample.improper.normal(D.time,1,1)$covar;
### this is the covariance of the improper prior
      Ht.wc <- sum(diag(D.time.pinv%*%ZZ))/(length(age.vec)*mean(t.list));    
      Ht.wc.list[[cntry.str]] <- Ht.wc;
      D.time.lst[[cntry.str]] <- D.time
      colt <-  colnames(D.time) 
      time.block <- lapply(age.char,function(ch, colt){
            ch  <- paste("^", ch, sep="")
            ind <- grep(ch, colt)
            return(ind)},colt)
      names(time.block) <- age.char
    }
### indices along columns or rows from D.age, D.age.time and D.time
### that select a given age group (i.e. names(d.block)) correlation
### blocks indices to be filled with matrices from who.C.cntry
    if( Ha.theta != 0)
      d.block <- age.block
    else if( Hat.theta != 0)
      d.block <- age.time.block
    else
      d.block <- time.block
  
   d.block <- sapply(names(d.block), function(x){
      a  <- NULL
      t  <- NULL
      at <- NULL 
     if(Ha.theta  != 0)
      a <- age.block[[x]]
     if(Ht.theta  != 0)
      t <- time.block[[x]]
      if(Hat.theta  != 0)
        at <- age.time.block[[x]]
      return <- unique.default(c(a,t,at))})
    
### COMPUTATION OF COEFFICIENTS
    D <- 0
   
    if(Ha.theta != 0)  
      D <- D+ Ha.theta*Ha.wc*D.age
    if( Hat.theta != 0)
      D <- D+ Hat.theta*Hat.wc*D.age.time
    if(Ht.theta != 0)
      D <- D+ Ht.theta*Ht.wc*D.time
       
    beta.hat <- NA*v
    S <- NULL
    QD <- Q + D
    if (!all(is.na(QD)) && svdonly <= 0)
      {
        beta.hat <- try(solve(QD,v),silent=T)
        if(inherits(beta.hat, "try-error"))
          beta.hat <- try(solve(QD,v,LINPACK=T),silent=T)
        
        if(inherits(beta.hat, "try-error")){
          messout(paste("Cannot invert matrix with solve code= ", cntry.str,
                        "\n", "Using singular value descomposition", sep=""), verbose)
        
          S <- svd.inv(QD)
          beta.hat <- S%*% v
        }
       
      
      }else if (!all(is.na(QD)) && svdonly <= 1){
        S <- svd.inv(QD)
       
        beta.hat <- S%*% v
      }else if(!all(is.na(QD)) && svdonly <= 2){
        message("Generalized inverse not implemented yet")
        
      }
    D.list[[cntry.str]] <- D;
    S.list[[cntry.str]] <- S;    
    beta.hat.list[[cntry.str]] <- beta.hat;
    
### now we need to partition beta.hat in blocks, one for each age group.
    beta.list <- split.matrix(beta.hat, beta.dim);
    names(beta.list) <- paste(cntry,age.char,sep="");
    nn <- names(beta.list)
    beta.list <- lapply(names(beta.list),FUN=function(x,beta.list){b <- beta.list[[x]];
                                                                   rownames(b) <- cov.name.list[[x]];
                                                                   b},beta.list)
    beta.list <- lapply(1:length(age.char),FUN=function(x,beta.list){b <- beta.list[[x]];
                                                                   colnames(b) <- age.char[x];
                                                                   b},beta.list)
    names(beta.list) <- nn
    beta[names(beta.list)] <- beta.list;
    coeff <- beta;
    if(Ha.theta != 0 || Ht.theta != 0 || Hat.theta != 0)
      beta.hat.csid.lst <- beta;
  }# end for loop over cntry
} ## if (any age, time, age-time Ha.theta, Ht.theta, Hat.theta !0)
  env.cxc <- environment();
  assign("env.cxc", env.cxc, envir=env.base)  
  
    D.cntry.lst <- as.list(1:length(age.vec));
    names(D.cntry.lst) <- age.char
    Hct.wc.lst <-  as.list(1:length(age.vec));
    names(Hct.wc.lst) <-  age.char;
## beta's are for each age since cntry are mixed up
    beta.hat.cntry.lst <-  as.list(1:length(age.vec));
    names(beta.hat.cntry.lst) <- age.char;
    csid.vec <- kronecker(as.character(cntry.vec), age.char, paste,sep="")
### beta.hat.csid.lst <- vector(mode="list",length=length(csid.vec))
### names(beta.hat.csid.lst) <- csid.vec
  
  if(length(cntry.vec) && Hct.theta !=0 ){

    for(i in 1:length(age.vec)){
      output   <- make.cntry.time.prior.matrix(age.vec[i],cntry.vec,W.cntry, who.C.cntry,env.base);
      Di <- output$D;
      isle.cntry <- output$isle.cntry
      D.cntry.lst[[age.char[i]]] <- Di;
      roli <-  rownames(Di) 
      cntry.block <- lapply(as.character(cntry.vec),function(ch, roli){
            ch  <- paste("^", ch, sep="")
            ind <- grep(ch, roli)
            return(ind)},roli)
      names(cntry.block) <- as.character(cntry.vec)
 
   Qcntry  <- age.Xprime.X(ag=age.vec[i],Xmat=X,isle.cntry =isle.cntry, ebase=env.base)
   vcntry  <- age.Xprime.Y(ag=age.vec[i],X=X,Y=Y,isle.cntry=isle.cntry, ebase=env.base)
   ZZcntry <- age.Xprime.X(ag=age.vec[i],Xmat=whocov,isle.cntry=isle.cntry, ebase=env.base);

### this is the covariance of the improper prior (the pseudoinverse of D)
   Di.pinv <- sample.improper.normal(Di,1,1)$covar;

   inc <- na.omit(match(isle.cntry, cntry.vec))
   cntry.related <- cntry.vec
   if (length(inc) > 0)
     cntry.related <- cntry.vec[-inc]
   island.char <- paste("^", isle.cntry,sep="")
      
   namc <- paste(as.character(cntry.related), age.char[i], sep="")

   cov.mat <- rmv.cntry(isle.cntry=isle.cntry, Xmat=whocov)

   tt.list  <- sapply(namc, function(n){nc <- nrow(cov.mat[[n]])});
     if ( class( try(Di.pinv%*%ZZcntry) )== "try-error"){
       messout(paste('Dimensions of Di.pinv\t ', dim(Di.pinv), sep=""), verbose)
       messout(paste('Dimensions of ZZcntry\t ', dim(ZZcntry), sep=""), verbose)
       messout(namc, verbose);
     }
       
    Hct.wc <- sum(diag(Di.pinv%*%ZZcntry))/(length(cntry.related)*mean(tt.list));
    Hct.wc.lst[[age.char[i]]] <- Hct.wc;
    Dcntry <-  Hct.theta * Hct.wc *Di
    ###############################################################
       beta.hat.cntry <- NA*vcntry
      if(!all(is.na(Qcntry + Dcntry))){
        
        if(class(try(solve(Qcntry + Dcntry) %*% vcntry,silent=T))=="try-error"){
          messout("try-error for cntry smoothing",verbose)
          fact <- try(svd.inv(Qcntry + Dcntry), silent=T)
          if(class(fact) != "try-error"){
 
            beta.hat.cntry <- fact %*% vcntry
            rownames(beta.hat.cntry) <- rownames(vcntry)
            colnames(beta.hat.cntry) <- age.char[i]
    
          }
        }else
        {
 
          beta.hat.cntry <- solve(Qcntry + Dcntry) %*% vcntry
          colnames(beta.hat.cntry) <- age.char[i]
                   
        }
      }
   
      beta.hat.cntry.lst[[age.char[i]]] <- beta.hat.cntry
     
###   print(beta.hat.cntry.lst[[age.char[i]]])
      
      limb <- make.beta.cn.lst(beta.hat.cntry, age.char[i],as.character(cntry.related))
      nmbeta <- rownames(beta.hat.cntry)
   
    
      beta.hat.related.csid <- lapply(limb, function(x, beta.hat.cntry, nmbeta){
        f <- x[1]
        s <- x[2]
        
        if(!is.na(f) && !is.na(s)){
          nms <- nmbeta[f:s]
         
          nms <- sapply(nms, function(nm){
            if(length(grep("\\.",nm)) >0)
               return(nm <-  strsplit(nm, "\\.")[[1]][2])
             
            return(nm <- substring(nm,who.cntry.digits+1))
            
          })
          nms <- unlist(nms)
          res <- as.matrix(beta.hat.cntry[f:s])
          rownames(res) <- nms
          return(res)
        }else
        return(as.matrix(beta.hat.cntry))
      }, beta.hat.cntry, nmbeta)

        
    nm <- names( beta.hat.related.csid)
    if(Hct.only)
    beta.hat.csid.lst[nm] <- beta.hat.related.csid
      
    }}
  coeff <- beta.hat.csid.lst
  ln <- 1:length(coeff)
  names(ln) <- names(coeff)
  ln <- sapply(ln, function(n,coeff) {
    x <- coeff[[n]]
    if (all(is.na(x))){
      messout(paste("Missing values for cntry code ", names(coeff)[n],sep=""), verbose)
      return(NULL)}
    else
      return(names(coeff)[n])}, coeff)
  ln <- unlist(ln)
  if(length(ln) <= 0)
    messout("Data set of countries coeff return NA's for MAP model...pls check your data",verbose)
  else if (length(ln) < length(coeff) ){
   
    mc <- setdiff(names(coeff), names(ln))
    
    stop(paste("Subset countries coeff all NA's ", names(mc), " pls check your data"))
  }
 
#################### DONE WITH LOOP OVER COUNTRIES #############################
###  return(list(coeff=beta,C.age=who.C.age,C.time=who.C.time,C.age.time=who.C.age.time,
###              Ha.wc.list=Ha.wc.list,Hat.wc.list=Hat.wc.list,Ht.wc.list=Ht.wc.list,D.list=D.list,beta.hat.list=beta.hat.list));
 
  return(env.cxc);}

 check.coeff.na <- function(coeff){
    ix <- 1:length(coeff)
    names(ix) <- names(coeff)
    ix <- sapply(ix, function(n, coeff)
                 {
                   vec <- coeff[[n]]
                   vec <- na.omit(vec)
                   if(length(vec) <=0)
                   return(list())
                   else
                     return(n)
                 }, coeff)
    return(ix)
  }
  
