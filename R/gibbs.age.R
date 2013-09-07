###################################################################
##
## FUNCTION NAME:  gibbs.age.cnst
##
## PLACE:    
##
## IMPORTED:    covariates matrices whocov from make.mortality.data()
##              environments with globals, including countries and age groups 
##
## DESCRIPTION: It implements the gibbs sampling algorithm for the betas and
##              time, age ,age-time priors. It calculates some of the building blocks
##              for the matrices Lambda. Whatever is in this function is constant
##              independent on the number of samples.Thus, you have to calculate it
##              only one and re-use it as you sample over the betas, sigmas and thetas.
##             
## INPUT:  environmnets with globals and the results from the cxc models for betas
##         and some ols results for sigma (or std), preprocessing matrices and lists,    
##         smoothing parameters and the C.matrices from posteriors.
##         Many building blocks from previous programs of prior and posterior calculations.   
##        
## OUTPUT: an environment with Ltheta.agprior.csid.lst,  Ltheta.tmprior.csid.lst, 
##         Ltheta.agtmprior.csid.lst.These are the Lamda matrices for age, time, age-time,
##         for every csid and without specific theta and sigma dependent contributions. 
##
## WRITTEN BY: Elena Villalon & Federico Girosi    
##             evillalon@latte.harvard.edu;
##             girosi@rand.org; 
##             CBRSS, Harvard University
## 
## Last modified: 05/01/2004
## 
## ************************************************************************
## *************************************************************************

gibbs.age.cnst <- function( ebase=env.base){
  ebase <- get("env.base", envir=parent.frame())
  env.base <- ebase
  ewho <- get("env.who", envir=ebase)
  ecxc <- get("env.cxc", envir=ebase)
## if the priors are not used then who.Ha.sigma = NA; same for all others
  who.Ha.sigma  <- get("who.Ha.sigma", envir=ewho)
### the average standard deviation of the prior. If NA the prior is not ### used
### (it is like having an infinite standard deviation )
  who.Hat.sigma <- get("who.Hat.sigma", envir=ewho)
  who.Ht.sigma  <- get("who.Ht.sigma", envir=ewho)
 
  
  age.prior <- ( !is.na(who.Ha.sigma) || !is.na(who.Ht.sigma) || !is.na(who.Hat.sigma))
  if ( age.prior == F )
    return (0);
                 
############################ PRIOR OVER AGE AND TIME GROUPS ###############################
######
  age.vec <- get("age.vec", envir=ewho)
  n.age <- length(age.vec)
  param <- list(Ha.sigma=who.Ha.sigma,Hat.sigma=who.Hat.sigma,Ht.sigma=who.Ht.sigma)
  whoinsampy <- get("whoinsampy", envir=ewho)
  whoinsampx <- get("whoinsampx", envir=ewho)
  cntry.vec  <- get("cntry.vec", envir=ewho)
  whocov     <- get("whocov", envir=ewho)
  who.age.digits <- get("who.age.digits", envir=ewho)
  who.cntry.digits  <- get("who.cntry.digits", envir=ewho)
  who.digit.first   <- get("who.digit.first", envir=ewho)
  digit.cntry.end   <- who.digit.first + who.cntry.digits; 
  digit.cntry.begin <- who.digit.first + 1
  digit.age.begin <- digit.cntry.end + 1
  digit.age.end   <- digit.cntry.end + who.age.digits
### cxc model produce basic matrices for computation and the coeffs
### now ols.result is also stored in ewho, besides being 
### stored also in the environmnet ecxc
### ols.result <- ols(env.base);
### cxc.result is an environmnet
### sigma.ols is stored in ewho and also in ecxc environments
### sigma.ols <- estimate.standard.deviations(ols.result,ebase= env.base)
### let us start with OLS estimation
### some storage lst for sample beta or coeff
### for who.C.time elements are cntry+a+a, only same age are correlated
###  env.gibbs.age <- environment();
###  assign("env.gibbs.age", env.gibbs.age, envir=ebase)

  if(!is.na(who.Ha.sigma)){
    W.age <- get("W.age", envir=ecxc)
    who.C.age <- get("who.C.age", envir=ecxc)
    D.age.lst <- get("D.age.lst", envir=ecxc)
    Ha.wc.list <- get("Ha.wc.list", envir=ecxc)}
    
  if(!is.na(who.Hat.sigma)){
    W.age.time <- get("W.age.time",envir=ecxc)
    who.C.age.time <- get("who.C.age.time", envir=ecxc)
    D.age.time.lst <- get("D.age.time.lst", envir=ecxc)
    Hat.wc.list <- get("Hat.wc.list", envir=ecxc) }
  
  if(!is.na(who.Ht.sigma)){
    W.time <- get("W.time", envir=ecxc)
    who.C.time <- get("who.C.time",envir=ecxc)
    D.time.lst <- get("D.time.lst", envir=ecxc)
    Ht.wc.list <- get("Ht.wc.list", envir=ecxc)}
### beta.hat.list covariates beta's classified for each cntry
### (with all ages for each list element)
  
### XX is a lst with as many elements as cntrys x age groups,
### where each element of lst is a matrix of dim= no.cov x no.cov; 
### no.cov = relevant covariates for csid identifier.  
### XX = X'%*%X, with X matrix of covariates as obtained in OLS
### Xy is also a list of csid's but now each elemnt is a one col
### matrix which result from X %*% y (y is death's)

### Here we compute basic quantities used by the priors
### Age prior: some relevant matrices derived from W.age
  if (!is.na(who.Ha.sigma)){
    omega.plus.age <- diag(diag(W.age)); ## diagonal of W.age
    omega.age <- omega.plus.age - W.age; ## diagonal = 0,
}
### Compute same basic quantities for time priors
### Time prior
  if (!is.na(who.Ht.sigma)){
    omega.plus.time <- diag(diag(W.time))
    omega.time <- omega.plus.time - W.time
}
### Compute same basic quantities for age time priors
### Age-time
  if (!is.na(who.Hat.sigma)){
    omega.plus.age.time <- diag(diag(W.age.time))
    omega.age.time <- omega.plus.age.time - W.age.time
}
### first beta.hat.csid estimates comes from OLS
### some storage for cntry weight or Ha.wc.lst,
### one number for each cntry
### estimated beta's are coeff or beta.hat.list
### clist with country unit
  clist <- vector(mode="list",length=length(cntry.vec))
  names(clist) <- as.character(cntry.vec);
  omega.Ca1a2.lst <- clist;
  omega.Caa.lst   <- clist;
  omega.Caa.time.lst <- clist
### note that omega.Ca1a2.time = 0 if a1 != a2
  omega.Ca1a2.age.time.lst <- clist;
  omega.Caa.age.time.lst   <- clist;
### Lambda^-1 withouth the theta for different priors and cntry unit
  Ltheta.age.lst <- clist
  Ltheta.time.lst <- clist
  Ltheta.age.time.lst <- clist 

###some list storage as fnction of csid = cntry+age
### unit is csid (cross sectional identifiers)
  names.ca <- names(whoinsampy)
  calst <- lapply(names.ca,function( x) x <- 0);
  names(calst) <- names.ca
### Ltheta= Lambda^{-1} for the priors (age, time, ag-time)
### but without the multiplier theta, which depends on sampling
### Ltheta is the constant part of Lambda.mn.csid (same for all samples)
  Ltheta.agprior.csid.lst  <- calst
  Ltheta.tmprior.csid.lst <- calst
  Ltheta.agtmprior.csid.lst <- calst
  beta.dim.lst <- calst
### our starting solution is the cxc.model (or ols model) of coeffs 
### beta.hat.lst  <- ols.beta;

  age.char <- formatC(age.vec, width=who.age.digits, format="d", flag="0")
  age.age  <- sapply(age.char,function(x) paste(x,x,sep=""));
#################### START LOOP OVER COUNTRIES #############################
 
  for (i in 1:length(cntry.vec)) {
    ttmp <- proc.time()
    cntry <-  cntry.vec[i]
    cntry.str <- as.character(cntry);
### print(length(cntry.vec)-i);
    csid <- paste(cntry.str,age.char,sep="");
### for the cntry's elements of who.C.age, and other lists of matrices 
### get cntry specific elements of lists
    ctr <- paste("^", cntry.str,sep="");
### get the likelihood matrices from the lst XX, Xy for given cntry
    whoiny <- whoinsampy[grep(ctr,names(whoinsampy))];
    whoinx <- whoinsampx[grep(ctr,names(whoinsampx))];
    cntryage  <-  kronecker(cntry,age.age,paste,sep="")
### t.list is a list like age.vec. Each element is the length of the time series
### of the covariates in whocov in the correspoding age group
### We are assuming that the time series in whocov have for each country,
### the same length independent of age groups.
### Otherwise, we need to take the average length of all age groups, 
### which affects the computation of wc.
    t.list  <- sapply(csid, function(n){nc <- nrow(whocov[[n]])});
    tt.lst <- t.list
    names(tt.lst) <- NULL
    chk <-  rep(t.list[1],length(t.list))
    names(chk) <- NULL
    try(if (!all.equal.numeric(tt.lst, chk))
      warning("Length of covariates differ in different age groups"))
### beta.dim is a list like age.vec. Each element is the number of covariates
### in the correspoding age group
    beta.dim  <- lapply(csid, function(n){nc <- ncol(whoinsampx[[n]])});
    beta.dim.lst[csid] <- beta.dim; 
### indx gives for each elemnt the first entry in a matrix (row or col)
### for corresponding first cov coeff (or beta) of any group of csid's
     ind.beta <- sapply(1:length(beta.dim), function(n) {
                   return (beta.dim[[n]] + n - 1)});    
### subset of list who.C.age (correlation matrices) for specific cntry and all ages
### *********** AGE PRIOR ***********************************
### age prior
    if (!is.na(who.Ha.sigma)){
      Ccntry.ages <- who.C.age[grep(ctr,names(who.C.age))];
      nm.C <- names(Ccntry.ages)
      sub.ind   <- match(cntryage, nm.C)
     ### correlation matrices for same age group: age, time and age-time priors
      Ccntry.ag.ag <- Ccntry.ages[sub.ind];
       
### diagonal elements of who.C.age for specific cntry and same age group
### autocorrelation, each multiply by the diagoanal elements of W.age
      omega.Caa <- lapply(1:length(Ccntry.ag.ag),function(x, Ccntry.ag.ag, W.age){
        Caa <- Ccntry.ag.ag[[x]];
        wa.plus <- W.age[x,x];
        return(wa.plus * Caa);}, Ccntry.ag.ag, W.age);
        names(omega.Caa) <- names(Ccntry.ag.ag)
   
      omega.D.age <- build.super.mat(omega.Caa,age.char,cntry.str)
         
### building D matrices compound of age groups correlation with W.age
### D.age <- make.age.time.prior.matrix(cntry,age.vec,W.age,who.C.age,env.base);
### Get them from envir=ecxc, where they are calculated for every cntry and
### stored in the lists: D.age.lst, D.time.lst, D,age.time.lst
     D.age <- D.age.lst[[cntry.str]]
     omega.Ca1a2 <- omega.D.age - D.age
     
### creates list for cntry dependent; age=a1, age=a2; and age, age
      omega.Ca1a2.lst[[cntry.str]] <- omega.Ca1a2;
      omega.Caa.lst[[cntry.str]] <- omega.Caa;
         
### the weigth for each cntry contribution to prior
### this is the covariance of the improper prior (the pseudoinverse of D), 
### which turns out the cntry weights, independent of age (w^{cntry}_{c})
### D.age.pinv <- sample.improper.normal(D.age,1,1)$covar;
### Ha.wc <- sum(diag(D.age.pinv%*%ZZ))/(length(age.vec)*mean(t.list));
      Ha.wc <- Ha.wc.list[[cntry.str]]
  
### Lambda.mn.prior is the super matrix of covariates correlation, which is 
### constructed with the sub-blocks matrices of every age group for cntry
### we indicate a super matrix for all age groups with the subscript ".cntry"
### If it is a list of age-cntry groups, one element for every age, then subscript is ".csid"
### For example, Lambda.mn.age (one large matrix with all age groups) for given cntry;
### and Lambda.mn.csid a list of elements matrices for each age group
      Ltheta   <-  Ha.wc * omega.D.age;
      Ltheta.age.lst[[cntry.str]] <- Ltheta
### Lambda.mn.age  <-  Ha.theta * Ha.wc * omega.D.age;
### we want to split Ltheta (or Lambda.mn.age) into a list, whose elements are the
### blocks sub-matrices corresponding to each age group
      lambda.indx <- make.Lambda.lst(Ltheta,age.char,cntry.str)
### Ltheta.csid is for each cntry+age-group; it does not depend on beta
### or the sampling oreder so it is a constant and equal for all samples
      Ltheta.csid <- lapply(lambda.indx, function(mm,Ltheta){
        fr <- mm[1,1]
        lr <- mm[1,2]
        fc <- mm[2,1]
        lc <- mm[2,2]
        as.matrix(Ltheta[fr:lr,fc:lc])}, Ltheta)
      la.nm <- names(Ltheta.csid)
### store in a list with csid identifiers
      Ltheta.agprior.csid.lst[la.nm] <- Ltheta.csid
        
    }
###
### **************TIME PRIOR ********************************
### time prior
    if (!is.na(who.Ht.sigma)){
     
      Ccntry.time <- who.C.time[grep(ctr,names(who.C.time))];
      nm.t <- names(Ccntry.time)
      sub.ind.t <- match(cntryage, nm.t)
      Ccntry.time.ag.ag <- Ccntry.time[sub.ind.t];
### autocorrelation for W.time
      omega.Caa.time <-  lapply(1:length(Ccntry.time.ag.ag),
                                function(x, Ccntry.time.ag.ag, W.time){
                                  Caa <- Ccntry.time.ag.ag[[x]];
                                  wa.plus <- W.time[x,x];
                                  return(wa.plus * Caa);}, Ccntry.time.ag.ag, W.time);
      names(omega.Caa.time) <- names(Ccntry.time.ag.ag)
      omega.D.time.ag.ag <- build.super.mat(omega.Caa.time,age.char,cntry.str)
      D.time <- D.time.lst[[cntry.str]]
      omega.Ca1a2.time <- omega.D.time.ag.ag - D.time  ### it should be = 0
### creates list for cntry dependent
      omega.Caa.time.lst[[cntry.str]] <- omega.Caa.time;
      Ht.wc <- Ht.wc.list[[cntry.str]]
      Ltheta.time <- Ht.wc * omega.D.time.ag.ag;
      Ltheta.time.lst[[cntry.str]] <- Ltheta.time
### Lambda.mn.time <- Ht.theta * Ht.wc * omega.D.time.ag.ag
### blocks for time prior
      lambda.indx <- make.Lambda.lst(Ltheta.time,age.char,cntry.str)
      Ltheta.time.csid <- lapply(lambda.indx, function(mm,Ltheta.time){
        fr <- mm[1,1]
        lr <- mm[1,2]
        fc <- mm[2,1]
        lc <- mm[2,2]
        as.matrix(Ltheta.time[fr:lr,fc:lc])}, Ltheta.time)
### Ltheta.tmprior.csid is for each cntry+age-group; it does not depend on beta
### or the sampling order so is constant and equal for all samples
### store in a list with csid identifiers
      lt.nm <- names(Ltheta.time.csid)
      Ltheta.tmprior.csid.lst[lt.nm] <-   Ltheta.time.csid;
    
    }
###
### ***********AGE-TIME PRIOR ********************  
### age-time prior
    if (!is.na(who.Hat.sigma)){
      Ccntry.age.time <- who.C.age.time[grep(ctr,names(who.C.age.time))]; 
      nm.at <- names(Ccntry.age.time)
### finding the elements of Ccntry.age for self-correlation with same age group
      sub.ind.at <- match(cntryage, nm.at)      
      Ccntry.age.time.ag.ag <- Ccntry.age.time[sub.ind.at];
### autocorrelation for W.age.time
      omega.Caa.age.time <-  lapply(1:length(Ccntry.age.time.ag.ag),
                                  function(x, Ccntry.age.time.ag.ag, W.age.time){
      Caa <- Ccntry.age.time.ag.ag[[x]];
      wa.plus <- W.age.time[x,x];
      return(wa.plus * Caa);}, Ccntry.age.time.ag.ag, W.age.time);
      names(omega.Caa.age.time) <- names(Ccntry.age.time.ag.ag)
      omega.D.age.time.ag.ag <- build.super.mat(omega.Caa.age.time,age.char,cntry.str)
      D.age.time <- D.age.time.lst[[cntry.str]]
      omega.Ca1a2.age.time <- omega.D.age.time.ag.ag - D.age.time
### creates list for cntry dependent
      omega.Ca1a2.age.time.lst[[cntry.str]] <- omega.Ca1a2.age.time;
      omega.Caa.age.time.lst[[cntry.str]] <- omega.Caa.age.time;
      Hat.wc <- Hat.wc.list[[cntry.str]]
      Ltheta.age.time <-  Hat.wc * omega.D.age.time.ag.ag
      Ltheta.age.time.lst[[cntry.str]] <- Ltheta.age.time; 
##  Lambda.mn.age.time <- Hat.theta * Hat.wc * omega.D.age.time.ag.ag
### blocks age.time prior
      lambda.indx <- make.Lambda.lst(Ltheta.age.time,age.char,cntry.str)
      Ltheta.age.time.csid <- lapply(lambda.indx, function(mm,Ltheta.age.time){
        fr <- mm[1,1]
        lr <- mm[1,2]
        fc <- mm[2,1]
        lc <- mm[2,2]
        as.matrix(Ltheta.age.time[fr:lr,fc:lc])}, Ltheta.age.time)
### Ltheta.agtime.csid is for each cntry+age-group; it does not depend on beta
### or the sampling order so is constant and equal for all samples
### store in a list with csid identifiers
      lat.nm <- names( Ltheta.age.time.csid)
      Ltheta.agtmprior.csid.lst[lat.nm] <-   Ltheta.age.time.csid;
     
    }
### this depends on the theta for each prior, age, time, age.time
### which values are taken for initial conditios as defined in
### who.Ha.sigma, who.Ht.sigma, who.Hat.sigma  
    
  }
  env.gibbs.age <- environment();
  assign("env.gibbs.age", env.gibbs.age, envir=ebase);
 
}

##################################################################
##
## FUNCTION NAME:  gibbs.age.beta
##
## PLACE:    
##
## IMPORTED:    covariates matrices whocov from make.mortality.data()
##              environments with globals, including countries and age groups 
##
## DESCRIPTION: It implements the gibbs sampling algorithm for the betas and
##              time, age ,age-time priors. The unit is csid but only ages and
##              time series within a cntry are smooth with the gibbs sampling.
##              Because we hav to integrate the smoothing over age-time with that of
##              cross-cntry, we need to define the unit csid for every cntry
##              and age combination.
##             
##         
## INPUT:  environmnets with globals and the results from the cxc models for betas
##         and some ols results for sigma (or std), preprocessing matrices and lists,    
##         smoothing parameters and the C.matrices from posteriors.
##         Many building blocks from previous programs  
##        
## OUTPUT:  the smooth betas one set of them for each csid and each covariates.
##          Thus if contributing covariates are say, i.e. gdp, tobacco, cnst, time
##          each cntry and age group will have one beta for each of the covariates;
##          for our example, each csid will have 4 betas for gdp, tobacco, cnst, time
##          lst <- list(btheta.bar.age.ca, btheta.bar.age.time.ca)
##          Only age and age-time priors contribute to betas, our resukt is lst
##          but we do not include any contributions specific to theta and sigma that
##          will be sample separately.  
##
## WRITTEN BY: Elena Villalon & Federico Girosi   
##             evillalon@latte.harvard.edu;
##             girosi@rand.org; 
##             CBRSS, Harvard University
## 
## Last modified: 03/01/2004
## 
## ************************************************************************
## *************************************************************************
gibbs.age.beta <- function(cntry, age.select, cxc.beta =NULL,sigma.ols =NULL,
                           Ha.theta=NULL,Ht.theta=NULL,Hat.theta=NULL) {   
#### **** beta calculations from OLS and CXC****************************
### OLS and CXC specific
### find cntry elements for cxc.model's coeff (or cxc.beta); 
### CXC model coefficients for covariates given cntry and all age groups
### ols.beta.csid <- ols.beta[grep(ctr,names(ols.beta))]
  ttmp <- proc.time()
  ebase <- get("env.base", envir=parent.frame())
  env.base <- ebase
  ewho <- get("env.who", envir=ebase)
  ecxc <- get("env.cxc", envir=ebase)
  egibbs  <- get("env.gibbs.age", envir=ebase)
  age.vec <- get("age.vec", envir=ewho)
  whoinsampx <- get("whoinsampx", envir=ewho)
  who.age.digits <- get("who.age.digits", envir=ewho)
  who.cntry.digits  <- get("who.cntry.digits", envir=ewho)
  who.digit.first   <- get("who.digit.first", envir=ewho)
  verbose <- get("verbose", envir=ebase)
  if(length(Ha.theta) <= 0)
    Ha.theta <- get("Ha.theta", envir= ecxc)
  if(length(Ht.theta) <= 0)
    Ht.theta <- get("Ht.theta", envir= ecxc)
  if(length(Hat.theta) <= 0)
    Hat.theta <-  get("Hat.theta", envir=ecxc)
  age.prior <- ( Ha.theta !=0 || Ht.theta !=0 || Hat.theta !=0)
  if (! age.prior )
    return(0); 
  digit.cntry.end   <- who.digit.first + who.cntry.digits; 
  digit.cntry.begin <- who.digit.first + 1
  digit.age.begin <- digit.cntry.end + 1
  digit.age.end   <- digit.cntry.end + who.age.digits
   
  cntry.str <- as.character(cntry) 
  if (Ha.theta != 0){
    omega.Ca1a2.lst <- get("omega.Ca1a2.lst", envir=egibbs)
    Ha.wc.list <- get("Ha.wc.list", envir=ecxc)
    Ha.wc  <- Ha.wc.list[[cntry.str]]
    omega.Ca1a2 <- omega.Ca1a2.lst[[cntry.str]]
  }
 
  if (Hat.theta != 0){
    omega.Ca1a2.age.time.lst <- get("omega.Ca1a2.age.time.lst", envir=egibbs)
    Hat.wc.list <- get("Hat.wc.list", envir=ecxc)
    Hat.wc <- Hat.wc.list[[cntry.str]]
    omega.Ca1a2.age.time <- omega.Ca1a2.age.time.lst[[cntry.str]]
}

  age.char <- formatC(age.vec, width=who.age.digits, format="d", flag="0")
  ctr <- paste("^", as.character(cntry), sep="")
  csid <- paste(cntry,age.char,sep="");
  beta.dim  <- sapply(csid, function(n){nc <- ncol(whoinsampx[[n]])});
  
  if(length(cxc.beta) <= 0)
    cxc.beta <- get("coeff", envir=ecxc)
  if(length(sigma.ols) <= 0)
    sigma.ols <- get("sigma.ols", envir=ecxc)
 
  cxc.beta.csid <- cxc.beta[grep(ctr,names(cxc.beta))]
### subscript ".csid" denotes a list of covariates matrices,
### each corresponding to an age group and specific cntry
### to first order approximation
### beta.hat.csid <- ols.beta.csid;
    beta.hat.csid <- cxc.beta.csid;
    bl <- names(cxc.beta.csid)
### some useful quantities    
    n.beta <- sum(sapply(cxc.beta.csid,nrow,simplify=T))
    if(n.beta != sum(unlist(beta.dim)))
       messout("Gibbs age:Dimensions for beta's do not match", verbose)
### building matrix of coefficients for given cntry,
### all age groups and their corresponding covariates
### find the names of a vector:cntry+ age group+ cov (say,245045gdp) 
    names.beta <- sapply(1:length(cxc.beta.csid), function(x, cxc.beta.csid) {
      res <- lapply(cxc.beta.csid[x], function(y){
        xx <- substring(names(cxc.beta.csid)[x],digit.age.begin, digit.age.end)
        nmx <- paste(xx,".",rownames(y),sep="")
        return(nmx) } )
      return(res)},cxc.beta.csid)
    names.beta <- unlist(names.beta, recursive=T)
### the .cntry denotes one matrix of 1 col and length =no covariates X age groups
### so all covariates for all age groups and specific cntry
    beta.hat.age <- matrix(unlist(beta.hat.csid, recursive = T, use.names=T));
    rownames(beta.hat.age) <- names.beta;
    colnames(beta.hat.age) <- cntry.str
### first estimation for beta.hat.cntry is from OLS
### calculate the parts contributing to the posteriors     
### age prior
     agend <- paste(age.select,"$", sep="")
     btheta.bar.age.ca <- 0
   if(Ha.theta !=0){
     beta.star.age <- omega.Ca1a2 %*% beta.hat.age;
     btheta <- Ha.wc * beta.star.age;
     colnames(btheta) <- cntry.str
     age.select  <- formatC(age.select, width=who.age.digits, format="d", flag="0")
     limb <- make.beta.lst(btheta,age.char, cntry.str)
     btheta.bar.age.csid <- lapply(limb, function(x, btheta){
       f <- x[1]
       s <- x[2]
       return(as.matrix(btheta[f:s]))}, btheta)
     nm <- names( btheta.bar.age.csid)
     indt <- grep(agend,nm)
     btheta.bar.age.ca <- btheta.bar.age.csid[[indt]]
     colnames( btheta.bar.age.ca) <- cntry.str
   }
### each elemnt of list is an age group for given cntry
### time prior gives beta.bar.time=0    
###************age-time prior**************************************
###
  btheta.bar.age.time.ca <- 0
  if(Hat.theta != 0){
      beta.star.age.time <- omega.Ca1a2.age.time %*% beta.hat.age
      btheta.age.time <-  Hat.wc * beta.star.age.time;
      limb <- make.beta.lst(btheta.age.time,age.char, cntry.str)
      btheta.bar.age.time.csid <- lapply(limb, function(x, btheta.age.time){
        f <- x[1]
        s <- x[2]
        return(as.matrix(btheta.age.time[f:s]))}, btheta.age.time)
      nm <- names( btheta.bar.age.time.csid)
      indt <- grep(agend,nm)
      btheta.bar.age.time.ca <- btheta.bar.age.time.csid[[indt]]
      colnames( btheta.bar.age.time.ca) <- cntry.str
    }
    
###
###   print(proc.time()-ttmp)
  
  lst <- list(btheta.bar.age.ca= btheta.bar.age.ca,
              btheta.bar.age.time.ca= btheta.bar.age.time.ca)
 return(lst)  }
 
###
##################################################################
##
## FUNCTION NAME:  check.mat
##
## PLACE:    
##
## IMPORTED:   matrices
##
## DESCRIPTION: to test some of the matrix calculations
##             
##
## WRITTEN BY: Elena Villalon    
##             evillalon@latte.harvard.edu;  
##             CBRSS, Harvard University
## 
## Last modified: 05/01/2004
## 
## ************************************************************************
check.mat <- function(Lambda, Lambda1, verbose=T){
  ebase <- try(get("env.base", envir=parent.frame()), silent=T)
  if(class(ebase)!="try-error")
  ###  verbose <- get("verbose", silent=T)
       lapply(1:length(Lambda), function(n, Lambda, Lambda1){
      x <- Lambda[[n]]
      y <- Lambda1[[n]]
      if( any(x != y)){
         messout("Gibss sampler; indexing from beta's and Lambda's do not agree", verbose)
         return (F)
      }else
        return(T)}, Lambda, Lambda1)}
###
##################################################################
##
## FUNCTION NAME:  build.super.mat
##
## PLACE:    
##
## IMPORTED:    globals and list of matrices omega
##
## DESCRIPTION: It takes the different elemnst of omega for each csid
##              and builds a matrix with them for each cntry
##             
##         
## INPUT:  environmnets with globals, omega list of matrices,
##          all age groups and one cntry
##
## OUTPUT:  the matrix with all age groups for each cntry           
##
## WRITTEN BY: Elena Villalon    
##             evillalon@latte.harvard.edu;  
##             CBRSS, Harvard University
## 
## Last modified: 03/01/2004
## 
## ************************************************************************
## *************************************************************************       
build.super.mat <- function(omega,agch,ctr){
  ebase <- get("env.base", envir=parent.frame())
  env.base <- ebase
  ewho <- get("env.who", envir=ebase)
  who.digit.first <- get("who.digit.first", envir=ewho)
  who.cntry.digits <- get("who.cntry.digits", envir=ewho)
  who.age.digits <- get("who.age.digits", envir=ewho)
  who.year.digits <- get("who.year.digits", envir=ewho)
  verbose <- get("verbose", envir=ebase)
### Structure of dataset cstsid: 
  digit.cntry.begin <- who.digit.first + 1 
  digit.cntry.end   <- who.digit.first + who.cntry.digits
  digit.age.begin   <- digit.cntry.end + 1
  digit.age.end     <- digit.cntry.end + who.age.digits
  digit.year.begin  <- digit.age.end   + 1 
  digit.year.end    <- digit.age.end   + who.year.digits
  if(length(agch) <= 0){
    age.vec <- get("age.vec", envir=ewho)
    agch <- formatC(age.vec, width=who.age.digits, format="d", flag="0")}
### omega is a list of matrices to be joined to form large mat
### agech a vector with age groups as chars
### id a vector with csid's identifiers= cntry+age.char (one cntry and all ages)
 each.col  <- sapply(omega, function(x) nc <- ncol(x) )
 each.row  <- sapply(omega, function(x) nc <- nrow(x) )
 s.col <- sum(each.col)
 s.row <- sum(each.row)
 itage <- sapply(agch,function(x){
   each <- paste(ctr, x, x,"$", sep="")
   dx <- grep(each, names(omega))
   return(dx) })
 nm.col <- sapply(itage, function(x) {
   nm2 <- substr(names(omega[x])[1],digit.age.begin,digit.age.end)
   nm2 <- paste(nm2,colnames(omega[[x]]),sep="") })
 nm.row <- sapply(itage, function(x){
   nm1 <- substring(names(omega[x])[1],digit.age.begin + who.age.digits )
   nm1 <- paste(nm1,rownames(omega[[x]]),sep="")})
### the block matrix 
### build the D matrices compound of age groups correlation with diagonal W.age
    omega.D <- matrix(0,nrow=s.row, ncol=s.col)
    rownames(omega.D) <- unlist(nm.row)
    colnames(omega.D)  <- unlist(nm.col)
    count <- 0; 
    for(r in 1:length(each.row)){
      count <- count + 1
      c = r;
       sc  <- ifelse(c > 1, sum(each.col[1:(c-1)]), 0)
        idc <- (sc + 1):(sc + each.col[c])
        sr  <- ifelse(r > 1, sum(each.row[1:(r-1)]), 0)
        idr <- (sr+1): (sr+ each.row[r])
        Drc <-  omega[[count]]
        omega.D[idr,idc] <- Drc}
    return(omega.D) }

##################################################################
##
## FUNCTION NAME:  make.beta.lst
##
## PLACE:    
##
## IMPORTED:   matrix beta (one column)  all ages
##
## DESCRIPTION: helper func 
##             
## OUTPUT : list for each csid and given cntry of betas 
##
## WRITTEN BY: Elena Villalon    
##             evillalon@latte.harvard.edu;  
##             CBRSS, Harvard University
## 
## Last modified: 05/01/2004
## 
## ************************************************************************
  make.beta.lst <- function(beta, agech,cntrych){
### find the indeces to split beta.bar.prior according to 
    finx <- lapply(1:length(agech), function(n, agech, beta){
    ch <- agech[n]
    ch <- paste("^", ch, sep="")
    mul <- grep(ch,rownames(beta))
    mulf <-  mul[1]
    mull <- mul[length(mul)]
    mult <- c(mulf,mull)
    return(mult)},agech,beta)
### apply the indeces to find elements of matrix
    bl <-  kronecker(cntrych,agech,paste,sep="")
    beta.indx <-as.list(finx)
    names(beta.indx) <- bl
   return(beta.indx)}
##################################################################
##
## FUNCTION NAME:  make.Lambda
##
## PLACE:    
##
## IMPORTED:   matrix Lambda, one cntry, all ages
##
## DESCRIPTION: helper func 
##             
## OUTPUT : list for each csid and given cntry of Lambda 
##
## WRITTEN BY: Elena Villalon & Federico Girosi    
##             evillalon@latte.harvard.edu;
##             girosi@rand.org; 
##             CBRSS, Harvard University
## 
## Last modified: 05/01/2004
## 
## ************************************************************************
  make.Lambda.lst <- function(Lambda, agech,cntrych, verbose=T){
    ebase <- try(get("env.base", envir=parent.frame()),silent=T)
    if(class(ebase)!="try-error")
      verbose <- get("verbose", envir=ebase)
      findx.c.r <- lapply(1:length(agech), function(x,Lambda){
      age <- agech[x]
      nc <- grep(age, colnames(Lambda))
      nr <- grep(age,rownames(Lambda))
      fc <- nc[1]
      lc <- nc[length(nc)]
      fr <- nr[1]
      lr <- nr[length(nr)]
      if(fc != fr || lc != lr)
        messout("Gibbs sampler error building matrices", verbose) 
      m <- matrix(c(fr,lr,fc,lc),nrow=2,byrow=T)
      rownames(m) <- c("rows", "cols")
      return(m)},Lambda)
      bl <-  kronecker(cntrych,agech,paste,sep="")
      names(findx.c.r) <- bl
      Lambda.mn.lst <- as.list(findx.c.r)
      names(Lambda.mn.lst) <- bl
      return(Lambda.mn.lst)}

###
##################################################################
##
## FUNCTION NAME:  build.Z
##
## PLACE:    
##
## IMPORTED:   matrix X with covariates (whoinsapmx) for each csid
##             and sigmas
##
## DESCRIPTION: helper func 
##             
## OUTPUT : X/sigma^2 
##
## WRITTEN BY: Elena Villalon    
##             evillalon@latte.harvard.edu;  
##             CBRSS, Harvard University
## 
## Last modified: 05/01/2004
## 
## ************************************************************************
 build.Z <- function(xx, ss, verbose=T){
   ebase <- try(get("env.base", envir=parent.frame()), silent=T)
   if(class(ebase)!="try-error")
     verbose <- get("verbose", envir=ebase)
     indx <- 1:length(xx)
     names(indx) <- names(xx) 
     zz <- lapply(indx, function(n){
       x <- xx[[n]]
       s <- ss[[n]]
       if (length(s) > 1)
         messout("Gibss sampler: we have more than one sigma for each csid", verbose)
       z <- x /s^2
       return(z)})
     return(zz)}
  
 ### beta= beta.hat.csid
 ### ctr = '^cntry'
 ### beta.dim = no of cols in whoinsampx
 ### names.beta as obtained from OLS
##################################################################
##
## FUNCTION NAME:  beta.by.cntry
##
## PLACE:    
##
## IMPORTED:   matrix beta (one column)  all ages and cntrys
##             listed for each csid
##
## DESCRIPTION: helper func
##
## INPUT: beta, cntry 
##             
## OUTPUT : one column matrix for each cntry with all age groups 
##
## WRITTEN BY: Elena Villalon    
##             evillalon@latte.harvard.edu;  
##             CBRSS, Harvard University
## 
## Last modified: 05/01/2004
## 
## ************************************************************************
beta.by.cntry <- function(ctr,beta,names.beta, beta.dim, verbose=T){
   ebase <- try(get("env.base", envir=parent.frame()), silent=T)
   if(class(ebase)!="try-error")
   verbose <- get("verbose", envir=ebase)
### find cntry elements for estimated model's coeff () and
### CXC model coefficients for covariates specific cntry and all age groups    
    beta.cntry <- beta[grep(ctr,names(beta))]
    n.beta <- sum(sapply(beta.cntry,nrow,simplify=T))
    if(n.beta != sum(unlist(beta.dim)))
       messout("Gibbs sampler:Dimensions for beta do not match", verbose)
### building matrix of coefficients for given cntry,
### all age groups and their corresponding covariates
### find the names of a vector:cntry+ age group+ cov (say,245045gdp) 
    beta.hat.cntry <- matrix(unlist(beta.cntry, recursive = T, use.names=T));
    rownames(beta.hat.cntry) <- names.beta;
    return(beta.hat.cntry)}
###
##################################################################
##
## FUNCTION NAME:  sigma.by.cntry
##
## PLACE:    
##
## IMPORTED:   matrix sigma (one column)  all ages and cntrys
##             listed for each csid
##
## DESCRIPTION: helper func
##
## INPUT: sigma, cntry 
##             
## OUTPUT : one column matrix for each cntry with all age groups 
##
## WRITTEN BY: Elena Villalon     
##             evillalon@latte.harvard.edu;  
##             CBRSS, Harvard University
## 
## Last modified: 05/01/2004
## 
sigma.by.cntry <- function(ctr,ss,names.ss){
### find cntry elements for estimated model's residuals 
### OLS model std for covariates given cntry and all age groups    
    ss.cntry <- ss[grep(ctr,names(ss))]
    n.ss <- sum(sapply(ss.cntry,nrow,simplify=T))
### building matrix of coefficients for given cntry,
### all age groups and their corresponding covariates
### find the names of a vector:cntry+ age group+ year (say,2450451999) 
    sigma.hat.cntry <- matrix(unlist(ss.cntry, recursive = T, use.names=T));
    rownames(sigma.hat.cntry) <- names.ss;
    return(sigma.hat.cntry)}
###
#### I'll need this later on to re-estimate my sigma###################  
 sigma.estimated.cntry <- function(ctr, beta.hat.csid){
     ebase <- get("env.base", envir=parent.frame())
     ewho <- get("env.who", envir= ebase)
     who.ols.sigma.param <- get("who.ols.sigma.param", envir=ewho)
     whoinsampy <- get("whoinsampy", envir=ewho)
     whoiny <- whoinsampy[grep(ctr, names(whoinsampy))]
     whoutsampy <- get("whoutsampy", envir=ewho)
     whouty <- whoutsampy[grep(ctr, names(whoutsampy))]
     whoinsampx <- get("whoinsampx", envir=ewho)
     whoutsampx <- get("whoutsampx", envir=ewho)
     whoinx <-    whoinsampx[grep(ctr, whoinsampx)]
     whoutx <-    whoutsampx[grep(ctr, whoutsampx)]
     whogender <- get("whogender", envir=ewho)
     whoyrest <- get("whoyrest", envir=ewho)
     age.vec <- get("age.vec", envir=ewho)
     ####
  cntry.names.lst <- get("cntry.names.lst", envir=ewho)
  cntry.names.lst <- cntry.names.lst[grep(ctr,names( cntry.names.lst))]
  gender.str <- ifelse(whogender==2, "m", "f")  
  clist <- vector(mode="list",length=length(whoiny));
  names(clist) <- names(whoiny)
  coeff <- clist; 
  std   <- clist;  
  sigma <- make.sigma(who.ols.sigma.param);
  sigma <- sigma[grep(ctr, names(sigma))]
    for (i in 1:length(whoiny)){
    y <- whoiny[[i]];
    x <- whoinx[[i]];
    beta <- beta.hat.csid[[i]]; 
### figure out the sigma[[i]]
    yw <- whoiny[[i]]/sigma[[i]];
    xw <- t(scale(t(whoinx[[i]]),center=FALSE,scale=sigma[[i]]));
    yxw <- na.omit(cbind(yw, xw));
    yw <- yxw[,1]
    xw <- yxw[,-1]
    n.obs <- length(yw);
    txw <- t(xw);
    invxw  <- svd.inv(txw%*% xw);
### beta <- invxw %*% txw %*% yw; beta is now given
    yin  <- x %*% beta; 
    std[[i]] <- sqrt(crossprod(na.omit(yin - y))/(n.obs-1));
    coeff[[i]] <- beta;
    rownames(coeff[[i]]) <- colnames(whoinx[[i]]);
  }
### now we make the forecasts of the dependent variable
  yhatin <- make.forecast(coeff,whoinx);
  yhatout <- make.forecast(coeff,whoutx);  
###  yhatin <- lapply(yhatin,FUN=y2logmortality,whotransform);
###  yhatout <- lapply(yhatout,FUN=y2logmortality,whotransform);
###  insampy <- lapply(whoiny,FUN=y2logmortality,whotransform);
###  outsampy <- lapply(whouty,FUN=y2logmortality,whotransform);      
  model <- model.string()
  lst <- list(coeff=coeff,yhatin=yhatin,yhatout=yhatout,std=std,insampy =whoinsampy,outsampy=whoutsampy)
 return(invisible(lst))}
