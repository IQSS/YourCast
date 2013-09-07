##################################################################
##
## FUNCTION NAME:  gibbs.cntry.cnst
##
## PLACE:    
##
## IMPORTED:    covariates matrices whocov from make.mortality.data()
##              environments with globals, including countries and age groups 
##
## DESCRIPTION: It calculates basic building matrices for
##              cross-cntry smoothing; Ltheta.ctprior.csid.lst
##              is a list which elemnts are csid with matrices Lambda but
##              with no contributions for theta and sigma priors.  Anything in
##              this function is constant independent on the sample number
##              It is calculated once for the sample algorithm and re-use for each
##              sample period.
##             
##         
## INPUT:  global environmnet
##          
##        
## OUTPUT: environmnet with calculations for  Lambda
##         for each csid,  Ltheta.ctprior.csid.lst
##           
##
## WRITTEN BY: Elena Villalon & Federico Girosi    
##             evillalon@latte.harvard.edu;
##             girosi@rand.org; 
##             CBRSS, Harvard University
## 
## Last modified: 05/01/2003
## 
## ************************************************************************
## *************************************************************************

gibbs.cntry.cnst <- function(ebase=env.base){
  ebase <- get("env.base", envir=parent.frame())
  env.base <- ebase
  ewho <- get("env.who", envir=ebase)
  ecxc <- get("env.cxc", envir=ebase)
  verbose <- get("verbose", envir=ebase)
### if the priors are not used then who.Hct.sigma = NA; same for all others
### the average standard deviation of the prior. If NA the prior is not ### used
### (it is like having an infinite standard deviation )
  who.Hct.sigma <- get("who.Hct.sigma", envir=ewho)
  who.Hct.t.deriv <- get("who.Hct.t.deriv", envir=ewho)
  who.Hct.time.weight <- get("who.Hct.time.weight", envir=ewho)
  if ( is.na(who.Hct.sigma))
    return(0);
  
############### PRIOR OVER CNTRY AND TIME GROUPS (for same age)################
######
  
  age.vec <- get("age.vec", envir=ewho)
  whoinsampy <- get("whoinsampy", envir=ewho)
  whoinsampx <- get("whoinsampx", envir=ewho)
  cntry.vec <- get("cntry.vec", envir=ewho)
  whocov <- get("whocov", envir=ewho)
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
###  ecxc <- cxc.model(param)
### some storage lst for sample beta or coeff
### dim(W.cntry) (or dim(omega.cntry)= length(cntry.vec) X length(cntry.vec)
### matrices for cntry's. Matrices of autocorrelation among cntrys,
### time and covariates, for given age
### only same age groups for same or different cntrys are correlated. 
### For who.C.cntry elements cntry1+cntry2+age correlation cntry1, cntry2;
### each elemnt is a matrix with multiplication among covariates
### 
  age.char <- formatC(age.vec, width=who.age.digits, format="d", flag="0")  
if(length(cntry.vec)>1){
  messout("Preparing for smoothing of time trend over countries", verbose);
  time.prior.param <- list(time.der=who.Hct.t.deriv,time.weight=who.Hct.time.weight)
  res <- build.C.cntry.time(cntry.vec,time.prior.param, env.base);
  who.C.cntry <- res$who.C.cntry
  W.cntry <- res$W.cntry
  }else
  return(0); 
###

    D.cntry.lst <- as.list(1:length(age.vec));
    names(D.cntry.lst) <- age.char;
    Hct.wc.lst  <-  as.list(1:length(age.vec));
    names(Hct.wc.lst)  <-  age.char;
  
  if (length(cntry.vec) > 1 && length(who.C.cntry)> 0){
  for(i in 1:length(age.vec)){
    namc <- paste(as.character(cntry.vec), age.char[i], sep="")
    output <- make.cntry.time.prior.matrix(age.vec[i],cntry.vec,W.cntry, who.C.cntry,env.base);
    Di <- output$D
### isolated cntry's not related to any in the data set
    isle.cntry <- output$isle.cntry
### cntry's isolated are not included in Di 
    D.cntry.lst[[i]] <- Di;
    Di.pinv <- sample.improper.normal(Di,1,1)$covar;
    inc <- NA; 
    if (length(isle.cntry) > 0)
      inc <- match(isle.cntry, cntry.vec)
    inc <- na.omit(inc)
    cntry.related <- cntry.vec
    if (length(inc) > 0)
      cntry.related <- cntry.vec[-inc]
    island.char <- paste("^", isle.cntry,sep="")  
    namc <- paste(as.character(cntry.related), age.char[i], sep="")
    cov.mat <- rmv.cntry(isle.cntry=isle.cntry, Xmat=whocov)
    tt.list  <- sapply(namc, function(n){nc <- nrow(cov.mat[[n]])});
    ZZcntry <- age.Xprime.X(ag=age.vec[i],Xmat=whocov,isle.cntry=isle.cntry, ebase=env.base);
    Hct.wc <- sum(diag(Di.pinv%*%ZZcntry))/(length(cntry.related)*mean(tt.list));
    Hct.wc.lst[[age.char[i]]] <- Hct.wc;
  
} }else
  return(0);
    
### Compute same basic quantities for cntry priors
### cntry smoothing
  omega.plus.cntry <- diag(diag(W.cntry))
  omega.cntry <- omega.plus.cntry - W.cntry
 
  age.age  <- sapply(age.char,function(x) paste(x,x,sep=""));
### include in cntry.char only correlated cntry's; isolated cntry's are excluded
  cntry.char <- formatC(cntry.related, width=who.cntry.digits,format="d", flag="0")
  ctr.ctr<- sapply(cntry.char, function(x) paste(x,x,sep=""))
### cntry matrices
  alist <- vector(mode="list",length=length(age.vec));
  names(alist) <- age.char
  omega.Cc1c2.cntry.lst <- alist;
  omega.Ccc.cntry.lst <- alist;
  Ltheta.cntry.lst <- alist
### our starting solution is the cxc.model (or ols model) of coeffs 
### beta.hat.lst  <- ols.beta; per csid
  names.ca <- names(cov.mat)
  calst <- lapply(names.ca,function( x) x <- 0);
  names(calst) <- names.ca
  Ltheta.ctprior.csid.lst <- calst;
  beta.dim.lst <- calst;
  
##   ttmp <- proc.time()
#################### START LOOP OVER AGES all countries#############################
### given a cntry and age add over all ages for cross-cntry smoothing 
  for (i in 1:length(age.vec)){
    age  <-  age.vec[i]
    age.str <- age.char[i];
### only correlated cntry's in cntry.char;
    csid <- paste(cntry.char,age.str,sep="");
### for the cntry's elements of who.C.age, and other lists of matrices 
### get cntry specific elements of lists
    atr <- paste(age.str,"$",sep="");
### get the likelihood matrices from the lst XX, Xy for given age
### and all cntry's in data set 
    whoiny  <- whoinsampy[grep(atr,names(whoinsampy))];
    whoinx  <-  whoinsampx[grep(atr,names(whoinsampx))];
     if (length(isle.cntry) > 0){
      whoiny  <- rmv.cntry(isle.cntry,whoiny)
      whoinx  <- rmv.cntry(isle.cntry, whoinx)}
    ctrage  <-  kronecker(ctr.ctr,age.str,paste,sep="")
### beta.dim is a list like age.vec. Each element is the number of covariates
### in the correspoding age group
    beta.dim  <- lapply(csid, function(n){nc <- ncol(whoinsampx[[n]])});
    names(beta.dim) <- csid
    beta.dim.lst[csid] <- beta.dim
### indx gives for each elemnt the first entry in a matrix (row or col)
### for corresponding first cov coeff (or beta) of any group of csid's
    ind.beta <- sapply(1:length(beta.dim), function(n) {
                   return (beta.dim[[n]] + n - 1)});
###
### *********** CNTRY PRIOR ***********************************
### subset of list who.C.cntry (correlation matrices) for specific age
### cntry prior
    Cage.cntrys <- who.C.cntry[grep(atr,names(who.C.cntry))]; 
    nm.Ca <- names(Cage.cntrys)
### country correlated with itself for given age=age.vec[i]
    sub.ind   <- match(ctrage, nm.Ca)
### correlation matrices for same cntry : cntry prior
    Cage.c.c <- Cage.cntrys[sub.ind];
### diagonal elements of who.C.cntry for specific age and same cntry for all cntry's
### autocorrelation, each multiply by the diagoanal elements of W.cntry
    omega.Ccc <- lapply(1:length(Cage.c.c),function(x, Cage.c.c, W.cntry){
      Ccc <- Cage.c.c[[x]];
      wc.plus <- W.cntry[x,x];
      return(wc.plus * Ccc);}, Cage.c.c, W.cntry);
    names(omega.Ccc) <- names(Cage.c.c)
### only for correlated cntry's in cntry.char; isolated cntry are excluded
    omega.D.cntry    <- build.super.cntry.mat(omega.Ccc,age.str,cntry.char)
### building D matrices compound of cntrys groups correlation with W.cntry
### stored in the lists: D.cntry.lst
    D.cntry <- D.cntry.lst[[age.str]]
    omega.Cc1c2 <- omega.D.cntry - D.cntry
### creates list for cntry dependent; cntry=c1, cntry=c2; and cntry, cntry
    omega.Cc1c2.cntry.lst[[age.str]] <- omega.Cc1c2;
    omega.Ccc.cntry.lst[[age.str]]   <- omega.Ccc;
### the weigth for each cntry contribution to prior
### this is the covariance of the improper prior (the pseudoinverse of D)
### D.cntry.pinv <- sample.improper.normal(D.cntry,1,1)$covar;
### Hct.wc <- sum(diag(D.cntry.pinv%*%ZZcntry))/(length(cntry.vec)*mean(t.list));
    Hct.wc <- Hct.wc.lst[[age.str]]
### (Ltheta) is the super matrix of covariates correlation, which is 
### constructed with the sub-blocks matrices of every  cntry for each age group
### we indicate a super matrix for all cntry with the subscript ".cntry"
### When it is a list of age-cntry groups, one element for every cntry, subscript is ".csid"
### For example, Ltheta.cntry.lst (one large matrix with all cntry's) for given age group;
### and Ltheta.csid a list of elements matrices for each cntry and age
    Ltheta <- Hct.wc * omega.D.cntry
    Ltheta.cntry.lst[[age.str]] <- Ltheta
### Lambda.mn.cntry  <- Hct.theta * Hct.wc * omega.D.cntry;
### we want to split Ltheta.cntry prior into a list, whose elements are the
### blocks sub-matrices corresponding to each age group and cntry (Lambda^{-1}_{ca})
    lambda.indx <- make.Lambda.age.lst(Ltheta,age.str,cntry.char)
### Ltheta.csid is for each cntry+age-group; it does not depend on beta
### or the sampling order so it is a constant and equal for all samples
    Ltheta.csid <- lapply(lambda.indx, function(mm,Ltheta){
      fr <- mm[1,1]
      lr <- mm[1,2]
      fc <- mm[2,1]
      lc <- mm[2,2]
    as.matrix(Ltheta[fr:lr,fc:lc])}, Ltheta)
    lamb.nm <- names(Ltheta.csid)
### store in a list with csid identifiers
    Ltheta.ctprior.csid.lst[lamb.nm]  <- Ltheta.csid       
  }
###
##  print(proc.time()-ttmp)   
  env.gibbs.cntry <- environment();
  assign("env.gibbs.cntry",  env.gibbs.cntry, envir=ebase); 
}


##################################################################
##
## FUNCTION NAME:  gibbs.cntry.beta
##
## PLACE:    
##
## IMPORTED:    covariates matrices whocov from make.mortality.data()
##              environments with globals, including countries and age groups 
##
## DESCRIPTION: It implements the gibbs sampling algorithm for the betas and
##              cntry priors. The unit is csid but only for the same age group
##              related cntry's are  smoothed-out. 
##              Because we have to integrate the smoothing over cntry with that of
##              age-time, we need to define the unit csid for every cntry
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
##          lst <- list(btheta.bar.cntry.ca) 
##          Only those betas contribute to betas, our resukt is lst
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
gibbs.cntry.beta <- function(age,cntry.select, cxc.beta =NULL, Hct.theta=NULL) {     
###
### **** beta calculations from OLS and CXC****************************
### OLS and CXC specific
### find cntry elements for cxc.model's coeff (or cxc.beta); 
### CXC model coefficients for covariates given age and all cntrys
### ols.beta.csid <- ols.beta[grep(ctr,names(ols.beta))]
##  ttmp <- proc.time()
  ebase <- get("env.base", envir=parent.frame())
  env.base <- ebase
  ewho <- get("env.who", envir=ebase)
  ecxc <- get("env.cxc", envir=ebase)
  egibbs  <- get("env.gibbs.cntry", envir=ebase)
  age.vec <- get("age.vec", envir=ewho)
  cntry.vec  <- get("cntry.vec", envir=ewho)
  whoinsampx <- get("whoinsampx", envir=ewho)
  who.age.digits <- get("who.age.digits", envir=ewho)
  who.cntry.digits  <- get("who.cntry.digits", envir=ewho)
  who.digit.first   <- get("who.digit.first", envir=ewho)
  who.Hct.sigma <- get("who.Hct.sigma", envir= ewho)
  verbose <- get("verbose", envir=ebase)
   if ( is.na(who.Hct.sigma))
     return (0); 
  digit.cntry.end   <- who.digit.first + who.cntry.digits; 
  digit.cntry.begin <- who.digit.first + 1
  digit.age.begin   <- digit.cntry.end + 1
  digit.age.end     <- digit.cntry.end + who.age.digits
  if (length(Hct.theta) <= 0)
    Hct.theta <- get("Hct.theta", envir=ecxc)
  if (Hct.theta == 0)
    return (0);
  
  Hct.wc.lst <- get("Hct.wc.lst", envir=egibbs)
  omega.Cc1c2.cntry.lst <- get("omega.Cc1c2.cntry.lst", envir=egibbs)
  beta.dim.lst <- get("beta.dim.lst", envir=egibbs)
  isle.cntry <- get("isle.cntry", envir=egibbs)
### the vector numeric form of cntry.char with only the subset of
### correlated cntry's of cntry.vec 
  cntry.related <- get("cntry.related", envir=egibbs)
### cntry.char is character (string) with correlated cntry's
  cntry.char <- get("cntry.char", envir=egibbs)
  
  age.str <- formatC(age, width=who.age.digits, format="d", flag="0")
  atr <- paste(age.str, "$", sep="")
  csid <- paste(cntry.char,age.str, sep="")
  Hct.wc  <- Hct.wc.lst[[age.str]]
  omega.Cc1c2   <- omega.Cc1c2.cntry.lst[[age.str]]
####
   
  if(length(cxc.beta) <= 0)
    cxc.beta <- get("coeff", envir=ecxc)
  
    cxc.beta.csid <- cxc.beta[grep(atr,names(cxc.beta))]
### subscript ".csid" denotes a list of covariates matrices,
### each corresponding to a cntry for specific age froup
### to first order approximation
### beta.hat.csid subset of cxc.beta.csid with correlated cntry's
    beta.hat.csid <- rmv.cntry(isle.cntry, cxc.beta.csid);
    bl <- names(beta.hat.csid)
### some useful quantities
    n.beta <- sum(sapply(beta.hat.csid,nrow,simplify=T))
    if(n.beta != sum(unlist(beta.dim.lst[csid]))){
       messout(paste("Gibbs cntry: Dimensions for beta's do not match", csid, sep=""), verbose)
     }
### building matrix of coefficients for given cntry
### among the list of cntry.related(cntry.char) ,or correlated cntrys
### all age groups and their corresponding covariates
### find the names of a vector:cntry+ age group+ cov (say,245045gdp) 
    names.beta <- sapply(1:length(beta.hat.csid), function(x, beta.hat.csid) {
      res <- lapply(beta.hat.csid[x], function(y){
        xx <- substring(names(beta.hat.csid)[x],digit.cntry.begin, digit.cntry.end)
        nmx <- paste(xx,".",rownames(y),sep="")
        return(nmx) } )
      return(res)},beta.hat.csid)
    names.beta <- unlist(names.beta, recursive=T)
### the .age denotes one matrix of 1 col and length =no covariates X cntry's
### so all covariates for specific age group and all cntry
    beta.hat.cntry <- as.matrix(unlist(beta.hat.csid, recursive = T, use.names=T));
    rownames(beta.hat.cntry) <- names.beta;
### first estimation for beta.hat.cntry is from OLS
### cntry prior 
    beta.star.cntry <- omega.Cc1c2 %*% beta.hat.cntry;
    colnames(beta.star.cntry) <- age.str
    btheta <-  Hct.wc * beta.star.cntry;
    colnames(btheta) <- age.str
    cntry.select  <- as.character(cntry.select)
   ctrstart <- paste("^", cntry.select,sep="")
   limb <- make.beta.cn.lst(btheta,age.str, cntry.char)
   btheta.bar.cntry.csid <- lapply(limb, function(x, btheta){
     f <- x[1]
     s <- x[2]
     return(as.matrix(btheta[f:s]))}, btheta)
   nm <- names( btheta.bar.cntry.csid)
   indt <- grep(ctrstart,nm)
   btheta.bar.cntry.ca <- btheta.bar.cntry.csid[[indt]]
   colnames( btheta.bar.cntry.ca) <- age.str
  
###   print(proc.time()-ttmp)
  lst = list(btheta.bar.cntry.ca=btheta.bar.cntry.ca)
  return(lst); 
     
}
  
##################################################################
##
## FUNCTION NAME:  build.super.cntry.mat
##
## PLACE:    
##
## IMPORTED:    globals and list of matrices omega
##
## DESCRIPTION: It takes the different elemnst of omega for each csid
##              and builds a matrix with them for each age group
##             
##         
## INPUT:  environmnets with globals, omega list of matrices,
##          all cntry's and specific age group
##
## OUTPUT:  the matrix with all cntry  for each age group           
##
## WRITTEN BY: Elena Villalon    
##             evillalon@latte.harvard.edu;  
##             CBRSS, Harvard University
## 
## Last modified: 03/01/2004
## 
## ************************************************************************
## *************************************************************************       

build.super.cntry.mat  <- function(omega,agch,ctr){
  ebase <- get("env.base", envir=parent.frame())
  env.base <- ebase
  ewho <- get("env.who", envir=ebase)
  who.digit.first <- get("who.digit.first", envir=ewho)
  who.cntry.digits <- get("who.cntry.digits", envir=ewho)
  who.age.digits <- get("who.age.digits", envir=ewho)
  who.year.digits <- get("who.year.digits", envir=ewho)
### Structure of dataset cstsid: 
  digit.cntry.begin <- who.digit.first + 1 
  digit.cntry.end   <- who.digit.first + who.cntry.digits
  digit.age.begin   <- digit.cntry.end + 1
  digit.age.end     <- digit.cntry.end + who.age.digits
  digit.year.begin  <- digit.age.end   + 1 
  digit.year.end    <- digit.age.end   + who.year.digits
  if(length(ctr) <= 0){
    cntry.vec <- get("cntry.vec", envir=ewho)
    ctr <- formatC(cntry.vec, width=who.cntry.digits, format="d", flag="0")}
### omega is a list of matrices to be joined to form large mat
### agech a vector with age groups as chars
### id a vector with csid's identifiers= cntry+age.char (one cntry and all ages)
 each.col  <- sapply(omega, function(x) nc <- ncol(x) )
 each.row  <- sapply(omega, function(x) nc <- nrow(x) )
  s.col <- sum(each.col)
 s.row <- sum(each.row)
 itage <- sapply(ctr,function(x){
   each <- paste("^",x, x,agch,"$", sep="")
   dx <- grep(each, names(omega))
   return(dx) })
 nm.col <- sapply(itage, function(x) {
   nm2 <- substr(names(omega[x])[1],digit.cntry.begin,digit.cntry.end)
   nm2 <- paste(nm2,colnames(omega[[x]]),sep="") })
 nm.row <- sapply(itage, function(x){
   nm1 <- substring(names(omega[x])[1],digit.cntry.begin + who.cntry.digits )
   nm1 <- paste(nm1,rownames(omega[[x]]),sep="")})
### the block matrix 
### build the D matrices compound of age groups correlation with diagonal W.age
    omega.D <- matrix(0,nrow=s.row, ncol=s.col)
    rownames(omega.D) <- unlist(nm.row)
    colnames(omega.D)  <- unlist(nm.col)
    totcount <- length(omega);
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
   if(count != totcount)
     stop("Wrong build.super.cntry.mat"); 
  return(omega.D) }

##################################################################
##
## FUNCTION NAME:  make.beta.cn.lst
##
## PLACE:    
##
## IMPORTED:   matrix beta (one column)  all ages
##
## DESCRIPTION: helper func 
##             
## OUTPUT : list for each csid and given agem for all cntrys of betas 
##
## WRITTEN BY: Elena Villalon    
##             evillalon@latte.harvard.edu;  
##             CBRSS, Harvard University
## 
## Last modified: 05/01/2004
## 
## ************************************************************************
  make.beta.cn.lst <- function(beta, agech,cntrych){
### find the indeces to split beta.bar.prior according to 
    finx <- lapply(1:length(cntrych), function(n, cntrych, beta){
    ch <- cntrych[n]
    ch <- paste("^", ch, sep="")
   
    mul <- try(grep(ch,rownames(beta)), silent=T)
    mulf <-  mul[1]
    mull <- mul[length(mul)]
    mult <- c(mulf,mull)
    return(mult)},cntrych,beta)
### apply the indeces to find elements of matrix
    bl <-  kronecker(cntrych,agech,paste,sep="")
    beta.indx <-as.list(finx)
    names(beta.indx) <- bl
   return(beta.indx)}

###
##################################################################
##
## FUNCTION NAME:   make.Lambda.age.lst 
##
## PLACE:    
##
## IMPORTED:   matrix Lambda, all cntry, one age group
##
## DESCRIPTION: helper func 
##             
## OUTPUT : list for each csid and cntry's of Lambda 
##
## WRITTEN BY: Elena Villalon    
##             evillalon@latte.harvard.edu;  
##             CBRSS, Harvard University
## 
## Last modified: 05/01/2004
## 
## ************************************************************************
  make.Lambda.age.lst <- function(Lambda, agech,cntrych, verbose=T){
    ebase <- try(get("env.base", envir=parent.frame()),silent=T)
    if(class(ebase)!="try-error")
    verbose <- get("verbose", envir=ebase)
      findx.c.r <- lapply(1:length(cntrych), function(x,Lambda){
      ctr <- cntrych[x]
      ctr <- paste("^", ctr, sep="")
      nc <- grep(ctr, colnames(Lambda))
      nr <- grep(ctr,rownames(Lambda))
      fc <- nc[1]
      lc <- nc[length(nc)]
      fr <- nr[1]
      lr <- nr[length(nr)]
      if(length(nc) == length(nr) && (fc != fr || lc != lr))
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


 ### beta= beta.hat.csid
 ### ctr = '^cntry'
 ### beta.dim = no of cols in whoinsampx
 ### names.beta as obtained from OLS




