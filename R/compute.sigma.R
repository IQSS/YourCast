##################################################################
##
## FUNCTION NAME:  compute.sigma 
##
## PLACE:    
##
## IMPORTED:  sigma.vc(), sumTc() function helpers and global environmnet
##
## DESCRIPTION: Computes for each sample period the standard deviation as
##              function of the current beta (depends on sample and the prior)
##              The output list elements are identified with csid (cntry+ age group)
##              Sigma is factorized as cntry contribution, cnstant independent of sample time
##              and age contribution, which depends on beta and on sample period
##             
## OUTPUT : a list with values of std sigma for each csid. After drawing a random
##          number from gamma distribution 
##
## WRITTEN BY: Elena Villalon    
##             evillalon@latte.harvard.edu;  
##             CBRSS, Harvard University
## 
## Last modified: 05/01/2004
## 
## ***************************************
compute.sigma <-  function(beta = NULL,ebase = env.base){
  ebase    <- get("env.base", envir=parent.frame())
  env.base <- ebase
  ewho <- get("env.who", envir=ebase)
  ecxc <- get("env.cxc", envir=ebase)
  whoinsampy <- get("whoinsampy", envir=ewho)
  cntry.vec <- get("cntry.vec", envir=ewho)
  age.vec <- get("age.vec", envir=ewho)
  who.age.digits <- get("who.age.digits", envir=ewho)
  who.cntry.digits <- get("who.cntry.digits", envir=ewho)
  verbose <- get("verbose", envir=ebase)
  if(length(beta) <=0 )
     beta  <- get("coeff", envir=ecxc)
  
  cntry.char <- as.character(cntry.vec)
  age.char  <- formatC(age.vec, width=who.age.digits, format="d", flag="0")
  sigma.ag.lst <- as.list(1:length(age.vec))
  names(sigma.ag.lst) <- age.char
  sigma.csid.lst <- as.list(1:length(whoinsampy))
  names(sigma.csid.lst) <- names(whoinsampy)
### a list of vector elements: vc one for each cntry,
### with the inverse sqrt of the average no. dth for each cntry.
### vc.csid one fore each csid with inverse square of
### average no. dth for each cntry age group
   sigma.vc <- sigma.vc();
   vc <- sigma.vc$vc  ##vector length of cntry.vec
   vc.csid <- sigma.vc$vc.csid ##lst length of whoinsampy
### a number with the sum for all cntry's of the time length
### insample period
   sumTc <- sumTc();  ##number
   for(a in 1:length(age.vec)){
      age.str <-  formatC(age.vec[a], width=who.age.digits, format="d", flag="0")
      calst <- kronecker(as.character(cntry.vec), age.str,paste,sep="")
      agend <- paste(age.str,"$", sep="")
## a number for each age group, after adding for all cntry's
      sse <- sigma.SSE(age=age.str, beta.ca=beta,vc=vc, vc.csid=vc.csid);
### lst with two numbers d, e

      param <- sigma.param()
      d <- param$d
      e <- param$e
      dTc <-  0.5 *(d + sumTc)## number for each age
      eSS <- 0.5 * (e + sse) ## number
      sa  <- rgamma(1,shape = dTc, scale = 1/eSS)
      sigma.a <- sqrt(1/sa) ## number for each age
      sigma.ag.lst[[age.str]] <- sigma.a
      invc <- grep(agend,names(vc.csid))
      vc.ca <- unlist(vc.csid[invc]) ##for the age group and all cntry's
      sigma.ca <- kronecker(vc.ca,sigma.a,FUN="*")
      sigma.csid.lst[calst] <- sigma.ca
       
## s <- rgamma(n,shape=d/2,scale=2/e)
## sigma <- 1/sqrt(s)
## m.sample <- mean(sigma)
## std.sample <- sd(sigma)
## print(paste("desired mean:",m))
## print(paste("sample mean:",m.sample))
## print(paste("desired standard deviation:",std))
## print(paste("sample standard deviation:",std.sample))
##
        }## for(a in 1:length(age.vec)) sigma calculations
  return(sigma.csid.lst)
}
  
#####################################################################################
##################################################################
##
## FUNCTION NAME:  sumTc() 
##
## PLACE:    
##
## IMPORTED:  globals and whoinsampy
##
## DESCRIPTION: Computes for each cntry the average number of death in the 
##              insample period, after eliminating NA's from whoinsampy
##
## OUTPUT : a vector with average observations for each cntry.
##
## WRITTEN BY: Elena Villalon    
##             evillalon@latte.harvard.edu;  
##             CBRSS, Harvard University
## 
## Last modified: 05/01/2004
## 
## ***************************************
sumTc <-  function(ebase=env.base){
  ebase    <- get("env.base", envir=parent.frame())
  env.base <- ebase
  ewho <- get("env.who", envir=ebase)
  whoinsampy <- get("whoinsampy", envir=ewho)
  cntry.vec <- get("cntry.vec", envir=ewho)
  
  csid <- as.list(1:length(whoinsampy))
  names(csid) <- names(whoinsampy)
  t.list  <- sapply(whoinsampy, function(x){
          x <- na.omit(x)
          nc <- nrow(x)});
  ct <- paste("^", as.character(cntry.vec), sep="")
  ct <- as.list(ct)
  names(ct) <- as.character(cntry.vec)
  tindx <- sapply(ct,function(x) {grep(x, names(t.list))})
  Tc <- sapply(tindx, function(vec) {
            tt     <- t.list[vec]
            n.tlst <- length(tt)
            Tc <- sum(tt)/n.tlst
            return(Tc)})
   sumTc <- sum(unlist(Tc));}
##################################################################
##
## FUNCTION NAME:  sigma.vc() 
##
## PLACE:    
##
## IMPORTED:  globals and whoinsampy
##
## DESCRIPTION: Computes for each csid  the average number of death in the 
##              insample period, after eliminating NA's from whoinsampy
##              Also calculates for each cntry the average no. dth over all csid
##              or age groups and given cntry
##
## OUTPUT : a list of three elements, the no. dth for each csid, the inverse
##          of sqrt of no. dth for each csid, and inverse of the sqrt of
##          the mean no.dth for each cntry
##
## WRITTEN BY: Elena Villalon    
##             evillalon@latte.harvard.edu;  
##             CBRSS, Harvard University
## 
## Last modified: 05/01/2004
## 
## ***************************************

    sigma.vc <- function(ebase=env.base){
             ebase    <- get("env.base", envir=parent.frame())
             env.base <- ebase
             ewho <- get("env.who", envir=ebase)
             whoinsampy <- get("whoinsampy", envir=ewho)
             cntry.vec  <- get("cntry.vec", envir=ewho)
             age.vec <- get("age.vec", envir=ewho)
             who.age.digits <- get("who.age.digits", envir=ewho)
             age.char <-  formatC(age.vec, width=who.age.digits, format="d", flag="0")
             ct <- paste("^", as.character(cntry.vec), sep="")
             ct <- as.list(ct)
             names(ct) <- as.character(cntry.vec)
             csid <- names(whoinsampy)
             names(csid) <- names(whoinsampy)
### number of death for each csid in the insample period
             no.dth <- lapply(csid, function(ca,whoinsampy){
                 x <- whoinsampy[[ca]]
                 x <- na.omit(x)
                 no.dth <- length(x)
                 return(no.dth)},whoinsampy)
### sqrt of the number of dth for each csid
             vc.csid <- 1./sqrt(unlist(no.dth))
             vc.csid <-  as.list(vc.csid)
             names(vc.csid) <- names(whoinsampy)
### each cntry get the indeces for no.dth
             vec.ind <- lapply(ct,function(x) {grep(x, names(no.dth))})
### calculate for every cntry the average number of dths insample period 
             meano.dth <- sapply(vec.ind,function(v,no.dth) {
               vv <- no.dth[v]
               vv <- unlist(vv)
               mn <- mean(vv)
               return(mn)},no.dth)
### sqrt of the average no dth for each cntry
             vc <- 1./sqrt(meano.dth)
             lst <- list(no.dth.csid = no.dth, vc.csid = vc.csid, vc=vc)
             return(lst)}

## ######################################################################################
## ######################################################################################
##
## FUNCTION NAME: sigma.param
##
## DESCRIPTION: this function is used to set the parameters e and d which appear in the prior for sigma
##
##
## FORMAT:    param <- sigma.param(m,std)
##
##
## INPUT:   m: positive scalar, the mean of the prior for the standard deviations
##
##        std: positive scalar, the standard deviation of the prior for the standard deviations
##
## OUTPUT: param: list with two elenents: param$d, and param$e. These are the parameters which
##                the prior for s = 1/sigma^2. More precisely, if s is sampled according to the
##                distribution Gamma(d/2,e/2) then sigma = 1/sqrt(s) will have mean and standard
##                deviation given by m and std respectively.
## IMPORTANT: the definition of the Gamma distribution used in R is different from the one used in the
##            Girosi/King book. When we say that the prior for s is Gamma(d/2,e/2) we mean
##            that
##                     P(s) = K s^{d/2-1} e^{-e/2*s}
##            In order to sample 100 times from such a distribution in R we would set
##                   s <- rgamma(100, shape=d/2,scale=2/e)
##
##             You can test the function sigma.param running the following code
## 
## n <- 100000
## m <- 0.3;
## std <- 0.1;
## 
## param <- sigma.param(m,std)
## d <- param$d
## e <- param$e
## 
## s <- rgamma(n,shape=d/2,scale=2/e)
## sigma <- 1/sqrt(s)
## m.sample <- mean(sigma)
## std.sample <- sd(sigma)
## print(paste("desired mean:",m))
## print(paste("sample mean:",m.sample))
## print(paste("desired standard deviation:",std))
## print(paste("sample standard deviation:",std.sample))
##
##
## WRITTEN BY: Federico Girosi 
##             girosi@rand.org
##             RAND, Santa Monica, CA
## 
## Last modified: 5/10/2004
## 
## ######################################################################################
## ######################################################################################


sigma.param <- function(m=NULL,std=NULL, ebase = env.base){
  ebase    <- get("env.base", envir=parent.frame())
  env.base <- ebase
  ewho <- get("env.who", envir=ebase)
  if (length(m) <= 0)
    m <- get("LI.sigma.mean", envir=ewho)
  if (length(std) <= 0)
    std <- get("LI.sigma.sd", envir=ewho)
##  cat(m,"\n")
##  cat(std, "\n")
  v <- std^2

fun <- function(x,a,v){r <- sqrt(x/2-1)*exp(lgamma(x/2-0.5) - lgamma(x/2)) - a/sqrt(a^2+v)};
aux <- m^2/(v+m^2);
interval <- c(2.0000001,1000)
x <- seq(from=interval[1],to=interval[2],length=100)
### plot(x,fun(x,m,v));
d  <- uniroot(fun,interval,a=m,v=v)$root
e <- (d - 2)*(v+m^2)
return(list(e=e,d=d))
}

##################################################################
##
## FUNCTION NAME:  sigma.SSE 
##
## INPUT: age group, beta for age cntry combo,
##        vc (or inverse sqrt of average no.dth each cntry) from sigma.vc
##        vc.csid (inverse sqrt no.dth for each csid) from sigma.vc 
##
## IMPORTED:  globals and sigma.vc
##
## DESCRIPTION: Computes the part of sigma prior with mortality and covariates  
##              contribution and beta that is unique for each age group and contains
##              all cntry's and add over all years.  
##
## OUTPUT : number for each age group to be input into
##          gamma function to obtain s_{a} of sigma prior.
##
## WRITTEN BY: Elena Villalon    
##             evillalon@latte.harvard.edu;  
##             CBRSS, Harvard University
## 
## Last modified: 05/01/2004
## 
## ***************************************
sigma.SSE <- function(age,beta.ca,vc,vc.csid) {
  ebase    <- get("env.base", envir=parent.frame())
  env.base <- ebase
  ewho <- get("env.who", envir=ebase)
  ecxc <- get("env.cxc", envir=ebase)
  age.vec <- get("age.vec", envir=ewho)
  who.age.digits <- get("who.age.digits", envir=ewho)
  cntry.vec <- get("cntry.vec", envir=ewho)
  who.age.digits <- get("who.age.digits", envir=ewho)
  whoinsampx <- get("whoinsampx", envir=ewho)
  whoinsampy <- get("whoinsampy", envir=ewho)
  who.cntry.digits <- get("who.cntry.digits", envir=ewho)
  who.digit.first  <- get("who.digit.first", envir=ewho)
  verbose <- get("verbose", envir=ebase)
  age.char <- formatC(age.vec, width=who.age.digits, format="d", flag="0")
  age.str <-  formatC(age, width=who.age.digits, format="d", flag="0")
  agend  <- paste(age.str,"$", sep="")
  whoiny <- whoinsampy[grep(agend,names(whoinsampy))]
  whoinx <- whoinsampx[grep(agend,names(whoinsampx))]                    
  beta   <- beta.ca[grep(agend,names(beta.ca))]
## a vector with all cntrys and one age group
  vc.ca  <- vc.csid[grep(agend,names(vc.csid))]
  calist <- as.list(1:length(whoinx))
  names(calist) <- names(whoinx)
  Zbeta <- sapply(calist, function(i,whoinx,beta) {
    xw  <- whoinx[[i]]
    yw  <- whoiny[[i]]
    yxw <- cbind(yw,xw)
    yxw <- na.omit(yxw)
    xw <- yxw[,-1]
    yw <- yxw[,1]
    b  <- beta[[i]]
    nmx <- substr(names(whoinx)[i],1,who.cntry.digits )
    nmb <- substr(names(beta)[i],1,who.cntry.digits)
    nmv <- substr(names(vc)[i],1, who.cntry.digits)
    nmva <- substr(names(vc.ca)[i],1,who.cntry.digits)
    if (nmx != nmb || nmv != nmb || nmx != nmv || nmx != nmva)
      messout("Inside sigma.SSE: Names do not agree", verbose)
    z  <- yw - xw %*% b
    sc <- (t(z) %*% z)/vc.ca[[i]]^2
### sc <- (t(z) %*% z)/vc[i]^2
    return(sc)},whoinx,beta)
  
    sse.a <- sum(unlist(Zbeta));
}
                      
##################################################################
##
## FUNCTION NAME:  theta.age.priors 
##
## INPUT: beta.csid for age cntry combo,
##        environmnets env.base, env.who, env.cxc, 
##        which contains the matrices and parameters for any of the priors
##        age, time, age-time
##
## IMPORTED:  globals and environmnets
##
## DESCRIPTION: Computes the part of theta priors with beta's 
##              contribution that is unique for each of the priors
##              all cntry's and add over all age groups.  
##
## OUTPUT : number for each prior to be input in gamma function.
##
## WRITTEN BY: Elena Villalon    
##             evillalon@latte.harvard.edu;  
##             CBRSS, Harvard University
## 
## Last modified: 06/08/2004
## 
## ***************************************                      
    
theta.age.priors <- function(beta.csid=NULL, ebase=env.base){
  ebase    <- get("env.base", envir=parent.frame())
  env.base <- ebase
  ewho <- get("env.who", envir=ebase)
  ecxc <- get("env.cxc", envir=ebase)
  age.vec <- get("age.vec", envir=ewho)
  cntry.vec <- get("cntry.vec", envir=ewho)
  if(length(beta.csid) <= 0)
    beta.csid <- get("coeff", envir=ecxc)
  who.Ha.sigma  <- get("who.Ha.sigma", envir=ecxc)
  who.Hat.sigma <- get("who.Hat.sigma", envir=ecxc)
  who.Ht.sigma  <- get("who.Ht.sigma", envir=ecxc) 
  Ha.wc.list  <-  get("Ha.wc.list", envir=ecxc)
  Hat.wc.list <-  get("Hat.wc.list", envir=ecxc)
  Ht.wc.list  <-  get("Ht.wc.list", envir=ecxc)
  D.age.lst <-  get("D.age.lst", envir=ecxc)
  D.age.time.lst <-  get("D.age.time.lst", envir=ecxc)
  D.time.lst <-  get("D.time.lst", envir=ecxc)
  beta.cntry <- lapply(as.character(cntry.vec), function(x,beta.csid){
         ctr <- paste("^",x,sep="")
         ind <- grep(ctr, names(beta.csid))
         beta.c <- beta.csid[ind]
         beta.c <-  matrix(unlist(beta.c, recursive = T, use.names=T));
         return(beta.c)}, beta.csid)
  names(beta.cntry) <- as.character(cntry.vec)
  
  th.Ha.sum  <- 0
  th.Ht.sum  <- 0
  th.Hat.sum <- 0
  rkag.lst <- 0
  rktm.lst <- 0
  rkagtm.lst <- 0
  
  if(!is.na(who.Ha.sigma)){
    th.Ha <- sapply(as.character(cntry.vec), function(x, Ha.wc.list,D.age.lst,beta.cntry){
            Ha.wc <- Ha.wc.list[[x]]
            D.age <- D.age.lst[[x]]
            beta  <- beta.cntry[[x]]
            sc <- t(beta) %*% D.age %*% beta
            sc <- 0.5 * Ha.wc * sc
            return(sc)}, Ha.wc.list,D.age.lst,beta.cntry)
    th.Ha.sum <- sum(unlist(th.Ha))
    if(th.Ha.sum <= 0) th.Ha.sum <- 0
    rkag.lst <- lapply(D.age.lst,function(x) {find.zero.eigen(x)$w.rank})
         }
  
    if(!is.na(who.Ht.sigma)){
    th.Ht <- sapply(as.character(cntry.vec), function(x, Ht.wc.list,D.time.lst,beta.cntry){
               Ht.wc  <- Ht.wc.list[[x]]
               D.time <- D.time.lst[[x]]
               beta <- beta.cntry[[x]]
               sc <- t(beta) %*% D.time %*% beta
               sc <- 0.5 * Ht.wc * sc
               return(sc)}, Ht.wc.list,D.time.lst,beta.cntry)
    
    th.Ht.sum <- sum(unlist(th.Ht))
    if(th.Ht.sum <= 0)  th.Ht.sum <- 0
       
    rktm.lst <- lapply(D.time.lst,function(x) {find.zero.eigen(x)$w.rank})
  }
  
    if(!is.na(who.Hat.sigma)){
    th.Hat <- sapply(as.character(cntry.vec), function(x, Hat.wc.list,D.age.time.lst,beta.cntry){
            Hat.wc <- Hat.wc.list[[x]]
            D.age.time <- D.age.time.lst[[x]]
            beta <- beta.cntry[[x]]
            sc <- t(beta) %*% D.age.time %*% beta
            sc <- 0.5 * Hat.wc * sc
            return(sc)}, Hat.wc.list,D.age.time.lst,beta.cntry)
    th.Hat.sum <- sum(unlist(th.Hat))
    if(th.Hat.sum <= 0) th.Hat.sum <- 0
    rkagtm.lst <- lapply(D.age.time.lst,function(x) {find.zero.eigen(x)$w.rank})
  }
  lst <- list(th.Ha.sum=th.Ha.sum, th.Ht.sum= th.Ht.sum, th.Hat.sum= th.Hat.sum,
              rkag.lst = rkag.lst, rktm.lst=rktm.lst,rkagtm.lst=rkagtm.lst )
  return(lst)
  
      
}
  ##################################################################
##
## FUNCTION NAME:  theta.cntry.prior 
##
## INPUT: beta.csid for age cntry combo,
##        environmnets env.base, env.who, env.cxc, 
##        which contains the matrices and parameters for cntry prior
##
## IMPORTED:  globals and environmnets
##
## DESCRIPTION: Computes the part of theta priors with beta's 
##              contribution that is unique for cntry prior
##               
##
## OUTPUT : number for each prior to be input in gamma function.
##
## WRITTEN BY: Elena Villalon    
##             evillalon@latte.harvard.edu;  
##             CBRSS, Harvard University
## 
## Last modified: 06/08/2004
## 
## ***************************************                      
    
theta.cntry.prior <- function(beta.csid=NULL, isle.cntry=NULL, ebase=env.base){
  ebase    <- get("env.base", envir=parent.frame())
  env.base <- ebase
  ewho <- get("env.who", envir=ebase)
  ecxc <- get("env.cxc", envir=ebase)
  age.vec <- get("age.vec", envir=ewho)
  cntry.vec <- get("cntry.vec", envir=ewho)
  who.age.digits <- get("who.age.digits", envir=ewho)
  who.Hct.sigma   <- get("who.Hct.sigma", envir=ecxc)
  Hct.wc.lst  <- get("Hct.wc.lst", envir=ecxc)
  D.cntry.lst <- get("D.cntry.lst", envir=ecxc)
  if(length(beta.csid) <= 0)
    beta.csid <- get("coeff", envir=ecxc)
  if(length(isle.cntry) <= 0)
    isle.cntry <- get("isle.cntry", envir=ecxc)
  age.char <- formatC(age.vec, width=who.age.digits, format="d", flag="0")
  if(length(isle.cntry) > 0){
     n.age <- length(age.vec)
    ln <- length(isle.cntry) * n.age
   
    ind.beta <- sapply(1:ln,function(x,age.char,isle.cntry,beta.csid){
     n.age <- length(age.char)
      xc <- x %/% n.age
      agc <- x %% n.age
      if(agc == 0 && xc !=0 ) agc <- agc + n.age
      if(x%%n.age != 0){
        xc <- xc + 1
      }else{
        agc <- agc+1}
      agc <- age.char[agc]
      agc[is.na(agc)] <- age.char[length(age.char)]
       is <- isle.cntry[xc]      
      isc <- paste(is,agc,sep="")
      ind <- grep(isc,names(beta.csid))
      return(ind)},age.char,isle.cntry,beta.csid ) 
    ind.beta <- unlist(ind.beta)
    beta.csid <- beta.csid[-ind.beta]}
  beta.age <- lapply(age.char, function(x,beta.csid){
         atr <- paste(x,"$",sep="")
         ind <- grep(atr, names(beta.csid))
         beta.a <- beta.csid[ind]
         beta.a <-  matrix(unlist(beta.a, recursive = T));
         return(beta.a)}, beta.csid)
  names(beta.age) <- age.char
##  beta <- matrix(unlist(beta.cntry,recursive=T, use.name=T))
  
  th.Hct.sum  <- 0
 
 
  
  if(!is.na(who.Hct.sigma)){
    th.Hct <- sapply(age.char, function(x, Hct.wc.lst,D.cntry.lst,beta.age){
 
            Hct.wc <- Hct.wc.lst[[x]]
            D.cntry <- D.cntry.lst[[x]]
            beta  <- beta.age[[x]]
            sc <- t(beta) %*% D.cntry %*% beta
            sc <- 0.5 * Hct.wc * sc
            return(sc)}, Hct.wc.lst,D.cntry.lst,beta.age)
    th.Hct.sum <- sum(unlist(th.Hct))
    if(th.Hct.sum <= 0) th.Hct.sum <- 0
  }
    rkct.lst <- lapply(D.cntry.lst,function(x) {find.zero.eigen(x)$w.rank})  
    lst <- list(th.Hct.sum=th.Hct.sum,rkct.lst=rkct.lst)
  return(lst)
  
      
  }
  

      
                      
