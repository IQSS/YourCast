##################################################################
##
## FUNCTION NAME:  gibbs.sampler
##
## PLACE:    
##
## IMPORTED:    gibbs.cntry and gibbs.age, for the cntry and age-time priors
##              environments with globals, including countries and age groups
##              compute.sigma 
##
## DESCRIPTION: The driver for the gibbs sampler; it loops through all cntry's in the
##              cntry set and age groups.  For every csid, it updates the coefficients
##              for as many samples periods as required until they converge.  The
##              starting input for the coefficients may be the solution from cxc or
##              the solution from ols. It invokes age-time priors and cntry prior
##              depending on values of who.Ha[Hat,Ht].sigma and who.Hct.sigma; 
##              the paramter sigma is also sample for age groups.  
##             
##         
## INPUT:    cxc solutions for beta, sigma ols, globals and outputs from priors. 
##          
##          
##        
## OUTPUT:  updated coefficients with cross cntry and/or age-time smoothing.
##          mortality prediction for outsample periods according to beta values. 
##           
##
## WRITTEN BY: Elena Villalon & Federico Girosi   
##             evillalon@latte.harvard.edu;  
##             CBRSS, Harvard University
##             girosi@rand.org
## 
## Last modified: 05/13/2004
## 
## ************************************************************************
## *************************************************************************
## param <- list(Ha.sigma=who.Ha.sigma,Hat.sigma=who.Hat.sigma,Ht.sigma=who.Ht.sigma)

gibbs.sampler <- function(nsample=NULL,cxc.beta=NULL,bool=TRUE ){
  
   ebase    <- get("env.base", envir=parent.frame())
   env.base <- ebase
   ewho <- get("env.who", envir=ebase)
   verbose <- get("verbose", envir=ebase)
   messout("Running Bayes model", verbose)
   lst <- cxc();
   
   
   ecxc <- get("env.cxc", envir=ebase)
   
   who.age.digits <- get("who.age.digits", envir=ewho)
   if(length(nsample) <= 0)
    nsample <- get("nsample", envir=ewho)
   XX <- get("XX", envir=ewho)
   ### the following is necessary because we need to subtract the mean of the prior from whoinsampy
   Xy <- Xy.only(env.base)

   
   sigma.hat.csid <- get("sigma.ols", envir=ecxc)
   who.Ha.sigma   <- get("who.Ha.sigma", envir=ewho)
   who.Ha.sigma.sd   <- get("who.Ha.sigma.sd", envir=ewho)
   who.Hat.sigma  <- get("who.Hat.sigma", envir=ewho)
   who.Hat.sigma.sd   <- get("who.Hat.sigma.sd", envir=ewho)
   who.Ht.sigma   <- get("who.Ht.sigma", envir=ewho)
   who.Ht.sigma.sd   <- get("who.Ht.sigma.sd", envir=ewho)
   who.Hct.sigma  <- get("who.Hct.sigma", envir=ewho)
   who.Hct.sigma.sd   <- get("who.Hct.sigma.sd", envir=ewho)
   who.zero.mean <- get("who.zero.mean", envir=ewho)
   whoinsampx <- get("whoinsampx", envir=ewho)
   whoutsampx <- get("whoutsampx", envir=ewho)
   whoinsampy <- get("whoinsampy",envir=ewho);
   whoutsampy <- get("whoutsampy", envir=ewho) 
   age.vec   <- get("age.vec", envir=ewho)
   cntry.vec <- get("cntry.vec", envir=ewho)
   whoyrest <- get("whoyrest", envir=ewho)
   cntry.names.lst <- get("cntry.names.lst", envir=ewho)
   tol <- get("solve.tol", envir=ewho)
   verbose <- get("verbose", envir=ebase)
   if (!is.logical(who.zero.mean)) {
     who.mean.age.profile <- who.zero.mean;
     func <- function(x,param){x-param};
     whoinsampy <- modify.age.profiles(whoinsampy,func,who.mean.age.profile)
   }else{
     if(!who.zero.mean){
       who.mean.age.profile <- make.average.age.profile(whoinsampy,ebase=env.base)
       func <- function(x,param){x-param};
       whoinsampy <- modify.age.profiles(whoinsampy,func,who.mean.age.profile)
     }
   }
  
   param <- list(Ha.sigma=who.Ha.sigma,Hat.sigma=who.Hat.sigma,Ht.sigma=who.Ht.sigma)
 
### theta parameter to multiply prior; if theta --> 0, prior --> 0
  Ha.sigma <- param$Ha.sigma;
  if (is.na(Ha.sigma)) {Ha.theta <-  0;
  } else { Ha.theta <- 1/Ha.sigma^2}
### theta parameter to multiply prior; if theta --> 0, prior --> 0
  Ht.sigma <- param$Ht.sigma;
  if (is.na(Ht.sigma)) {Ht.theta <-  0;
  } else { Ht.theta <- 1/Ht.sigma^2}
### theta parameter to multiply prior; if theta --> 0, prior --> 0
  Hat.sigma <- param$Hat.sigma;
  if (is.na(Hat.sigma)) {Hat.theta <-  0;
  } else { Hat.theta <- 1/Hat.sigma^2}
###
   Hct.sigma <- who.Hct.sigma;
  if (is.na(Hct.sigma)) {Hct.theta <-  0;
  } else{ Hct.theta <- 1/Hct.sigma^2}

### get the constants for gibbs model with cntry and age priors
### age, time and age-time priors
  age.prior <- (!is.na(who.Ha.sigma) || !is.na(who.Ht.sigma) || !is.na(who.Hat.sigma) )

   if (age.prior == T){
    evage <- try(get("env.gibbs.age", envir= ebase), silent=T)
  
    if (class(evage) == "try-error" || bool == F )
      evage <- try(gibbs.age.cnst(ebase))
  
  }
     
### 
### cntry prior
 
  if (!is.na(who.Hct.sigma)){
    evcntry <-  try(get("env.gibbs.cntry", envir= ebase), silent=T)
    if (class(evcntry) == "try-error" ||  bool== F)
      evcntry <- gibbs.cntry.cnst(ebase)
    isle.cntry    <- get("isle.cntry", envir= evcntry)
    messout("Isolated cntry's are ", verbose, obj=isle.cntry)
    cntry.related <- get("cntry.related", envir=evcntry)
    messout("Related cntry's are ", verbose, obj=cntry.related)
   
     }
###
  cntry.char <- as.character(cntry.vec)
  
  age.char  <- formatC(age.vec, width=who.age.digits, format="d", flag="0")
  if ( age.prior == T){
    
    Ltheta.agprior.csid.lst <- get("Ltheta.agprior.csid.lst", envir=evage)
      
    Ltheta.tmprior.csid.lst <- get("Ltheta.tmprior.csid.lst", envir=evage)
    Ltheta.agtmprior.csid.lst <- get("Ltheta.agtmprior.csid.lst", envir=evage)
  }
    
  if (!is.na(who.Hct.sigma))
    Ltheta.ctprior.csid.lst <- get("Ltheta.ctprior.csid.lst", envir=evcntry)
  
  if(length(cxc.beta) <=0 )
     beta.hat.csid  <- get("coeff", envir=ecxc)
  else
     beta.hat.csid <- cxc.beta

   beta.hat.csid.old <- beta.hat.csid
### to obtain the average over the number of samples
  beta.hat.csid.total  <- lapply(beta.hat.csid, function(x) x/nsample) 
  sigma.hat.csid.total <- lapply(sigma.hat.csid, function(x) x/nsample);
### for the standard errors
  S.csid <- vector(mode="list", length=length(whoinsampy))
  names(S.csid) <-  names(whoinsampy)
   
  ttmp <- proc.time(); 
  for (n in 1:nsample){
   messout(n,verbose)
    for (c in 1:length(cntry.vec)) {
      
      if (!is.na(who.Hct.sigma) && age.prior == F && length(isle.cntry) > 0 ){
             
          is.isle <- grep(as.character(cntry.vec[c]), isle.cntry)
           if(length(is.isle) > 0)
            next; }
     for(a in 1:length(age.vec)){
     if (!is.na(who.Ha.sigma) || !is.na(who.Hat.sigma))  
      betas <- gibbs.age.beta(cntry =cntry.char[c], age.select= age.char[a],
                              cxc.beta=beta.hat.csid, sigma.ols= sigma.hat.csid,
                              Ha.theta=Ha.theta,Ht.theta=Ht.theta,Hat.theta=Hat.theta)    
      beta.bar.age.ca <-  0
     if (Ha.theta != 0)
      beta.bar.age.ca   <- (betas$btheta.bar.age.ca) * Ha.theta
     
      beta.bar.age.time.ca <- 0 
     if (Hat.theta != 0)
        beta.bar.age.time.ca <- (betas$btheta.bar.age.time.ca) * Hat.theta
      
     beta.bar.cntry.ca <- 0
      
     if (Hct.theta != 0){
       is.isle <- NA
       if (length(isle.cntry) > 0 ){
         is.isle <- match(cntry.vec[c], isle.cntry)
     
        if (length(na.omit(is.isle)) > 0 && n==1)
          messout(paste("Isolated cntry ",cntry.vec[c], " and age group ", age.vec[a], sep=""),verbose)
       }
       if( is.na(is.isle)){
         betac <- gibbs.cntry.beta(age=age.char[a], cntry.select=cntry.char[c],
                                   cxc.beta=beta.hat.csid, Hct.theta=Hct.theta)
         beta.bar.cntry.ca <- (betac$btheta.bar.cntry.ca) * Hct.theta}
       }
     
      row.vec <- 0
      n <-  0
      csid  <- paste(cntry.char[c],age.char[a],sep="")     
    if (Ha.theta != 0){
      row.vec1 <- rownames(beta.bar.age.ca)
      n <- length(row.vec1)
      row.vec <- row.vec1
    }else if (Hat.theta != 0){
      row.vec2 <-  rownames(beta.bar.age.time.ca)
      n2 <- length(row.vec2)
      if( n2 > n ) row.vec <- row.vec2
    }else if (Hct.theta != 0 && is.na(is.isle)){
      row.vec3 <-  rownames(beta.bar.cntry.ca)
      n3 <- length(row.vec3)
      if(n3 > n) row.vec <- row.vec3
    }else{
      coeff <- get("coeff", envir=ecxc)
      inx <- grep(csid,names(coeff))
      row.vecn <- rownames(coeff[[inx]])
      nn <- length(row.vecn)
      row.vec <- row.vecn}                 
      nm.sg <- names(sigma.hat.csid)
      nm.X  <- names(XX)
      nm.Xy <- names(Xy)
      inds  <- grep(csid,nm.sg)
      indX  <- grep(csid,nm.X)
      indY  <- grep(csid,nm.Xy)
      sigma.hat.ca <- unique.default(sigma.hat.csid[[inds]])
     
      XX.ca <- XX[[indX]]/(sigma.hat.ca)^2
      Xy.ca <- Xy[[indY]]/(sigma.hat.ca)^2

     beta <- 0 
     beta.bar.ca <- Xy.ca 
     if (Ha.theta != 0 || Hat.theta != 0 ){
     beta   <- as.matrix( beta.bar.age.ca + beta.bar.age.time.ca)
      colnames(beta) <- csid}
     if(Hct.theta != 0 && is.na(is.isle))
       beta <- as.matrix(beta + beta.bar.cntry.ca)
     
     beta.bar.ca <- beta.bar.ca  + beta 
        
      Lambda.mn.ca  <- XX.ca

     if (Ha.theta != 0){
      Lambda.mn.age.ca <- Ha.theta * Ltheta.agprior.csid.lst[[csid]]
      Lambda.mn.ca <- Lambda.mn.ca + Lambda.mn.age.ca

    }
     if (Hat.theta != 0){
      Lambda.mn.age.time.ca <-  Hat.theta * Ltheta.agtmprior.csid.lst[[csid]]
      Lambda.mn.ca <-  Lambda.mn.ca +  Lambda.mn.age.time.ca

    }
     if (Ht.theta != 0){
      Lambda.mn.time.ca <-  Ht.theta * Ltheta.tmprior.csid.lst[[csid]]
      Lambda.mn.ca <-  Lambda.mn.ca + Lambda.mn.time.ca

    }
     if (Hct.theta != 0 && is.na(is.isle) ){
      Lambda.mn.cntry.ca <- Hct.theta * Ltheta.ctprior.csid.lst[[csid]]
      if (length(XX.ca) != length(Ltheta.ctprior.csid.lst[[csid]])){
        messout(length(XX.ca), verbose)
        messout(length( Ltheta.ctprior.csid.lst[[csid]]), verbose)
          stop("Error in matrices of gibbs sampler");}
      Lambda.mn.ca <- Lambda.mn.ca + Lambda.mn.cntry.ca
 
    }
    
     if(!class(try(svd.inv(Lambda.mn.ca),silent=F))=="try-error"){
        
        Lambda.ca <- svd.inv(Lambda.mn.ca)
     }else if ( !class(try(solve(Lambda.mn.ca, tol=tol),silent=F))=="try-error"){
         
         Lambda.ca <- solve(Lambda.mn.ca,tol=tol)
    }else if(all(is.na(Lambda.mn.ca)))
         Lambda.ca <- Lambda.mn.ca
     else
             Lambda.ca <- NULL
   
             
      if(!any(is.na(Lambda.ca))){
        s <- svd(Lambda.ca)
        Lambda.ca.root <- s$u %*% diag(sqrt(s$d)) %*% s$v;
      }else{
        messout("No inverting matrices at gibss.model...check data", verbose)
        Lambda.ca.root <- Lambda.ca
      }
      beta.bar.ca <- Lambda.ca %*% beta.bar.ca
      len <- length(beta.bar.ca)
### for checking purposes only
### wn <- matrix(0, nrow=len);
### end of testing
### white noise
     wn <- rnorm(len,mean=0,sd=1);
     beta.hat.csid[[csid]] <-  beta.bar.ca + Lambda.ca.root %*% matrix(wn);
     rownames(beta.hat.csid[[csid]]) <- rownames(beta.hat.csid.old[[csid]]) 
     beta.hat.csid.total[[csid]] <- beta.hat.csid[[csid]]/nsample  +
                                    beta.hat.csid.total[[csid]]

} ## length(age.vec)
   } ## length(cntry.vec)
### end of beta prior for each sample
### here comes the sigma prior, which is the same for all beta priors

  sigma.hat.csid <- compute.sigma(beta=beta.hat.csid)
    
  sigma.hat.csid.total <- lapply(1:length(sigma.hat.csid), function(x) {
    sigma <- sigma.hat.csid[[x]]/nsample + sigma.hat.csid.total[[x]];
    return(sigma); })
      
  names(sigma.hat.csid.total) <- names(sigma.hat.csid)
## standard errors for death predictions insample and outsample   
  S.csid <-  standard.errors(n,beta.hat.csid,sigma.hat.csid,whoinsampx,whoutsampx,S.csid)
### end of sigma prior
### here comes theta priors, one for each of the beta priors

 if(age.prior){
   that <- theta.age.priors(beta.csid=beta.hat.csid)
   Ha.sum  <- that$th.Ha.sum
   rkag.lst <- that$rkag.lst
   rkag <- sum(unlist(rkag.lst))
   Ht.sum  <- that$th.Ht.sum
   rktm.lst <- that$rktm.lst
   rktm <- sum(unlist(rktm.lst))
   Hat.sum <- that$th.Hat.sum
   rkagtm.lst <- that$rkagtm.lst
   rkagtm <- sum(unlist(rkagtm.lst))
 }
  
  if(!is.na(who.Hct.sigma)){
    thct <- theta.cntry.prior(beta.csid=beta.hat.csid, isle.cntry=isle.cntry)
    Hct.sum <- thct$th.Hct.sum
    rkct.lst <- thct$rkct.lst
    rkct <- sum(unlist(rkct.lst))
  }

  if(!is.na(who.Ha.sigma)){
     if (is.na(who.Ha.sigma.sd) || who.Ha.sigma.sd == 0){
      Ha.theta <- 1/who.Ha.sigma^2
    } else {
      sp <- sigma.param(m=who.Ha.sigma,std=who.Ha.sigma.sd)
      g <- sp$e
      f <- sp$d
      
      Ha.theta <- rgamma(1,shape=0.5*(f +rkag), scale=1/(0.5*g + Ha.sum))
    }
  }
 
  
    if(!is.na(who.Ht.sigma)){
      if (is.na(who.Ht.sigma.sd) || who.Ht.sigma.sd == 0){
        Ht.theta <- 1/who.Ht.sigma^2
      } else {   
        sp <- sigma.param(m=who.Ht.sigma,std=who.Ht.sigma.sd)
        g <- sp$e
        f <- sp$d
   
        Ht.theta <- rgamma(1,shape=0.5*(f +rktm), scale=1/(0.5*g + Ht.sum))
      
      }
      
    }
    if(!is.na(who.Hat.sigma)){
      if (is.na(who.Hat.sigma.sd) || who.Hat.sigma.sd == 0){
      Hat.theta <- 1/who.Hat.sigma^2
    } else {         
      sp <- sigma.param(m=who.Hat.sigma,std=who.Hat.sigma.sd)
      g <- sp$e
      f <- sp$d
      Hat.theta <- rgamma(1,shape=0.5*(f +rkagtm), scale=1/(0.5*g + Hat.sum))
    }
  }   
    if(!is.na(who.Hct.sigma)){
       if (is.na(who.Hct.sigma.sd) || who.Hct.sigma.sd == 0){
         Hct.theta <- 1/who.Hct.sigma^2
       } else {             
         sp <- sigma.param(m=who.Hct.sigma,std=who.Hct.sigma.sd)
         g <- sp$e
         f <- sp$d
         Hct.theta <- rgamma(1,shape=0.5*(f +rkct), scale=1/(0.5*g + Hct.sum))
       }
    
     }
    
          
   }## n=1:nsample
  ### now we make the forecasts of the dependent variable
###  print(proc.time()-ttmp)  
  coeff <- beta.hat.csid.total
  sigma <- sigma.hat.csid.total
   yhatin <- make.forecast(coeff,whoinsampx);
   yhatout <- make.forecast(coeff,whoutsampx);
     
 if (!is.logical(who.zero.mean)) {
    who.mean.age.profile <- who.zero.mean;
    func <- function(x,param){x+param};
    yhatin <- modify.age.profiles(yhatin,func,who.mean.age.profile)
    whoinsampy <- modify.age.profiles(whoinsampy,func,who.mean.age.profile)    
    yhatout <- modify.age.profiles(yhatout,func,who.mean.age.profile)  
   
  }else{
    if(!who.zero.mean){
      who.mean.age.profile <- make.average.age.profile(whoinsampy)
      func <- function(x,param){x+param};
      yhatin <- modify.age.profiles(yhatin,func,who.mean.age.profile)
      whoinsampy <- modify.age.profiles(whoinsampy,func,who.mean.age.profile)          
      yhatout <- modify.age.profiles(yhatout,func,who.mean.age.profile)  
         }
  }
      
  
### standard errors
  
   V.csid <- variance.csid(S=S.csid,yin=yhatin,yout=yhatout)
 

   model <- model.string()
   
 lst <- list(yrest=whoyrest,model=model,age.vec=age.vec, cntry.lst=cntry.names.lst,
              coeff=coeff,yhatin=yhatin,yhatout=yhatout,std=V.csid,
              insampy =whoinsampy,outsampy=whoutsampy)
  assign("lst.output", lst, envir=ewho)
 return(invisible(lst))
###  return(beta.hat.csid);
## end function gibbs.sampler
}
#####################################################################################
##
## FUNCTION NAME:  standard.errors
##
## PLACE:    
##
## IMPORTED:   gibbs.sampler function 
##
## DESCRIPTION: For every iteration of the Gibbs we obtain beta, sigma for
##              cntry and ages (i.e.csid).  Those currenta values are used to
##              obtain S.csid part of the variance for standar errors.
##              See Federico's manual for the theory behins and the Gauss code. 
##              Death is first obtained for insample and outsample according to
##              values of beta and whoinsampx, whoutsampx and with the noise of
##              gaussian distribution of 0 mean and sd=sigma.  S.csid
##              is the part of variance proportional to y^2
##              (or square of dth for each csid and year)
##
## INPUT:    current estimated values of beta=beta.hat.csid, sigma=sigma.hat.csid
##           whoinsampx, whoutsampx covarites for insample and outsample periods.
##           n, number of iterations in Gibbs, last value of S.csid (part of variance)
##                    
## OUTPUT:  S.csid updated value with new iteration of Gibbs, S.csid ~ y^2
##           
##
## WRITTEN BY: Elena Villalon & Federico Girosi    
##             evillalon@latte.harvard.edu;
##             girosi@rand.org; 
##             CBRSS, Harvard University
## 
## Last modified: 06/15/2004
## 

standard.errors <- function(n,beta.hat.csid,sigma.hat.csid,whoinsampx,whoutsampx,S.csid){
  
   yhatin <- make.forecast(beta.hat.csid,whoinsampx);
   yhatout <- make.forecast(beta.hat.csid,whoutsampx);
  
   yhat.csid <- lapply(names(whoinsampx),function(ch,yhatin,yhatout,sigma.hat.csid){
     y1 <- na.omit(yhatin[[ch]])
     y2 <- na.omit(yhatout[[ch]])
     y <- rbind(as.matrix(y1), as.matrix(y2))
     sigma <- na.omit(sigma.hat.csid[[ch]])
     ln <- nrow(y)
     if(length(sigma) <= 0) sigma <- 1
     wn= rnorm(1,mean=0,sd=sigma)
### wn <- as.matrix(rnorm(ln,mean=0,sd=sigma));
     y <- y + wn
     return(y)},yhatin,yhatout,sigma.hat.csid)
   
   yhat2.csid <- lapply(yhat.csid,function(mat) mat^2)
   names(yhat2.csid) <- names(whoinsampx)
   
   if ( length(unlist(S.csid)) <= 0 && n <= 1){
     S.csid <- yhat2.csid
    
     return(S.csid)
   }
   if (length(unlist(S.csid)) <= 0){
     S.csid <- lapply(yhat2.csid, function(x) x /n )
     
     return(S.csid)
   }
     ind <- 1:length(S.csid)
     names(ind) <- names(whoinsampx)
     S.csid <- lapply(ind, function(n,S.csid,yhat2.csid) {
       div <- n + 1.e-10
       y <- S.csid[[n]]*(1.0 - 1.0/div) + yhat2.csid[[n]]/div
       
       return(y)},S.csid, yhat2.csid)
   
  
     return(S.csid)   
}
## part of standar errors calculation after completing all iterations
## obtain variance of errors for every csid and years insample and outsample.

 variance.csid <- function(S,yin,yout){
      
      yhat2.csid <- lapply(names(yin),function(ch,yin,yout){
                  y1 <- na.omit(yin[[ch]])
    
                  y2 <- na.omit(yout[[ch]])
    
                  y <- rbind(as.matrix(y1), as.matrix(y2))
                  y <- y^2
                  return(y)},yin,yout)
      names(yhat2.csid) <- names(S)
     
      V <- lapply(names(S),function(ch,S,yhat2.csid)
                  {return(abs(S[[ch]] - yhat2.csid[[ch]]))},S,yhat2.csid)
      names(V) <- names(S)
      return(V)}
#####################################################################################
##
## FUNCTION NAME:  gibbs.fft
##
## PLACE:    
##
## IMPORTED:   gibbs.sampler function 
##
## DESCRIPTION: The gibbs.sampler produces a list with the coeff and death predictions
##              for the insample and outsample period for the number of sample, nsample
##              We take the fft for every csid elemnt of the coefficients and
##              compare two estimates due to different sample numbers. After taking
##              the Modulus of each beta we compare the spectral density for the
##              two samples coefficients.  
##         
## INPUT:    two set of coefficients, coeff1 and coeff2 for nsample=n1, n2
##           each set is a list of elemnts identifiers csid 
##                    
## OUTPUT:  Spectral densities for each csid of the coeff's; takes the ratio
##          for acuracy of sampling predictions. 
##           
##
## WRITTEN BY: Elena Villalon    
##             evillalon@latte.harvard.edu;  
##             CBRSS, Harvard University
## 
## Last modified: 06/15/2004
## 
## ************************************************************************
gibbs.fft <- function(coeff1, coeff2){
### Actually
### coeff1 <- gibbs.sampler(nsample=n1)$coeff
### coeff2 <- gibbs.sampler(nsample=n2)$coeff
  coeff1.fft <- lapply(coeff1, fft)
  coeff2.fft <- lapply(coeff2, fft)
  spec1 <- lapply(coeff1.fft, Mod)
  spec2 <- lapply(coeff2.fft, Mod)
  ind <- 1:length(spec1)
  names(ind) <- names(spec1)
  rate.fft <- lapply(ind,function(x) {
              spec1[[x]]/spec2[[x]]})
  return(rate.fft)}

### comparing death predictions for two samples and insample period
gibbs.hatin <- function(yhatin1, yhatin2){
### Actually
### yhatin1 <- gibbs.sampler(nsample=n1)$yhatin
### yhatin2 <- gibbs.sampler(nsample=n2)$yhatin
  ind <- 1:length(yhatin1)
  names(ind) <- names(yhatin1)
  rate.hatin <- lapply(ind,function(x) {
              yhatin1[[x]]/yhatin2[[x]]})
  return(rate.hatin)}

### comparing death predictions for two samples and insample period
gibbs.hatout <- function(yhatout1, yhatout2){
### Actually
### yhatout1 <- gibbs.sampler(nsample=n1)$yhatout
### yhatout2 <- gibbs.sampler(nsample=n2)$yhatout
  ind <- 1:length(yhatout1)
  names(ind) <- names(yhatout1)
  rate.hatout <- lapply(ind,function(x) {
              yhatout1[[x]]/yhatout2[[x]]})
  return(rate.hatout)}
