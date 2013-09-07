

cs.mean <- function(x,extended=FALSE){
  func <- function(z,extended){
    m <-  mean(as.data.frame(z),na.rm=TRUE);
    if (extended == TRUE) {
      z[T,T] <- 0; ### this makes a matrix of 0s preserving the names of z
      m <- scale(z, center=-m,scale=rep(1,ncol(z)));
    }
      return(m);
  }
  y <- lapply(x,FUN="func",extended);
  return(invisible(y))
}


## ************************************************************************
## ************************************************************************
##
## FUNCTION NAME:  impute.ts
##
## DESCRIPTION: it imputes missing values in a multivariate time series
##
## IMPORTED FUNCTIONS: 
##
## FORMAT: x.imp <- impute.ts(x)
##
## INPUT:        x:  (array) a multivariate time series with missing values 
##
## 
## OUTPUT:   x.imp: (array) <- the multivariate time series x with the
##                  missing values imputed by local linear interpolation.
##                  names of x are preserved.
##
##     
##
## WRITTEN BY: Federico Girosi 
##             fgirosi@latte.harvard.edu
##             CBRSS, Harvard University
## 
## Last modified: 08/14/2003
## 
## ************************************************************************
## ************************************************************************

impute.ts <- function(x){
  nr <- nrow(x);
  nc <- ncol(x);
  ind <- 1:nc;
  z <- apply(as.array(ind),MARGIN=1,FUN=function(n,x){s <- x[,n];
                                                      if (any(is.na(s))){
                                            s <- approx(1:length(s),s,xout=1:length(s),rule=2)$y;}
                                            return(s)},x)
  colnames(z) <-  colnames(x);
  rownames(z) <-  rownames(x);
  return(z);
}



## ************************************************************************
## ************************************************************************
##
## FUNCTION NAME:  impute.csts
##
## DESCRIPTION: it imputes missing values in a cross-sectional multivariate time series
##
## IMPORTED FUNCTIONS: 
##
## FORMAT: x.imp <- impute.csts(x)
##
## INPUT:        x:  (list) a cross-sectional multivariate time series with missing values 
##
## 
## OUTPUT:   x.imp: (array) <- the cross-sectional multivariate time series x with the
##                  missing values imputed by local linear interpolation.
##                  names of x are preserved.    
##
## WRITTEN BY: Federico Girosi 
##             fgirosi@latte.harvard.edu
##             CBRSS, Harvard University
## 
## Last modified: 08/14/2003
## 
## ************************************************************************
## ************************************************************************

impute.csts <- function(y){
  s <- lapply(y,FUN=impute.ts);
  return(invisible(s))
}


prior.variance <- function(ecxc,N=100,graphics=FALSE){
  ebase <- get("env.base", envir=parent.frame())
  env.base <- ebase
  ewho <- get("env.who", envir=ebase)
  whocov <- get("whocov", envir=ewho)
  verbose <- get("verbose", envir=ebase)
  whoinsampy <- get("whoinsampy", envir=ewho)
  who.age.digits <- get("who.age.digits", envir=ewho)
  ecxc <- get("env.cxc", envir=ebase)
  t.list  <- sapply(whoinsampy, function(x){
          x <- na.omit(x)
          nc <- nrow(x)});
  beta.dim  <- sapply(whocov, function(x){nc <- ncol(x)});
  cntry.vec <- get("cntry.vec", envir=ewho)
  age.vec <- get("age.vec", envir=ewho)
  age.char <- formatC(age.vec, width=who.age.digits, format="d", flag="0")
  who.mean.age.profile <- get("who.mean.age.profile", envir=ecxc)
  messout("Running prior.variance..",verbose)
  vlist <- ls(ecxc)
  for (v in vlist){
    assign(v,get(v,ecxc))
}
 
  S <- sample.improper.normal(D,1,N)
  D.pinv <- S$covar;
  beta <- S$x

## make covariates  matrix 
 
  Z <- matrix(0,nrow=sum(unlist(t.list)),ncol=sum(unlist(beta.dim)))

  row.max <- cumsum(unlist(t.list))
  
  row.min <- row.max-unlist(t.list)+1

  col.max <- cumsum(unlist(beta.dim))

  col.min <- col.max-unlist(beta.dim)+1

  Z.rownames <- rep(NA,nrow(Z))

  
  for (i in 1:length(whocov)){
    Z[row.min[i]:row.max[i],col.min[i]:col.max[i]] <- whocov[[i]]
    Z.rownames[row.min[i]:row.max[i]] <- rownames(whocov[[i]])
  }
  rownames(Z) <- Z.rownames

  var.vec <- diag(Z%*%D.pinv%*%t(Z))
  names(var.vec) <- rownames(Z)

  mu <- Z%*%beta
  
  aux <- split.matrix(mu,t.list)
  
  mu.time <- NULL
  for (i in 1:length(aux)){
    mu.time <- cbind(mu.time,aux[[i]]+who.mean.age.profile[i])
  }
  Ha.sigma <- get("who.Ha.sigma", envir=ewho)
  Ht.sigma <- get("who.Ht.sigma", envir=ewho)
  Hat.sigma <- get("who.Ha.sigma", envir=ewho)
  
  if(graphics){
    op <- par(cex.lab=1.4,omi=c(0,0.1,0,0),cex.main=1,bg="lavenderblush1",font=2,font.axis=2);
    matplot(mu.time[,1:N],type="l",ylim=c(-13,-10),main=paste("Ha.sigma=",round(Ha.sigma,2),
                                                     ", Ht.sigma=",round(Ht.sigma,2),
                                                     ", Hat.sigma",round(Hat.sigma,2)))
  }
  mu.age <- NULL
  for (i in 1:t.list[1]){
    mu.age  <- cbind(mu.age,mu[row.min+i-1,])
  }

  mu.age <- mu.age + who.mean.age.profile
  
  mu.age.time <- lapply(1:ncol(mu),
                        FUN=function(k,mu,t.list,who.mean.age.profile){
                          m <- matrix(mu[,k],nrow=t.list[1]);
                          m <-  scale(m,center=-who.mean.age.profile,scale=FALSE)},mu,t.list,who.mean.age.profile)

  
  return(list(mu.age.time=mu.age.time,var.vec=var.vec,mu.time=mu.time,mu.age=mu.age))
}

######################################################################

long.to.square.format <- function(mat){
  if (ncol(mat) != 3) stop("long.to.square.format argument must have 3 columns")
  rn <- paste(mat[,1],mat[,2],sep="_")
  if(length(rn) != length(unique.default(rn))){print(length(rn))
             print(length(unique.default(rn)));stop("long.to.square.format argument contains duplicates")}
  rownames(mat) <- rn
  x <- sort(unique.default(mat[,1]))
  y <- sort(unique.default(mat[,2]))
  z <- matrix(0,nrow=length(x),ncol=length(y))
  rownames(z) <- paste(x)
  colnames(z) <- paste(y)
  aux <- NULL
  for (i in 1:length(x)){
    for(j in 1:length(y)){
      z[i,j] <- mat[paste(x[i],y[j],sep="_"),3]
      aux <- rbind(aux,c(x[i],y[j],z[i,j]))
    }
  }

### sanity check
  rownames(aux) <- paste(aux[,1],aux[,2],sep="_")
  if(length(intersect(rownames(aux),rownames(mat))) != nrow(mat)) stop("mismatch")
  aux <- aux[rownames(mat),]
#  if (sum((aux-mat)^2) > 0)  stop("mismatch")
  return(list(x=x,y=y,z=z))
}

######################################################################


## ************************************************************************
## ************************************************************************
##
## FUNCTION NAME:      dt.da.stat
##
##
## DESCRIPTION: Use this to compute a summary statistics: absolute change in time trend
##              of log-mortality from one age group to the next:
##              [mu(a,t) - mu(a, t-1)] - [mu(a-1,t) - mu(a-1,t-1)]
##              
##              
##
## INPUT:   y: a list of lof-mortality matrices (T X A), one for each country
##
## OUTPUT: a vector of summary statistics, ready to be histogrammed
##
## WRITTEN BY: Federico Girosi, 7/19/2005
## 
## ********************************************************************
## ********************************************************************


dt.da.stat <- function(mu.age,mu.time,mu.age.time,mat=NA){
  f <- function(mu,mat){
    dmu <- mat %*% mu
    dmu.mean <- mean(apply(dmu,1,FUN=function(x){mean(abs(diff(x)))}))}
  return(mean(sapply(mu.age.time,FUN=f,mat)))
}




######################################################################

## ************************************************************************
## ************************************************************************
##
## FUNCTION NAME:      sd.stat
##
##
## DESCRIPTION: Use this to compute a summary statistics: age-specific standard deviation (over time)
##              
##
## INPUT:   y: a list of lof-mortality matrices (T X A), one for each country
##
## OUTPUT: a vector of summary statistics, ready to be histogrammed
##
## WRITTEN BY: Federico Girosi, 11/5/2005
## 
## ************************************************************************
## ************************************************************************

sd.stat <- function(mu.age.time){
  SD <- NULL
  for (i in 1:length(mu.age.time)){
    m <-  mu.age.time[[i]]
    SD <- c(SD,as.vector(apply(m,2,sd)))    
  }
  return(SD)
}


## ************************************************************************
## ************************************************************************
##
## FUNCTION NAME:      d1.a.stat
##
##
## DESCRIPTION: Use this to compute a summary statistics: change in log-mortality
##              from one age group to the next (averaged over time)
##              
##
## 
## ************************************************************************
## ************************************************************************

d1.a.stat <- function(mu.age,mu.time,mu.age.time,mat=NA){
  return(mean(apply(mu.age,2,FUN=function(x){mean(abs(diff(x)))})))
}




## ************************************************************************
## ************************************************************************
##
## FUNCTION NAME:      d1.t.stat
##
##
## DESCRIPTION: Use this to compute a summary statistics: change in log-mortality
##              from one year to the next (averaged over age groups and countries)


d1.t.stat <- function(mu.age,mu.time,mu.age.time,mat=NA){
  return(mean(apply(mu.time,2,FUN=function(x){mean(abs(diff(x)))})))
}


find.best.sigma <- function(G,d1.a,d1.t,dtda,SD,
                            n.row=25,summary.measures=c("SD","d1.a","d1.t","dtda"), verbose=T){
   verb <- try(get("verbose", envir=parent.frame()),silent=T)
   if(class(verb)!="try-error")
     verbose <- verb
  rownames(G) <- 1:nrow(G)
  G <- as.data.frame(G)
  G <- subset(G, SD < 10)
  X <-  G[,c("SD","d1.a","d1.t","dtda")]
  sd <- sqrt(diag(var(X)))
  mean.x <- colSums(X)/nrow(X)
  w <- sd/mean.x
  messout(summary.measures, verbose)

  X <-  scale(X,center=F,scale=c(SD,d1.a,d1.t,dtda))
  X <-  scale(X,center=c(1,1,1,1),scale=F)
  excluded <- setdiff(c("SD","d1.a","d1.t","dtda"),summary.measures)
 if (length(excluded) > 0) X[,excluded] <- 0
#  w[excluded] <- 0
#  w["SD"] <- 0
#  print(X)
#  print(w)
 dist <- apply(X,MARGIN=1,FUN=function(x,w){sum(abs(x)^2*w)},w)
  z <-  sort(dist)
  if(n.row > nrow(G)) n.row <- nrow(G)
  G.small <- as.matrix(G[names(z)[1:n.row],])
  rownames(G.small) <- NULL
###  print(G.small)
  return(G.small)  
}


######################################################################


## ************************************************************************
## ************************************************************************
##
## FUNCTION NAME:      dt.da.stat
##
##
## DESCRIPTION: Use this to compute a summary statistics: absolute change in time trend
##              of log-mortality from one age group to the next:
##              [mu(a,t) - mu(a, t-1)] - [mu(a-1,t) - mu(a-1,t-1)]
##              
##
## INPUT:   y: a list of lof-mortality matrices (T X A), one for each country
##
## OUTPUT: a vector of summary statistics, ready to be histogrammed
##
## WRITTEN BY: Federico Girosi, 7/19/2005
## 
## ************************************************************************
## ************************************************************************

dt.da.stat.old <- function(mu.age.time){
  dtda <- NULL
  for (i in 1:length(mu.age.time)){
    y <-  mu.age.time[[i]]
    mat <- matrix(0,nrow=nrow(y)-1,ncol=nrow(y))
    diag(mat) <-  -1
    for (k in 1:nrow(mat)){
      mat[k,k+1] <- 1
    }
    dmu <- mat %*% y
    dtda <- c(dtda,apply(dmu,1,FUN=function(x){mean(abs(diff(x)))}))
  }
  return(dtda)
}

######################################################################

## ************************************************************************
## ************************************************************************
##
## FUNCTION NAME:      sd.stat.old
##
##
## DESCRIPTION: Use this to compute a summary statistics: age-specific standard deviation (over time)
##              
##
## INPUT:   y: a list of lof-mortality matrices (T X A), one for each country
##
## OUTPUT: a vector of summary statistics, ready to be histogrammed
##
## WRITTEN BY: Federico Girosi, 11/5/2005
## 
## ************************************************************************
## ************************************************************************

sd.stat.old <- function(mu.age.time){
  SD <- NULL
  for (i in 1:length(mu.age.time)){
    m <-  mu.age.time[[i]]
    SD <- c(SD,as.vector(apply(m,2,sd)))    
  }
  return(SD)
}


## ************************************************************************
## ************************************************************************
##
## FUNCTION NAME:      d1.a.stat.old
##
##
## DESCRIPTION: Use this to compute a summary statistics: change in log-mortality
##              from one age group to the next (averaged over time)
##              
##
## INPUT:   y: a list of lof-mortality matrices (T X A), one for each country
##
## OUTPUT: a vector of summary statistics, ready to be histogrammed
##
## WRITTEN BY: Federico Girosi, 11/5/2005
## 
## ************************************************************************
## ************************************************************************

d1.a.stat.old <- function(mu.age.time){
  d1.a <- NULL
  for (i in 1:length(mu.age.time)){
    m <-  mu.age.time[[i]]
    d1.a  <- c(d1.a,as.vector(apply(m,1,FUN=function(x){mean(abs(diff(x)))})))
  }
  return(d1.a)
}



## ************************************************************************
## ************************************************************************
##
## FUNCTION NAME:      d1.t.stat
##
##
## DESCRIPTION: Use this to compute a summary statistics: change in log-mortality
##              from one year to the next (averaged over age groups and countries)
##              
##
## INPUT:   y: a cross-sectional list of log-mortality time series
##
## OUTPUT: a vector of summary statistics, ready to be histogrammed
##
## WRITTEN BY: Federico Girosi, 11/5/2005
## 
## ************************************************************************
## ************************************************************************

d1.t.stat.old <- function(y.csid){
  d1.t <- NULL
  for (i in 1:length(y.csid)){
    m <-  y.csid[[i]]
    d1.t  <- c(d1.t, mean(abs(diff(m))))
  }
  return(d1.t)
}




Xy.only <- function(ebase=env.base){
   ebase <- get("env.base", envir=parent.frame())
   env.base <- ebase
   ewho <- get("env.who", envir=ebase)
  whoinsampy <- get("whoinsampy", envir=ewho)
  whoutsampy <- get("whoutsampy", envir=ewho)
  whoinsampx <- get("whoinsampx", envir=ewho)
  whoutsampx <- get("whoutsampx", envir=ewho)
  who.zero.mean <- get("who.zero.mean",envir=ewho)   
  who.ols.sigma.param <-  get("who.ols.sigma.param", envir=ewho)
  clist <- vector(mode="list",length=length(whoinsampy));
  names(clist) <- names(whoinsampy)
  coeff <- clist; 
  std   <- clist;  
  sigma <- make.sigma(who.ols.sigma.param);
### I need the following list for the computations of Gibbs' sample
   Xy  <- clist;

  
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
    
   
  for (i in 1:length(whoinsampy)){
    y <- whoinsampy[[i]];
    x <- whoinsampx[[i]];    
    yw <- whoinsampy[[i]]/sigma[[i]];
    xw <- t(scale(t(whoinsampx[[i]]),center=FALSE,scale=sigma[[i]]));
    yxw <- na.omit(cbind(yw, xw));
    yw <- yxw[,1]
    xw <- yxw[,-1]
    txw <- t(xw);
    Xy[[i]] <- txw %*% yw
  }
  return(Xy);       
}
## ************************************************************************
##
## FUNCTION NAME:      prior.param
##
## IMPORTED functions: 
##
##
## DESCRIPTION: it returns the smoothness parameters for smoothing over age,
##              time and age/time for one particular country and one particular
##              set of covariates. 
##              
##
## FORMAT:  G <- prior.param(Ha.sigma.vec=seq(0.1,1.5,length=5),
##                        Ht.sigma.vec=seq(0.1,1.5,length=5),Hat.sigma.vec=seq(0.05,1.5,length=5),
##                        stats=c(d1.a=d1.a.stat,d1.t=d1.t.stat,dtda=dt.da.stat),sims=20,
##                        autoset = c(d1.a=0.43,d1.t=0.022,dtda=0.004,SD=0.2))
##
## INPUT: 
##          Ha.sigma.vec: vector of values to try for the standard deviation of
##                        the prior over ages
##          Ht.sigma.vec: vector of values to try for the standard deviation of
##                        the prior over time
##          Hat.sigma.vec: vector of values to try for the standard deviation of
##                         the prior over age and time
##          stats: vector of functions which compute summary statistics
##                 at the moment only the defaults are available
##          sims: number of samples to draw from the prior in order to compute the standard deviation
##          autoset vector with the following elements
##                  d1.a: for the summary measure d1.a.stat (difference in log-mortality
##                        between adjacent age groups)
##                  d1.t: target for the summary measure d1.t.stat (difference in log-mortality
##                        between adjacent years)
##                  dtda: target for the summary measure dt.da.stat (difference in log-mortality
##                        between difference in adjacent age groups between adjacent years)
##
## OUTPUT: G: matrix with n.row rows. In each row we report the values of
##            the smoothness parameters and the corresponding summary measures.
##            The rows are ordered so that the first row contains the "best" set
##            of smoothness parameters (that is the set of smoothness parameters
##            which best reproduces the target values)
##
## WRITTEN BY: Federico Girosi
##             girosi@rand.org
##             (EV updated it with the current version of yourcast)
##             
## 
## Last modified: 2/8/2006
## 
## ************************************************************************
prior.param <- function(Ha.sigma.vec=seq(0.1,1.5,length=5), Ht.sigma.vec=seq(0.1,1.5,length=5),
                        Hat.sigma.vec=seq(0.05,1.5,length=5),sims=20, 
                        stats=c(d1.a=d1.a.stat,d1.t=d1.t.stat,dtda=dt.da.stat), 
                        autoset=c(d1.a=0.43,d1.t=0.022,dtda=0.004,SD=0.2))
  
{
 
 ebase <- get("env.base", envir=parent.frame())
 ewho <- get("env.who",ebase)
 env.base <- ebase
 env.who <- ewho
 verbose <- get("verbose", envir=ebase)
 summary.measures <- get("summary.measures", envir=env.base)
 if(is.list(autoset) && length(autoset) > 0){
   names(autoset) <- NULL
   autoset <- unlist(autoset)
 }
 SD.target <- autoset["SD"]
 d1.a.target <- autoset["d1.a"]

 d1.t.target <- autoset["d1.t"]
 dtda.target <- autoset["dtda"]
 
 G <- expand.grid(Ha.sigma=Ha.sigma.vec,Ht.sigma=Ht.sigma.vec,Hat.sigma=Hat.sigma.vec)
 aux <- matrix(NA,nrow=nrow(G),ncol=length(stats))
 colnames(aux) <- names(stats)
 G <- cbind(G,SD=rep(NA,nrow(G)),aux)
 ###print(G)
 mat <-  NA
 messout(paste("Number of samples per prior ", sims,sep=""),verbose) 
 for (i in 1:nrow(G)){
   Ha.sigma <- G[i,"Ha.sigma"]
   Ht.sigma <- G[i,"Ht.sigma"]
   Hat.sigma <- G[i,"Hat.sigma"]
  messout(paste("Ha.sigma = ",Ha.sigma, " Ht.sigma = ",Ht.sigma," Hat.sigma = ", Hat.sigma, sep=""),verbose)
   
   assign("who.Ha.sigma",Ha.sigma,envir=env.who)
   assign("who.Ht.sigma",Ht.sigma,envir=env.who)
   assign("who.Hat.sigma",Hat.sigma,envir=env.who)    
   
   m <- cxc()
   ecxc <- m$ecxc
  
### m contains obj= yhatin, yhatout, insampy, outsampy but listed with csid
### conversion.cntry.mat uses list.by.cntry to convert into country list matrices
### so it is a more useful format for the graphics output.  
   m <- conversion.cntry.mat(m); 
   P <- prior.variance(ecxc,N=sims,graphics=F)

   V <- P$var.vec
   mu.age <- P$mu.age
   mu.time <- P$mu.time
   mu.age.time <- P$mu.age.time
   
    if (identical(mat,NA)){
      mat <- matrix(0,nrow=nrow(mu.time)-1,ncol=nrow(mu.time))
      diag(mat) <-  -1
      for (k in 1:nrow(mat)){
        mat[k,k+1] <- 1
      }
    }

    G[i,"SD"] <- sqrt(mean(V))    
    for (j in 1:length(stats)){
      f <- stats[[j]]
      G[i,names(stats)[j]] <- f(mu.age,mu.time,mu.age.time,mat=mat)
    }
 }

   n.row    <- get("n.row", envir=ewho)
   digits   <- get("digits", envir=ewho)
   filename <- get("filename", envir=ewho)
 
   gsmall <- find.best.sigma(G,d1.a.target,d1.t.target,dtda.target,SD.target,
                          n.row=n.row,summary.measures=summary.measures)
 ###  if (!is.na(filename))
 ###    try(latex.default(round(gsmall,digits),file=paste(whooutpath,filename,sep="")), silent=T)
  return(gsmall);
}

######################################################################
### DESCRIPTION consistency checks to run model ebayes
consistent.checks.ebayes <- function(autoset0, summary.measures0,
                                     env.who=parent.frame())
    {
      autoset <- autoset0
      summary.measures <- summary.measures0
      
      if(class(autoset0 <- try(get("autoset", envir=env.who), silent=T)) !="try-error")
        autoset <- autoset0
      else
        assign("autoset", autoset0, envir=env.who)
  
      if(class(summary.measures0 <- try(get("summary.measures", envir=env.who), silent=T)) !="try-error")
        summary.measures <- summary.measures0
      else
        assign("summary.measures", summary.measures0, envir=env.who)
  
      if(length(na.omit(autoset)) != length(summary.measures))
        stop("You  must provide target values for summary measures")
       
      if(length(grep("SD", names(autoset))) <= 0 && length(autoset) > 0)
        if( "SD" %in% summary.measures )
          stop("you must pass target for SD")
        else 
          autoset <- c(autoset, SD=NA)     
        
      
      lst <- list(summary.measures=summary.measures, autoset=autoset)
      return(lst)
    }

### DESCRIPTION generates histograms for the input vectors d1.a, d1.t, dt.da, SD
###
#### INPUT: vectors d1.a ,d1.t, dtda , SD, with the statistics  
###         depvar, chararacter string with the name of dependent variable
###         graphics.file, output.file graphics parameters
###
###  OUTPUT the hsitograms plots for d1.a, d1.t, dt.da and SD
### AUTHOR: Federico Girosi
###         girosi@rand.org
###################################################################
histograph <- function(d1.a, d1.t, dt.da, SD, depvar=" ",
                       model="ebayes", graphics.file=NA)
  {
    whooutpath <- getwd()
    main.title <- paste(depvar,sep=" ")
    oop <- par(cex.lab=1.4,omi=c(0,0.1,0,0),cex.main=1.4,bg="lavenderblush1",font=2,font.axis=2)
    op  <- par(mfrow=c(2,2))

    z <- hist(d1.a,col="firebrick2",xlab="delta age",main=main.title)
    mean.text <- paste("Mean:",round(mean(d1.a),digits=2))
    median.text <- paste("Median:",round(median(d1.a),digits=2))
    text(x=max(z$breaks), y=max(z$counts),mean.text,pos=2)

    z <- hist(d1.t,col="firebrick2",xlab="delta time",main=main.title)
    mean.text <- paste("Mean:",round(mean(d1.t),digits=3))
    median.text <- paste("Median:",round(median(d1.t),digits=3))
    text(x=max(z$breaks), y=max(z$counts),mean.text,pos=2)
  
    z <- hist(dt.da,col="firebrick2",xlab="delta age/time",main=main.title)
    mean.text <- paste("Mean:",round(mean(dt.da),digits=3))
    median.text <- paste("Median:",round(median(dt.da),digits=3))
    text(x=max(z$breaks), y=max(z$counts),mean.text,pos=2)
  
    z <- hist(SD,col="firebrick2",xlab="Standard Deviation",main=main.title)
    mean.text <- paste("Mean:",round(mean(SD),digits=2))
    median.text <- paste("Median:",round(median(SD),digits=3))
    text(x=max(z$breaks), y=max(z$counts),mean.text,pos=2)
    strata <- NULL
    if(length(grep("[0-9]+", depvar)) > 0)
      {
        vesp <- strsplit(depvar, NULL)[[1]]
        ix <- grep("[0-9]+", vesp)
        ix <- max(ix)
        strata <- vesp[ix]
      }
    par(op)
    if(identical(graphics.file,NA)) graphics.file <- paste("empirical",depvar,strata,model,sep = "_")
    if (is.character(graphics.file)){
      try(savePlot(filename=paste(whooutpath,graphics.file,sep=""),type="pdf"), silent=T)
    }
  }

### DESCRIPTION calculates the vectors d1.a.vec , d1.t.vec, dtda.vec, SD.vec
### and their mean values from the results outputs of model=cxc() in yourcast
### It uses the functions: conversion.cntry.mat, d1.a.stat.old,
### d1.t.stat.old, dt.da.stat.old, left.side.formula
###
### FORMAT:  G <- prior.param(m, formula, graphics.file=NA,output.file=NA)
###
### INPUT m the outputlist of running cxc(_) function in yourcast
###       formula the input formula in yourcasst
###       parameters for graphics files
###
### OUTPUT
###        G: list with the following elements:
##                   m: the output of yourcast (same as input)
##                   summary: vector with the mean of the 4 summary measures.
##                            The summary measures are delta age, delta tme, delta age/time, satandard deviation.
##                   d1.a: vector with the values of delta age
##                   d1.t: vector with the values of delta time
##                   d1.t: vector with the values of delta age/time
##                   SD: vector with the values of the standard deviation (one for each age/country combination)
##
## WRITTEN BY: Federico Girosi 
##             girosi@rand.org,
##             (updated by EV for the latest version of yourcast) 
##             
## Last modified: 2/9/2006
## 
## ************************************************************************
## ************************************************************************


empirical.bayes <- function(m, formula, depvar, graphics.file=NA,verbose=T)
  {
   
    env.base <- get("env.base", envir=parent.frame())
    ebase <- env.base
    verbose <- get("verbose", envir=ebase)
    messout("Running empirical bayes...",verbose)
    m <- conversion.cntry.mat(m);
    insampy <- m$insampy
    yhatin <- m$yhatin
    yhatin.csid <- list.by.csid(yhatin)
            
    SD.vec <- sd.stat(yhatin)
    d1.a.vec <- d1.a.stat.old(yhatin)
    d1.t.vec <- d1.t.stat.old(yhatin.csid)
    dtda.vec  <- dt.da.stat.old(yhatin)
 
    mean.vec <- c(SD=mean(SD.vec),d1.a=mean(d1.a.vec),d1.t =mean(d1.t.vec),dtda=mean(dtda.vec))
    return(list(model.output=m, summary=mean.vec,d1.a.vec=d1.a.vec,
                d1.t.vec=d1.t.vec,dtda.vec=dtda.vec,SD.vec=SD.vec));
  }
### DESCRIPTION: Model ebayes calculations of relevant parameters and histogram plots
###              call the functions:
###              consistent.checks.ebayes, empirical.bayes
###
### INPUT: autoset, a vector with default values for SD, d1.a, d1.t, dtda
###                 d1.a: vector with the values of delta age
###                 d1.t: vector with the values of delta time
###                 dtda: vector with the values of delta age/time
###                 SD: vector with the values of the standard deviation
###                     (one for each age/country combination)
###        formula, the formula for regression
###        m, a list aoutput of running model MAP: the output of cxc()
###        summary.measures a vector of characters with elemnts
###                         the summary of statistics from the model
###        env.who, the environment wher all parameters are stored.
###
### OUPUT the list with the stats calculations for the elements of autoset
###       mean values and distributions
###       summary.mean vector with means d1.a, d1.t, dtda, SD
###       summary.vec vector with distributions for d1.a, d1.t, dtda, SD
###
### AUTHOR Federico Girosi
###        Rand, CA
###         girosi@rand.org
###       (modified by EV to update with latest version of yourcast)
###########################################################################################

  model.ebayes <- function(autoset,formula,m, summary.measures,depvar,
                           graphics.file,env.who, show.histo=FALSE)
                             
  {
    
    ebase <- get("env.base", envir=parent.frame())
    env.who <- get("env.who", envir=ebase)
    env.base <- ebase
    model <- get("model", envir=ebase)
    
 ###initialization
    summary <- c(SD=autoset$SD,d1.a=autoset$d1.a,d1.t=autoset$d1.t,dtda=autoset$dtda)
    dist.vec <- NULL
    lst.ret <- list(summary=summary)
  
    z <- try(empirical.bayes(m, formula, depvar, graphics.file), silent=T)
 
 
    if(class(z) != "try-error")
      {
        summary <- z$summary
###    print(summary)
        d1.a  <- summary["d1.a"]
        d1.t  <- summary["d1.t"]
        dtda <- summary["dtda"]
        SD <- summary["SD"]
        summary <- c(list(SD=SD),list(d1.a=d1.a),list(d1.t=d1.t),list(dtda=dtda))
        dist.vec <- c(list(d1.a =z$d1.a.vec),list(d1.t= z$d1.t.vec),
                      list(dtda=  z$dtda.vec), list(SD= z$SD.vec))
        lst.ret <- list(summary=summary, summary.vec = dist.vec)

        if(show.histo)
          histograph(z$d1.a.vec, z$d1.t.vec, z$dtda.vec,z$SD.vec, depvar,model, graphics.file)
      }
    
    return(lst.ret)
  }

### DESCRIPTION: helper function to yourcast; it calls prior.param that Federico wrote
###              and performs some sanity checks to make sure that alll smoothness
###              parameters either they have a value or NA
### CALL prior.param
### INPUT smooth, vector with the values for smoothness parameters, Ha.sigma,Ht.sigma,Hatsigma
###       count.cntry, number of countries in data set
###       Ha.sigma.vec, Ht.sigma.vec, Hat.sigma.vec, vectors of possible values for smoothness paramters
###       sims, integer
###       stats, vector with function names
###       autoset, vector with the values for SD, d1.a, d1.t, dtda
###       env.who, an environmnet
### OUTPT matrix G with values for smoothness parameters; the output of p[rior.param
### 05/04/2006

    
build.prior.param <- function(smooth,count.cntry, Ha.sigma.vec,
                              Ht.sigma.vec, Hat.sigma.vec, sims,
                              stats, autoset, model, env.who)
  {
    env.base <- get("env.base", envir=parent.frame())
    verbose <- get("verbose", envir=env.base)
    G <- NULL
    vsigma <- c("Ha.sigma", "Ht.sigma", "Hat.sigma")
    ln.smooth <- length(smooth)
    
### ln.smooth < 3 means that one of Ha.sigma, Ht.sigma or Hat.sigma is NULL
### then if model map and only one country in dataset runs prior.param
### otherwise fill the NULL value with NA and print a message
    if(ln.smooth < 3 ) 
      if(length(count.cntry) <= 1 && identical(model, "MAP"))
        G <- prior.param(Ha.sigma.vec, Ht.sigma.vec,Hat.sigma.vec,sims,stats, autoset)
      else {
        messout("Setting smoothness parameters that are NULL equal to NA",verbose)
        ix <- match(names(smooth), vsigma)
        nona <- 3 - length(na.omit(ix))
        smooth <- c(smooth, rep(NA, nona))
        names(smooth) <- c(names(smooth)[1:ln.smooth], vsigma[-ix])
        for(n in 1:length(smooth)){
          ass <- paste("who.",names(smooth)[n], sep="")
          if(n > ln.smooth)
            assign(ass,smooth[n], envir=env.who)
        }
      }
    if(length(G) <= 0)
      G <- rbind(G, smooth) 
    
    return(G)
  }
### DESCRIPTION helper function to yourcast; it is called when only one country, model=map
###             Loops over the rows of matrix G extracting, for each row, the values
###             of the smoothing parameters, Ha.sigma, Ht.sigma, Hat.sigma and runs model cxc()
### USES cxc()
###
### OUTPUT: The first result of model cxc, which is not a try-error and correspond to lowest row of G.
### run with Federico Girosi code
### (Elena Villalon, 05/04/2006)
##############################################################
onecntry.bayes.empirical <- function(G)
  {
    ebase <- get("env.base", envir=parent.frame())
    env.who <- get("env.who", envir=ebase)
    env.base <- ebase
    m <- list()
    
    maxprint <- get("maxprint", envir=env.who)
    flag <- 0 
    for(i in   1:(min(nrow(G),maxprint)))
      {
       
        Ha.sigma <- try(G[i,"Ha.sigma"],silent=T)
        Ht.sigma <- try(G[i,"Ht.sigma"], silent=T)
        Hat.sigma <- try(G[i,"Hat.sigma"], silent=T)
        smooth <- NULL
        if(class(Ha.sigma) != "try-error"){
          assign("who.Ha.sigma", Ha.sigma, envir=env.who)
          smooth <- c(smooth, Ha.sigma=Ha.sigma)
        }
        if(class(Ht.sigma) != "try-error"){
          assign("who.Ht.sigma", Ht.sigma, envir=env.who)
          smooth <- c(smooth, Ht.sigma=Ht.sigma)
        }
        if(class(Hat.sigma) != "try-error"){
          assign("who.Hat.sigma", Hat.sigma, envir=env.who)
          smooth <- c(smooth, Hat.sigma=Hat.sigma)
        }
       
        ml <- try(cxc(ebase=env.base),silent=T)
        
        if(class(ml) != "try-error" && flag <= 0){
          m <- ml
          lst <- c(list(m=m), list(params=smooth), list(G=G))
          flag <- 1
        }
        
      }
    return(lst)
  }



######################################################################
###
### AUTHOR Federico Girosi
###        Rand, CA
###         girosi@rand.org
### Federico wrote the utility functions to join matrices and list of matrices
### for calculating the std of MAP model
### this function concatenates 2 matrices diagonally
concat.mat <-  function(mat1,mat2){
  newmat <-  matrix(0,nrow=nrow(mat1)+nrow(mat2),ncol=ncol(mat1)+ncol(mat2))
  rownames(newmat) <- 1:nrow(newmat)
  colnames(newmat) <- 1:ncol(newmat)
  
  newmat[1:nrow(mat1),1:ncol(mat1)] <- mat1
  rownames(newmat)[1:nrow(mat1)] <- rownames(mat1)
  colnames(newmat)[1:ncol(mat1)] <- colnames(mat1)
  
  newmat[(nrow(mat1)+1):nrow(newmat),(ncol(mat1)+1):ncol(newmat)] <- mat2
  rownames(newmat)[(nrow(mat1)+1):nrow(newmat)] <- rownames(mat2)
  colnames(newmat)[(ncol(mat1)+1):ncol(newmat)] <- colnames(mat2)
  return(newmat)
}

######################################################################

### this function concatenates a list of  matrices diagonally
diag.mat <- function(matlist){
  if (length(matlist) == 1) return(matlist[[1]])
  newmat <- concat.mat(matlist[[1]],matlist[[2]])
  if (length(matlist) == 2) return(newmat)
  for (i in 3:length(matlist)){
    newmat <- concat.mat(newmat,matlist[[i]])
  }
  return(newmat)
}

######################################################################

### Module for ols/gaussian linear models:
### ols: normal linear regression, driver ols()
### ols.functional: same as ols, but uses glm( and slower); driver ols.functional()
### Contains following functions:
### covdel(xmat, tol=tol) eliminates colinearities of xmat for tolerance=tol
### eliminate.colinear(xw,nmx, delta.tol, tol) call covdel recursively
### inverse.trans(yin,yout,y,whotransform) inverse to trans (make.mortality.data)
### ols() calls covdel,eliminate.colinear and inverse.trans
### ols.functional() calls covdel, glm.colinear and inverse.trans
### glm.colinear(xw, yw, nmx, delta.tol, tol) uses glm (built in function)
### and checks for colinearities with covdel 
###
## ************************************************************************
## ************************************************************************
##
## FUNCTION NAME: covdel
##
##
## INPUT:     elemnt of whoinsampx, matrix of covariates for given csid (cntry+age group)
##            rows = years of observations insample matrix and cols=covariates; 
##            tol, tolerance factor for linear dependencies among cols whoinsampx
## OUTPUT:    matrix whoinsampx modified after eliminating linear dependencies, 
##            or those cols of covariates which are correlated with tolerance >= tol; 
##            also the index of those cols eliminated by colinearity;
##            and the names of the eliminated cols 
##
##  WRITTEN BY: Elena Villalon
##              evillalon@latte.harvard.edu
##              CBRSS, Harvard University
##
## Last modified: 05/12/2003
## 
## ***********************************************************************
## ***********************************************************************

covdel <- function(xmat, tol=tol, ebase=env.base) {
  ebase <- get("env.base", envir=parent.frame());
  env.base <- ebase;
  ewho <- get("env.who", envir=ebase)
  who.digit.first  <- get("who.digit.first", envir=ewho) 
  who.cntry.digits <- get("who.cntry.digits", envir=ewho)
  who.age.digits   <- get("who.age.digits",envir=ewho)
  who.year.digits  <- get("who.year.digits", envir=ewho)
if (dim(xmat)[2] <= 1)
  return(list(covx=xmat,covind=NULL,covnam = NULL))
  
x  <- data.frame( xmat)
vn <- colnames(x)
bf <- length(x)
### find cols with constant elements of sd =0
ind0 <- seq(along=x)[sd(x,na.rm=T)==0]
### index for cols other than "cnst" ( last?)
if(length(ind0) > 0){
cn <- grep("cnst", vn[ind0])
cn <- cn[length(cn)]
cti  <- ind0[-cn]
x    <- x[-ind0]}

mat <- cor(x, use="pairwise")
r <- row(mat)[row(mat) < col(mat) & abs(mat) >= tol]
r <- na.omit(r)
if(length(r) > 0)  {
  x <- x[-r]
  ri <- c(r, cti)
  rn <- vn[ri]
} else{
  ri <- cti 
  rn <- vn[cti]}

x  <- cbind(x, "cnst"=1)
af <-  length(x) 
csid <-  as.numeric(row.names(x)[1])%/%10^(who.year.digits)
### Because of the I/O we only print once the message below
if( abs(bf -af) >0 & csid %% 10^(who.age.digits) == 45) 
print(paste("Colinears, ncols deleted = ", bf - af, " for tol= ", tol,
            "; and cntry= ",csid%/%10^(who.age.digits)))
if (length(ri) > 0) xmat <- xmat[,-ri]
### returns covariate matrix with colinearities removed if appropiate
### index of columns being removed, and their names
return(list(covx=xmat,covind=ri,covnam = rn))
}  
