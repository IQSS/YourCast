lc <- function(ebase){
 
  ebase <- get("env.base", envir=parent.frame())
  env.base <- ebase
  ewho <- get("env.who", envir=ebase);
  verbose <- get("verbose", envir=ebase)
 
  age.vec <- get("age.vec", envir=ewho)

  who.cntry.digits <- get("who.cntry.digits", envir=ewho)
  who.age.digits <- get("who.age.digits", envir=ewho)
  who.year.digits  <- get("who.year.digits", envir=ewho)
  who.digit.first  <- get("who.digit.first", envir=ewho)
  whoyrest <- try(get("whoyrest", envir=ewho), silent=T)
  if(class(whoyrest)=="try-error"){
    yrest <- get("yrest", envir=ebase)
    whoyrest <- yrest
  }
  
  whoinsampy <- get("whoinsampy", envir= ewho)
  whoutsampy <- get("whoutsampy", envir=ewho)
 
### to otput with model 
 
  cntry.vec <- get("cntry.vec", envir=ewho)
  n.cntry <- length(cntry.vec)
  age.vec <- get("age.vec", envir=ewho)
  cntry.names.lst <- get("cntry.names.lst", envir=ewho)
  
### the structure of dataset cstsid: 
  digit.cntry.begin <- who.digit.first + 1 
  digit.cntry.end   <- who.digit.first + who.cntry.digits
 
  whofore <- max(unlist(lapply(whoutsampy , function(mat,whoyrest ) {nr <- nrow(mat) + whoyrest}, whoyrest)))
          
      
  messout("Using LC model", verbose);
  cs.vec <- as.numeric(names(whoinsampy));
  cs.cntry.vec <- digitpull(cs.vec,digit.cntry.begin,digit.cntry.end);
  nfore <- whofore - whoyrest;
  insampy.c <- list.by.cntry(whoinsampy);
  yhatin.c <- vector(mode="list",length=length(insampy.c));
  yhatout.c <- vector(mode="list",length=length(insampy.c));
  gamma.c <- vector(mode="list",length=length(insampy.c));
  beta.c <- vector(mode="list",length=length(insampy.c));
  names(yhatin.c) <- names(insampy.c);
  names(yhatout.c) <- names(insampy.c);  
  for (i in 1:n.cntry){
    l <- lc.model(insampy.c[[i]],nfore);
    yhatin.c[[i]] <- l$yhatin;
    yhatout.c[[i]] <- l$yhatout
    gamma.c[[i]] <- l$gamma;
    beta.c[[i]] <- l$beta;
  }
    yhatin <- list.by.csid(yhatin.c);
    yhatout <- list.by.csid(yhatout.c);
    coeff <- list(gamma=gamma, beta=beta.c);
    model <- model.string()
    lst <- list(yrest=whoyrest,model=model,age.vec=age.vec, cntry.lst=cntry.names.lst,
                coeff=coeff,yhatin=yhatin,yhatout=yhatout, insampy =whoinsampy, outsampy=whoutsampy)
    assign("lst.output", lst, envir=ewho)
    return(lst);
}



## ######################################################################################
## ######################################################################################
##
## FUNCTION NAME: lc.model
##
## DESCRIPTION: it fits a LC model to a log-mortality matrix and forecasts up to a
##              a given number of years ahead
##
## IMPORTED FUNCTIONS: none
##
## GLOBALS:   none
##
##
## FORMAT: lst <- lc.model(mat,horizon)
##
## INPUT:      mat: (matrix) a T x A matrix of log-mortality rates. The rows must be have names
##                   corresponding to the years. MIssing values are allowed
##
##         horizon: (scalar) how many years we want to predict 
##
##
## 
## OUTPUT:      lst: (list) an object of type "LC forecast" with the
##                   following non-empty fields:
##
##                  m: (matrix) the original data matrix
##         m.centered: (matrix) the centered data matrix
##     y.insample.hat: (matrix) the LC prediction for the insample period
##             yhatin: (matrix) the fit of the first principal component to the data,
##                     with the missing values filled in using the LC insample
##                     prediction
##            yhatout: (matrix) the LC predictions for the outsample period,
##                     which start one year after the end of the insample period
##                     and ends horizon years later
##               beta: the first principal component, normalized to have length 1,
##                     which represents the rate of decrease (increase) of log-mortality
##               gamma: the projection of the age profiles on the first
##                      principal component (the "force" of mortality)
##  gamma.insample.fit: same sa gamma, but with missing values replcaed by the LC
##                      insample prediction of gamma
##     gamma.outsample: the forecast of gamma, given by the random walk with drift model
##               theta: the drift parameter in the random walk for gamma
##                mbar: the empirical average age profile
##               drift: the drift vector of the age profiles, obtained
##                      multiplying theta by beta
##
##                   
##
## WRITTEN BY: Federico Girosi 
##             fgirosi@latte.harvard.edu
##             CBRSS, Harvard University
## 
## Last modified: 09/8/2003
## 
## ######################################################################################
## ######################################################################################


lc.model <-  function(mat,horizon){
   years <- as.numeric(rownames(mat));
###
### for the computation of the  parameters we use the log mortality matrix 
### with the NAs removed
###
  m <- na.omit(mat);
  n <- nrow(m);
###
### mb is the empirical average age profile
### matbar is the centered log mortality matrix , and mbar is the centered
### log mortality matrix with  NA removed
###
  mb <- colMeans(as.data.frame(m),na.rm=TRUE);
  mbar <- scale(m,center = mb,scale=FALSE);
  matbar <-  scale(mat,center = mb,scale=FALSE);
###
### cmat is the correlation matrix of the centered data
###
  cmat <- t(mbar)%*%mbar;
###
### beta is the first eigenvalue of cmat, normalized to have length 1
###
  e <- eigen(cmat,symmetric=TRUE);
  beta.lc <- e$vectors[,1];
  beta.lc <-  beta.lc/sqrt(crossprod(beta.lc));
###
### then gamma at time t is obtained by projecting the centered age profile at time t
### onto the first principal component. notice that gamma.insample has missing
### values where data are missing
###
  gamma.insample <-  matbar%*% beta.lc;
###
### now we can compute the drift parameter from the first and last non-missing
### elements of gamma
###
  gamma.lc <- na.omit(gamma.insample);
  nr <- nrow(gamma.lc);
  tvec <- as.numeric(rownames(gamma.lc));
  theta <- (gamma.lc[nr,1] - gamma.lc[1,1])/(tvec[nr] - tvec[1]);
###
### gamma.insample.fit is the prediction of LC applied to the insample period.
### We are allowing the insample data to begin with an arbitrary long series
### of missing values (useful for drawing some figures)
###
  gamma.insample.fit <- gamma.insample[as.character(tvec[1]),1] + (years-tvec[1])*theta;
  names(gamma.insample.fit) <- names(gamma.insample);
###
### now we fill the missing values of gamma with the values predicted by the model
### 
  isna <- is.na(gamma.insample);
  gamma.insample[isna] <- gamma.insample.fit[isna];
###
### in lc.insample we store the fit of the first principal component to the data
### filling in the missing values with the LC prediction (this is not very interesting)
### lc.insample.hat is the LC prediction for the insample period. there has
### to be a better way to compute this, without a for loop.
###
  lc.insample <- 0*mat;
  lc.insample.hat <- 0*mat;  
  for (i in 1:nrow(mat)){
    lc.insample[i,] <- gamma.insample[i]* beta.lc + mb;
    lc.insample.hat[i,] <- gamma.insample.fit[i]* beta.lc + mb;    
  };
###
### lc.outsample is the log mortality matrix with the LC predictions
### gamma.outsample is the prediction of gamma for the out of sample period
###
  lc.outsample <- matrix(0,nrow=horizon,ncol=ncol(mat));
  rownames(lc.outsample) <- years[length(years)] + 1:horizon;
  colnames(lc.outsample) <- colnames(lc.insample);
  tout <- as.numeric(rownames(lc.outsample));
  gamma.outsample <- 0*lc.outsample[,1];
  for (i in  1:nrow(lc.outsample)){
    gamma.outsample[i] <-  gamma.lc[nr] + (tout[i]-tvec[nr])*theta;
    lc.outsample[i,] <- gamma.outsample[i]*beta.lc + mb;
  };
  
  l <- list(m=mat,m.centered=matbar,y.insample.hat=lc.insample.hat,yhatin=lc.insample,
            yhatout=lc.outsample,beta=beta.lc,
            gamma=gamma.lc,gamma.insample.fit=gamma.insample.fit,gamma.outsample=gamma.outsample,
            theta=theta,mbar=mb,drift=theta*beta.lc)
};




## ######################################################################################
## ######################################################################################
##
## FUNCTION NAME: plot.lc.model
##
## DESCRIPTION: it plots an object of type "LC forecast", the result of lc.model(...)
##
## IMPORTED FUNCTIONS: long.causes
##
## GLOBALS:   none
##
##
## FORMAT: plot.lc.model(lc.result,cause,gender.str,cntry.name,outpath){
##
## INPUT:  lc.result: (list) an object of type "LC forecast", typically the
##                    result of a call to lc.model
##             cause: (string) 4 letters identifiers of a cause of death
##        gender.str: (string or numeric) a string identifying the gender,
##                    usually "m" or "f", to be used in plot legends and graph
##                    filenames. If numeric 2 then it is converted to "m" and
##                    if numeric 3 it is converted to "f" (for compatibility
##                    with other functions)
##        cntry.name: (string) a tring identifying the country, to be used in
##                    plot legends and graph filenames.
##           outpath: (string) pathname of the directory where we want to save
##                    the graphs
##
## 
## OUTPUT: 3 graphs are produced.
##
##         graph 1: Time series of the data and the forecasts, one for each age group,
##                  all in one graph. Each time series is labeled according to the
##                  corresponding age group. Labels are put near the end point of
##                  the time series. The file name begins with "lc.forecasts_"
##                  followed by the cause of death, the gender and the country.
##
##         graph 2: The betas are plotted as a function of age. The file name
##                  begins with "beta_", followed by the cause of death,
##                  the gender and the country.
##
##         graph 3: Age profiles, both data and the forecasts, one for each year,
##                  all in one graph. The graphs are color coded according rainbow
##                  color, from the first year up. A legend on the top left corner
##                  shows the color palette and the years corresponding to the
##                  first and last color. The file name begins with "lc.age.profiles_"
##                  followed by the cause of death, the gender and the country.
##                   
##
## WRITTEN BY: Federico Girosi 
##             fgirosi@latte.harvard.edu
##             CBRSS, Harvard University
## 
## Last modified: 09/10/2003
## 
## ######################################################################################
## ######################################################################################



plot.lc.model <- function(lc.result,cause,gender.str,cntry.name,outpath){
  wide <-  7.5;
  high <-  7.5;
  lname <- long.causes(cause);
  if (is.numeric(gender.str)){
    if (gender.str==2) {gender.str <- "m"
    } else {
      gender.str <- "f"}
  }
  gender.str.title <- paste("(",gender.str,")",sep="");
  titlestr <- paste(lname$disnam1,lname$disnam2,lname$disnam3,gender.str.title,cntry.name)  
  xl <- "Time";
  yl <- "Data and Forecasts";
  y.c <- lc.result$m;
  colo <- rainbow(ncol(y.c));
  colo[3] <- colo[2];
  colo[5] <- colo[4];
  age.vec <- colnames(y.c);

###
### First we plot the forecasts and the data, all in one graph
###
  
  y.plot <- rbind(y.c,lc.result$yhatout);
  y.age.labels <- y.plot[nrow(y.plot),];
  ###
  ### I add 5 years of NAs to I have space for the age labels
  ###
  aux <- matrix(NA,nrow=5,ncol=ncol(y.plot));
  rownames(aux) <- seq(as.numeric(rownames(y.plot)[nrow(y.plot)])+1,length=5);
  y.plot <- rbind(y.plot,aux);
  
  outfile <- paste(outpath,"lc.forecasts_",cause,"_",gender.str,"_",gsub(" ","_",cntry.name),".eps",sep="")
  postscript(outfile,paper="letter",pointsize = 14,width=wide,height=high,horizontal=FALSE);
  par(cex.lab=1.5,cex.main=1.4,bg="tan1");                          
  matplot(rownames(y.plot),y.plot,type="l",lty=1,col=colo,main=titlestr,xlab=xl,ylab=yl);
  last.year <- as.numeric(rownames(y.plot)[nrow(y.plot)]);
  text(x=matrix(c(last.year,last.year-4),nrow=length(age.vec),ncol=1),y=y.age.labels,age.vec,cex=0.7);
  dev.off();

###
### Then we plot the beta's
###

  beta <- lc.result$beta;
  ###
  ### The sign of beta is arbitrary
  ### I adopt the convention that the coefficient beta for the first age group is negative.
  ### 
  ###
  if (beta[1] > 0) beta <- - beta;
  outfile <- paste(outpath,"lc.beta_",cause,"_",gender.str,"_",gsub(" ","_",cntry.name),".eps",sep="")
  postscript(outfile,paper="letter",pointsize = 14,width=wide,height=high,horizontal = FALSE);
  par(cex.lab=1.5,cex.main=1.4);
  xl <- "Age";
  yl <- "Beta (decay rate)";  
  plot(age.vec,beta,type="l",col="red",main=titlestr,xlab=xl,ylab=yl);  
  dev.off();                        


###
### And finally we plot the age profiles
###
  y.plot <- na.omit(y.plot)
  colo <- rainbow(nrow(y.plot));
  outfile <- paste(outpath,"lc.age.profiles_",cause,"_",gender.str,"_",gsub(" ","_",cntry.name),".eps",sep="")
  postscript(outfile,paper="letter",pointsize = 14,width=wide,height=high,horizontal = FALSE);
  par(cex.lab=1.5,cex.main=1.4,bg="tan1");
  xl <- "Age";
  yl <- "Data and Forecasts";
  r <- range(na.omit(y.plot));
  a.vec <- as.numeric(age.vec);
  xp <- seq(from=a.vec[2],to=a.vec[round(length(a.vec)/3)],length=nrow(y.plot));
  yp <- rep(r[2]-0.1*(r[2]-r[1]),nrow(y.plot));
  matplot(age.vec,t(y.plot),type="l",lty=1,col=colo,main=titlestr,xlab=xl,ylab=yl);
  points(xp,yp,col=colo);
  text(a.vec[1],yp[1],rownames(y.plot)[1],cex=0.8);
  text(a.vec[round(length(a.vec)/3)+1],yp[1],rownames(y.plot)[nrow(y.plot)],cex=0.8);
  dev.off();                        
  
  
}
