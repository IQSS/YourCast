
## ************************************************************************
## ************************************************************************
##
## FUNCTION NAME:  plot.sigmas
##
##
## IMPORTED:    long.causes,cntry.names,strpad,sigma.profile.from.ols
##              death.from.logmortality,poisson.sigma.profile,
##              compare.model.based,
##
## GLOBALS: none
##
##
##
## DESCRIPTION: This function is used to compare the standard deviations derived from      
##              equation by equation OLS and those implied by the Poisson                  
##              specification. It takes as input ols.result, that is an object which       
##              is the output of the function ols, and several control parameters.
##              It produces 3 kinds of plots, which are either shown on the screen
##              or printed to files.
##              Plots of type 1: one plot for each country. We plot the age profiles
##                               of the standard deviations from OLS (continous line)
##                               and from the Poisson specification (dotted line).
##                               If show.adjusted=TRUE the age profile of the Poisson
##                               specification is shifted upward to fit the OLS profile.
##           
##               Plot of type 2: a histogram of the distribution of the insample relative
##                               error over the countries. The insample relative error
##                               for each country is defined as the age-average of the absolute 
##                               difference between the OLS and the adjusted Poisson age 
##                               profile, divided by the age-average of the OLS age profile.
##               Plot of type 3: for every age group and country we compute the insample
##                               relative error as defined above. Then we average this over
##                               countries and plot the result as a function of age.
##
##
## FORMAT: aux <- plot.sigmas(ols.result,show.adjusted=FALSE,graphdir=NA,compact=TRUE)
##         
## INPUT:       ols.result: (list) an ols object, the output of the function ols.
##
##          show.adjusted: (logical) if TRUE the  age profile of the Poisson
##                         specification is shifted upward to fit the OLS profile.
##                         This is the smoothed estimated as defined in the manual, with
##                         a smoothing parameter theta equal to infinity.
##
##               graphdir: (string or NA) if NA plots are showed on the user screen. The user can hit
##                         return to show the next graph, or any other key to skip
##                         the next graph.
##                         if string it must be the path of a directory, where the graphs
##                         will be printed to files, rather than being shown on screen.
##                         Plots of type 1 are saved with file name given by
##                         sigma_disease/gender_countryname (eg: sigma_allc_3_Italy.ps)
##                         The plot of type 2 is saved with file name given by
##                         sigma_hist_disease/gender (e.g. sigma_hist_allcg2.ps)
##                         The plot of type 3 is saved  with file name given by
##                         sigma_hist_disease/gender (e.g. sigma_ages_allcg2.ps)
##
##                compact: (logical) if TRUE we collapse all the printed files in one,
##                          multipage ps file, with 8 ps files in each page and delete
##                          the single files.
##
## WRITTEN BY: Federico Girosi
##             fgirosi@latte.harvard.edu
##             CBRSS, Harvard University
## 
## Last modified: 08/11/2003
## 
## ************************************************************************


plot.sigmas <- function(ols.result,show.adjusted=FALSE,graphdir=NA,compact=TRUE, ebase=env.base){
  ebase <- get("env.base", envir=parent.frame())
  env.base <- ebase
  ewho <- get("env.who", envir=ebase)
  popul <- get("whopopul", envir=ewho)
  ages <- get("age.vec", envir=ewho)
  disease <- get("whodisease", envir=ewho)
  gender <- get("whogender", envir=ewho)
                
  std <- ols.result$std;       ### standard deviations from OLS (cross-sectional list)
  yhatin <- ols.result$yhatin; ### insample prediction from ols (cross-sectional time series)
###  popul <- ols.result$popul;   ### insample population (cross-sectional time series)
###  ages <- ols.result$age.vec;  
###  disease <- ols.result$disease;
###  gender <- ols.result$gender;
  
###
### First we define quantities used in plots and file naming
###

  disease.long <- long.causes(disease)$name;
  if (gender == 2) {gender.str <- "(males)"
                  } else {
                    gender.str <- "(females)"
                  }
  c.names <- cntry.names(ebase=ebase);

  c.names.padded <- strpad(c.names); ### Padded to have all same length

  ###
  ### We store in sigma.ols.profile the age profiles of standard deviations
  ### obtained by equation by equation ols. The elements of the list are called
  ### $\sigma_c^\text{OLS}$ in the manual. This is a cross-sectional list
  ### of age profiles, indexed by countries
  ###
  
  sigma.ols.profile <- sigma.profile.from.ols(std,bycntry=TRUE)

  ###
  ### We estimate the expected number of deaths using the predicted values of
  ### equation by equation OLS (stored in yhatin). This is called $\lambda_{cat}$
  ### in the manual and it is a cross-sectional time series.
  ###
  
  death.cs <- death.from.logmortality(yhatin,list.by.cntry(popul));

  ###
  ### From death.cs we compute the age profiles of standard deviations
  ### implied by the Poisson specification and store it in sigma.poisson.profile
  ### This is called $s_{ca}$ in the code, and it is a cross-sectional list
  ### of age profiles, indexed by countries.
  ###
  
  sigma.poisson.profile <- poisson.sigma.profile(list.by.csid(death.cs),bycntry=TRUE);

  ###
  ### Since we want to use lapply over the countries in the lists defined above
  ### we need a way to access a number of quantities on a country by country
  ### basis. This is achieved by performing lapply on ind, which is a list
  ### of the names of the cross-sections (country codes), with the names themselves
  ### as names of the elements. For example ind[["2450"]] = "2450".
  ###
  
  ind <- as.list(names(sigma.poisson.profile));
  names(ind) <-  names(sigma.poisson.profile);
  
  if (dev.cur() != 1) dev.off() ### this is just to get a new window

  ###
  ### This call produces the plots of the sigmas from OLS and the sigmas implied by
  ### the Poisson specification, one for each country. It also returns a list
  ### indexed by countries. Each elements of the list is a matrix with A rows
  ### and 4 columns.
  ### Column 1: age groups
  ### Column 2: sigma profile from OLS
  ### Column 3: sigma profile adjusted
  ### Column 4: sigma profile from Poisson specification
  ### the adjustment parameter is any element of (column 3 - column 4)

  
  mat <- lapply(ind,FUN=compare.model.based,sigma.ols.profile,sigma.poisson.profile,
                ages,c.names,disease,gender,graphdir,show.adjusted);

  ###
  ### Now we wish to produce an histogram of the adjustment parameter
  ### However the adjustment parameters are not comparable across countries,
  ### and we need some relative measure. We choose the adjustment parameter
  ### divided by the average sigma from OLS. In formulas this is
  ### $\alpha_c \over \bar \sigma_c^\text{OLS}$. 
  ###

  alpha.perc <-  lapply(mat,FUN=function(x){adj <- abs(x[,3] - x[,4])[1];
                                            bar.sigma.ols <- mean(x[,2]);
                                            return(adj/bar.sigma.ols)})

  ###
  ### Now we wish to produce an histogram of the insample error, that is
  ### the difference between the OLS estimates (data) and the smoothed estimates
  ### To this end, first we store in err the absolute value of the difference
  ### between OLS estimates and smoothed estimates, divided by the
  ### average OLS estimate. err is a cross-sectional list of age profiles
  ### of relative error measures.
  ### The age profiles have names given by the age groups
  ###
  
  err <- lapply(mat,FUN=function(x){z <- abs(x[,2]-x[,3])/mean(x[,2]);
                                    names(z) <- conv.char(x[,1]);
                                    return(z)});

  ###
  ### Now we average the age profiles of relative error over age groups,
  ### obtaining err.by.cntry, a cross-sectional list, indexed by country,
  ### of average, relative error measures.
  ###
  
  err.by.cntry <- lapply(err,FUN=function(x){m <-  mean(x)});

  ###
  ### If graphdir is missing we plot on the screen, otherwise we save to a file
  ### named sigma_hist_disease_gender.ps (Example: sigma_hist_allcg2.ps)
  ### I hate to have to write twice the command hist, in the if and in the else,
  ### but this is so much simpler
  ###
  
  if (is.na(graphdir)){
    a <- readline("Press return to plot the distribution over countries of relative \n  insample errors between the OLS estimates and the smoothed estimates: ");
    if (a=="") {      
      hist(unlist(err.by.cntry),col="green",xlab="Relative Insample Error",
           main=paste(disease.long,gender.str,"\n","Distribution of Relative Insample Errors Over Countries"))
    }
  } else {
    outfile <- paste(graphdir,"sigma_hist_",disease,"g",gender,".ps",sep="");
    postscript(outfile,paper="letter",pointsize = 14,width=7,height=7, horizontal = TRUE, onefile = TRUE);
    hist(unlist(err.by.cntry),col="green",xlab="Relative Insample Error",
         main=paste(disease.long,gender.str,"\n","Distribution of Relative Insample Errors Over Countries"));
    dev.off();
  }


  ###
  ### In addition to the histogram we also print a list of the relative insample
  ### errors by country. ind plays the same role as above
  ### We pass the country names padded with blanks to have all the same length
  ### so that they print nicely.
  ###
  
  ind <- as.list(names(err.by.cntry));
  names(ind) <-  names(err.by.cntry);

  aux <- lapply(ind,FUN=function(n,err.by.cntry,c.names.padded){
    cntry <- c.names.padded[as.numeric(n)];
    e <- err.by.cntry[[n]];
    print(paste(cntry,round(e,3)))
  },err.by.cntry,c.names.padded);


  ###
  ### Now we want to plot the relative insample error by age groups, averaging over countries
  ### This information is contained in the list err, which we need to convert to a data frame.
  ### as.data.frame(err) has age groups on the rows and countries on the columns. We needed
  ### the opposite so we can use the function mean, which explains the following contorsion:
  ###

  err <- colMeans(as.data.frame(t(as.data.frame(err))));

  ###
  ### Now err is simply an age profile, with names given by age groups.
  ### If graphdir is missing we plot on the screen, otherwise we save to a file
  ### named sigma_ages_disease_gender.ps (Example: sigma_ages_allcg2.ps)
  ### I hate to have to write twice the commands plot and title, in the if and in the else,
  ### but this is so much simpler
  ###
  
  if (is.na(graphdir)){    
    a <- readline("Press return to plot relative insample errors by age (averaged over countries)");
    if (a=="") {      
      plot(names(err),err,type="l",col="red",xlab="Age",ylab="Relative Insample Error")
      title(main=paste(disease.long,gender.str,"\n","Relative Insample Errors, Country Averaged"));
    }
  } else {
    outfile2 <- paste(graphdir,"sigma_ages_",disease,"g",gender,".ps",sep="");
    postscript(outfile2,paper="letter",pointsize = 14,width=7,height=7, horizontal = TRUE, onefile = TRUE)
    plot(names(err),err,type="l",col="red",xlab="Age",ylab="Relative Insample Error");
    title(main=paste(disease.long,gender.str,"\n","Relative Insample Errors, Country Averaged"));
    dev.off();
  }

  ###
  ### If we saved to ps files and the option compact is TRUE we condense all the printed files
  ### in one, multipage ps file, with 8 ps files in each page and delete the single files
  ###
  
  if (!is.na(graphdir) && compact==TRUE){
    outps <- paste(graphdir,"sigma_",disease,"g",gender,".ps",sep="");
    mpage.str <- paste("mpage -c8",outfile,outfile2,paste(graphdir,"sigma_",disease,"g",gender,"_* >",sep=""),outps);
    system(mpage.str);
    system(paste("/bin/rm",outfile,outfile2,paste(graphdir,"sigma_",disease,"g",gender,"_*",sep="")))
  }
  
  return(invisible(mat))
}


## ************************************************************************
## ************************************************************************
##
## FUNCTION NAME:  estimate.standard.deviations
##
##
## IMPORTED:    sigma.profile.from.ols,death.from.logmortality,poisson.sigma.profile,
##
## GLOBALS: none
##
## DESCRIPTION: this is the functions which estimates the standard deviations
###             as explained in the manual, using equation by equation OLS
##
## FORMAT: std  <- estimate.standard.deviations(ols.result)
##         
## INPUT:   ols.result: (list) an ols object, the output of the function ols.
##
## OUTPUT:         std: (list) a cross-sectional list, indexed by countries and ages,
##                      with the values of the estimated standard deviations
##
## WRITTEN BY: Federico Girosi
##             fgirosi@latte.harvard.edu
##             CBRSS, Harvard University
## 
## Last modified: 08/11/2003
## 
## ************************************************************************

estimate.standard.deviations <- function(ols.result, ctr = NULL,ebase=env.base){
  ebase <- get("env.base", envir=parent.frame())
  env.base <- ebase
  ewho <- get("env.who", envir=ebase)
  popul <- get("whopopul", envir=ewho)
  if(length(ctr) > 0)
    popul <- popul[grep(ctr,names(popul))]
  insampy <- get("whoinsampy", envir=ewho)
  if(length(ctr) > 0)
    insampy <- insampy[grep(ctr,names(insampy))]               
  age.vec <- get("age.vec", envir=ewho)
  who.cntry.digits <- get("who.cntry.digits", envir=ewho)
  who.age.digits <- get("who.age.digits", envir=ewho)
  who.year.digits  <- get("who.year.digits", envir=ewho)
  who.digit.first  <- get("who.digit.first", envir=ewho)
  digit.year.begin <- who.digit.first + who.cntry.digits + who.age.digits + 1;  
  ages <- age.vec 
  std <- ols.result$std;       ### standard deviations from OLS (cross-sectional list)
  yhatin <- ols.result$yhatin; ### insample prediction from ols (cross-sectional time series)
###  popul <- ols.result$popul;   ### insample population (cross-sectional time series)  
###  insampy <- ols.result$insampy;
###  ages <- ols.result$age.vec;

  ###
  ### We store in sigma.ols.profile the age profiles of standard deviations
  ### obtained by equation by equation ols. The elements of the list are called
  ### $\sigma_c^\text{OLS}$ in the manual. This is a cross-sectional list
  ### of age profiles, indexed by countries
  ###
  sigma.ols.profile <- sigma.profile.from.ols(std,bycntry=TRUE)
 
  ###
  ### We estimate the expected number of deaths using the predicted values of
  ### equation by equation OLS (stored in yhatin). This is called $\lambda_{cat}$
  ### in the manual and it is a cross-sectional time series.
  ###
  
  death.cs <- death.from.logmortality(yhatin,popul);
 
  ###
  ### From death.cs we compute the age profiles of standard deviations
  ### implied by the Poisson specification and store it in sigma.poisson.profile
  ### This is called $s_{ca}$ in the code, and it is a cross-sectional list
  ### of age profiles, indexed by countries.
  ###
  
  sigma.poisson.profile <- poisson.sigma.profile(death.cs,bycntry=TRUE);
 
  ###
  ### Since we want to use lapply over the countries in the lists defined above
  ### we need a way to access a number of quantities on a country by country
  ### basis. This is achieved by performing lapply on ind, which is a list
  ### of the names of the cross-sections (country codes), with the names themselves
  ### as names of the elements. For example ind[["2450"]] = "2450".
  ###
  
  ind <- as.list(names(sigma.poisson.profile));
  names(ind) <-  names(sigma.poisson.profile);
 
  func <- function(n,sigma.ols.profile,sigma.poisson.profile){
    sigma.ols.profile.c <- sigma.ols.profile[[n]];
    sigma.poisson.profile.c <- sigma.poisson.profile[[n]];
    adj <- mean(sigma.ols.profile.c) - mean(sigma.poisson.profile.c);
    adj <- max(c(adj,0));
    estimated.sigma <- sigma.poisson.profile.c + adj;
    return(estimated.sigma);
  }

  ###
  ### I store in mat, a cross-sectional list indexed by countries, the age profiles
  ### of the estimated standard deviations
  ###
 
  mat <- lapply(ind,FUN=func,sigma.ols.profile,sigma.poisson.profile);
  age.char <- conv.char(age.vec);  
  mat <-  lapply(mat,FUN=function(x,age.char){names(x) <- age.char;return(x)},age.char)
 
  ###
  ### I have to store the value in mat as a cross-sectional list indexed
  ### by countreis and age groups. I use the fact unlist creates names
  ### by joining string with a dot (.). There has to be a better way of doing this.
  ###
  
  sigmas <- unlist(mat);
  names(sigmas) <- sapply(names(sigmas),FUN=function(x){gsub("[.]","",x)});

  ###
  ### I need to return the standard deviations as a cross-sectional time series
  ### constant over time, and I store it in res
  ###
  
  ind <- as.list(names(sigmas));
  names(ind) <-  names(sigmas);
  res <- lapply(ind,FUN=function(n,sigmas,insampy){
    s <- sigmas[[n]];
    s <- rep(s,length(insampy[[n]]))
    names(s) <- substring(rownames(insampy[[n]]),digit.year.begin)
    s <- as.matrix(s)},sigmas,insampy)
  return(res)
}



## ************************************************************************
## ************************************************************************
##
## FUNCTION NAME:  poisson.sigma.profile 
##
##
## IMPORTED:    make.average.age.profile
##
## GLOBALS: none
##
## DESCRIPTION: given a cross-sectional time series of expected number of deaths
##              it returns the corresponding age profiles of the standard deviations
##              according to the Poisson specification.
##
## FORMAT: std <- poisson.sigma.profile(death.cs,bycntry=FALSE)
##         
## INPUT:   death.cs: (list) a cross-sectional time series of expected number of deaths
##
##           bycntry: if TRUE the standard deviations are averaged over time only, and
##                    and an age profile is returned for every country. Otherwise the
##                    standard deviations are averaged over time and countries, and a single
##                    age profile is returned
##
## OUTPUT: std: vector (age profile) or list (of age profiles),
##              depending on the value of bycntry
##
## WRITTEN BY: Federico Girosi
##             fgirosi@latte.harvard.edu
##             CBRSS, Harvard University
## 
## Last modified: 08/11/2003
## 
## ************************************************************************

poisson.sigma.profile <- function(death.cs,bycntry=FALSE, ebase=env.base){
  ebase <- get("env.base", envir=parent.frame())
  env.base <- ebase
  ewho <- get("env.who", envir=ebase)
  sigma.cs <- lapply(death.cs,FUN=function(x){sqrt(log(1+1/x))});
  sigma.average.profile <- make.average.age.profile(sigma.cs,bycntry=bycntry);
  return(sigma.average.profile)
}



######################################################################
######################################################################
######################################################################

sigma.profile.from.ols <- function(std,bycntry=FALSE, ebase=parent.frame()){
  ebase <- get("env.base", envir=parent.frame())
    ewho <- get("env.who", envir=ebase)
  who.digit.first  <- get("who.digit.first", envir=ewho)
  who.cntry.digits <- get("who.cntry.digits", envir=ewho)
  who.age.digits   <- get("who.age.digits", envir=ewho)
  who.year.digits  <- get("who.year.digits", envir=ewho)  
    digit.cntry.begin <- who.digit.first + 1 
    digit.cntry.end   <- who.digit.first + who.cntry.digits
    digit.age.begin   <- digit.cntry.end + 1
    digit.age.end     <- digit.cntry.end + who.age.digits
    
    
  sigma <- unlist(std)
  n <- as.numeric(names(std))
  age.f <- factor(digitpull(n,digit.age.begin,digit.age.end))
  cntry.f <- factor(digitpull(n,digit.cntry.begin,digit.cntry.end))
  if (bycntry==FALSE){
    mat <- tapply(sigma,list(cntry.f,age.f),FUN=function(x){x})
    sigma.profile <- colSums(mat)/nrow(mat);
  } else {
    sigma.profile <- tapply(sigma,list(cntry.f),FUN=function(x){x})
  }
  return(sigma.profile)
}



######################################################################
######################################################################
######################################################################


compare.model.based <- function(n,sigma.profile,sigma.poisson.profile,ages,c.names,
                                disease,gender,graphdir,show.adjusted){
  sigma.profile.c <- sigma.profile[[n]];
  sigma.poisson.profile.c <- sigma.poisson.profile[[n]];
  disease.long <- long.causes(disease)$name;
  cntry <- c.names[as.numeric(n)];
  if (gender == 2) {gender.str <- "(males)"
                  } else {
                    gender.str <- "(females)"
                  }      
  adj <- mean(sigma.profile.c) - mean(sigma.poisson.profile.c);
  adj <- max(c(adj,0)); ### if adj is negative we do not adjust
  sigma.poisson.profile.c.adj <- sigma.poisson.profile.c + adj
  if (show.adjusted == TRUE) {
    sigmas <- cbind(sigma.profile.c,sigma.poisson.profile.c.adj);
  } else {
    sigmas <- cbind(sigma.profile.c,sigma.poisson.profile.c);
  }
  mat <- cbind(ages,sigma.profile.c,sigma.poisson.profile.c.adj,sigma.poisson.profile.c);
  if (is.na(graphdir)){
    dev.set(dev.prev());    
    a <- readline("Press return for plot ");
    if (a=="") {      
      matplot(ages,sigmas,type="l",xlab="Age",ylab="Standard Deviation of Log Mortality");
      title(main=paste(disease.long,gender.str,"\n",cntry));}
  } else {
    outfile <- paste(graphdir,"sigma_",disease,"g",gender,"_",gsub(" ","_",cntry),".ps",sep="");
    postscript(outfile,paper="letter",pointsize = 14,width=7,height=7, horizontal = TRUE, onefile = TRUE);
    matplot(ages,sigmas,type="l",xlab="Age",ylab="Standard Deviation of Log Mortality");
    title(main=paste(disease.long,gender.str,"\n",cntry));
    dev.off();
  }
  return(mat);
}



######################################################################
######################################################################
######################################################################

   compare.model.based.alternative <- function(n,sigma.profile,sigma.poisson.profile,ages,c.names,
                                               disease,gender,graphdir,show.adjusted){
      sigma.profile.c <- sigma.profile[[n]];
      sigma.poisson.profile.c <- sigma.poisson.profile[[n]];
      disease.long <- long.causes(disease)$name;
      if (gender == 2) {gender.str <- "(males)"
                      } else {
                        gender.str <- "(females)"
                      }      
      if (show.adjusted==TRUE){
        adj <- mean(sigma.profile.c^2) - mean(sigma.poisson.profile.c^2);
        adj <- max(c(adj,0));
        if (adj > 1/3) cat(c.names[as.numeric(n)], "\n");
        sigma.poisson.profile.c <- sqrt(sigma.poisson.profile.c^2 + adj)
      }
      cntry <- c.names[as.numeric(n)];
      sigmas <- cbind(sigma.profile.c,sigma.poisson.profile.c);      
      mat <- cbind(ages,sigmas);
      if (is.na(graphdir)){
        dev.set(dev.prev());    
        a <- readline("Press return for plot ");
        if (a=="") {      
          matplot(ages,sigmas,type="l",xlab="Age",ylab="Standard Deviation of Log Mortality");
          title(main=paste(disease.long,gender.str,"\n",cntry));}
      } else {
        outfile <- paste(graphdir,"sigma_",disease,"g",gender,"_",gsub(" ","_",cntry),".ps",sep="");
        postscript(outfile,paper="letter",pointsize = 14,width=7,height=7, horizontal = TRUE, onefile = TRUE);
        matplot(ages,sigmas,type="l",xlab="Age",ylab="Standard Deviation of Log Mortality");
        title(main=paste(disease.long,gender.str,"\n",cntry));
        dev.off();
      }
      return(mat);
    }
 
##**************************************************************************
##**************************************************************************

### for display in the graphs
### Elena Villalon: Because with the revised version of yourcast the data part of dataobj contains
### columns such as allc2, allc3, cvds2, cvds3 for male (2) and females we need to change
### the input x so that it recognizes the disease
###
long.causes <- function(x,dispgender=T){
  nc <- nchar(x)
  xs <- strsplit(x, NULL)[[1]]
  ix2 <- grep("2", xs)
  ix3 <- grep("3", xs)
  ix <- c(ix2, ix3)
  str <- NULL
  if(length(ix2) > 0)
    str <- "m"
  if(length(ix3) > 0)
    str <- "f"
  
  if(length(ix) > 0){
    xs <- xs[-ix]
    x <- paste(xs, collapse="")
  }
  disnam1 <- switch(paste(x),
                    allc =    "All",
                    malr =    "Malaria", 
                    aids =    "AIDS", 
                    tubr =    "Tuberculosis", 
                    otin =    "Other",
                    lung =    "Lung", 
                    molp =    "Cancer of",
                    livc =    "Liver", 
                    stom =    "Stomach",
                    brst =    "Breast",
                    cerv =    "Cervix",
                    omal =    "Other",
                    rspi =    "Respiratory",
                    rspc =    "Respiratory",
                    cvds =    "Cardiovascular",
                    dgst =    "Digestive",
                    matc =    "Maternal",
                    pern =    "Perinatal",
                    allo =    "All",
                    trns =    "Transportation",
                    unin =    "Other",
                    suic =    "Suicide", 
                    homi =    "Homicide",   
                    ward =    "War");

  disnam2 <- switch(paste(x),
                    allc =    "Causes", 
                    malr =    "", 
                    aids =    "", 
                    tubr =    "", 
                    otin =    "Infectious",
                    lung =    "Cancer", 
                    molp =    "Mouth and",
                    livc =    "Cancer", 
                    stom =    "Cancer", 
                    brst =    "Cancer", 
                    cerv =    "Cancer", 
                    omal =    "Malignant", 
                    rspi =    "Disease,", 
                    rspc =    "Disease,", 
                    cvds =    "Disease", 
                    dgst =    "Disease", 
                    matc =    "Conditions",
                    pern =    "Conditions", 
                    allo =    "Other", 
                    trns =    "Accidents", 
                    unin =    "Unintentional", 
                    suic =    "", 
                    homi =    "",   
                    ward =    "");

  disnam3 <- switch(paste(x),
                    allc =    "", 
                    malr =    "", 
                    aids =    "", 
                    tubr =    "", 
                    otin =    "Diseases", 
                    lung =    "", 
                    molp =    "Esophagus", 
                    livc =    "", 
                    stom =    "", 
                    brst =    "",
                    cerv =    "",
                    omal =    "Neoplasms", 
                    rspi =    "Infectious", 
                    rspc =    "Chronic", 
                    cvds =    "",
                    dgst =    "",
                    matc =    "",
                    pern =    "",
                    allo =    "Diseases", 
                    trns =    "",
                    unin =    "Injuries", 
                    suic =    "",
                    homi =    "",
                    ward =    "");
  
  name <- paste(disnam1,disnam2,disnam3);
  if(length(str) >0 && dispgender)
    name <- paste(name,"(",str,")", sep="")
  
  return(list(name,disnam1,disnam2,disnam3));
}
