### DESCRIPTION: It applies a formula to  dataobj$data. 
###            
### USES: dthpop(formula) for left-hand-side of formula
###       convert.timeseries(datamat= dataobj$data)
###       to get the ts values of elements of datamat
###       lag.datamat(...), with 4 arguments
###       that aplies the formula to the elements of datamat=dataobj$data
###       timecolumn that completes the "time" column of elements datamat
###       if an expansion in time has been achieved with the lag.
###       data.after.first.obvs to eliminate rows of missing depvar
###       before first actual observations of death
###
### INPUT: formula to apply and dataobj$data list of csid elements
###        
### OUTPUT: transform datamat= dataobj$data with lag in time if any
###         Example, formula can be
###         ff <- lag(log(allc.txt), 3) ~ lag(tobacco.txt^4, -7) +
###               fat.txt + lag(gdp.txt, -2) + time
###
### AUTHOR: Elena Villalon
###         IQS, Harvard Univ
##         evillalon@iq.harvard.edu
###         1737 Cambridge St
###         Cambridge MA
###
### Modified Dec 5, 2005
###
########################################################################
lagdataobj <- function(formula, dataobj)
{
  fc <- as.character(formula)
  ls <- fc[2]
  depvar <- findregex(ls)
  depvar <- split.ops(depvar)
  if(length(depvar) > 1){
    dth <- depvar[1] ## numerator
    pop <- depvar[2] ## denominator
  }else{
    dth <- depvar[1] 
    pop <- NA
  }
  dthlst <- list(dth = dth, pop=pop)
  pop <- na.omit(pop)
  dthcov <- grep(dth, as.character(formula)[3]) ## any depvar as covariate: numerator?
  popcov <-  NULL
  if(length(pop) > 0) ## denominator of depvar as covariate ? 
    popcov <- grep(pop, as.character(formula)[3])
  vec <- NULL
  ln  <- 0
  
  if(length(dthcov) > 0){
    vec <- c(vec, numerator=dth)
    ln <- 1}
  
  if(length(popcov) > 0){
    vec <- c(vec, denominator=pop)
    ln <- ifelse(ln > 0, 2, 1)}
  
  vec <- unique.default(vec)
 
  datamat <- dataobj$data
  if(length(vec) > 0){  ## denominator or denominator aso depvar as covariate ?  
    nmvec <- names(vec)
    for(n in 1:ln){
      dcomp <- vec[n]     
      nm <- nmvec[n]
      addlst  <- catdepvar(formula, datamat,dcomp)
      datamat    <- addlst$datamat   
      formula <- addlst$ff
    }
 }
  
  dataobj$data <- datamat
  ix <- grep("lag", formula)
  data <- datamat
  nm <- names(data)
 
  yearlst <- lapply(data,FUN="first.nonay", dth, pop)

  firstyear <- lapply(yearlst,function(lst) lst$yearobs)
  foy <- lapply(yearlst, function(lst) lst$ind)

  firstyear <- lapply(firstyear,as.numeric)
 
  datamatseries <- convert.timeseries(datamat= dataobj$data);
 
  covnames <- lapply(dataobj$data, colnames)
     
  datamatseries <- lag.datamat(ff=formula, datamatseries, dthlst, covnames, foy)
   
  itime <- grep("time", colnames(datamatseries[[1]]))
   
  if(length(itime) >0)
    datamatseries <- lapply(datamatseries, FUN="timecolumn")
     
  datamat <- lapply(datamatseries, as.matrix)

  datamat <- data.after.first.obv(datamat, formula, noadjust=F)
    
  names(datamat) <- nm
  ind <- as.list(1:length(datamat))
  names(ind) <- nm
  datamat <- lapply(ind,FUN="rowyears", datamat, firstyear)
  lst <- list(datamat = datamat, ff=formula)

  return(lst)
}


### DESCRIPTION it converts matrices into time-series objects.
### INPUT : a list, datamat, of cross-sectional time series and
###        
### OUTPUT: Every element of datamat is converted into an object ts
###
###         Elena Villalon
###         IQS, Harvard Univ
###         evillalon@iq.harvard.edu
###         1737 Cambridge St
###         Cambridge MA
### November 2005
#######################################################

  convert.timeseries <- function(datamat)
  {
    
    datamat <- lapply(datamat, function(mat){
      if (class(mat) == "mts" || class(mat) == "ts")
        return(mat)
      
      rec    <- nrow(mat)
      serial <-  as.numeric(rownames(mat))
      covnames <- colnames(mat)
      freq  <- serial[2] -serial[1]
      noyr  <- serial[length(serial)] - serial[1] + 1
      covyr <- rec - noyr
      st <- ifelse(covyr >= 0, (serial[1] - covyr), serial[1]) 
      tmat <- ts(mat, start=st, frequency=freq, names=covnames)
      lst <- list(tsmat = tmat, covs=covnames)
      return(lst$tsmat)})
    
   return(datamat)
  }
### DESCRIPTION: It actually applies the formula to  datamat,
###              cross-sectional-time-series matrices
###            
### USES: lag.depvar, lag.cov which applies the lag to elements
###       of datamat for depvar; 
###       left-hand-side of formula, and/or covaraites
###       right-hand-sides of formula. find.lag to get years of lag 
###
### INPUT: formula, datamat,  dlst with depvar string values
###        cvname with names for covariates for every elemnt of datamat
###        
### OUTPUT: datamat modified with years lag. 
###         transform datamat with lag in time if any
###         Example, formula can be
###         ff <- lag(log(allc.txt), 3) ~ lag(tobacco.txt^4, -7) + fat.txt + lag(gdp.txt, -2) + time
###         
###
### AUTHOR: Elena Villalon
###         IQS, Harvard Univ
###         evillalon@iq.harvard.edu
###         1737 Cambridge St
###         Cambridge MA
###
### Modified Dec 5, 2005
###
########################################################################
lag.datamat <- function(ff, datamat, dlst, cvnm,fix, verbose=T)
{
###  cat("Inside lag.datamat...","\n")
 verb <- try(get("verbose", envir=get("env.base",envir=parent.frame())),silent=T)
    if(class(verb)!="try-error")
      verbose <- verb
  lst <- findlag(ff)
  lg  <- lst$lag ### lagging number of years & 
  ind <- lst$ind ### if 2 coming from depvar; if 3 coming from covariates
### print(lst)
  dth <- dlst$dth
  pop <- dlst$pop
  
  if(length(pop) > 0)
    pop <- na.omit(pop)
  
  if(length(ind) <= 0)
    return(datamat);
   
  ixmat <- as.list(1:length(datamat))
  messout(paste("Number of csid is  ", length(ixmat), sep=""),verbose)
 
  datamat <- lapply(ixmat, function(mm){
       
    mat <- datamat[[mm]]
    cnames <- cvnm[[mm]]
    cnames <- unlist(sapply(cnames, FUN="trim.blanks"))
    dth0 <- paste("^",dth,"$", sep="")
    indth  <- grep(dth0, cnames)
    indpop <- NULL

    if(length(pop) > 0){
      pop0 <- paste("^", pop, "$", sep="")
      indpop <- grep(pop0, cnames)
    }
    
    indth <- c(indth, indpop)
 
###    print(indth)
    ixl   <- (ind == 0)
    ixr   <- (ind != 0)
 ###   print(ind)
    if(all(ind == 0 )){ ##only depvar is lag
    
      mat <- lag.depvar(mat,indth, lg[1], dth, pop)  
      return(mat)
    }
      
    if(all(ind != 0)) { ##only lag the covariates or rhs of formula
      for(n in 1:length(lg))##l
        {
          nm <- names(lg)[n]        
          mat <- lag.cov(lg[n], nm, mat, fix[[mm]])
        }
      
      return(mat)
    }
    onlyonce <- F

    for(n in 1:length(lg))##lag apply to depvar and covariates
      {
    
        if(n <= 1 && !onlyonce){
            
          mat <- lag.depvar(mat,indth, lg[1], dth, pop)
          onlyonce  <- T
        }else{
         
          nm <- names(lg)[n]
         
          mat <- lag.cov(lg[n], nm, mat, fix[[mm]])
       
        }
      }
    
    return(mat)})
  return(datamat)
        
}
   

### DESCRIPTION takes a matrix and finds the number of years
###             and add the column "time" to the matrix with
###             years of observations (or insample) and outsample years
###
### Elena Villalon
### December 2005
  
 timecolumn <- function(mat)
   {
     itime <- grep("time", colnames(mat))
     if(length(itime) <= 0)
       return(mat)
     cnm <- colnames(mat)
     indtm <- grep("time", cnm)
     nr <- nrow(mat)
     ncl <- ncol(mat)
     vst <- start(mat)[1]
     ved <- end(mat)[1]
###     print(vst)
###     print(ved)
     nrw <- ved[1] - vst[1] + 1 ##number of rows
     timeseries <- vst:ved
     if(itime <= 1){
       mat <- cbind(time=timeseries,mat[,2:ncl])
       colnames(mat) <- c("time",cnm[-1]) 
       return(mat)
     }
     if(itime >= ncl){
       mat <- cbind(mat[,1:(itime-1)], time=timeseries)
       colnames(mat) <- c(cnm[1:(ncl-1)], "time")
       return(mat)
     }
     mat <- cbind(mat[,1:(itime-1)], time=timeseries, mat[,(itime+1):ncl])
     colnames(mat) <- c(cnm[1:(itime-1)], "time", cnm[(itime+1):ncl])
     return(mat)
   }

### DESCRIPTION lag the dependent variable, numerator and denominator
###             Instead of applying the lag to the left-hand-side of
###             the formula, it applies it to the covariates or right-hand-side
###             the years to lag becomes -lag, when lagging covariates.
###             All covariates in the rhs of equation are lagged the
###             same numbers of years, -lg. 
###
###
### Elena Villalon
### evillalon@iq.harvard.edu
##################################################################

 lag.depvar <- function(mat,indth, lg1, dth, pop){
  
   if(length(indth) <= 0)
     stop("Inconsistency in lag.datamat")
   
   indth <- unique.default(indth)
   nm <- colnames(mat)
   ip <-  NULL
   dth <- paste("^", dth, sep="")
   pop <- paste("^", pop, sep="")
   if(length(pop) > 0)
     ip <- grep(pop, nm)
   nmdth <- nm[indth]
   nm <- nm[-indth]
   if(length(pop) > 0)
     indth <- unique.default(c(indth, ip))
 
   mat <- cbind(mat[,indth], lag(mat[,-indth], -lg1))
   colnames(mat) <- c(nmdth, nm)
   return(mat)
 }

### DESCRIPTION lag the covariates, with number of years = lag,
###             covariate name = cov for matrix element mat, 
###             Instead of applying the lag to the left-hand-side of
###             the formula, it applies it to the covariates or right-hand-side
###             the years lag becomes -lag, when lagging covariates.
###             It lags only one covariate at a time, and it calls for several
###             covariates the lag or number of years may be different.
###             The covariates are lagged independently of each other. 
###
### Elena 3Villalon
### evillalon@iq.harvard.edu
##################################################################
lag.cov <- function(lag, cov, mat, fix, verbose=T){   
 verb <- try(get("verbose", envir=get("env.base",envir=parent.frame())),silent=T)
    if(class(verb)!="try-error")
      verbose <- verb
  nr  <- nrow(mat)
  ncl <- ncol(mat)
  nm  <- colnames(mat)
###   print(cov)
###  cov <- paste("^", cov, "$", sep="")
  
  ind <- grep(cov, nm)
###  print(ind)
  if(length(ind) <= 0)
    return(mat)

  covf <-first.nonay(mat, cov, pop=NULL, fromcov=TRUE)$ind 
  permit <- covf - fix
  
  if (permit > 0 && lag > permit) {
    lag <- permit
    str <- paste("Wasting too much s=data; setting lag to ", permit)
    messout(str,verbose)
  }
    
  if(ind == 1){
    
    mat <- cbind(lag(mat[,1], lag), mat[,-ind])
    colnames(mat) <- nm
    return(mat)
  }
  if(ind == ncl){
   
    inc <- ncl - 1; 
    mat <- cbind(mat[,1:inc],lag(mat[,ncl], lag))
    colnames(mat) <- nm
    return(mat)
  }
###  print(colnames(mat))  
  icov <- grep(cov, colnames(mat))
  if(length(icov) > 1)
    print("complaint")
  mat <- cbind(mat[,1:(icov-1)], lag(mat[,icov], lag),mat[,(icov+1):ncl]) 
  colnames(mat) <- nm
  
  return(mat)
}

### DESCRIPTION given the formula ff, it finds 
###             the lag years for each of the elements of ff; 
###             includes the rhs and lhs.
### INPUT    the formula: lag(log(allc.txt), 3) ~
###          lag(tobacco.txt^4, -7) + fat.txt + lag(gdp.txt, -2) + time
### OUTPUT a list with two elements: lag with the number of years to lag, named
###        with the covariates names.
###        lag = 3, 2 10. and name(lag) = "gdp", "tobacco", "fat"
###        The index ind is a flag; if ind=0 then the corresponding entry in ff
###        is in the lhs of the formula; if ind > 0 it is on the rhs of the formula.
###
### USES findlagyear and find.regx
###
### Elena Villalon
### evillalon@iq.harvard.edu
####################################################################
                                       
findlag <- function(ff)
  {
    ind <- grep("lag", ff)
    if(length(ind) <= 0 || ind == 1)
      return(ff)
  
    fc <- as.character(ff)
    ls <-  fc[2]
    rs <- fc[3]
    laglst <-  NULL
    ix  <- NULL
    lagg <- NULL
    nm <- NULL
    
    if(length(grep("lag",ls)) > 0) ##left-hand-side or depvar 
      {
### lagging in depvar: ls
        lagg <- findlagyear(ls)
        ix <- 0
        names(lagg) <- ls
        laglst <- list(lag=lagg, ind=ix)
        if (length(grep("lag", rs)) <= 0) 
          return(laglst)
          
      }
    if (length(grep("lag", rs)) > 0 ) ## lag the covariates data or rhs
      {
        vec <- strsplit(rs, "\\+")[[1]] ## get all covariates
        vec <- unlist(sapply(vec, FUN="trim.blanks"))
        indlst <- as.list(1:length(vec))
        
        indlst <- sapply(indlst, function(n, vec){
          comp <- vec[n]
          ind <- grep("lag", comp)
          ret <- NULL
          if(length(ind) > 0)
           ret <- n
          return(ret)}, vec)
       
        loc <- unlist(indlst)
        if(length(loc) <= 0)
          laglst <-  laglst ##nothing to lag 
        else if (length(loc) <= 1){ 
          lagadd <- findlagyear(vec[loc])
          lagg <- c(lagg, lagadd)
        }else{
          
          lagadd  <- sapply(as.list(vec[loc]), FUN="findlagyear")
      
          lagg <- c(lagg, unlist(lagadd))
      
        }
        
        if(length(ix) > 0 && ix ==  0)
          nm <- ls
        ix <- c(ix, loc)
       
        if(length(nm) >0 )
          rside <- vec[ix[-1]]
        else
          rside <- vec[ix]
       
        nmrs <- sapply(names(rside), FUN="find.regex")
        nrms <- unlist(nmrs)
           
        nm <- c(nm, nmrs)
        names(lagg) <- nm
        laglst <- list(lag=lagg, ind=ix)
    
        return(laglst)
      }

  }
### DESCRIPTION Takes a string of chars and numbers what and parses it
###             It looks for lag of years such as lag(gdp, -5),finds -5
### OUTPUT      A numeric with the years lag
###
### Elena Villalon
### evillalon@iq.harvard.edu
####################################################
findlagyear <- function(what)
  {
          
    resto <- strsplit(what, ",")[[1]][2]
###    print(resto)
      
    lag <- sub("([- 0-9]*)([\\)]+)([\\^ . 0-9]*)([a-z A-Z / \\* \\)]*)", "\\1", resto)
    
    
    lag <- as.numeric(lag)
         
    if (is.na(lag))
      return(NULL)
    else
      return(lag)
  }


### DESCRIPTION Takes two strings with names for depvar
###             numerator, denominator; and a matrix with data
###             It looks for first year obv (not NA) for depvar 
### OUTPUT      A numeric with the first observation year
###
### Elena Villalon
### evillalon@iq.harvard.edu
####################################################
       
first.nonay <- function(mat, dth, pop, fromcov=FALSE)
  {
   
    mat <- as.matrix(mat)
    rnm <- rownames(mat)
    ncl <- ncol(mat)


    md <- mat[,dth]
    dix <-  seq(along= md)[!is.na(md)]
    do <- ifelse(length(dix) > 0, min(dix), NA)
    mp <- NULL
    po <- NA
   
    if(length(pop) > 0){
      mp <- mat[,pop]
      pix <- seq(along=mp)[!is.na(mp)]
      po <- ifelse(length(pix) > 0, min(pix), NA)
    }
    if (!is.na(do) && !is.na(po))
      foy <- max(do, po)
    else if(!is.na(do))
      foy <- do
    else if(!is.na(po))
      foy <- po
    else
      foy <- NA
     
    if(length(rownames(mat)) >0 && !is.na(foy)){
      yearobs <- rownames(mat)[foy]
      yearlst <- list(yearobs=yearobs, ind = foy)

      return(yearlst)
    }
    nr <- nrow(mat)
    ncl <- ncol(mat)
    vst <- start(mat)[1]
    ved <- end(mat)[1]
###     print(vst)
###     print(ved)
    nrw <- ved[1] - vst[1] + 1 ##number of rows
    timeseries <- vst:ved
    rownames(mat) <- timeseries
       
    yearobs <- vst[1] 
    yearlst <- list(yearobs=yearobs, ind = foy)

    return(yearlst)
  }
### DESCRIPTION Takes an index, and two lists with datamatrices and first year 
###             give names to the rows of the corresponding matrix of index n 
### OUTPUT      matrix with named rows according to years from first observed. 
###
### Elena Villalon
### evillalon@iq.harvard.edu
####################################################
rowyears <- function(n, datamat, firstyear){
   mat <- datamat[[n]]
   fy <- firstyear[[n]]
   nr <- nrow(mat)
   vec <- fy:(fy+nr-1)
   rownames(mat) <- vec
   return(mat)
 }

catdepvar <- function(ff,datamat,dth){
  fc <-  as.character(ff)
  rs <- fc[3] ###right-hand-side
  ls <- fc[2] ###left-hand-side
  symb <- fc[1] ### "~"

  resto <- strsplit(rs, "\\+")[[1]]
  resto <- sapply(resto, FUN="trim.blanks")
  splrs <- resto
    
  ind <- grep(dth, resto)
    
  for(n in ind){
    
    morte <- resto[n] ###"lag(allc.txt^4, -2)"
    nudemorte <- trim.blanks(find.regex(morte)) ###"allc.txt"
    parts <- strsplit(morte, nudemorte)[[1]]
    copyofmorte <- paste("copyof_", nudemorte, sep="") ###"copyof_allc.txt"
    putback <- paste(parts[1], copyofmorte, parts[2], sep="") ##"lag(allc.txt.cov^4, -2)"

    res1 <- strsplit(resto[n], ",")[[1]][2]
    lag <- sub("([- 0-9]*)([\\)]+)([\\^ . 0-9]*)([a-z A-Z _ / \\*]*)", "\\1", res1)
    lag <- as.numeric(lag)
 #####   ****
  
    datamat <- lapply(datamat, function(mat){
           nm <- colnames(mat)
           nm <- c(nm, copyofmorte)
###           print(nudemorte)
           ind <- grep(nudemorte, colnames(mat))
           if(length(ind) > 0 ){
             mat <- cbind(mat, mat[, nudemorte])
             colnames(mat) <- nm
           }
           return(mat)})
    splrs[n] <- putback
  }
 
rs <- paste(splrs, collapse=" + ")
ffnew <- paste(ls, " ~ ", rs)
ffnew <- as.formula(ffnew)
lst <- list(datamat=datamat,ff= ffnew) 
return(lst)          
}
  
### helper function check.ops(c, cols, opvec)
### for any element, c, of a formula and for corresponding columns of
### given matrix of covariates and a vector of operations. It checks if
### c contains the operation in opvec (basically, multiplication and division), 
### and after anayzing its componets, if
### all components are included as elemenst of cols.
### For example gdp * hc, will spot the two elements gdp and hc and check
### if they are included in cols. Same for dth/pop, which checks if the
### division operation may be performed between dth and pop after checking
### if dth and pop are columns of cols.  It returns the index of cols
### to find the relevant components or NULL if they do not exist
##
### AUTHOR: Elena Villalon
###         IQS, Harvard Univ
###         evillalon@iq.harvard.edu
###         1737 Cambridge St
###         Cambridge MA
###
### Modified Sept 21, 2005
########################################################################

check.ops <- function(c, cols, opvec){
  
  ix <- unlist(sapply(opvec, grep,c))
  
  if(length(ix) <=0 ) ### no operation found; return the index
    return(grep(c,cols))

  op <- opvec[grep(names(ix), opvec)]

  dos <- strsplit(c,op)[[1]]
  dos <- unlist(sapply(dos, FUN="trim.blanks"))
  nix <-  unlist(sapply(dos, grep, cols))
  
  if(length(nix) <= 1) ### we need two cols entry for an operation
    return(NULL)
  else
    return(grep(dos[1],cols)) ### if both include return one of the index as OK 
}
