### Methods to generate the matrices of covaraites and dths
### according to a formula supplied with the arguments of yourcast
###
### driver for all functions
###
### AUTHOR: Elena Villalon
###         IQS, Harvard Univ
###         evillalon@iq.harvard.edu
###         1737 Cambridge St
###         Cambridge MA
###
### Modified Sept 27, 2005
###
########################################################################
input.to.model <- function(datamat,ff, all.pow=T, sample.frame=c(1950, 2000, 2001, 2030),
                           index.code=NULL, Gnames= NULL,  
                           standard=T, elim.collinear=T,
                           delta.tol=.Machine$double.eps,
                           tol=0.9999, solve.tol=1.e-10, log.poiss=F,verb=T)
  {
    

### argument log.poiss is just passed directly to output
    verbose <- verb
    yrest <- sample.frame[2]
    
    lst <- build.covs.depvar.lst(datamat, ff, all.pow, sample.frame, standard,verbose=verb)
 
    list.ff <- lst$list.ff
  

    whocov <- lst$whocov
    whoinsampx <- lst$whoinsampx
    whoutsampx <- lst$whoutsampx
    whoinsampy <- lst$whoinsampy
    whoutsampy <- lst$whoutsampy
    whopopul <- lst$whopopul
    whopopulos <- lst$whopopulos
    cov.lst0 <- list.covariates(whocov)
    
    ix.code <- parse.index.code(icode=index.code,verbose)
   
    cdgts <- ix.code["cntry.digits"]
    adgts <- ix.code["age.digits"]
    ydgts <- ix.code["year.digits"]
    fdgts <- ix.code["digit.first"]
 
    if(elim.collinear)
      {
    
        messout("Checking for collinear covariates...",verbose)
        lst <- find.collinear.covs(whoinsampy, whoinsampx, whoutsampx,
                                   whocov,delta.tol,tol, solve.tol,
                                   age.digits= adgts,year.digits=ydgts,
                                   cntry.digits=cdgts,verbose)
        
    
        whoinsampx <- lst$whoinsampx
        whoutsampx <- lst$whoutsampx
        whocov <- lst$whocov
        cov.lst <- list.covariates(whoinsampx)
  
        cov.omit <- find.omitted.cov(cov.lst, cov.lst0)
        cov.omit <- unlist(cov.omit)
       
        if(length(cov.omit) > 0){
         
          messout("Covariates deleted due to collinearities:", verbose,obj=cov.omit)}
          
      }else
    cov.lst <- list.covariates(whoinsampx)
    
    age.vec <- age.lst(tag = names(whoinsampy), cdigits=cdgts, adigits=adgts)
    cntry.vec <- cntry.lst(tag = names(whoinsampy), cdigits=cdgts, adigits=adgts)
   
    geoindex <- from.Gnames.to.geoindx(Gnames, cntry.vec)
###    print(geoindex)
### we need to check for    
### dataobj,sample.frame, formula,standardization, elim.collinear,
### tol?, and sove.tol?
        
    dataobj <- list(data=datamat, index.code=index.code, G.names=Gnames)
###
    lst <- list(dataobj=dataobj, sample.frame=sample.frame, formula=ff,
                standardization=standard, elim.collinear= elim.collinear,
                tol=tol, solve.tol=solve.tol,
### generated from the previous lists parameters,                  
                whocov=whocov, whoinsampx=whoinsampx, whoutsampx=whoutsampx,
                whoinsampy=whoinsampy, whoutsampy=whoutsampy,
                whopopul=whopopul, whopopulos=whopopulos,
                cov.lst= cov.lst, age.vec = age.vec,
                cntry.vec=cntry.vec, geo.index= geoindex, log.poiss=log.poiss)
    
###   print(names(lst))
    
    return(lst)
      
  }


  
### DESCRIPTION: Takes two list of covariates names and compare them
###              element by element, which every element is a csid unit
###
### INPUT: Two lists of covariates names generated before and after
###        elimination of linear colinearities.
### OUTPUT: the differences between the two input lists term by term
###         
###
### AUTHOR: Elena Villalon
###         IQS, Harvard Univ
###         evillalon@iq.harvard.edu
###         1737 Cambridge St
###         Cambridge MA
###
### Modified Sept 23, 2005
###
########################################################################
find.omitted.cov <- function(cov, cov0){
  ind <- 1:length(cov)
  names(ind) <- names(cov)
  
  omissions <- sapply(ind, function(n, cov,cov0){
    covmat  <- sort(cov[[n]])
    covmat0 <- sort(cov0[[n]])
    dff <- setdiff(covmat, covmat0)
    return(dff)}, cov, cov0)
  return(omissions)
}

### DESCRIPTION: Takes  list of datamat (with all the data for covs and depvars)
###              and the formula. Apply all functions to generate all the 
###              the necessary input lists for the yourcast's models.  
###
### INPUT: A list of csid units with cov and depvar data and the formula
###       
### OUTPUT: the insample, outsample lists for depvar, population, and covarites.
###         The list of formulas for every csid unit. 
###         
###
### AUTHOR: Elena Villalon
###         IQS, Harvard Univ
###         evillalon@iq.harvard.edu
###         1737 Cambridge St
###         Cambridge MA
###
### Modified Sept 27, 2005
###
########################################################################
build.covs.depvar.lst <- function(datamat, ff, all.pow=T,sample.frame,
                                  standard=T,verbose=T )
  {
    
### increasing years per row and for each matrix elements; consitency checks
   
    yrest <- sample.frame[2]
    datamat <- lapply(datamat, FUN="order.by.increasing.year", yrest)
   
### check for first observation in variable on the l-h-s of the formula
### eliminate the part with NA's before first observation of depvar or pop
    datamat <- data.after.first.obv(datamat, ff)
   
### generates a list of formulas one for every csid unit.
    messout("Creating formulas for all cross-sections...", verbose)
    
    
    lst.ff <- makeformulas(datamat, ff,verbose)
   
### parses lst.ff and builds the interpreted formulas
### one for every csid, to be used with model.matrix and model.frame
    messout("Applying formulas to data in each cross-section...", verbose)
    lst.ff <- lapply(lst.ff, FUN="parse.formula",all.pow)
  
### creates model.frames to be used to build the covs and depvars matrices
### according to the entries in datamat and the formulas for each csid
    messout("Running model frames for data matrices...", verbose)
    modfrm.lst <-   build.model.frames(datamat, lst.ff)
     
## takes modfrm.lst with the list of covs and depvars, 
## and for every csid unit extract the matrices of covariates
## according what is in the corresponding element of list formulas, lst.ff
## if standardization is required standard=T:
    messout("Building the covariates list...", verbose)   
    whocov <- build.covs(modfrm.lst,lst.ff, standard,verbose)
###    print("The covariates are: ");
    cov.lst <- list.covariates(whocov)
##    print(cov.lst)

## takes mod.frm.lst with list of depvar and population, 
## and for every csid unit extract the depvar and population
## Builds elements of two columns matrices one with the depvar as in the
## left-hand-side of the formula and the second column is the denominator
## term in the depvar, if any, otherwise creates
## a column of 1 for the population
       
    messout("Constructing the list of dependent variables...", verbose)
    dth <- build.depvar(modfrm.lst, datamat, lst.ff)

### We need to divide between the insample and outsample periods
### for the covariates according to the year yrest
    messout("Creating the in-sample and out-sample periods for covariates...", verbose)
    whox <- split.list(whocov, sample.frame)
    whoinsampx <- build.insampx(whox)
    whoutsampx <- build.outsampx(whox)

### split insample and outsample periods for depvar and population
    
    whoy <- split.list(dth, sample.frame,verbose)
### the following two lists' elements are  matrices of two columns: 
### the first is the depvar as specified in the r-h-s of the formula, 
### the second column is population. We have the insample and outsample
### periods depending on the year=yrest, which cretaes the T and F components
### for every list element of whoy
    messout("Creating lists for dependent variable,in-sample and out-sample...", verbose)
    whoinsampdepvar <- build.insampx(whoy)
    whoutsampdepvar <- build.outsampx(whoy)
 
### depvar is first column; population the second
  
    whoinsampy <- lapply(whoinsampdepvar,FUN="select.column",1, "depvar")
    whoutsampy <- lapply(whoutsampdepvar,FUN="select.column",1, "depvar")
   
  
    whopopul <-  lapply(whoinsampdepvar,FUN="select.column", 2, "popus")
    whopopulos <-  lapply(whoutsampdepvar,FUN="select.column", 2, "popus")

       
    lst <- list(list.ff = lst.ff, whocov=whocov, whoinsampx=whoinsampx, whoutsampx=whoutsampx,
                whoinsampy=whoinsampy, whoutsampy=whoutsampy, whopopul=whopopul, whopopulos=whopopulos)
     
     return(lst)
   }


### DESCRIPTION: takes a matrix and a year number representative of time series
###              Checks if time series is one of the columns of mat, if its
###              the case reorder the rows of mat with increasing years provided
###              the parameter choice is set to true (which is not the default). 
###              If time series is not a column of mat or choice=F, then checks if
###              it can be extracted from the rownames(mat), and if so
###              sort rows for increasing years.
###              If not a column or not in rownames it stops the simulation.
###
### INPUT: matrix mat and yrest to count the digits or characters of
###        numbers in time series. Latter is used in case that rownames
###        are the entire cstid (cntry+age+time) identifiers
###
### OUTPUT: The modified matrix mat order by rows names. 
###
### AUTHOR: Elena Villalon
###         IQS, Harvard Univ
###         evillalon@iq.harvard.edu
###         1737 Cambridge St
###         Cambridge MA
###
### Modified Sept 23, 2005
###
########################################################################
order.by.increasing.year <- function(mat, yrest, choice=F)
{
  cnm <- colnames(mat)
  cnm <- sapply(cnm, FUN="trim.blanks")
  ix <- ifelse( length(cnm) > 0, grep("time", cnm), NULL)
  if(length(ix) > 0)
    ix <- na.omit(ix)

  if(length(ix) > 0 && choice)
    {
      ord <- order(mat[,"time"])
      mat <- mat[ord,]
      return(as.data.frame(mat) )
    }
  
  ny <- nchar(yrest)
  if(length(rownames(mat)) <= 0)
    stop("Input datalist must have time series")
  
  ord <- order(as.numeric(rownames(mat)))
  mat <- mat[ord, ]
 
  return(as.data.frame(mat))
  
}

### DESCRIPTION: it takes the formula, ff, and split it into three elements
###              the left, "~", and rigt sides.  Splits the right-side into
###              components according to the separation sign="+"
###              If it finds the carat "^" sign, then it reads ^ as a power, 
###              and convert it into I(cov^pow) to be interpreted correctly
###              by the R-functions model.frame and model.matrix.
###
### INPUT : the formula ff <- log(dth) ~ gdp^3 + log(tobacco) + gdp * tobacco + hc
###         and a boolean to indicate if lower powers
###         should be included in the final formula
###
### OUTPUT: The final formula to be interprated with model.matrix and model.frame
###         for the example above with low.pow=T:
###         ff <- log(dth) ~ I(gdp^3) + I(gdp^2) + I(gdp^1) + log(tobacco) +
###               I(gdp* tobacco) + hc
###         If low.pow=F, then ff <- log(dth) ~ I(gdp^3) + log(tobacco) +
###               I(gdp* tobacco) + hc
###
### AUTHOR: Elena Villalon
###         IQS, Harvard Univ
###         evillalon@iq.harvard.edu
###         1737 Cambridge St
###         Cambridge MA
###
### Modified Sept 21, 2005
###
########################################################################
parse.formula <- function(ff, low.pow=T)
{
  fc  <- as.character(ff)
  rs0 <- fc[3]
  ls  <- fc[2]
 
  if(!identical(trim.blanks(fc[1]), "~") )
    stop("Bad parsing in formula")
  
  vec <- strsplit(rs0,"\\+")[[1]]
  ind <- grep("\\^", vec)
  
  if(length(vec[ind]) <= 0 )
    return(ff)
 
  carats <- vec[ind]### power of covs "^"
  gpars  <- vec[-ind] ### the rest
  
  if(length(gpars) > 0)
     rfx <- paste(gpars, collapse="+")
   else
     rfx <- ""
  
  for(ch in carats)
   {
     
     velem <- strsplit(ch, "\\^")[[1]]
    
###     print(velem)
     cov <- trim.blanks(velem[1])
     v2 <- velem[2]
     ind <- grep("\\)", v2)
     
     if(length(ind) > 0){

       add <- sub("([0-9.]+)([- a-z A-Z _ . 0-9 \\)]*)", "\\2", v2)
       cov <- paste(cov, add)
       v2  <- strsplit(v2, add)[[1]]
       
     }

     pow <- as.numeric(v2)
     if(is.na(pow))
       stop("wrong parsing")
     
 ###    print(pow)
 ###    print(cov)
     devp <- T
     fch <- paste("1 ~", ch)
     ndvar <- all.vars(as.formula(fch))
     op1 <- paste("log(", ndvar,"^",as.character(v2),")", sep="")
     op0 <- strsplit(op1, NULL)[[1]]
     ch0 <- strsplit(trim.blanks(ch), NULL)[[1]] 
     
     if(identical(intersect(op0, ch0),unique.default(op0))){
      
       devp <- F
     }
     ln <- length(ch0)
     rvar <- paste(ch0[1:ln],collapse="")
   
     if(identical(trim.blanks(op1), trim.blanks(rvar)))
        devp <- F
     
     if(pow > 1 && low.pow && devp)
       rs <- paste(" I(", cov,"^",paste(1:pow, ")",sep=""),sep="", collapse="+")
     else{
       if (!devp){
         rs <- paste("I(",ch, ")", sep="")
        
       }else 
       rs <- paste(" I(", cov, "^", pow, ")", sep="")
     }
     if(!identical(rfx,""))
       rfx <- paste(rfx, " + ", rs, sep="")
     else
       rfx <- rs
   
   } 
      
  ret <- as.formula(paste(ls, "~", rfx))
}


### 
### DESCRIPTION: it takes the formula, ff, and checks if the rigth-hand-side 
###              contains elements that are also columns of the matrices of datamat
###              It creates a list of formulas one for every element of datamat
###              Each of the element formula contains all or part of the right-hand-side
###              of ff, depending on whether the elements of ff are included in the
###              corresponding columns of datamat element's matrices.
###
### INPUT : the formula ff <- log(dth) ~ gdp^3 + log(tobacco) + gdp * tobacco + hc
###         and the list datamat of matrices of covariates.  Each element of datamat
###         is identified by csid (or cntry+age combo) and is a matrix with as many columns
###         as covariates that are used for the specified csid unit.
###
### OUTPUT: The list of formulas one for each csid units with
###          covariates that are included in both the formula input ff and the
###          matrix of covariates for the corresponding element of datamat
##
### AUTHOR: Elena Villalon
###         IQS, Harvard Univ
###         evillalon@iq.harvard.edu
###         1737 Cambridge St
###         Cambridge MA
###
### Modified Sept 21, 2005
########################################################################
  
makeformulas <- function(datamat, ff, choice=F,verbose=T)
{
  csid <- names(datamat)
  
  lst <- 1:length(csid)
  names(lst) <- csid
  
  fc  <- as.character(ff)
  rs  <- fc[3]
  ls  <- fc[2]
  rsvec <- strsplit(rs,"\\+")[[1]]
  ind <- grep("time", rsvec)
 
  if(length(ind) <= 0 && choice)
    {
      messout("Adding time to the covariate list...", verbose)
      rsvec <- c(rsvec, "time")
    }
  
  rsmod <- sapply(rsvec, FUN="findregex")
  
### I do not think I want to check operations yet
###  rsmod <- sapply(rsmod, FUN="split.ops")
  rsmod <- unlist(rsmod)
   
  lst.ff <- lapply(lst,function(n, datamat, rsmod){
    frm <- datamat[[n]]
    cols <- colnames(frm)

    cols <- sapply(cols, FUN="trim.blanks")
    
    vv <- unlist(sapply(rsmod, FUN="check.ops", cols, c("/", "\\*", ":")))
    
    vnm <- names(vv)
    
    vnm <- sapply(vnm, FUN="trim.blanks")

    rs <- paste(vnm, collapse=" + ")
    
    if(!identical(trim.blanks(rs), "") )
      ff <- as.formula(paste(ls, "~", rs))
    else
      ff <- NULL
    
    return(ff)
  }, datamat, rsmod)
    
return(lst.ff)
}
### DESCRIPTION: It takes a char or string to find the arguments of a function
###              or combinations of functions with possibles powers and other
###              operations among arguments. It returns its atomic components
###              ch <- "lag(log(gdp * time)^4, -5)"
###              find.regex(ch) returns "gdp *time"
###              ch <- sqrt(lag(gdp^0.5, -10)
###              returns "gdp"
###              ch <- lag(log(gdp/time)^4, -5)
###              returns "gdp/"time"
###
### Elena Villalon
### November 2005

find.regex <- function(ch)
  {
   findregex(ch)
 }


### helper function to buil.covs and buil.depvar
### for any element of datamat, and of formula list, lst.ff
### Given matrix of covs and depvars, and the corresponding 
### formula, it construct the model.frame
###
### INPUT: datamat a list with the data (covs and depvars)
###         for every csid, and the list of formulas to apply to datamat matrices.
### 
### OUTPUT: the list of model frames,
###         after applying model.frame to datamat and lst.ff
###
### AUTHOR: Elena Villalon
###         IQS, Harvard Univ
###         evillalon@iq.harvard.edu
###         1737 Cambridge St
###         Cambridge MA
###
### Modified Sept 23, 2005
########################################################################

 build.model.frames <- function(datamat, lst.ff)
     {
       indx <- 1:length(datamat)
       names(indx) <- names(datamat)
       m.lst <- lapply(indx, function(n,datamat,lst.ff)
                       {
                         mat <- datamat[[n]]
                         ffn <- lst.ff[[n]]
                         fcn  <- as.character(ffn)
                         lsn  <- fcn[2]
                
                         
             ###            mat[is.na(mat)] <- 9999
                         stress <- as.data.frame(mat)
                       
                         m <- model.frame(ffn, stress, na.action=NULL)
                         return(m)}, datamat, lst.ff)
       return(m.lst)
     }
### 
### DESCRIPTION: it takes the list of formulas, lst.ff, and the list of model.frames build
###              with function build.model.frame and 
###              for every cross sectional unit, and builds the expanded
###              matrices of covaraites depending on what is in the list of
###              formulas for every csid. If stand = T, standardize covs
###
### INPUT : the list of formulas of components of the form,
###         ff <- log(dth) ~ gdp^3 + log(tobacco) + gdp * tobacco + hc,
###         one for every csid unit with names as in datamat.
###         Any covaraites included in each element of lst.ff is also a column
###         of the corresponding unit of datamat or modfrm, which is the model frame we have built.  
###         The list modfrm contains matrices of covariates and depvars.
###         Each element of modfrm (and lst.ff), 
###         is identified by csid (or cntry+age combo) and is a matrix with as many columns
###         as covariates and depvars that are used for the specified csid unit.
###
### OUTPUT: A list of matrices similar to modfrm one for evry csid unit but
###         any matrix is the result of applying model.matrix to  modfrm
###         and with the formula for the csid unit in the input list of formulas.
###         They are extended matrices
###         with the operations prescribed in the list of formulas. 
###         The matrices of covaraites may be standardized if stand=T.
###         Collinearities may be eliminated is elim.coll = T
###        
### AUTHOR: Elena Villalon
###         IQS, Harvard Univ
###         evillalon@iq.harvard.edu
###         1737 Cambridge St
###         Cambridge MA
###
### Modified Sept 23, 2005
########################################################################
  
build.covs <- function(modfrm,lst.ff, stand,verbose=T)
  {
    indx <- 1:length(modfrm)
    names(indx) <-  names(modfrm)
    
    covlist <- lapply(indx, FUN="extract.covs", modfrm, lst.ff)
    
     if(stand){
      messout("Standardizing covariates...",verbose)     
      covlist <- lapply(covlist, FUN="standardX")
  
       } 
        return(covlist)
  
     
  }
### given the list with covariates, whose elements are identified with
### csid tags (cntry+age); the matrix elements contain the covarites specific
### to the csid (say, "245045" USA + age=45).  The function select the names
### of the covaraites that are included in every csid unit, and return that list
### which has as many elemnts as the list cov
###
### AUTHOR: Elena Villalon
###         IQS, Harvard Univ
###         evillalon@iq.harvard.edu
###         1737 Cambridge St
###         Cambridge MA
###
### Modified Sept 27, 2005
###
########################################################################
list.covariates <- function(cov){
  lst <- sapply(cov, colnames)
 ## names(lst) <-  NULL
 return(lst)
}
### 
### DESCRIPTION: it takes the list of formulas, lst.ff, and the list of model.frames build
###              with function build.model.frame and 
###              for every cross sectional unit, it builds the depvar matrix
###              Each matrix has two columns one with the left-hand-side
###              of the corresponding formula in lst.ff, for every cross-sectional unit, 
###              and the second column is population (popus)
###              If population is supplied in datamat, then it fills the second
###              column of depvar matrices with those values, 
###              otherwise it fills the popus column with 1's
###
### INPUT : the list of formulas of components of the form,
###         ff <- log(dth) ~ gdp^3 + log(tobacco) + gdp * tobacco + hc,
###         one for every csid unit with names as in datamat.
###         modfrm or model frames list build with build.model.frame function
###         that applies R-function model.frame to elements of datamat;  
###         datamat to extract the columns of population or the
###         denominator use for the depvar if any. 
###         The list modfrm contains matrices of covariates and depvars.
###         Each element of modfrm, datamat and lst.ff are  
###         identified by csid (or cntry+age combo) and is a matrix with as many columns
###         as covariates and depvars that are used for the specified csid unit.
###
### OUTPUT: A list of matrices similar to modfrm one for evry csid unit but
###         any matrix is the result of applying model.matrix to  modfrm
###         and with the formula for the csid unit in the input list of formulas.
###         They are extended tow columns matrices with depvar and population
###         or denominator. 
###        
### AUTHOR: Elena Villalon
###         IQS, Harvard Univ
###         evillalon@iq.harvard.edu
###         1737 Cambridge St
###         Cambridge MA
###
### Modified Sept 23, 2005
########################################################################
  
build.depvar <- function(modfrm,datamat, lst.ff)
  {
    
    ind <- 1:length(modfrm)
     
    names(ind) <- names(modfrm)
  
    dthlist <- lapply(ind, FUN="extract.depvar", modfrm,datamat,lst.ff)

   
    return(dthlist)
  }
### helper function to build.covs
### for given matrix applies model.matrix to the matrix
### modfrm[[n]] and lst.ff[[n]], n is integer
###
### INPUT : list of model frames, modfrm,  and list of
###         formulas lst.ff, index n to extract corresponding matrices
###
###OUTPUT: The resulting matrix after applying model.matrix
###        
### AUTHOR: Elena Villalon
###         IQS, Harvard Univ
###         evillalon@iq.harvard.edu
###         1737 Cambridge St
###         Cambridge MA
###
### Modified Sept 23, 2005
########################################################################

extract.covs <- function(n, modfrm,lst.ff)
  {
   
      m <- modfrm[[n]]
      ffn <- lst.ff[[n]]
      mat <- model.matrix(ffn, m)
      
      nm <- sapply(colnames(mat),FUN="colnames.cov")
###      print(length(m))
      colnames(mat) <- nm
      
      return(mat)
    }

### helper function to build.depvar
### for given matrix applies model.matrix to the matrix
### modfrm[[n]] and lst.ff[[n]], n is integer.
### Also extract population matrix from datamat[[n]]
### if included as one of its columns, otherwise
### it fills population column of depvar matrix with 1's.
###
### INPUT : list of model frames, modfrm,  datamat (original list of data),
###         and list of formulas lst.ff,
###         the index n to extract corresponding matrices
###
###OUTPUT: The resulting matrix after applying model.matrix
###        and extracting the depvar
###        and adding population as a column
###        
### AUTHOR: Elena Villalon
###         IQS, Harvard Univ
###         evillalon@iq.harvard.edu
###         1737 Cambridge St
###         Cambridge MA
###
### Modified Sept 23, 2005
########################################################################
    
extract.depvar <- function(n, modfrm,datamat,lst.ff)
  {

    m <- modfrm[[n]]
    ffn <- lst.ff[[n]]
    mat <- as.data.frame(datamat[[n]])
    
    fc  <- as.character(ffn)
    rs  <- fc[3]
    ls  <- fc[2]
    
    morta  <- m[1]
    morta1 <- m[ls]
                        
    if(any(na.omit(morta) != na.omit(morta1)))
      stop("Error: Applying formula")
  
    morta <- as.matrix(morta)
    nr <- nrow(morta)
 
      
    ixdiv <- grep("/", colnames(morta))
    
    if(length(ixdiv) <= 0){
      popus <- as.matrix(rep(1, nr))
      colnames(popus) <- "popus"
       
      depvar <- cbind(depvar=morta, popus = popus)
      return(depvar)
    }
    
    
    lst <- split.ops(findregex(ls))
    dth <- lst[1]
    pop <- lst[2]
###    print(pop)
    popus <- mat[pop]
### Fixing zeros
    
###   print(nm)
   
 
    return(cbind(depvar=morta, popus=popus))
     
  }
### DESCRIPTION: Takes  list of datamat (with all the data for covs and depvars)
###              and the formula. Find the l-h-s of formula with the depvar 
###              and it search for one or two components (or numerator and denominator)
###              Call function first.obvy for each element of datamat and
###              changes the matrices by substracting those rows before first observations
###              of depvar and population. 
###
### INPUT: dtamat list of csid units and formula
###       
### OUTPUT: A list of csid units with modified datamat if applicable. 
###         
###
### AUTHOR: Elena Villalon
###         IQS, Harvard Univ
###         evillalon@iq.harvard.edu
###         1737 Cambridge St
###         Cambridge MA
###
### Modified Sept 27, 2005
###
########################################################################


data.after.first.obv <- function(datamat=NULL, ff, noadjust=T)
  {
    nm <- names(datamat)    
    fc   <- as.character(ff)
    ls   <- fc[2]
    indl <- grep("lag", ff)
    if(length(indl) > 0 && noadjust)
      return(datamat)
      
###    lst <- dthpop(ff)
###    dth <- lst$dth
###    pop <- lst$pop
    lst <- split.ops(findregex(ls))
    if(length(lst) <= 1){
          dth <- lst
          pop <- NULL
        }else{
          dth <- lst[1]
          pop <- lst[2]
        }
### make sure you do not have values of death below 0.5: user should know better
###    print("Setting 0.5 any death value lower")
###    datamat <- lapply(datamat, function(mat){
###      ix <- grep(dth, colnames(mat))
###      if(length(ix) > 0){
###        dd <- mat[, ix]
###       dd[ dd <= 0.5 ] <- 0.5
###        mat[, ix] <- dd}
###        return(mat)})

  
### Examples: from allc/popu + 1 to pop= popu, dth=allc
### from allc/popu^4 + 5 to pop=popu, dth=allc    
  
   datamat <- lapply(datamat, FUN= "first.obvy", dth, pop)
      
   names(datamat) <- nm
  

  return(datamat)
  }



### DESCRIPTION: Takes a matrix mat (with all the data for covs and depvars)
###              for given csid unit, and the names of two columns of the matrix mat.
###              Find first observation for any of those two columns and delete those rows
###              of mat that are before first observed along dth and pop.
###
### INPUT: mat matrix, and two columns
###       
### OUTPUT: The modified mat by deleting rows before first observed along the columns dth and pop. 
###         
###
### AUTHOR: Elena Villalon
###         IQS, Harvard Univ
###         evillalon@iq.harvard.edu
###         1737 Cambridge St
###         Cambridge MA
###
### Modified Sept 27, 2005
###
########################################################################
first.obvy <- function(mat, dth, pop)
  {
    mat <- as.matrix(mat)
    rnm <- rownames(mat)
    
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
    if(!is.na(foy) && foy > 1)
      {  
        mat <- mat[-(1:(foy-1)), ]
        if(length(rnm) >0){
          rnm <- rnm[-(1:foy-1)]
          rownames(mat) <- rnm
        }else{
          nr <- nrow(mat)
          ncl <- ncol(mat)
          vst <- start(mat)[1]
          ved <- end(mat)[1]
###     print(vst)
###     print(ved)
          nrw <- ved[1] - vst[1] + 1 ##number of rows
          timeseries <- vst:ved
          rownames(mat) <- timeseries
        }
        
      }
    nr <- nrow(mat)
  
   
    return(mat)
  }
#### helper function to standardize covariates
standardX <- function(ch)
{
 
  ch <- as.data.frame(ch)
  ch <- scale(ch, center=T, scale=T)
  ch[is.nan(ch)] <- 1
  return(ch)
}
### helper function to name the columns of the matrices of covaraites
### with meningful names that include any interpretaion of formulas.
###
### AUTHOR: Elena Villalon
###         IQS, Harvard Univ
###         evillalon@iq.harvard.edu
###         1737 Cambridge St
###         Cambridge MA
###
### Modified Sept 23, 2005
###
########################################################################
colnames.cov <- function(ch)
{
  if(identical(ch, "(Intercept)"))
    return("cnst")
           
 ### all omnipotent regular expresssions
 ### from I(gdp^4) to gdp^4
 ### from I(gdp^0.4) to gdp^0.4
 ### from I(gdp^3 * hc) to gdp^3 * hc
 ### from I(gdp/hc^4) to gdp/hc^4
  
  ch <- sub("I\\(([a-z ^ 0-9 . * _ / \\( \\)]*)+\\)",'\\1', trim.blanks(ch))
 
  ch <- trim.blanks(ch)
}
###
### DESCRIPTION: build the insample and outsample lists of covaraites
### INPUT : datamat is a list of covariates matrices
###         yrest is the year to spli between insample and outsample periods.
###
### OUTPUT: the list whose elements, mat, are composed of two matrices corresponding  
###         to two values mat$T and mat$F for insample and outsample
###
### AUTHOR: Elena Villalon
###         IQS, Harvard Univ
###         evillalon@iq.harvard.edu
###         1737 Cambridge St
###         Cambridge MA
###
### Modified Sept 22, 2005
########################################################################
split.list <- function(datamat, sample.frame,verbose=T)
  {
    yrest <- sample.frame[2]
     
    whox <- lapply(datamat, FUN="split.whocov",sample.frame,verbose)
    return(whox)
   }

build.insampx <- function(whox)
     whoinsampx <- lapply(whox, function(mat) mat$T)
build.outsampx <- function(whox)
     whoutsampx <- lapply(whox, function(mat) mat$F)
   
### Same as build.insamp for the outsample period



### helper function to build.insamp(outsamp) and split.list
### It takes any given matrix in the covariates list whocov or in dth from datamat
### and according to specifications of the yrest split between
### insample and outsample periods. If time series is present
### such as for the covariates list, the it takes time series and divide
### according if years <= yrest or years > yrest.  If time series
### is not present such as for the dependent variable list,
### then it looks at the row names of the matrix and split according
### to years.
###
### AUTHOR: Elena Villalon
###         IQS, Harvard Univ
###         evillalon@iq.harvard.edu
###         1737 Cambridge St
###         Cambridge MA
###
### Modified Sept 27, 2005
########################################################################
 adjust.digits <- function(y, nch) {
   nchy <- nchar(y)
          
   if(nchy < nch){
     pad <- rep(0,nch-nchy)
     pad <- paste(pad, collapse="")
     y <- paste(pad, y, sep="")

   }
   return(y)
 }

split.whocov <- function(mat,sample.frame,verbose)
  {

    yrest <- sample.frame[2]
    ylast <- sample.frame[4]
    
    mat  <- as.matrix(mat)
    rwnm <- rownames(mat)
    nca  <- nchar(rwnm[1])
    ncy  <- nchar(ylast)
    anos  <- substring(rwnm[length(rwnm)],nca- ncy +1)
    lanos <- trim.blanks(anos) 
    
### otherwise used the tags names for rows. 
    yrest <- trim.blanks(as.character(yrest))
    ylast <- trim.blanks(as.character(ylast))
   
    if(!identical(lanos, ylast)){
      messout("Expanding matrices to total number of years",verbose)  
      tags <- sapply(rwnm, substr, 1, nca - ncy)
      ind  <- (as.numeric(lanos)+1):as.numeric(ylast)
      matexpand <- matrix(data=NA, nrow=length(ind), ncol=ncol(mat))
      yy   <- sapply(1:length(ind), function(n){
        ret <- paste(trim.blanks(tags[n]), trim.blanks(ind[n]), sep="")
        return(ret)
      })
      rownames(matexpand) <- yy
      mat <- rbind(mat, matexpand)
    }

    ncy <- nchar(ylast) - nchar(yrest)
    if(ncy > 0){
      pad <- paste(rep(0,ncy), collapse="")
      yrest <- paste(pad,yrest,sep="")
    }
   
    nch <- nchar(yrest)
    tags  <- rownames(mat)
  ###  print(rownames(mat))
    
    nc <- nchar(tags[length(tags)])
    lb <- nc - nchar(yrest)
    csid0 <- NULL
    years <- sapply(tags, substring,lb+1)
    ultimo <- years[length(years)]
    
    if (lb > 0 && length(lb) > 0)
      {
    
        years1 <- sort(years)
        
        if(any(years != years1))
           messout("Matrices should be sorted by year",verbose)

        ln <- length(years1)
        ultimo <- years1[ln]
        
        csid <- sapply(tags, substr,1,lb)
        nch  <- nchar(as.character(yrest))
        years <- sapply(years, FUN="adjust.digits",nch) 
     
        ix <- grep(as.character(yrest),as.character(years))
        
        if(length(ix) <= 0){
          yrest <- min(as.numeric(years))
        
        messout(paste("yrest for insample-outsample no in years range; adjusted yrest is " ,yrest,sep=""),verbose)
        }
        
        csid0 <- trim.blanks(csid[ix])
      
        todiv <- paste(csid0, yrest, sep="") 
        mat <- split.data.frame(mat,rownames(mat) <= todiv)
       
      }else{
        
        mrw <- sapply(rownames(mat), FUN="adjust.digits",nch)
       rownames(mat) <- mrw
        mat <- split.data.frame(mat,rownames(mat) <= yrest)
      }
    
    if(ultimo > ylast)
      {
        mxy <- paste(csid0,ylast, sep="")
        matout <- split.data.frame(mat$F, rownames(mat$F) <= mxy)$T
        matin  <- mat$T
        mat <- c(list(matin), list(matout))
        names(mat) <- c(TRUE, FALSE)
      }
     
    return(mat)
  }
                            
                            
### given a matrix x, a column number and a tag
### Select the column from the matrix and give it a name.

 select.column <- function(x, col=1, name="depvar")
      {
       
        x <- try(as.matrix(x[,col]), silent=T)
        if(class(x) == "try-error"){
      
          stop("Selecting the wrong column")
        }
        colnames(x) <- name
        return(x)
      }
###
### INPUT :the list names for either whoinsampy, whoutsampy or
###        any of covs, such as whoinsampx
###        the number of digits for country codes, cdigits and/or age difits, adigits
###
### OUTPUT: vector with country codes, cntry.lst,
###         and with age group numbers, age.lst
###
### AUTHOR: Elena Villalon
###         IQS, Harvard Univ
###         evillalon@iq.harvard.edu
###         1737 Cambridge St
###         Cambridge MA
###
### Modified Sept 27, 2005
########################################################################

cntry.lst <- function(tag, cdigits, adigits)
      {
        cc <- nchar(tag[1])
        if(cc != (cdigits + adigits))
          stop("tags names in list are wrong")
        
        cntry.vec <- sapply(tag, substring, 1, cdigits)
        cntry.vec <- unique.default(unlist(cntry.vec))
        
        return(cntry.vec)
      }

age.lst <- function(tag =NULL, cdigits, adigits)
  {
    
    cc <- nchar(tag[1])
    if(cc != (cdigits + adigits))
      stop("tags names in list are wrong")
        
    st <- cdigits + 1 + adigits
    age.vec  <- sapply(tag, substring, cdigits+1, st)
 
    age.vec0 <- sapply(tag, substring, cdigits+1)
    
    if(any(age.vec0 != age.vec))
      stop("tags names in list are wrong")
    
    age.vec <- unique.default(unlist(age.vec))       
    return(age.vec)
  }
    
  ## ************************************************************************
## ************************************************************************
##
## FUNCTION NAME:      find.colinear.covs
##
## 
## IMPORTED functions: elim.colinear (also in this module)
##
## USED GLOBALS:We use whoinsampy, whoinsampx, whoutsampx, and whocov,
##             and a bunch of tolerance parameters. Number of digits
##             in age groups and in country codes. 
##
## DESCRIPTION: It eliminates some of the covariates for csid identifiers
##              of list elements whoinsampx,whoutsampx, whocov.
##              Criteria for colinearities with tolerance values, delta.tol & tol
##              as defined with elim.colinear.
##
##
## VALUE:  Returns modified values of whoinsampx, whoinsampy,whocov and
##
## WRITTEN BY: Elena Villalon 
##             evillalon@iq.harvard.edu, 
##             CBRSS, Harvard University
## 
## Modified Sept 27, 2005
##
## ************************************************************************
## ************************************************************************

  
find.collinear.covs <- function(whoinsampy, whoinsampx, whoutsampx,
                                whocov,delta.tol=.Machine$double.eps,
                                tol=0.9999, solve.tol=1.e-10,
                                age.digits, year.digits, cntry.digits,verbose=T)
{
    
  csid <- 0 *vector(, length=length(whoinsampy))
  indx <- 1:length(whoinsampy)
  age.vec <- age.lst(names(whoinsampy), cntry.digits, age.digits)
  age.vec <- unlist(age.vec)
  n.age <- length(age.vec)
  
  for(i in indx){
###    print(names(whoinsampy)[i])
    y <- whoinsampy[[i]];
    n.obs <- nrow(na.omit(y));
   
    x <- whoinsampx[[i]];
    csid[i] <- names(whoinsampy[i])
    nmx <-  csid[i]
    w <- as.numeric(is.na(y) == FALSE) ### 0 if NA, 1 otherwise ###
    xout <- whoutsampx[[i]];
    yw  <- y;
    xw  <-  x;
    yxw <- cbind(yw, xw)
    yxw <- na.omit(yxw)
    yw  <- yxw[,1]
    xw  <- yxw[,-1]
    
    colin  <- elim.colinear(xw,nmx, delta.tol, tol, solve.tol,age.digits, year.digits,verbose)
    xw     <- colin$xw
    invxw  <- colin$invxw
    tol    <- colin$tol
    ri     <- colin$covind
   
 ### Just not to have a huge I/O we only print for one age group, since that
 ### except for tobacco they behave very similarly
    if (length(ri) >0  ) {
### printing for anormal conditions 
###      if(i%% n.age <= 0 || i%%n.age == floor(n.age/2))
###        print(paste("Covariate(s)",whocovariates[ri],"eliminated from cs",
###                    substr(nmx, 1, cntry.digits)));
      
      xout <- xout[, -ri]
      x <- x[,-ri]
      whocov[[i]] <- whocov[[i]][,-ri]
      whoinsampx[[i]] <- x
      whoutsampx[[i]] <- xout
      
    }
  
  }
 lst <- list(whocov=whocov, whoinsampx=whoinsampx, whoutsampx=whoutsampx)
 return(lst)
}
## ************************************************************************
## ************************************************************************
##
## FUNCTION NAME: elim.colinear 
##
## INPUT:     element of whoinsampx or a matrix of covariates for given csid 
##            rows = years of observations insample matrix & cols= covariates; 
##            tol, tolerance for linear dependencies among cols of whoinsampx
##
## DESCRIPTION: if matrix t(xw) %*% xw (where xw= whoinsampx) cannot be inverted, 
##              it calls covariate.delete to eliminate the linear dependencies with 
##              tol that decreases from call to call by delta.tol
##              until built-in function solve() inverts t(xw) %*% xw   
##
## OUTPUT:    matrix element of whoinsampx modified after eliminating lineariarities 
##            the inverted matrix of (t(xw) %*% xw); and the tol that remains after 
##            eliminating correlations in decrement steps of delta.tol
##
##  WRITTEN BY: Elena Villalon
##              evillalon@latte.harvard.edu
##              CBRSS, Harvard University
##
## Last modified: 10/12/2003
## 
## ************************************************************************
## ************************************************************************
 elim.colinear <- function(xw,nmx, delta.tol, tol, stol, age.digits, year.digits, paso=10,verbose=T)
{
   
   ri <- vector(,length=0)
   xw <- as.matrix(xw)
   xww <- xw
   options(show.error.messages=F)
   
### tol=1.e-10 may be changed as part of the options of solve
    invxw <- try(solve(t(xw) %*% xw, tol=stol))
   
   if( !inherits(invxw, "try-error"))
      return(list(xw =xw, covind=ri,invxw=invxw, tol= tol))
  
     n <- floor( tol/delta.tol)
  
     while(n > 0 ){
      
       tol <- tol - delta.tol
       options(show.error.messages=T)
    
       ret <- covariate.delete(xw, tol, ydigits=year.digits, adigits=age.digits,verbose)
       options(show.error.messages=F)               
       xww  <- as.matrix(ret$covx)
       ri   <- ret$covind
       invxw <- try(solve(t(xww) %*% xww, tol=stol))
       
       if(inherits(invxw, "try-error") && tol > delta.tol){
         n   <-  n - 1
         delta.tol <- delta.tol * paso
       }else{
         n <- 0
         xw <- xww
         delta.tol <- tol
       }
     
     }
   if(inherits(invxw, "try-error")) 
     messout(paste("Collinearity not eliminated for csid= ", nmx, sep=""),verbose)
  
   options(show.error.messages=T)        
  return(list(xw =xw, covind=ri,invxw=invxw, tol= tol))
 }


covariate.delete <- function(xmat, tol, ydigits, adigits,verbose=T) {
 

  if (dim(xmat)[2] <= 1)
    return(list(covx=xmat,covind=NULL,covnam = NULL))

  x  <- data.frame( xmat)
  vn <- colnames(x)
  bf <- length(x)

### find cols with constant elements of sd =0

  if(length(sd(x, na.rm=T) <= 0))
        return(list(covx=xmat,covind=NULL,covnam = NULL))
     
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
  csid <-  as.numeric(row.names(x)[1])%/%10^(ydigits)
### Because of the I/O we only print once the message below
  if( abs(bf -af) >0 & csid %% 10^(adigits) == 45) 
    messout(paste("Colinears, ncols deleted = ", bf - af, " for tol= ", tol,
            "; and cntry= ",csid%/%10^(adigits), sep=""),verbose)
  if (length(ri) > 0) xmat <- xmat[,-ri]
### returns covariate matrix with colinearities removed if appropiate
### index of columns being removed, and their names
  return(list(covx=xmat,covind=ri,covnam = rn))
}  



## ************************************************************************
## ************************************************************************
##
## FUNCTION NAME: parse.index.code 
##
## INPUT:     A string, icode, with the coding for input data cross-sectional time series id.
##
## DESCRIPTION: Parses the string icode using extending regular expressions 
##              It look for symbols, "-", "g G", "a A", "t T", and counts them. 
##
## OUTPUT:    The number of entries for each of the symbols, 
##            nodiscard for "-" at begin of string,
##            nogeo for "g" and "G"; noage for "a" and "A",
##            notime for "t" or "T", and noends for "-" at end of string
##
##  WRITTEN BY: Elena Villalon
##              evillalon@iq.harvard.edu
##              IQSS, Harvard University
##
## Last modified: 09/29/2005
## 
## ************************************************************************
## ************************************************************************

parse.index.code <- function(icode, verbose=T,yrest=NULL){
### example icode = "--ggggaatttt--"
  
  icode <- trim.blanks(icode)
  nodiscard <- 0 ### number of "-" at beginning of string
  nogeo  <- 0 ### number of "g" 
  noage  <- 0 ### number of "a"
  notime <- 0 ### number of "t"
  noends <- 0 ### number of "-" end of string
  
### parsing any number of "-", at beginning only  
  ix <- grep("-", icode)
  ln.char <- nchar(icode) 
  if(length(ix) > 0)
    {
      icode <- trim.blanks(sub("-*","",icode))
      nodiscard <- ln.char - nchar(icode) 
      ln.char <- nchar(icode)
    }
### parsing geo index  either G or g: unespecified number
  ix <- grep("[G g]", icode)
  if(length(ix) > 0){
    
    icode <- trim.blanks(sub("([g G]*)", "",icode))    
    nogeo <- ln.char - nchar(icode)
    ln.char <- nchar(icode)
    
  }
 ### parsing age index  either A or a
  ix <- grep("[A a]", icode)
  if(length(ix) > 0)
    {
      icode <- trim.blanks(sub("([a A]*)", "", icode))
      noage <- ln.char  - nchar(icode)
      ln.char <- nchar(icode)
    }
### parsing time index  either T or t
  ix <- grep("[T t]", icode)
  if(length(ix) > 0)
    {
      icode <- trim.blanks(sub("([T t]*)","", icode))
      notime <- ln.char - nchar(icode)
    }
  noends <- nchar(icode)
  digit.first <- nodiscard
  cntry.digits <- nogeo
  age.digits <- noage
  year.digits <- notime

  if(length(yrest) > 0 && year.digits != nchar(as.character(yrest)))
      messout("sample.frame and index.code different digits",verbose)
  
  lst <- c(digit.first=nodiscard, cntry.digits=nogeo, age.digits=noage,
           year.digits=notime, noends=noends)
  return(lst)
}

from.Gnames.to.geoindx <- function(Gnames=NULL, cvec)
{
  

  if(length(Gnames) > 0)
    {
     
      geo.indx <- as.vector(Gnames[,2])
      names(geo.indx) <- as.vector(Gnames[,1])
      return(geo.indx)
    }
  
  geo.indx <- rep(" ", length(cvec))
  names(geo.indx) <- cvec

  return(geo.indx)
}


### DESCRIPTION: test for the equality of two formulas
###              If either the left-hand-sides or right-hand-sides
###              are different the formulas are differents
### INPUT two formulas, ff0, ff1
### OUTPUT a boolean either True or False depending if ff0 == ff1 or not
###
### AUTHOR Elena Villalon
###        IQSS, Harvard Univ
###        evillalon@iq.harvard.edu
###
###  DATE October 7 2005
#################################################################

equal.formulas <- function(ff0,ff1)
{
  fc0  <- as.character(ff0)
  rs0  <- fc0[3]
  ls0  <- trim.blanks(fc0[2])
  
  fc1  <- as.character(ff1)
  rs1 <- fc1[3]
  ls1  <- trim.blanks(fc1[2])
  
  if(!identical(trim.blanks(fc0[1]), "~") && !identical(trim.blanks(fc1[1]), "~") )
    stop("Bad parsing in formulas")
  
  vec0 <- strsplit(rs0,"\\+")[[1]]
  vec1 <- strsplit(rs1,"\\+")[[1]]
  vec0 <- sapply(vec0, FUN="trim.blanks")
  vec1 <- sapply(vec1, FUN="trim.blanks")
  vec0 <- sort(vec0)
  vec1 <- sort(vec1)

  equality <- ifelse(length(vec0) == length(vec1), T, F)
  if(!equality)
    return(equality)
  
  equality <- ifelse(identical(ls0, ls1), T, F)
  if(!equality)
    return(equality)

  if(any(vec0 != vec1))
    equality <- F
  else
    equality <- T
  return(equality)
}
