## ************************************************************************
## ************************************************************************
##
##   THIS DOCUMENTATION HAS TO BE CHANGED
##
## FUNCTION NAME:  build.C.age.time (driver) & cov.age.cor (executed with build.C.age.time)
##
## PLACE:    
##
## IMPORTED:    build.laplacian and covariates matrices whocov from
##              make.mortality.data()
##
## DESCRIPTION: It produces correlation matrices among all covariates for all age groups
##              and each contry. Matrices are named with the code of cntry plus
##              the two age groups (say "24503555"), and are products of covariates.
##              It calls cov.age.cor to construct for every cntry the
##              correlations matrices among covariates and ages
##
## FORMAT: res <- build.C.age.time(wcntry), which requires
##                the function cov.age.cor(cntry, W, age.char,age.comb)
##         
## INPUT:   whocov, cntry.vec, age.vec (globals),and wcntry
##        
## OUTPUT:  global, who.C.age,  list with  cntrys, and for each cntry,
##          correlation of covariates for age groups 
##
## WRITTEN BY: Elena Villalon & Federico Girosi  
##             evillalon@latte.harvard.edu; fgirosi@latte.harvard.edu
##             CBRSS, Harvard University
## 
## Last modified: 07/01/2003
## 
## ************************************************************************
## ************************************************************************

build.C.age.time <- function(W.age,time.prior.param, env.base) {
  ebase <- get("env.base", envir=parent.frame())
  ewho <- get("env.who", envir=ebase)
  age.vec <- get("age.vec", envir=ewho)
  cntry.vec <- get("cntry.vec", envir=ewho)
  verbose <- get("verbose", envir=ewho)
  n.cntry <- length(cntry.vec)
  n.age <- length(age.vec)
 
  who.digit.first  <-  get("who.digit.first", envir=ewho)
  who.cntry.digits <-  get("who.cntry.digits", envir=ewho)
  who.age.digits   <- get("who.age.digits", envir=ewho)
  whocov <- get("whocov", envir=ewho)
 
### structure of dataset cstsid: 
  digit.cntry.begin <- who.digit.first + 1 
  digit.cntry.end   <- who.digit.first + who.cntry.digits
  digit.age.begin   <- digit.cntry.end + 1
  digit.age.end     <- digit.cntry.end + who.age.digits

  age.char <- formatC(age.vec,width= who.age.digits, format="d", flag="0")
  paste.char <- function(x, y){ paste(x, y, sep="")}
  age.comb <- kronecker(age.char, age.char, FUN=paste.char)
### check Wij to store only those elements of age.comb, for which Wij != 0
### it helps to give names

  if(length(rownames(W.age)) <= 0 || length(colnames(W.age)) <= 0){
    rownames(W.age) <- age.char
    colnames(W.age) <- age.char}
 
  who.C.age <- list()
  for(i in 1:n.cntry){
   
    Ci <- cov.age.time.cor(cntry.vec[i],W.age,time.prior.param,age.char,age.comb,env.base)
    who.C.age <- c(who.C.age, Ci)}
### just for consistency checking; testing but not needed in code 
  g1 <- grep(paste("^",as.character(cntry.vec[1]),sep="") ,names(who.C.age))
  n.g1 <- length(g1) 
  n.g <- sapply(1:n.cntry, function(i){
    gi <- grep(paste("^",as.character(cntry.vec[i]),sep="") ,names(who.C.age))
    n.gi <- length(gi)
    return(n.gi)})
  if(any(n.g != n.g1))
    messout("Error in build.Caa, number of matrices for countries not agree", verbose)
### end testing
  
  return(invisible(who.C.age))
}

## FUNCTION NAME:  cov.age.time.cor; (executed with build.C.age.time)
##
## PLACE:    
##
## IMPORTED:    covariates matrices whocov from make.mortality.data()
##              environments with globals, including countries and age groups 
##
## DESCRIPTION: It produces a correlation matrix among contributing covariates for all age groups
##              and chosen contry. Matrices are named with the code of cntry plus
##              the two age groups that are correlatd (say "24503555"), and are products of covariates.
##              Because of symmetry we assume first age group <= second group;  
##              
##
## FORMAT: res <- cov.age.time.cor(cntry, W,time.prior.param, age.char,age.comb, ebase)
##         
## INPUT:   one cntry from cntry.vec; W for correlation among ages;
##          time.prior.param (time correlation); age group as a char age.char;
##          and age.comb for age combinations to be correlated;
#3          the environment, ebase, for globals
##        
## OUTPUT:  one element of list who.C.age,  which is a matrix of 
##          correlation among covariates and all contributing age groups 
##
## WRITTEN BY: Elena Villalon & Federico Girosi  
##             evillalon@latte.harvard.edu; fgirosi@latte.harvard.edu
##             CBRSS, Harvard University
## 
## Last modified: 07/01/2003
## 
## ************************************************************************
## ************************************************************************

cov.age.time.cor <- function(cntry, W,time.prior.param, age.char,age.comb, ebase){
  ebase <- get("env.base", envir=parent.frame())
  ewho <- get("env.who", envir=ebase)
  verbose <- get("verbose", envir=ewho)
  env.base <- ebase
  
  age.vec <- get("age.vec", envir=ewho)
  who.digit.first  <-  get("who.digit.first", envir=ewho)
  who.cntry.digits <-  get("who.cntry.digits", envir=ewho)
  who.age.digits   <- get("who.age.digits", envir=ewho)
  whocov <- get("whocov", envir=ewho)
### structure of dataset cstsid: 
  digit.cntry.begin <- who.digit.first + 1 
  digit.cntry.end   <- who.digit.first + who.cntry.digits
  digit.age.begin   <- digit.cntry.end + 1
  digit.age.end     <- digit.cntry.end + who.age.digits
  
  ctr <- paste("^", as.character(cntry),sep="")

### can only grep one value or element at a time

  indx <- grep(ctr, names(whocov))
  if(length(indx) <= 0) {
  messout(paste("Warning wrong set of countries and names...check preprocessing", ctr, sep=""),verbose)
  return(NULL); }
  
  n.ind <- length(indx)
  n.age <- length(age.char)
  vnm <- paste(as.character(cntry),age.comb,sep="")
  whocov[indx] <- lapply(whocov[indx], na.omit)
### to pass as a reference ord in/out, create the nev to save it
  etrial <- environment(); 
  assign("indx0",vector( ,length=0), envir=etrial)
  ages <- substr(names(whocov[indx]),digit.age.begin,digit.age.end)
### f.age <- sort(ages)[1])
  f.age <- min(as.numeric(ages))
  f.age <- as.character(f.age)
  nd <- who.age.digits - nchar(f.age)  
  if(nd > 0 ){
   vec00 <- rep("0", nd)
   app <- paste("",vec00,collapse="", sep="")
   f.age <- paste(app, f.age, sep="")
  }
  #### better way is
  ages.sort <- sort(ages)
  f.age <- ages.sort[1]
  

  C.list <- sapply(indx, function(y,time.prior.param,etrial) { ##fun1(sapply) loops indx
### y points to column indeces; ncol = n.ind
  
    ny <- names(whocov[y])
    ny <- substr(ny,digit.age.begin, digit.age.end)
    y1 <- y %% n.age
    if( y1 == 0) y1 <- n.age
### for symmetry defines inds
    inds <- indx[1]:(y-1)
    lapply(indx, function(x,time.prior.param, etrial){  ##fun2(lapply) loops indx
### x points to row indeces; nrow= n.ind
### elements for every pair (x, y) are matrices
  
      nx <- names(whocov[x])
      nx <- substr(nx,digit.age.begin,digit.age.end)
      x1 <- x %% n.age
      if( x %% n.age == 0) x1 <- x1 + n.age 
### matrices indexing  dx; just counts element of whocov with cntry 
      dx <- x1 + n.ind * (y1 -1)
      indx0 <- get("indx0", envir=etrial, inherits =T)
### if symmetry, do not calculate matrices that belongs to inds
      if (is.element(x,inds) && ny != f.age){  ##if1 symmetry
        cxy <- NA
        assign("indx0",c(indx0, dx), envir=etrial)
      }else{  ##if1 symmetry
###  do not store matrices for W=0; just for the extra space    
        if (W[dx] == 0) {  ##if2 W=0
          assign("indx0",c(indx0, dx), envir=etrial)
          cxy <- NA  
        } else  {cxy <- C.matrix(whocov[[x]],whocov[[y]],time.prior.param)} ##if2 W != 0 
      } ## end if1 symmetry
### this part of code use for consistency checking of names
      wc <- colnames(W)[y1]
      wr <- rownames(W)[x1]
      wp <- paste(wr,wc,sep="")
###   if(wp != paste(nx,ny,sep="")) print("check names")
      names(cxy) <- paste(as.character(cntry),ny,nx, sep="")
      return(cxy)
    },time.prior.param,etrial) ## closing fun2
  },time.prior.param,etrial) ##closing fun1
  
   names(C.list) <- vnm
   indx0 <- get("indx0", envir=etrial, inherits=T)
   rm(etrial)
### check the right names
  if( any(names(C.list) != attr(C.list,"names"[1])))
    messout("Wrong naming", verbose)
### remove unnecessary matrices that are either repetition (for symmetry)
### or do not contribute because W==0
  if (length(indx0) > 0 ) C.list <- C.list[-indx0]
### do some cleaning
   rm(indx0)
### typing C.list is invisible;index C.list for matrix   
### you want to see results without indexing
  C.list <-  lapply(C.list,eval)
return(invisible(C.list))}



### utility function convert ages from int to char with two digits;
###  conv.char <- function(agx){
###    age.char <- as.character(agx)
###    sapply(age.char, function(x){
###                    if(nchar(x) <= 1) x <- paste("0", x, sep="")
###                     return(x)})}
### this function is not actually needed, instead
conv.char <- function(agx, dg=NULL,ebase=env.base){
             ebase <- get("env.base", envir=parent.frame())
             env.base <- ebase
          
             if (length(dg) <= 0)
               dg <- get("who.age.digits", envir=get("env.who", envir=ebase))
             return (formatC(agx,width=dg,format="d", flag="0"))}


## ************************************************************************
## ************************************************************************
##
## FUNCTION NAME:  cntry.Xprime.X 
##
## PLACE:    
##
## IMPORTED:    none

## DESCRIPTION: For any given country and a matrix of (possibly
##              weighted) covariates it returns a block diagonal matrix, with as
##              many blocks as number of age groups. The blocks on the diagonal are
##              the age specific correlation matrices of the covariates.
##              Matrix rows and columns are named with a combination of age
##              and covariate name (i.e. 45cnst).
##
## FORMAT: XX <- cntry.Xprime.X(cntry,Xmat)
##         
## INPUT:   cntry: scalar, a country code
##
##           Xmat: list, a cross-sectional time series of covariates,
##                 with cross-sections ranging over countries and age groups.
##                 It cannot contain missing values.
##        
## OUTPUT:  X: square diagonal block matrix. The number of blocks is
##             the number of age groups. If X_ca is the covariate matrix
##             for country cntry and age group a, the a-th block on the diagonal
##             is t(X_ca) %*% X_ca
##
## WRITTEN BY: Elena Villalon & Federico Girosi
##             evillalon@latte.harvard.edu;  fgirosi@latte.harvard.edu
##             CBRSS, Harvard University
## 
## Last modified: 07/26/2003
## 
## ************************************************************************
    
  
cntry.Xprime.X <- function(cntry,Xmat, ebase=env.base){
  ebase <- get("env.base", envir=parent.frame())
  env.base <- ebase
  ewho <- get("env.who", envir=ebase)
  who.digit.first  <-  get("who.digit.first", envir=ewho)
  who.cntry.digits <-  get("who.cntry.digits", envir=ewho)
  who.age.digits   <- get("who.age.digits", envir=ewho)
  who.year.digits <- get("who.year.digits", envir=ewho)
### Structure of dataset cstsid: 
  digit.cntry.begin <- who.digit.first + 1 
  digit.cntry.end   <- who.digit.first + who.cntry.digits
  digit.age.begin   <- digit.cntry.end + 1
  digit.age.end     <- digit.cntry.end + who.age.digits
  digit.year.begin  <- digit.age.end   + 1 
  digit.year.end    <- digit.age.end   + who.year.digits
  
   ctr  <- paste("^", cntry, sep="")
   indx <- grep(ctr, names(Xmat))
   xcntry   <- lapply(Xmat[indx], function(x) t(x) %*% x)
   each.col <- sapply(xcntry, function(x) nc <- ncol(x))
   s.col <- sum(each.col)
   each.row <- sapply(xcntry, function(x) nr <- nrow(x))
   s.row <- sum(each.row)
   nm.col <- sapply(1:length(indx), function(x) {
                            nm2 <- substr(names(xcntry[x])[1],digit.age.begin,digit.age.end)
                            nm2 <- paste(nm2,colnames(xcntry[[x]]),sep="") })
   nm.row <- sapply(1:length(indx), function(x) {
                            nm2 <- substr(names(xcntry[x])[1],digit.age.begin,digit.age.end) 
                            nm2 <- paste(nm2,rownames(xcntry[[x]]),sep="") })
   X <- matrix(0, nrow=s.row, ncol=s.col)
   rownames(X) <- unlist(nm.row)
   colnames(X) <- unlist(nm.col )
   count <- 0
    for(c in 1:length(each.col)){ ##for1 (cols)
            count <-  count + 1
            r <- c
            sc  <- ifelse(c > 1, sum(each.col[1:(c-1)]), 0)
            idc <- (sc + 1):(sc + each.col[c])
            sr  <- ifelse(r > 1, sum(each.row[1:(r-1)]), 0)
            idr <- (sr+1): (sr+ each.row[r])
            Xrc <-  xcntry[[count]]
            X[idr,idc] <-  Xrc} ##for1 (cols)
   return(X)
 }


## ************************************************************************
## ************************************************************************
##
## FUNCTION NAME:  cntry.Xprime.Y 
##
## PLACE:    
##
## IMPORTED:    none
##
## DESCRIPTION: it produces the vector v in eq. \ref{eq:qv} in the manual.
##              For given cntry code, a matrix is constructed with the product of
##              t(whoinsampy) %*% whoinsampx, for every age group 
##              Then, a matrix is build stucking all ages in one sincle column 
##
## FORMAT: v <- cntry.Xprime.Y(cntry,X,Y)
##         
## INPUT:   cntry: scalar, a country code
##
##          X: list, a cross-sectional time series of (possibly weighted) covariates
##
##          Y: list, a cross-sectional time series (possibly weighted)
##        
## OUTPUT:  v:  a block matrix with 1 column, with as many blocks as age groups.
##              if X_ca and Y_ca are the covariate matrix and the data for  country c
##              and age group a, the a-th block is t(X_ca)%*% Y_ca. The rows are named
##              using a combination of age group and covariate name (eg: 80.lngdp)
##              If Y_ca has missing values in certain years, those years are
##              removed from X_ca and Y_ca.
##
## WRITTEN BY: Elena Villalon & Federico Girosi
##             evillalon@latte.harvard.edu; fgirosi@latte.harvard.edu
##             CBRSS, Harvard University
## 
## Last modified: 07/26/2003
## 
## ************************************************************************
    
    
cntry.Xprime.Y <- function(cntry,X,Y, ebase=env.base){
   ebase <- get("env.base", envir=parent.frame())
   env.base <- ebase
   ewho <- get("env.who", envir=ebase)
   age.vec <- get("age.vec", envir=ewho)
   who.age.digits <- get("who.age.digits", envir=ewho)
   ctr  <- paste("^", cntry, sep="")
   indx <- grep(ctr, names(X))
   indy <- grep(ctr, names(Y))
   if(length(indx) != length(indy))  
      stop("build.Y: X(y) no match")
   whoyg <- Y[indy]
   whoxg <- X[indy]
   age.char <- formatC(age.vec, width=who.age.digits, format="d", flag="0")
   names(whoyg) <- age.char
   names(whoxg) <- age.char
### create etrial so that whoyg and whoxg are stored
   etrial <- environment()
   assign("whoyg", whoyg, envir=etrial)
   assign("whoxg", whoxg, envir=etrial)
### assuming that whoyg may have years with NA's; remove them
### for dth (whoyg) and covariates (whoxg) 
   isna <- lapply(1:length(indy), function(y, etrial){
                mat  <- Y[[indy[y]]]
                isna <- unique.default(row(mat)[is.na(mat)])
                if(length(isna) > 0){
                  whoyg <- get("whoyg", envir=etrial, inherits=T)
                  whoxg <- get("whoxg", envir=etrial, inherits=T)
                whoyg[[y]] <- whoyg[[y]][-isna,] 
                whoxg[[y]] <- whoxg[[y]][-isna,]} ##end ifisna
                assign("whoyg", whoyg, envir=etrial)
                assign("whoxg", whoxg, envir=etrial)
                bool <- nrow(whoxg[[y]]) != length(whoyg[[y]]) &&
                length(whoxg[[y]]) > 0 && length(whoyg[[y]]) > 0
                  bool <- na.omit(bool)
                if (length(bool) >0 && bool)
                    stop("build.Y:elimenation of na's not correct...data inconsistent", "\n")
                if (length(whoxg[[y]]) > 0 && length(whoyg[[y]]) > 0)
                  whoxg[[y]] <- data.frame(t(whoyg[[y]]) %*% whoxg[[y]])
                
                assign("whoxg", whoxg, envir=etrial)
                
                },etrial)
   yb <- get("whoxg", envir=etrial) 
   yb <- unlist(yb,recursive=T, use.names=T)
   vnm <- names(yb)
   yb <- as.matrix(yb)
### some cleaning; in other languages, like C, we would have passed
### whoxg, whoyg as references (or pointers),  in/out for Fortran 90
   rm(whoxg, whoyg,etrial)
   return(y=yb)
}

###
### A bunch of tests to check that the constructions of D, X, and Y are correct
### these tests are of no concern to others than software developers
###
### Driver test0 
all.test0.D <- function(wmat=NULL,ebase=env.base){
ebase <- get("env.base", envir=parent.frame())
env.base <- ebase
  ewho <- get("env.who", envir=ebase)
               cntry.vec <- get("cntry.vec", envir=ewho)
               n.cntry <- length(cntry.vec)
               who.C.age <- get("who.C.age", envir=ewho)
               age.vec <- get("age.vec", envir=ewho)
               for(i in 1:n.cntry){
 
               Di <- make.age.time.prior.matrix(cntry.vec[i],age.vec,wmat,who.C.age, env.base)
               test0.D(Di)}}
###    
test0.D <- function(D=NULL, tol=1.e-7){
eg <- eigen(D)$values
indx <- seq(along=eg)[eg <0]
eg0 <- abs(eg[indx])
###  .Machine$double.eps^0.5=  1.490116e-08
if (any(eg0 > tol ))
  cat("Test for eigenvalues of D fails", "\n")
else
  cat("Test Pass", "\n")
}
###
### Driver test1
### good only if wmat comes from 1st or 2nd derivatives; i.e.
### wmat <- t(M) %*% M, where M <- first(17) or M <- n.diff(2, 17)
### see last test for more details on how to set only 1st and 2nd
all.test1.D <- function(wmat=NULL, ebase= env.base){
              ebase <- get("env.base", envir=parent.frame())
              env.base <- ebase
              ewho <- get("env.who", envir=ebase)
               cntry.vec <- get("cntry.vec", envir=ewho)
               n.cntry <- length(cntry.vec)
                who.C.age <- get("who.C.age", envir=ewho)
               age.vec <- get("age.vec", envir=ewho)
               for(i in 1:n.cntry){
 
               Di <- make.age.time.prior.matrix(cntry.vec[i],age.vec,wmat,who.C.age, env.base)
             
               test1.D(Di, cntry.vec[i])}}

test1.D <- function(D=NULL, cntry=NULL,ebase=parent.frame()){
   ebase <- get("env.base", envir=parent.frame())
   env.base <- ebase
  ###identical(.GlobalEnv, environment())
  print(environment())
  ewho <- get("env.who", envir=ebase)
  cntry.vec <- get("cntry.vec", envir=ewho)
  age.vec <- get("age.vec", envir=ewho)
   who.C.age <- try(get("who.C.age", envir=ewho))
   if(class(who.C.age) == "try-error")
     who.C.age <- get("who.C.age", envir=ebase)
 if(class(who.C.age) == "try-error")
   {
     print("Missing who.C.age")
     return(list())
   }
   n.cntry <- length(cntry.vec)
  ctr <- paste("^", cntry, sep="")
  indx <- grep(ctr, names(who.C.age))
  age.char <- conv.char(age.vec)
### find indeces for submatrices correlating same age groups, i.e. 
### itage are indeces of who.C.age whose corresponding matrices
### are the diagonal sub-matrices of the block matrix
  itage <- sapply(age.char,function(x){
                    each <- paste(ctr, x, x,"$", sep="")
                    dx <- grep(each, names(who.C.age))
                    return(dx) })
  p1 <- runif(1,min=-5, max=5)
  betai <- lapply(itage, function(i) {nr <- nrow(who.C.age[[i]])
                                      beta <- c(rep(0, nr - 1), p1)})
  beT  <- matrix(unlist(betai))
 test <- abs(D %*% beT)
###  print(test)
  n.test <- length(test)
  z <- matrix(rep(0, n.test))
  ifelse(identical(all.equal.numeric(test,z), T),  
          print("Test for D (null space) matrix pass"),
          print("Test fail"))
  
  if (any(abs(test- z) > 1.e-8))
     print("Test fail")
  else
      print("Test for D (null space) matrix pass")
  return(test)
   }

###
### Driver for test2.C
all.test2.C <- function(ebase=env.base){
              ebase <- get("env.base", envir=parent.frame())
              env.base <- ebase
               ewho <- get("env.who", envir=ebase)
               cntry.vec <- get("cntry.vec", envir=ewho)
               n.cntry <- length(cntry.vec)
               print("test done with W matrix all 1; for all ages combo")
               for(i in 1:n.cntry){
               test2.C(cntry.vec[i])}
               print("if no messages previous to this you did well")}

test2.C <- function(cntry=NULL, ebase=parent.frame()){
  ebase <- get("env.base", envir=parent.frame())
  env.base <- ebase
  ewho <- get("env.who", envir=ebase)
  cntry.vec <- get("cntry.vec", envir=ewho)
  n.cntry <- length(cntry.vec)
  age.vec <- get("age.vec", envir=ewho)
  whocov <- get("whocov", envir=ewho)
  ctr <- paste("^", as.character(cntry),sep="")
### can only grep one value or element at a time 
  indx <- grep(ctr, names(whocov))
  n.ind <- length(indx)
  whocov[indx] <- lapply(whocov[indx], na.omit)
  Qcov <- lapply(whocov[indx], function(x) {
    res <- qr.Q(qr(x))
    if(ncol(res) >= ncol(x))
      res <- res[,1:ncol(x)]
    return(res)
                              })
  age.char <- conv.char(age.vec)
  paste.char <- function(x, y){ paste(x, y, sep="")}
  age.comb <- kronecker(age.char, age.char, FUN=paste.char )
### I need to store whocov and then redefine it with Qcov;
### cov.age.cor uses modified values of whocov; these changes 
### need to be global, otherwise cov.age.cor do not see them
  whostore <- whocov
  whocov[indx] <- Qcov
  assign("whocov", whocov, envir=ewho)
  wmat <- matrix(1, nrow=length(age.char), ncol=length(age.char))
  res <- cov.age.time.cor(cntry, wmat, age.char,age.comb)
  res <- lapply(res, function(x){
              nt <- min(ncol(x),nrow(x))
              v1 <- rep(x[1,1], nt)
              mm <- diag(v1)
              if(any(abs(mm- x) > 1.e-7))
                print("test2.D fails")
                
              pt <- ifelse(identical(all.equal.numeric(mm,x),T),0 ,
                     print("test2.D fails; identical")) 
             
              return(x)})
### go back to the good values
  whocov <- whostore
  assign("whocov", whocov, envir=ewho)
  return(invisible(res))
                            }
### Driver for test3.X
all.test3.X <- function(ebase=env.base){
              ebase <- get("env.base", envir=parent.frame())
              env.base <- ebase
               ewho <- get("env.who", envir=ebase)
               cntry.vec <- get("cntry.vec", envir=ewho)
               n.cntry <- length(cntry.vec)
               for(i in 1:n.cntry){
               test3.X(cntry.vec[i])}
              print("if no messages previous to this you did well")}

test3.X <- function(cntry=NULL, ebase=env.base){
   ebase <- get("env.base", envir=parent.frame())
              env.base <- ebase
  ewho <- get("env.who", envir=ebase)
  whoinsampx <- get("whoinsampx", envir=ewho)
  age.vec <- get("age.vec", envir=ewho)
 ctr <- paste("^", as.character(cntry),sep="")
### can only grep one value or element at a time 
  indx <- grep(ctr, names(whoinsampx))
  n.ind <- length(indx)
  Qcov <- whoinsampx[indx]
  Qcov <- lapply(Qcov, function(x) {
                              res <- qr.Q(qr(x))
                              if(ncol(res) > ncol(x))
                              res <- res[,1:ncol(x)]
                              return(res)
                              })
  age.char <- conv.char(age.vec)
### I'm just copying here the the code from build.X, but taking away 
### the name assignments for rows and cols, which do not change the testing
  xcntry   <- lapply(Qcov, function(x) t(x) %*% x)
  each.col <- sapply(xcntry, function(x) nc <- ncol(x))
  s.col <- sum(each.col)
  each.row <- sapply(xcntry, function(x) nr <- nrow(x))
  s.row <- sum(each.row)
  X <- matrix(0, nrow=s.row, ncol=s.col)
  count <- 0
    for(c in 1:length(each.col)){ ##for1 (cols)
            count <-  count + 1
            r <- c
            sc  <- ifelse(c > 1, sum(each.col[1:(c-1)]), 0)
            idc <- (sc + 1):(sc + each.col[c])
            sr  <- ifelse(r > 1, sum(each.row[1:(r-1)]), 0)
            idr <- (sr+1): (sr+ each.row[r])
            Xrc <-  xcntry[[count]]
            X[idr,idc] <-  Xrc
### added extra for testing
            mm  <- 0* Xrc
            mm[col(mm)==row(mm)] <- 1
            ifelse (identical(all.equal(Xrc,mm),T), 0,
                    print("identical:test fails for X"))
            if( any(abs(mm - Xrc) > 1.e-7))
              print("test fails for X")
            } ##for1 (cols)

   return(invisible(X))}
### Driver test4.D
### this test as test1.D, includes all cntrys and only works with
### 1st and 2nd derivatives, which are explicitly defined here
###
all.test4.D <- function(derivative=NULL, ebase=env.base){
   ebase <- get("env.base", envir=parent.frame())
              env.base <- ebase
                 ewho <- get("env.who", envir=ebase)
                 cntry.vec <- get("cntry.vec", envir=ewho)
                 who.C.age <- get("who.C.age", envir=ewho)
                 n.cntry <- length(cntry.vec)
                 if(derivative==2){
                 M  <- n.diff(2,17)
                 w2 <- t(M) %*% M
                 if(length(who.C.age)>0)  rm(who.C.age, inherits=T)
               ##  build.Caa(w2)
                 wmat <- w2
               }else if(derivative==1) {
                 M  <- first.diff(17)
                 w1 <- t(M) %*% M
                 if(length(who.C.age)>0)  rm(who.C.age, inherits=T)
         ###        build.Caa(w1)
                 wmat <- w1
               }else{
                 cat("derivative must be 1 or 2; you gave me:  ")
                 print(derivative)
                 stop("enter correct derivative")}
   age.vec <- get("age.vec", envir=ewho)
               for(i in 1:n.cntry){
                   Di <- make.age.time.prior.matrix(cntry.vec[i],age.vec,wmat,who.C.age, env.base)
           
               test4.D(Di, cntry.vec[i],derivative=derivative)}}

test4.D <- function(D=NULL, cntry=NULL,derivative=NULL, ebase=env.base){
   ebase <- get("env.base", envir=parent.frame())
              env.base <- ebase
  ewho <- get("env.who", envir=ebase)
  age.vec <- get("age.vec", envir=ewho)
  identical(.GlobalEnv, environment())
  print(environment())
  ctr  <- paste("^", cntry, sep="")
  indx <- grep(ctr, names(who.C.age))
  age.char <- conv.char(age.vec)
   who.C.age <- try(get("who.C.age", envir=ewho))
   if(class(who.C.age) == "try-error")
     who.C.age <- get("who.C.age", envir=ebase)
 if(class(who.C.age) == "try-error")
   {
     print("Missing who.C.age")
     return(list())
   }
### find indeces for submatrices correlating same age groups, i.e. 
### itage are indeces of who.C.age whose corresponding matrices
### are the diagonal sub-matrices of the block matrix
  itage <- sapply(age.char,function(x){
                    each <- paste(ctr, x, x,"$", sep="")
                    dx <- grep(each, names(who.C.age))
                    return(dx) })
  p1 <- runif(1,min=-1,max=1)
  p2 <- runif(1,min=-1,max=1)
  betai <- lapply(itage, function(i) {nr <- nrow(who.C.age[[i]])
                                      beta <- rep(0, nr - 2)
                                      a  <- as.numeric(substr(names(who.C.age[i]), 5,6) )
                                      a <- sort(unique.default(a))
                                      tm <- p1 * a + p2
                                      beta <- c(beta, tm,0)})
  beT  <- matrix(unlist(betai))
  test <- D %*% beT
###  print(test)
### to give an idea what is going on with deriv=1 and first & last ages
### they exactly match each other with opposite signs;
  if (derivative == 1) {
    itn <- length(itage)
    nr1 <- nrow(who.C.age[[(itage[1])]])
    nrn <- nrow(who.C.age[[(itage[itn])]])
    t1 <- 1:nr1
    tn <- (length(beT)-nrn + 1):length(beT)
    tt <- c(t1,tn)
    test1 <- test[t1]
    testn <- test[tn]
    if(any( abs(test1 + testn) > 1.e-7))
       print("first:derivative fails")
    else
      test <- matrix(test[-tt])}
###      
  z  <- matrix(rep(0, length(test))) ##second derivative
  if (identical(all.equal.numeric(test,z),T)){
    print(paste("derivative good: ", derivative))
  }else{
    print(paste("derivative FAILS order= ", derivative))}
  
  if(any(abs(test) > 1.e-8)){
    print(paste("derivative FAILS order= ", derivative))
  }else{
     print(paste("derivative good: ", derivative))
   }
      
    }

        
######################################################################


make.sigma <- function(param=list(sigma.bar=1,use.deaths=FALSE,average=TRUE,model.based=FALSE), ebase=env.base){
  ebase <- get("env.base", envir=parent.frame())
  env.base <- ebase
  ewho <- get("env.who", envir=ebase)
  sigma.bar <- param$sigma.bar;
  use.deaths <- param$use.deaths;
  average <- param$average;
  model.based <- param$model.based;
  ### We allow for some of the parameters to be NULL, so the user does
  ### not have to specify the entire list of parameters, but we need
  ### to convert them to the default so that they can appear in if statements
  ### without causing problems due to 0 length
  if (is.null(model.based)) model.based <- FALSE;
  if (is.null(use.deaths)) use.deaths <- FALSE;
  if (is.null(average)) average <- TRUE;
  if (is.null(sigma.bar)) sigma.bar <-  1;
   whoinsampy <- get("whoinsampy", envir=ewho)
   whopopul <- get("whopopul", envir=ewho)
  ### First some sanity checks
  ### We have to be clear about which model we use
  if (use.deaths == TRUE && model.based == TRUE)
    stop("Cannot have both use.deaths and model.based = TRUE in make.sigma");

  ### use.deaths == TRUE: this is the model in which we use the observed number of
  ### deaths, possibly averaged over time, to specify sigma ( a la Wilmoth)
  if (use.deaths == TRUE) {
    d <-  death.from.logmortality(whoinsampy,whopopul);
    if (average == TRUE){d <- cs.mean(d,extended=TRUE)};
    sigma <- lapply(d,FUN=function(x,sigma.bar){sigma.bar/sqrt(x)},sigma.bar);
    return(sigma);
  }
  ### model.based == TRUE: in this model the variance within an age group is inversely
  ### proportional to the the expected number of deaths plus a country-specific intercept.
  ### The expected number of deaths is modeled using an equation by equation OLS.
  if (model.based == TRUE) {
    message("Option model.base = true is not available, setting it to false") 
 ###   ols.result <- ols();
 ###   sigma <- make.model.based.sigma(ols.result);
 ###   return(sigma);
  }
  sigma <- lapply(whoinsampy,FUN=function(x,sigma.bar){rep(sigma.bar,length(x))},sigma.bar);  
  return(sigma);  
}



## ************************************************************************
## ************************************************************************
##
## FUNCTION NAME:  make.age.time.prior.matrix 
##
## PLACE:    
##
## GLOBALS: none
##
## DESCRIPTION: For a given cntry code, the correlation matrices among
##              all covariates for all age groups (i.e. who.C.age) are stuck together
##              Matrix rows are named with combo of age and cov (i.e. 45cnst);
##              and say for cols; thus each value is correlation among ages and covariates
##
## FORMAT: res <- make.age.prior.matrix(cntry,wmat,C.mat)
##         
## INPUT:   cntry: (scalar) a country code
##
##           wmat: (matrix) an A x A matrix of prior weights
##
##          C.mat: (list) a list of C matrices (only those corresponding to
##                  non-zero elements of wmat
##        
## OUTPUT:  res: ((matrix) a square symmetric matrix, of the type
##               defined in eq. \ref{eq:da} in the manual.
##
## WRITTEN BY: Elena Villalon & Federico Girosi
##             evillalon@latte.harvard.edu;  fgirosi@latte.harvard.edu
##             CBRSS, Harvard University
## 
## Last modified: 07/08/2003
## 
## ************************************************************************


### notice that this return Cji, not Cij as it may look
C.matrix <- function(Zi,Zj,time.prior.param){
  di <- dim(Zi)[1]
  dj <- dim(Zj)[1]
### when including the depvar or any of dth varibles, such as cvds,
### countries may start with data at year =1950 or a later, say year = 1951 
### which case we need to take the largest of (1951) to multiply matrices
  df <- di - dj
  if (df >= 1)
   Zi <- Zi[-(1:df), ]
  else if (df <= -1)
    Zj <- Zj[-(1:(-df)), ]
    
  tl <- nrow(Zi);
  time.der <- time.prior.param$time.der;
  time.weight <- time.prior.param$time.weight;
  W.time <- derivative.prior(tl,der.v=time.der,weight=time.weight);

  cxy <- t(Zj) %*% W.time %*% Zi;
  if (identical(Zj,Zi)){cxy <- 0.5*(cxy+t(cxy))}
  return(cxy)
}


## ************************************************************************
## ************************************************************************
##
## FUNCTION NAME:  make.age.time.prior.matrix 
##
## PLACE:    
##
## GLOBALS: none
##
## DESCRIPTION: For a given cntry code, the correlation matrices among
##              all covariates for all age groups (i.e. who.C.age) are stuck together
##              Matrix rows are named with combo of age and cov (i.e. 45cnst);
##              and say for cols; thus each value is correlation among ages and covariates
##              The name of the cntry is stored in first element of names(D)
##
## FORMAT: res <- make.age.prior.matrix(cntry,wmat,C.mat)
##         
## INPUT:   cntry: (scalar) a country code
##
##           wmat: (matrix) an A x A matrix of prior weights
##
##          C.mat: (list) a list of C matrices (only those corresponding to
##                  non-zero elements of wmat
##        
## OUTPUT:  res: ((matrix) a square symmetric matrix, of the type
##               defined in eq. \ref{eq:da} in the manual.
##
## WRITTEN BY: Elena Villalon & Federico Girosi
##             evillalon@latte.harvard.edu;  fgirosi@latte.harvard.edu
##             CBRSS, Harvard University
## 
## Last modified: 07/08/2003
## 
## ************************************************************************

make.age.time.prior.matrix <- function(cntry,age.vec=NULL, wmat,C.mat, ebase=env.base){
  ebase <- get("env.base", envir=parent.frame())
  env.base <- ebase
  ewho <- get("env.who", envir=ebase)
  if(length(age.vec) <= 0)
    age.vec <- get("age.vec", envir=ewho)
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
  
### assume that who.C.age is a global
### who.C.age <- build.Caa(wmat)
### select only one country from cntry.vec
    ctr  <- paste("^", cntry, sep="")
    indx <- grep(ctr, names(C.mat))
    age.char <- formatC(age.vec,width=who.age.digits,format="d", flag="0")
### find indeces of submatrices correlation for same age groups, i.e. 
### itage are indeces of C.mat whose corresponding matrices
### are the diagonal sub-matrices of the block matrix
    itage <- sapply(age.char,function(x){
                    each <- paste(ctr, x, x,"$", sep="")
                    
                   
                    dx <- grep(each, names(C.mat))
                   
                    return(dx) })
  itage <- unlist(itage)
### no.of cols of each block matrix at the diagonal
  
    each.col  <- sapply(itage, function(x){
                        
                        nc <- ncol(C.mat[[x]])
                      } )
    s.col  <- sum(each.col)
    nm.col <- sapply(itage, function(x) {
                            nm2 <- substr(names(C.mat[x])[1],digit.age.begin,digit.age.end)
                            nm2 <- paste(nm2,colnames(C.mat[[x]]),sep="") })
### no.of rows of each block matrix at the diagonal
    each.row  <- sapply(itage, function(x) nr <- nrow(C.mat[[x]]))
    s.row <- sum(each.row)
    nm.row <- sapply(itage, function(x){
                            nm1 <- substring(names(C.mat[x])[1],digit.age.begin + who.age.digits )
                            nm1 <- paste(nm1,rownames(C.mat[[x]]),sep="")})
### the block matrix 
    D <- matrix(0, nrow=s.row, ncol=s.col )
    rownames(D) <- unlist(nm.row)
    colnames(D)  <- unlist(nm.col)
### fill the block matrix with elements of C.mat for cntry
    count <- 0
    for(r in 1:length(each.row)){ ##for1 (cols)
      for(c in r:length(each.col)){  ##for2 (rows)
        if(wmat[r, c] != 0){ ##if1 (wmat != 0)
            count <-  count + 1
            sc  <- ifelse(c > 1, sum(each.col[1:(c-1)]), 0)
            idc <- (sc + 1):(sc + each.col[c])
            sr  <- ifelse(r > 1, sum(each.row[1:(r-1)]), 0)
            idr <- (sr+1): (sr+ each.row[r])
            Drc <-  C.mat[[indx[count]]]
            D[idr,idc] <-  wmat[r,c] * Drc
            D[idc,idr] <-  t( D[idr,idc] )
            
### just for the sake of finding if name rows & cols of D are correct
            nm1 <- substr(names(Drc)[1],digit.age.begin, digit.age.end)
            nm1 <- paste(nm1,rownames(Drc),sep="")
            nm2 <- substring(names(Drc)[1],digit.age.begin + who.age.digits )
            nm2 <- paste(nm2,colnames(Drc),sep="")
           if(any(nm2 != colnames(D[idr,idc])) ||  any(nm1 != rownames(D[idr,idc])))
              { ##print("make.block.matrix:names do not match")
               messout("make.block.matrix: col and row naming wrong", verbose)}              
         } ##end if1 (wmat != 0)
      } ## for2 (rows)
    } ## for1 (cols)
    if(any(D != t(D)))
      messout(paste("Complaint you did not build a symmetric D ", cntry, sep=""),verbose)
    if(count != length(C.mat[indx]))
       messout("make.block.matrix: counting wrong", verbose)
    names(D) <- paste("COUNTRY = ", cntry, sep="")
    return(invisible(D))}


##
## FUNCTION NAME:  omega.cntry.all, and build.adjacency 
##
## PLACE:    
##
## IMPORTED:    build.priors, 
##              
##
## DESCRIPTION: It builds the matrix of correlation among different countries. 
##              It uses build.priors to read the file "adjacency.txt"
##              with the correlation weight among countries, which can be 0, 1, 2. 
##              omega.cntry.all is square symmetric matrix whose dimension= 191x191,
##              which rows and cols are named with the country codes in order
##              of increasing value of the code, and the diagonal with 0.
##              build.adjacency is a subst of omega.cntry.all with only rows and cols
##              corresponding to countries in the group cntry.vec;
##              the dimension is length(cntry.vec) x length(cntry.vec)
##
## FORMAT: res <- omega.cntry.all(), build.adjacency(cntry.vec)
##         
## INPUT:  
##        
## OUTPUT:  Two Square, symmetric matrices of elements = 0,1,2; dimensions are
##          dim(omega.cntry.all) = 191 x 191, dim(build.adjacency) =length(cntry.vec)^2
##
## WRITTEN BY: Elena Villalon 
##             evillalon@latte.harvard.edu; 
##             CBRSS, Harvard University
## 
## Last modified: 01/07/2004
## 
## ************************************************************************
## ************************************************************************

  omega.cntry.all <- function(){
    ebase <- get("env.base", envir=parent.frame())
    env.base <- ebase
    verbose <- get("verbose", envir=ebase)
     ewho <- get("env.who", envir=ebase)
     Hct.c.deriv <- get("who.Hct.c.deriv", envir=ewho)
     c.vec  <- get("cntry.vec", envir=ewho)
  
###    print(dim(Hct.c.deriv))
     cntry.weight <- (build.priors(Hct.c.deriv,c.vec))$cntry.weight
     
     min.cntry <- min(cntry.weight[,1])
     max.cntry <- max(cntry.weight[,1])
     min.cntry2 <- min(cntry.weight[,2])
     max.cntry2 <- max(cntry.weight[,2])
     if (min.cntry != min.cntry2 || max.cntry != max.cntry2 )
       stop("Matrix omega.cntry.all  with all country correlations wrong")

     ord <- order(as.numeric(rownames(cntry.weight)))
     cntry.weight <- cntry.weight[ord, ]
     third <- cntry.weight[,3]
     l.cntry <- length(third)
     omega.cntry <- matrix(third, nrow=sqrt(l.cntry), byrow=T)
     nm <- unique.default(cntry.weight[,1])
     rownames(omega.cntry) <- nm 
     colnames(omega.cntry) <- nm
     if (!identical((omega.cntry),omega.cntry))
       messout("Wrong building omega.cntry", verbose)
     return(omega.cntry); 
    
  }
##
## We select a group of cntry's from the data set, and
## return the matrix with correlations among the selected group; 
## it is a different version of function build.adjacency.mat

build.adjacency <- function(cntry.vec=NULL, Wcntry=NULL,verbose=T){
 env.base <- try(get("env.base", envir=parent.frame()),silent=T)
 if(class(env.base)!="try-error")
   verbose <- get("verbose", envir=env.base)
  if (length(Wcntry) <= 0)
    Wcntry <- omega.cntry.all()
 
   
     if (length(cntry.vec) <= 0 && identical(.GlobalEnv, parent.frame())){
      messout("All countries in data set are included in correlation matrix.", verbose)
      return(Wcntry)
    }else if(length(cntry.vec) <= 0 ){
      ebase <- get("env.base", envir=parent.frame())
      env.base <- ebase
      ewho <- get("env.who", envir=ebase)
      cntry.vec <- get("cntry.vec", envir=ewho)}
 ind.row <- match(cntry.vec, rownames(Wcntry))
 ind.col <- match(cntry.vec, colnames(Wcntry))
 Wcntry <- Wcntry[ind.row, ind.col]
  return(Wcntry)}
##
##
## FUNCTION NAME:  cntry.lst.correlation 
##
## PLACE:    
##
## IMPORTED:    omega.cntry.all, build.adjacency
##              
##
## DESCRIPTION: It takes the matrices of correlation among cntry's either
##              the complete data set with omega.cntry.all or a subset of countries 
##              with build.adjacency and finds the group of countries that are 
##              correlated (values=1,2) with the input cntry, a code for country
##
## FORMAT: res <- cntry.lst.correlation(cntry.vec, cntry);
##         res <- cntry.lst.correlation(cntry);
##         
## INPUT:   cntry.vec ( group of selected cntrys) or NULL if all are included
##          cntry the country code to search correlation for
##        
## OUTPUT:  data frame one column only, named with the cntry code and   
##          as many elements of cntry's correlated to cntry.
##
## WRITTEN BY: Elena Villalon 
##             evillalon@latte.harvard.edu; 
##             CBRSS, Harvard University
## 
## Last modified: 01/07/2004
## 
## ************************************************************************
## ************************************************************************

cntry.lst.correlation <- function(cntry,cntry.vec=NULL,verbose=T){
  ebase <- try(get("env.base", envir=parent.frame()),silent=T)
  if(class(ebase)!="try-error")
    verbose <- get("verbose", envir=ebase)
  if(length(cntry.vec) <=0){
    m <- paste("Correlation for cntry-code=", cntry, "and all available cntrys in data set.", sep="")
    messout(m, verbose)
    Wcntry <- omega.cntry.all();
  }else if(length(cntry.vec) > 0){
    m <- paste("Correlations for cntry-code=", cntry, "and selected cntrys from data set.", sep="")
    messout(m,verbose)
    Wcntry <- build.adjacency(cntry.vec); 
  }
    Wcntry <- as.data.frame(Wcntry)
    Vcntry <- as.vector(Wcntry[as.character(cntry)])
    Vcntry <- Vcntry[Vcntry > 0]
    Vcntry <- data.frame(Vcntry)
    names(Vcntry) <- names(Wcntry[as.character(cntry)])
    return(Vcntry)}


    
## FUNCTION NAME:  build.C.cntry.time
##
## PLACE:    to be used with the cxx model
##
## IMPORTED:    cntry.vec, cov.cntry.age.cor, time.prior.param, 
##              the environments where the global variables are stored
##              laplacian.fom.adjacency(cntry.vec)
##
## DESCRIPTION: It produces correlation matrices among selected covariates for each age group
##              and all countries, from the list whocov that outputs make.mortality.data.
##              Matrices are products of covariates and are named with the codes
##              of correlated cntrys plus the age group.
##              Example "2450428045" for USA, France and age 45. 
##              The list elements of who.C.cntry are countries, which were selected at
##              the start of the simulation, and whose corresponding elements of omega.cntry
##              have either values 1, 2 contribute.  In addition because of symetry
##              the code of the first country is smaller than the value of the second.  
##              It calls cov.cntry.cor to construct for every contributing country the
##              correlations matrices among covariates, countries for every age group.  
##
## FORMAT: res <- build.C.cntry.time(omega.cntry, time.prior.param), which requires
##                cov.cntry.cor(agech, W.time,time.prior.param, cntry.char,cntry.comb)
##         
## INPUT:   whocov, cntry.vec, age.vec (globals),and omega.cntry;
##          plus other globals from environmnets
##        
## OUTPUT:  global, who.C.cntry,  list with  cntrys, for each pair of cntrys,
##          correlation among covariates for each age group.   
##
## WRITTEN BY: Elena Villalon
##             evillalon@latte.harvard.edu; 
##             CBRSS, Harvard University
## 
## Last modified: 01/08/2004
## 
## ************************************************************************
## ************************************************************************
    
build.C.cntry.time <- function(cntry.vec,time.prior.param, ebase=env.base) {
  
  ebase <- get("env.base", envir=parent.frame())
   env.base <- ebase
   ewho <- get("env.who", envir=ebase)
   age.vec <- get("age.vec", envir=ewho)
   cntry.vec <- get("cntry.vec", envir=ewho)
   who.age.digits <- get("who.age.digits", envir=ewho)
   who.cntry.digits <- get("who.cntry.digits", envir=ewho)
   verbose <- get("verbose", envir=ebase)
    omega.cntry <- build.adjacency(cntry.vec)
 
  if(identical(omega.cntry, 0*omega.cntry)){
    messout("Correlation no existent for this set of countries", verbose)
    messout("Setting who.Hct.sigma = NA", verbose)
    assign("who.Hct.sigma", NA, envir=get("env.who", envir=env.base))
    return (NULL)}
    
  ind.row <- match(cntry.vec, rownames(omega.cntry))
  ind.col <- match(cntry.vec, colnames(omega.cntry))
### the wrong matrix to use we need the laplaian, not the adjacency 
### omega.cntry <- (build.adjacency(cntry.vec))[ind.row, ind.col]
  W.cntry <- laplacian.from.adjacency(c.vec=cntry.vec, Wcntry=omega.cntry)
   if (nrow(W.cntry) != length(ind.row) || ncol(W.cntry) != length(ind.col))
     stop("Wrong building laplacian and adjacency matrices for cross-cntry correlation");
  n.age <- length(age.vec)
  n.cntry <- length(cntry.vec)
 
  age.char <- formatC(age.vec, width=who.age.digits, flag="0", format = "d")
  paste.char <- function(x, y){ paste(x, y, sep="")}
  age.comb <- kronecker(age.char, age.char, FUN=paste.char)
  cntry.vec <- sort(cntry.vec)
  cntry.char <- formatC(cntry.vec, width=who.cntry.digits, flag="0", format = "d")
  cntry.comb <- kronecker(cntry.char, cntry.char, FUN=paste.char) 
### check Wij to store only those elements of age.comb, for which Wij != 0
  who.C.cntry <- list()

  for(i in 1:n.age){
    Ci <- cov.cntry.age.cor(age.char[i],W.cntry,time.prior.param,cntry.char, cntry.comb, env.base)
   who.C.cntry <- c(who.C.cntry, Ci)}
   who.C.cntry <- who.C.cntry 
### just for consistency checking
 
  g1 <- grep(paste(age.char[1],"$",sep="") ,names(who.C.cntry))
  n.g1 <- length(g1)
   n.g <- sapply(1:length(age.vec), function(i){
    gi <- grep(paste(as.character(age.char[i]),"$",sep="") ,names(who.C.cntry))
    n.gi <- length(gi)
    return(n.gi)})

  who.C.cntry <- who.C.cntry
###  assign("who.C.cntry",who.C.cntry,envir=get("env.who",envir=env.base))
  if(any(n.g != n.g1))
      messout("Error in build.Caa, number of matrices for countries not agree",verbose)
  lst <- list(omega.cntry=omega.cntry, W.cntry=W.cntry, who.C.cntry= who.C.cntry)
  return(invisible(lst))
}

## FUNCTION NAME:  cov.cntry.time.cor; (executed with build.C.cntry.time)
##
## PLACE:    
##
## IMPORTED:    covariates matrices whocov from make.mortality.data()
##              environments with globals, including countries and age groups 
##
## DESCRIPTION: It produces a correlation matrix among contributing covariates for each age groups
##              and chosen countrys. Matrices are named with the code of correlated cntrys plus
##              the age group that (say "24504028035":USA+FRance+age=35), and are products of covariates.
##              Because of symmetry we assume first cntry code <= second cntry code;  
##              
## FORMAT: res <- cov.cntry.age.cor(agech, W,time.prior.param,cntry.char,cntry.comb, ebase)
##         
## INPUT:   one age group from age.vec; W for correlation among countries;
##          time.prior.param (time correlation); cntry.char contributing countries;
##          and cntry.comb for country combinations to be correlated;
#3          the environment, ebase, for globals
##        
## OUTPUT:  one element of list who.C.cntry,  which is a matrix of 
##          correlation among covariates and one contributing age,
##          for all countries.  
##
## WRITTEN BY: Elena Villalon  
##             evillalon@latte.harvard.edu; 
##             CBRSS, Harvard University
## 
## Last modified: 01/09/2004
## 
## ************************************************************************
## ************************************************************************

cov.cntry.age.cor <- function(agech, W,time.prior.param, cntry.char,cntry.comb, ebase=env.base){
  ebase <- get("env.base", envir=parent.frame())
  env.base <- ebase
  ewho <- get("env.who", envir=ebase)
  age.vec <- get("age.vec", envir=ewho)
  who.digit.first <- get("who.digit.first", envir=ewho)
  who.cntry.digits <- get("who.cntry.digits", envir=ewho)
  who.age.digits <- get("who.age.digits", envir=ewho)
  who.year.digits <- get("who.year.digits", envir=ewho)
  whocov <- get("whocov", envir=ewho)
  verbose <- get("verbose", envir=ebase)
### structure of dataset cstsid: 
  digit.cntry.begin <- who.digit.first + 1 
  digit.cntry.end   <- who.digit.first + who.cntry.digits
  digit.age.begin   <- digit.cntry.end + 1
  digit.age.end     <- digit.cntry.end + who.age.digits
  n.age <- length(age.vec)
  whocov <- cntry.to.age.order(whocov, who.age.digits, who.cntry.digits) 
 
  ache <- paste(agech,"$",sep="")
### can only grep one value or element at a time 
  indx <- grep(ache,names(whocov))
  n.ind <- length(indx)
  n.cntry <- length(cntry.char)
  vnm <- paste(cntry.comb,agech, sep="")
  whocov[indx] <- lapply(whocov[indx], na.omit)
  etrial <- new.env(TRUE,parent.frame())
  assign("indx0",vector( ,length=0), envir=etrial)
  cntrys <- substr(names(whocov[indx]),digit.cntry.begin,digit.cntry.end)
  f.cntry <- min(as.numeric(cntrys))
  f.cntry <- as.character(f.cntry)

  C.list <- sapply(indx, function(y,time.prior.param, etrial) { ##fun1(sapply) loops indx
### y points to column indeces; ncol = n.ind 
    ny <- names(whocov[y])
    ny <- substr(ny,digit.cntry.begin, digit.cntry.end)
    y1 <- y %% n.cntry 
    if( y1 == 0) y1 <- n.cntry 
### for symmetry defines inds; thus we only have
### matrices for combinations such as 24504080age = USA+France+age
### but not 40802450+ age= France+USA+age.  First cntry code <= second.
    inds <- indx[1]:(y - 1)
    lapply(indx, function(x,time.prior.param, etrial){  ##fun2(lapply) loops indx
### x points to row indeces; nrow= n.ind
### elements for every pair (x, y) are matrices
      nx <- names(whocov[x])
      nx <- substr(nx,digit.cntry.begin,digit.cntry.end)
      x1 <- x %% n.cntry 
      if( x1  == 0) x1 <- x1 + n.cntry 
### matrices indexing  dx; just counts element of whocov with cntry 
      dx <- x1 + n.ind * (y1 - 1)
### if symmetry, do not calculate matrices that belongs to inds
      if (is.element(x,inds) && ny != f.cntry){  ##if1 symmetry
        cxy <- NA
        indx0 <- get("indx0", envir=etrial)
        assign("indx0", c(indx0, dx), envir=etrial)
      }else{  ##if1 symmetry
###  do not store matrices for W=0; just for the extra space    
        if (W[dx] == 0) {  ##if2 W=0
            indx0 <- get("indx0", envir=etrial)
            assign("indx0", c(indx0, dx), envir=etrial)
            cxy <- NA  
        } else  {
    
          cxy <- C.matrix(whocov[[x]],whocov[[y]],time.prior.param)} ##if2 W != 0 
      } ## end if1 symmetry
### this part of code use for consistency checking of names
      wc <- colnames(W)[y1]
      wr <- rownames(W)[x1]
      wp <- paste(wr,wc,sep="")
###   if(wp != paste(nx,ny,sep="")) print("check names")
      names(cxy) <- paste(ny,nx,agech, sep="")
      return(cxy)
    },time.prior.param, etrial) ## closing fun2
  },time.prior.param, etrial) ##closing fun1 
  
   names(C.list) <- vnm 
### check the right names
  if( any(names(C.list) != attr(C.list,"names"[1])))
    messout("Wrong naming", verbose)
### remove unnecessary matrices that are either repetition (for symmetry)
### or do not contribute because W==0
  indx0 <- get("indx0", envir= etrial)
  if (length(indx0) > 0 ) C.list <- C.list[-indx0]
### We have no pointers or references to work with,do some cleaning
  rm(indx0)
### typing C.list is invisible;index C.list for matrix   
### you want to see results without indexing
  C.list <-  lapply(C.list,eval)
return(invisible(C.list))}
###
### helper function to reorder whocov according to age groups
### useful only for cross country smoothing
cntry.to.age.order <-  function(whocov,who.age.digits,who.cntry.digits){
### We need to rename elemnts of list whocov, so that first comes age group
### and second= cntry, i.e. we flip the country with age group around.
### For example instead of 245045=USA+45 age, become 452450=age 45 + USA
### Then we order elemnts of whocov according to increasing values of the names tags.
### Thus, first age group 0 with all cntry's; second age group=5 with all cntrys.
### Every age group comes first with all countries, and same age group
### has consecutive indeces to cover all group of countries.
  rindx <- (as.numeric(names(whocov))%%10^who.age.digits )*10^who.cntry.digits +
            as.numeric(names(whocov))%/%10^who.age.digits
  whocovr <- whocov
  names(whocovr) <- rindx
  ord <- order(rindx)
  whocovr <- whocovr[ord]
  names(whocovr) <- (as.numeric(names(whocovr))%%10^who.cntry.digits )*10^who.age.digits +
            as.numeric(names(whocovr))%/%10^who.cntry.digits
  whocov <-  whocovr
  return(whocov)}
#######################################################################################
## ************************************************************************
## ************************************************************************
##
## FUNCTION NAME:  make.age.time.prior.matrix 
##
## PLACE:    
##
## GLOBALS: cntry.vec, and others globals stored in the environmnets. 
##
## DESCRIPTION: For a given age group, the correlation matrices among
##              the country's that were selected for the simulation
##              (cntry.vec) and for all contributing covariates, 
##              (i.e. who.C.age) are stuck together to form matrix square D;
##              Some of the cntry's need to be eliminated from final list 
##              ctry.subset, because they are isolated and correlation matrix =0
##              the dimensions(D) is = number of covariates X length(cntry.subset)
##              Matrix rows are named with combo of cntry and cov (i.e. 2450cnst);
##              and same for cols; thus each value is the correlation
##              among different countries and contributing covariates,
##              in every age group.  
##
## FORMAT: res <- make.age.prior.matrix(age,cntry.vec,wmat,C.mat)
##         
## INPUT:   age: (scalar) for age group in char
##          cntry.vec = country vector with all countries
##          wmat: correlation matrix among countries as obtaiend from omega.cntry
##
##          C.mat: (list) a list of C matrices (only those corresponding to
##                  non-zero elements of wmat
##        
## OUTPUT:  res: ((matrix) a square symmetric matrix, of the type
##               defined in eq. \ref{eq:da} in the manual, whose
##               number of rows or cols = (no. of countries) X (no. covariates)
##
## WRITTEN BY: Elena Villalon
##             evillalon@latte.harvard.edu; 
##             CBRSS, Harvard University
## 
## Last modified: 01/12/2004
## 
## ************************************************************************

make.cntry.time.prior.matrix <- function(age,cntry.vec=NULL, wmat,C.mat, ebase=env.base){
  ebase <- get("env.base", envir=parent.frame())
  env.base <- ebase
  ewho <- get("env.who", envir=ebase)
  if(length(cntry.vec) <= 0)
    cntry.vec <- get("cntry.vec", envir=ewho)
  age.vec <- get("age.vec", envir=ewho)
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
  
### the countries in the data set that need to be correlated  
  cntry.vec <- sort(cntry.vec)
  assign("cntry.vec", cntry.vec, envir=ewho)
  n.cntry <- length(cntry.vec)
  cntry.ind <- 1:n.cntry
  names(cntry.ind) <- cntry.vec
  cntry.char <- formatC(cntry.vec,width=who.cntry.digits,format="d", flag="0")
### assume that who.C.cntry is a global
### who.C.cntry <- build.C.cntry.time(wmat)
### select only one age from age.vec and from who.C.cntry
  if(nchar(age) < 2 )
    age <- formatC(age, width=who.age.digits, flag="0", format = "d")  
  agech  <- paste(age,"$", sep="")
  indx <- grep(agech, names(C.mat))
    
### find indeces of submatrices correlation, 
### for first country of non-zero correlation with other cntrys. 
### itage are indeces of C.mat whose corresponding matrices
### are  sub-matrices of the block matrix; choose smallest
### cntry code with correlation for next cntry code.
### Submatrix of C.mat naming after code1, code2
### and age=code1code2age, where code1 < code2.
### We need two vectors that will stored indeces of submatrices.
### Vector vec contains those submatrices, which code2 is the largest for all
### the correlated cntrys of its group, and then it will always appear  
### second in the names of who.C.cntry after another code.
### Vector torem contains those countries that are isolated (no correlated
### to any of those in cntry.vec)
  
  e <- environment()
  vec <- vector(,length=0)
  torem <- vector(,length=0)  
  itage <- sapply(cntry.ind,function(x, cntry.char){
      for(i in x:n.cntry){
        each <- paste(cntry.char[x], cntry.char[i],agech, sep="")
        dx <- grep(each, names(C.mat))
        if (length(dx) > 0 )
          break;}
      if(length(dx) <= 0){
        for(i in 1:n.cntry){
        each <- paste(cntry.char[i], cntry.char[x],agech, sep="")
        dx <- grep(each, names(C.mat))
        if (length(dx) > 0 ){
           vec <- get("vec", envir=e)
           assign("vec", c(vec,x), envir=e)
          break;}
      }
      }
      if (length(dx) <= 0 ){
        torem <- get("torem", envir=e)
        assign("torem", c(torem,x), envir=e)
     ###   if (age == age.vec[1])
       ###   cat("Country", cntry.char[x],"no correlated with any of country selection.", "\n")
      }
        return(dx) },cntry.char)
  
### after unlist(itage), it will only contain cntry's which are correlated
### If some of the cntry's in cntry.vec are isolated, length(torem) >0,
### then redifining cntry.ind along itage contains only the non-isolated cntry's
### vec must also be shifted according to the values in torem
  
  itage <- unlist(itage)
  cntry.ind <- 1:length(itage)
  names(cntry.ind) <- names(itage)
  totorem <- torem
  if(length(torem) > 0){
    for(i in 1:length(torem)){
      r <- split(vec, vec > totorem[i])
      vec <- c(r$F, r$T - 1)
      totorem <- totorem -1 }
  }
     
### end of itage or indeces of correlated matrices.
### no.of rows of each block matrix 
  each.col <- sapply(cntry.ind, function(i){
    x <- itage[i]
    if(!is.element(i, vec))
      nc <- ncol(t(C.mat[[x]]))
    else
      nc <- ncol(C.mat[[x]])
  })
  s.col <- sum(each.col)
### no.of cols of each block matrix 
  each.row <- sapply(cntry.ind, function(i){
    x <- itage[i]
    if(!is.element(i, vec))
      nr <- nrow(C.mat[[x]])
    else
      nr <- nrow(t(C.mat[[x]]))})
  s.row <- sum(each.row)
  names.C <- names(C.mat[itage])
### these are the codes for the cntry's of cntry.vec that are correlated
### cntry.subset is a group within cntry.vec, which migth be all of them
  
  cntry.subset <- sapply(cntry.ind, function(x){
    xc <- names.C[x]
    if(identical(names.C[length(itage)], xc) || is.element(x,vec) ){
       nm1 <- substr(xc,1+who.cntry.digits, 2*who.cntry.digits)
       }else 
      nm1 <- substr(xc,1,who.cntry.digits)
      return(nm1)},simplify=T)
  
  uniq.cntry <- unique.default(cntry.subset)
  ind <- match(uniq.cntry, cntry.subset)
  if(length(ind) < length(cntry.subset) )
    messout("Error building D cntry smoothing: make.cntry.time.prior.matrix", verbose)
###  cat("For age group",age,", the cntry subset after eliminating isolated cntry's is", "\n")
###  cat (cntry.subset,"\n")
  assign("cntry.subset", cntry.subset,envir=ewho)
### for naming the rows of D 
    nm.row <- sapply(cntry.ind, function(x,C.mat,cntry.subset){
      if(!is.element(x, vec))
        nm1 <- paste(cntry.subset[x],rownames(C.mat[[itage[x]]]), sep="")
      else
         nm1 <- paste(cntry.subset[x],colnames(C.mat[[itage[x]]]), sep="")
      return(nm1)}, C.mat, cntry.subset)
    nm.col <- nm.row
  
### the block matrix, which needs to be filled
### with the submatrix of who.C.cntry  
    D <- matrix(0, nrow=s.row, ncol=s.col )
    rownames(D) <- unlist(nm.row)
    colnames(D)  <- unlist(nm.col)
### remove unwanted country with no correlation
### (or corresponding entry rows/cols all zeros)
### from correlation matrix wmat 
  if(length(torem) > 0){
    cntry.del <- as.character(cntry.vec[torem])
    row.to.del <- match(cntry.del, rownames(wmat))
    col.to.del <- match(cntry.del, colnames(wmat))
### the same as wmat <- wmat[-torem, -torem]
    wmat <- wmat[-row.to.del,-col.to.del]}
### fill the block matrix with elements of C.mat for cntry
    count <- 0
    for(r in 1:length(each.row)){ ##for1 (cols)
      for(c in r:length(each.col)){  ##for2 (rows)
       
        if(wmat[r, c] != 0){ ##if1 (wmat != 0)
            count <-  count + 1
            sc  <- ifelse(c > 1, sum(each.col[1:(c-1)]), 0)
            idc <- (sc + 1):(sc + each.col[c])
            sr  <- ifelse(r > 1, sum(each.row[1:(r-1)]), 0)
            idr <- (sr+1): (sr+ each.row[r])
            Drc <-  C.mat[[indx[count]]]
            D[idr,idc] <-  wmat[r,c] * Drc
            D[idc,idr] <-  t( D[idr,idc] )
            
         } ##end if1 (wmat != 0)
      } ## for2 (rows)
    } ## for1 (cols)
    if(any(D != t(D)))
      messout(paste("Matrix.block D is not symmetric for cross-cntry smoothing"), verbose)
    if(count != length(C.mat[indx]))
       messout("make.block.matrix, counting wrong", verbose)
### let us store the age group. 
    names(D) <- paste("age.group = ", age, sep="")
    names(D) <- na.omit(names(D))
   if(length(torem) > 0)
     isle.cntry  <- cntry.del
   else
     isle.cntry  <- NULL
    assign("isle.cntry", isle.cntry, envir=ewho)
    lst <- list(D = D, isle.cntry=isle.cntry)
    return(invisible(lst))}

rmv.cntry <- function(isle.cntry, Xmat){
  if (length(isle.cntry) > 0){
    isolated <- paste("^", isle.cntry, sep="")
    nmXmat <- names(Xmat)
    indX   <- lapply(isolated, function(x){
              ind <- grep(x,nmXmat)})
    indX <- unlist(indX)
  if (length(indX) > 0 )
    Xmat <- Xmat[-indX]}
 return(invisible(Xmat)) }
    
### I need this for ZZ that goes into the calculation of Hc.wc
### age is agechar; Xmat is whocov
### age version 
age.Xprime.X <- function(ag,Xmat, isle.cntry =NULL, ebase=env.base){
  ebase <- get("env.base", envir=parent.frame())
  env.base <- ebase
  ewho <- get("env.who", envir=ebase)
  who.digit.first  <-  get("who.digit.first", envir=ewho)
  who.cntry.digits <-  get("who.cntry.digits", envir=ewho)
  who.age.digits   <- get("who.age.digits", envir=ewho)
  who.year.digits  <- get("who.year.digits", envir=ewho)
  verbose <- get("verbose", envir=ebase)
  if(length(isle.cntry) <= 0)
    isle.cntry <- get("isle.cntry", envir=ewho)
  cntry.vec <- get("cntry.vec", envir=ewho)
   
### Structure of dataset cstsid: 
  digit.cntry.begin <- who.digit.first + 1 
  digit.cntry.end   <- who.digit.first + who.cntry.digits
  digit.age.begin   <- digit.cntry.end + 1
  digit.age.end     <- digit.cntry.end + who.age.digits
  digit.year.begin  <- digit.age.end   + 1 
  digit.year.end    <- digit.age.end   + who.year.digits
### re-order Xmat (=whocov) according to age groups so a given age follows
### with all the countries group together for each age group
  Xmat <- rmv.cntry(isle.cntry, Xmat)
###
  Xmat <- cntry.to.age.order(Xmat, who.age.digits, who.cntry.digits)
  agech <- formatC(ag, width=who.age.digits, format="d", flag="0")
  agech  <- paste(agech,"$", sep="")
   indx <- grep(agech, names(Xmat))
   xage   <- lapply(Xmat[indx], function(x) t(x) %*% x)
   each.col <- sapply(xage, function(x) nc <- ncol(x))
   s.col <- sum(each.col)
   each.row <- sapply(xage, function(x) nr <- nrow(x))
   s.row <- sum(each.row)
   nm.col <- sapply(1:length(indx), function(x) {
                            nm2 <- substr(names(xage[x])[1],digit.cntry.begin,digit.cntry.end)
                            nm2 <- paste(nm2,colnames(xage[[x]]),sep="") })
   nm.row <- sapply(1:length(indx), function(x) {
                            nm2 <- substr(names(xage[x])[1],digit.cntry.begin,digit.cntry.end) 
                            nm2 <- paste(nm2,rownames(xage[[x]]),sep="") })
   X <- matrix(0, nrow=s.row, ncol=s.col)
   rownames(X) <- unlist(nm.row)
   colnames(X) <- unlist(nm.col )
   count <- 0
    for(c in 1:length(each.col)){ ##for1 (cols)
            count <-  count + 1
            r <- c
            sc  <- ifelse(c > 1, sum(each.col[1:(c-1)]), 0)
            idc <- (sc + 1):(sc + each.col[c])
            sr  <- ifelse(r > 1, sum(each.row[1:(r-1)]), 0)
            idr <- (sr+1): (sr+ each.row[r])
            Xrc <-  xage[[count]]
            X[idr,idc] <-  Xrc} ##for1 (cols)
   return(X)
 }
## ************************************************************************
## ************************************************************************
##
## FUNCTION NAME:  age.Xprime.Y 
##
## PLACE:    
##
## IMPORTED:    none
##
## DESCRIPTION: it produces the vector v in eq. \ref{eq:qv} in the manual.
##              For given age group, a matrix is constructed with the product of
##              t(whoinsampy) %*% whoinsampx, for every country
##              Then, a matrix is build stucking all countries in one single column 
##
## FORMAT: v <- age.Xprime.Y(age,X,Y)
##         
## INPUT:   age: scalar, a country code
##
##          X: list, a cross-sectional time series of (possibly weighted) covariates
##
##          Y: list, a cross-sectional time series (possibly weighted)
##        
## OUTPUT:  v:  a block matrix with 1 column, with as many blocks as countries
##              if X_ca and Y_ca are the covariate matrix and the data for  country c
##              and age group a, the a-th block is t(X_ca)%*% Y_ca. The rows are named
##              using a combination of country and covariate name (eg: 2450.lngdp)
##              If Y_ca has missing values in certain years, those years are
##              removed from X_ca and Y_ca.
##
## WRITTEN BY: Elena Villalon & Federico Girosi
##             evillalon@latte.harvard.edu; fgirosi@latte.harvard.edu
##             CBRSS, Harvard University
## 
## Last modified: 07/26/2003
## 
## ************************************************************************
    
    
age.Xprime.Y <- function(ag,X,Y, isle.cntry=NULL,ebase=env.base){
   ebase <- get("env.base", envir=parent.frame())
   env.base <- ebase
   ewho <- get("env.who", envir=ebase)
   cntry.vec <- get("cntry.vec", envir=ewho)
  who.digit.first  <-  get("who.digit.first", envir=ewho)
  who.cntry.digits <-  get("who.cntry.digits", envir=ewho)
  who.age.digits   <- get("who.age.digits", envir=ewho)
  who.year.digits  <- get("who.year.digits", envir=ewho)
  ebase <- get("verbose", envir=ebase)
   if (length(isle.cntry) <= 0)
     isle.cntry <- get("isle.cntry", envir=ewho)
### Structure of dataset cstsid: 
  digit.cntry.begin <- who.digit.first + 1 
  digit.cntry.end   <- who.digit.first + who.cntry.digits
  digit.age.begin   <- digit.cntry.end + 1
  digit.age.end     <- digit.cntry.end + who.age.digits
  digit.year.begin  <- digit.age.end   + 1 
  digit.year.end    <- digit.age.end   + who.year.digits
### re-order X (=whoinsampx *weight) according to age groups so a given age follows
### with all the countries group together for each age group
   X <- rmv.cntry(isle.cntry=isle.cntry, Xmat = X)
   Y <- rmv.cntry(isle.cntry=isle.cntry, Xmat = Y)
   X <- cntry.to.age.order(X, who.age.digits, who.cntry.digits)
### re-order Y (=whoinsampy) according to age groups so a given age follows
### with all the countries group together for each age group
  Y <- cntry.to.age.order(Y, who.age.digits, who.cntry.digits) 
   agech <- formatC(ag, width=who.age.digits, format="d", flag="0")
   agech  <- paste(agech,"$", sep="")
   indx <- grep(agech, names(X))
   indy <- grep(agech, names(Y))
   if(length(indx) != length(indy))  
      stop("build.Y: X(y) no match")
   whoyg <- Y[indy]
   whoxg <- X[indy]
   indisle <- na.omit(match(isle.cntry, cntry.vec))
   cntry.related <- cntry.vec
   if(length(indisle) > 0)
     cntry.related <- cntry.vec[-indisle]
   cntry.char <- formatC(cntry.related, width=who.cntry.digits, format="d", flag="0")
   names(whoyg) <- cntry.char
   names(whoxg) <- cntry.char
### create etrial so that whoyg and whoxg are stored
   etrial <- environment()
   assign("whoyg", whoyg, envir=etrial)
   assign("whoxg", whoxg, envir=etrial)
### assuming that whoyg may have years with NA's; remove them
### for dth (whoyg) and covariates (whoxg) 
   isna <- lapply(1:length(indy), function(y, etrial){
                mat  <- Y[[indy[y]]]
                isna <- unique.default(row(mat)[is.na(mat)])
                if(length(isna) > 0){
                  whoyg <- get("whoyg", envir=etrial, inherits=T)
                  whoxg <- get("whoxg", envir=etrial, inherits=T)
                whoyg[[y]] <- whoyg[[y]][-isna,] 
                whoxg[[y]] <- whoxg[[y]][-isna,]} ##end ifisna
                assign("whoyg", whoyg, envir=etrial)
                assign("whoxg", whoxg, envir=etrial)
                bool <- nrow(whoxg[[y]]) != length(whoyg[[y]])&&
                length(whoxg[[y]]) > 0 && length(whoyg[[y]])
                bool <- na.omit(bool)
                if (length(bool) > 0 && bool)
                  stop("build.Y:elimenation of na's not correct...datamay be inconsistent","\n")
                 
                if(length(whoxg[[y]]) > 0 && length(whoyg[[y]]) > 0)
                    whoxg[[y]] <- data.frame(t(whoyg[[y]]) %*% whoxg[[y]])
                assign("whoxg", whoxg, envir=etrial)
                
                },etrial)
   yb <- get("whoxg", envir=etrial) 

   if(length(yb) > 0){
     yb <- unlist(yb,recursive=T, use.names=T)
     vnm <- names(yb)
     yb <- as.matrix(yb)
   }
### some cleaning; in other languages, like C, we would have passed
### whoxg, whoyg as references (or pointers),  in/out for Fortran 90
   rm(whoxg, whoyg,etrial)
   return(y=yb)
}

  
