### DESCRIPTION It takes a matrix of country correlations and optionally a country vector
###              The country correlation matrix is called proxac and has three columns
###              First and second columns are country codes,and third column 
###              is the correlation level (0, 1, 2), where 0 means no correlation
###              Since correlation is symmetric it builds a "symmetric" matrix
##               (with two colums for countries and one for level) and 
###              as many countries as in data set and if, for example, one row is
###              2450 2460 1 then the row 2460 2450 1 is included.
###              If cvec has counties not included in dataset then they are added
###              to matrix cntry.weight with correlations to all other countries and
###              third coulmn value equal to 0.
###
### OUTPUT       Returns two matrices of three columns each (code code level);
###              cntry.cor contains those countries that are correlated, i.e. level >0
###              cntry.weight contains all countries even if the correlation is 0
###
### AUTHOR Elena Villalon
###        IQS, Harvard Univ
###        evillalon@iq.harvard.edu
###
########################################################################################


build.priors <- function(proxac, cvec=NULL){
### given a matrix of correlations proxac.
### Assume that only entries that are different from 0,
### are given and that they are indicated with first column cntry,
### interacting with second column cntry and third column being
### the level of interaction
  
###     print(dim(proxac))   
    ind <- as.list(1:length(proxac[,1]))
    rnm <- sapply(ind, function(n) 
                  paste(as.character(proxac[n,1]), as.character(proxac[n,2]), sep=""))
    rownames(proxac) <- rnm
 
    lst <- expand.matrix(proxac)
   
    cntry.cor <- lst$mat
    proxac <- lst$proxac
    cntry.vec  <- sort(proxac[, 1])
    ucntry.vec <- unique.default(cntry.vec)
    lncntry <- length(ucntry.vec)
   
    cntry.weight <- matrix(0, nrow=lncntry*lncntry, ncol=3)                        
    allcntry <- sapply(ucntry.vec, paste, ucntry.vec, sep="")
    rownames(cntry.weight) <- allcntry 
    nn <- round(nchar(allcntry[1])/2)
    cntry.weight[,1] <- as.numeric(substring(allcntry,1, nn))
    cntry.weight[,2] <- as.numeric(substring(allcntry,nn+1))
    ind <- match(rownames(cntry.cor), rownames(cntry.weight))
    cntry.weight[ind,3] <- cntry.cor[,3]
  
  if(length(cvec) > 0)
   
      cntry.weight <- check.expansion(cntry.weight, cvec)
   
##asume we are given a cntry vec to extract countries
  c.weight <- NULL
  c.cor <- NULL

  
   if(length(c.weight) > 0 && length(c.cor) > 0)
     lst <- list(cntry.weight=c.weight, cntry.cor = c.cor)
  else
    lst <- list(cntry.weight=cntry.weight, cntry.cor = cntry.cor)
  
  return(invisible(lst))
    }
### DESCRIPTION: It takes a matrix proxac and includes all reciprocal
###              correlations between countroies, say if 2450 2090 1 then
##               includes 2090 2450 1 if not in dataset. Eliminates
##               duplicates and order the matix, first by first column
##               and second given first column by second column
##
## OUTPUT The input matrix proxac moidified in the function, and
##        the matrix of correlations
##
## AUTHOR Elena Villalon
##        evillalon@iq.harvard.edu
##        11/31/2005
########################################################################
 
expand.matrix <- function(proxac){
  tmp <- proxac
  tmp[,1] <- proxac[, 2]
  tmp[,2] <- proxac[,1]
    
  ind <- as.list(1:length(proxac[,1]))
  rnm <- sapply(ind, function(n) 
                paste(as.character(tmp[n,1]), as.character(tmp[n,2]), sep=""))
  rownames(tmp) <- rnm
### contains all cross correlation, say 5198 5107 1,
### then we added 5107 5198 1
  proxac <- rbind(proxac, tmp)
  ord <- order(as.numeric(proxac[, 1]))
  proxac <- proxac[ord, ]  
  ccorr <- split.data.frame(proxac, proxac[,1])
### given a list element order by second column
### also eliminate duplicates  
  ccorr <- lapply(ccorr, function(mat){
    rname <- rownames(mat)
    urname <- unique.default(rname)
    if(length(rname) > length(urname)){
      ind <- match(urname, rname)
      mat <- mat[ind, ]
    }
    ord <- order(as.numeric(mat[,2]))
    mat <- mat[ord, ]
    return(mat)})
  
  ind <- as.list(1:length(ccorr))
  names(ind) <- names(ccorr)
   
  exp.mat <- NULL
  ev <- environment()
  
  ##cntry.cor becomes a pointer
  ccorr <- lapply(ccorr, function(mat, ev){
    exp.mat <- get("exp.mat", envir=ev)
    exp.mat <- rbind(exp.mat,mat)
    assign("exp.mat", exp.mat, envir=ev)
    return(mat)}, ev)
  lst <- list(mat=exp.mat, proxac=proxac)
  return(lst)
}
### DESCRIPTION Takes the cntry.weight matrix of three columns: code, code, level
###             and a cntry.vec of country codes, cvec.  It checks if the codes
###             in cvec are included in the matrix cntry.weight.
###             If they are not, it adds them and correlates them with all
###             countries in data set of cntry.weigth and among themselves with
###             correlation level 0.
###
### OUTPUT The matrix of cntry.weight with all the countries and their correlation
###        The matrix is symmetric col1 col2 level and col2 col1 level for all
###        codes and all levels of correlation.
###
### AUTHOR Elena Villalon
###        IQSS, Harvard Univ
###        evillalon@iq.harvard.edu
###        11/1/20005
##############################################################
check.expansion <- function(cntry.weight, cvec)
  {
    ind <- match(cvec, cntry.weight[,1])
    ind <- is.na(ind)
###    print(ind)
    if(all(ind ==F))
      return(cntry.weight)
   
    no.include <- cvec[ind]
    ucntry <- unique.default(cntry.weight[,1])
    ucntry <- c(ucntry, no.include)
    nn <- length(ucntry) * length(no.include)
    matext <- matrix(0, nrow= nn, ncol=3)
    ccvec <- sort(rep(no.include, length(ucntry)))
    recntry <- rep(ucntry, length(no.include))
    matext[, 1]  <- ccvec
    matext[, 2]  <- recntry
    nm <- sapply(as.list(1:length(matext[,1])), function(n){
      f <- matext[n,1]
      s <- matext[n,2]
      nm <- paste(f,s,sep="")
    })
        rownames(matext) <- nm
        cntry.weight <- rbind(cntry.weight,matext)
        cntry.weight <- expand.matrix(cntry.weight)$mat
        return(cntry.weight)
      }
##########################################################################################
##########################################################################################
##
## FUNCTION NAME:      build.adjacency.mat 
##
## PLACE:    
##
## IMPORTED:    build.priors()
##
## DESCRIPTION: from weight matrix, cntry.weight, builds a submatrix with  
##              only cntry codes in c.vec and no.rows=no.cols= length(c.vec)
##              matrix elements from correlation col of cntry.weight(0, 1, 2)
##              There is another version of this function called simply build.adjacency
##
## FORMAT: res <- buil.adjacency.mat(c.vec)
##         res <- build.adjacency.mat()
##
## INPUT:   cntry.weight (matrix no.row=191x191, no.col=3) and cntry.vec (global)
##        
## OUTPUT:   adjacency matrix s.cntry = cntry.vec.weight, symmetric square matrix
##           no.row=length(c.vec) elements (0,1,2),depend on correlation 
##
## WRITTEN BY: Elena Villalon & Federico Girosi 
##             fgirosi@latte.harvard.edu, evillalon@latte.harvard.edu
##             CBRSS, Harvard University
## 
## Last modified: 06/10/2003
## 
##########################################################################################
##########################################################################################
 

build.adjacency.mat <- function( c.vec=NULL, ebase){
### given a set of cntry's from whousercntrylist or cntry.vec  
### search for  cntry.vec in data cntry.weight and get submatrix
  ebase <- get("env.base", envir=parent.frame());
  env.base <- ebase; 
ewho <- get("env.who", envir=ebase)
if (length(c.vec) <= 0)
  c.vec <- get("cntry.vec", envir=ewho)

c.weight <-build.priors(proxac=get("Hct.c.deriv", envir=ewho), c.vec)$cntry.weight
ind1 <- is.element(c.weight[,1],c.vec)
ind2 <- is.element(c.weight[,2],c.vec)
ind  <- ind1 & ind2
### cntry.wec with all countries from c.vec;
### no.row = n.cntry x n.cntry (n.cntry = length(cntry.vec)),
### no.col = 3; (3rd col = 0, 1, 2) 
cntry.wec <- c.weight[ind, ]
### cntry.vec.weight is square matrix with no.row= n.cntry
### elements of matrix the correlation among cntry's 
cntry.vec.weight <- matrix(cntry.wec[,3], nrow=length(c.vec), byrow =T)
rownames(cntry.vec.weight) <- c.vec
colnames(cntry.vec.weight) <- c.vec
if( any(cntry.vec.weight != t(cntry.vec.weight)))
  stop("matrix adjacency is not symmetric")
return(invisible(cntry.vec.weight))}

##########################################################################################
##########################################################################################
##
## FUNCTION NAME:      laplacian.from.adjacency 
##
## PLACE:    
##
## IMPORTED:    build.adjacency.mat(c.vec)
##
## DESCRIPTION: takes adjacency, s.cntry, find vector of sum of cols for every row  
##              build splus (equal size of s.cntry) with sum of cols in diagonal
##              matrix laplacian = splus - s.cntry 
##
## FORMAT: s.cntry <- laplacian.from.adjacency(c.vec)
##         s.cntry <- laplacian.from.adjacency()
##
## INPUT:   adjacency, s.cntry, square matrix and cntry.vec (global)
##        
## OUTPUT:   matrix laplacian, W.cntry=lapl.cntry, same size as s.cntry
##           no.rows= no.cols= length(cntry.vec); elements <0 and >0 
##
## WRITTEN BY: Elena Villalon & Federico Girosi 
##             fgirosi@latte.harvard.edu, evillalon@latte.harvard.edu
##             CBRSS, Harvard University
## 
## Last modified: 06/9/2003
## 
##########################################################################################
##########################################################################################
 
laplacian.from.adjacency <- function(c.vec=NULL,Wcntry=NULL){
### calculate sum of cols elements for every row
### turn vector sum.col of length=nrow(s.cntry)
  ebase <- get("env.base", envir=parent.frame())
  env.base <- ebase
  ewho <- get("env.who", envir=env.base)
  if(length(c.vec) <= 0)
    c.vec <- get("cntry.vec", envir=ewho)
  if(length(Wcntry) <= 0)
    s.cntry <- build.adjacency.mat( c.vec)
  else
    s.cntry <- Wcntry
sum.col <- rowsum(t(s.cntry), rep(1,ncol(s.cntry)))
splus <- 0 * s.cntry
### fill diagonal with sum.col 
splus[col(splus) == row(splus)] <- sum.col
### build Laplacian
lapl.cntry <- splus - s.cntry
### check if
sum.lapl <- rowsum(t(lapl.cntry), rep(1,ncol(lapl.cntry)))
if (any(sum.lapl != 0)) stop("Building wrong laplacian & adjancy")
if( any(lapl.cntry != t(lapl.cntry)))
  stop("matrix laplacian is not symmetric")
return(invisible(lapl.cntry))}

##########################################################################################
##########################################################################################
##
## FUNCTION NAME:      adjacency.from.laplacian 
##
## PLACE:    
##
## IMPORTED:    none
##
## DESCRIPTION: from square, symmetric matrix wcntry,find diagonals elements and   
##              assigned them to the diagonal matrix splus, then define
##              matrix adjacency = splus - wcntry 
##
## FORMAT: res <- adjacency.from.laplacian(wcntry)
##
## INPUT:   wcntry, square matrix 
##        
## OUTPUT:   adjacency, s.cntry = adj.cntry same size as wcntry
##
## WRITTEN BY: Elena Villalon & Federico Girosi 
##             fgirosi@latte.harvard.edu, evillalon@latte.harvard.edu
##             CBRSS, Harvard University
## 
## Last modified: 06/9/2003
## 
##########################################################################################
##########################################################################################
 
adjacency.from.laplacian <- function(wcntry) {
### wcntry must be symmetric:
  if( any(wcntry != t(wcntry)))
  stop("matrix input laplacian not symmetric")
  
### find the diagonal elements of matrix w.cntry
splus <- 0 * wcntry
diag.wcntry <- wcntry[col(wcntry) == row(wcntry)]
splus[col(splus) == row(splus)] <- diag.wcntry
### build adjacency: square matrix with nrow = n.cntry
### and 0 for all diagonal elements
adj.cntry <- splus - wcntry
### check
diag.adj.cntry <- adj.cntry[col(adj.cntry) == row(adj.cntry)]
if(any(abs(diag.adj.cntry) > 1.e-10))
  stop("building wrong adjacency.from.laplacian")
return(invisible(adj.cntry)) }

##########################################################################################
##########################################################################################
## FUNCTION  NAME:      find.correlation  
##
## PLACE:    
##
## IMPORTED:    build.priors()
##
## DESCRIPTION: from correlation matrix, cntry.cor, finds rows for cntry.code 
##              obtain neighbors or correlation  (3rd col =1 or 2) for cntry.code
##
## FORMAT: res <- find.correlation(cntry.code)
##
## INPUT:   cntry.cor (matrix with 191 cntrys, no.col=3), and cntry.code 
##        
## OUTPUT:   cntry.code.cor, 1 col matrix with no.rows those 
##           countries that are neighbours or correlated with cntry.code 
##
## WRITTEN BY: Elena Villalon & Federico Girosi 
##             fgirosi@latte.harvard.edu, evillalon@latte.harvard.edu
##             CBRSS, Harvard University
## 
## Last modified: 06/9/2003
## 
##########################################################################################
##########################################################################################
 
find.correlation <- function(cntry.code, ebase){
### finding correlated cntry's (3rd col=1,2) of cntry.code
### given cntry.code find its correlated neighbors
  ewho <- get("env.who", envir=ebase)
  who.cntry.digits <- get("who.cntry.digits", envir=ewho)
  c.vec <- get("cntry.vec", envir=ewho)
 
  cntry.cor <- build.priors(proxac=get("Hct.c.deriv", envir=ewho), c.vec)$cntry.cor
c1 <- paste("^", as.character(cntry.code), sep="")
### ind1 for cntry.code at beginning of string 
ind1 <- grep(c1, rownames(cntry.cor)) 
c2   <-  paste(as.character(cntry.code),"$", sep="")
### ind2 for cntry.code at end of string 
ind2 <- grep(c2, rownames(cntry.cor)) 
indx.rows <- sort(c(ind1, ind2))
### build subset of cntry.cor, which 
### contains appearances of cntry.code either col 1 or 2 and any row 
cntry.code.mat <- cntry.cor[indx.rows, ]
### because of symmetry counts correlations once 
cntry.code.cor <- cntry.cor[ind1,3]
cntry.code.cor <- matrix(cntry.code.cor)
vr <- rownames(cntry.cor)[ind1]
rownames(cntry.code.cor) <- substr(vr, who.cntry.digits +1, 2 *who.cntry.digits)
colnames(cntry.code.cor) <- cntry.code
### cntry.code.r with ncols = 1 and nrows = cntry's correlated 
return(invisible(cntry.code.cor))}

##########################################################################################
##########################################################################################
##
## FUNCTION NAME:     qr.rank
##
## DESCRIPTION: calculate qr descomposition of square symmetric matrix wcntry,
##              
## FORMAT: res <- qr.rank(wcntry)
##          res <- qr.rank(wcntry, eps)
##                
## INPUT:   wcntry, square symmetric matrix,
##          eps (or tolerance); may use default tol=1.e-7  
##        
## OUTPUT:  rank of the matrix wcntry
##
## WRITTEN BY: Elena Villalon  
##             evillalon@latte.harvard.edu
##             CBRSS, Harvard University
## 
## Last modified: 06/20/2003
## 
##########################################################################################
##########################################################################################
qr.rank <- function(wcntry, eps=NA){
  if(is.na(eps))
    return(qr(wcntry)$rank)
  else
    return(qr(wcntry,tol=eps)$rank)}

##########################################################################################
##########################################################################################
##
## FUNCTION NAME:     find.zero.eigen 
##
## PLACE:   this is an alternative to R-function qr(x)$rank 
##
## IMPORTED:    none
##
## DESCRIPTION: from square symmetric matrix wcntry, calculate eigenvalues   
##              (and possibly eigenvectors for nu & nv > 0) with La.svd(x)  
##              1)standardize eigenvalues and find those which are zeros; 
##              2)calculate log10(eigenvalues) and find first derivative 
##      
##    
## FORMAT: res <- find.zero.eigen(wcntry,only.values,eps)
##                find.zero.eigen(wcntry)
##
## INPUT:   wcntry, square symmetric matrix;
##          only.values=T, only eigenvalues, and F for eigenvectors
##          eps < 1 (small parameter) to find identical eigenvalues; 
##          eps=NA use default 
##        
## OUTPUT:  list of four: w.eigen,w.eigv, w.rank and no.zeros
##          1) & 2) eigenvalues and eigenvectors (for only.values=F)
##          3) rank of input matrix wcntry 
##          2) n.zeros length of zeros eigenvalues of matrix wcntry 
##
## WRITTEN BY: Elena Villalon  
##             evillalon@latte.harvard.edu
##             CBRSS, Harvard University
## 
## Last modified: 06/23/2003
## 
##########################################################################################
##########################################################################################
find.zero.eigen <- function(wcntry,only.values=T, eps= .Machine$double.eps){
### wcntry must be symmetric:
  if( any(wcntry != t(wcntry)) )
    stop("find.zero.eigen: incorrect input parameters")
### La.svd returns a list of three elements:
### where La.svd()$d (eigenvalues descending order)
### and two other matrices with eigenvectors, which are NULL if nu=nv = 0
  if(only.values==T){
    if(class(try( La.svd(wcntry,nu=0, nv=0)))=="try-error")
          res <- svd(wcntry,nu=0, nv=0)
    else
         res <- La.svd(wcntry,nu=0, nv=0)
       
    eig.w <- res$d
    eig.v <- NA
  }else{
    
###    res <- La.svd(wcntry,nv=0) 
###     eig.w <- res$d
###     eig.v <- res$u}
       res <- eigen(wcntry,symmetric=TRUE)
       eig.w <- res$values;
       eig.w[eig.w <= 0] <-  .Machine$double.eps;
       eig.v <- res$vectors;
     }
       
  n.eig <- length(eig.w)
### standardize eigenvalues, find mean, std and minimum values  
  centered.z <- scale(eig.w, center=T, scale=T)
  mean.z <- attr(centered.z,"scaled:center")
  std.z  <- attr(centered.z,"scaled:scale")
  norm.z <- - mean.z /std.z
  zmin   <- min(centered.z)
  zmax   <- max(centered.z)
  seig.w <- centered.z - zmin
  
  if(!identical(all.equal(zmin, norm.z), T)){ ##if1: there exists no 0 eigenvalue  
     rank <- length(centered.z)
     rw <- rank;
     n.zeros <- 0
  }else{  ##if1: unless one 0 eigenvalue exists 
     ind0 <- seq(along=seig.w)[seig.w > 0]
     ind0.l <- ind0[length(ind0)]
### calculate log of eigenvalues
     leig.w <- log10(eig.w)
     deig.w <- diff(leig.w, 1,1)
     indl <- seq(along=deig.w)[deig.w == min(deig.w)]
     z.jump <- indl[length(indl)]
     if(ind0.l > z.jump && !is.na(z.jump) && !is.na(ind0.l)){ ##if2: jump 
       rw <- z.jump 
     }else{  ##if2:jump
         if(is.na(eps)){ ##if3: eps
           rw <- ind0.l
         }else{  ##if3: eps
           inds <- seq(along=seig.w)[seig.w > eps]
           inds.l <- inds[length(inds)]
           del <- (ind0.l - inds.l)/ 2 
           rw <- inds.l + trunc(del) 
         
          } ##if3: eps
       }  ##if2:jump
         n.zeros <- n.eig - rw
   } ## if1:  !identical(all.equal(zmin, norm.z)) 
   return(list(w.eigv= eig.v, w.eigen = eig.w, w.rank= rw, n.zeros= n.zeros)) }

##########################################################################################
##########################################################################################
## FUNCTION NAME:    sample.improper.normal 
##
## PLACE:    
##
## IMPORTED:    find.zero.eigen 
##
## DESCRIPTION: after the gauss code written by Federico
##              find range of input matrix, eigenvalues (possibly some = 0)
##              and eigenvectors; using svd descomposition to invert matrix
##              for improper prior, but if no null space used choleski descomposition  
##    
## FORMAT: res <- sample.improper.normal(wmat,theta, n.sample, dim.null=NA)
##          
## INPUT:   wmat (prior) square symmetric matrix;
##          either positive-(or semi)definite (improper) or proper
##          theta > 0; n.sample integer for random number generators
##          dim.null (dimension null space or zero eigenvalues)
##        
## OUTPUT:  list of two: matrix x, which is the inverted of wmat (prior)
##          and covar the covariance matrix  for improper (or proper) prior
##
## WRITTEN BY: Elena Villalon & Federico Girosi 
##             evillalon@latte.harvard.edu
##             CBRSS, Harvard University
## 
## Last modified: 06/16/2003
## 
##########################################################################################
##########################################################################################
 
sample.improper.normal <- function(wmat, theta, n.sample, dim.null=NA,verbose=T) {
  ebase <- try(get("env.base", envir=parent.frame()),silent=T)
  if(class(ebase)!="try-error")
    verbose <- get("verbose", envir=ebase)
  
  dim <- nrow(wmat)
  if (theta <= 0)
    stop("sample.improper.normal:incorrect input parameters")
  reg <- find.zero.eigen(wmat, only.values=F)
  if(is.na(dim.null)){ ## if1 (dim.null)
    dim.null <- reg$n.zeros 
    rw <- reg$w.rank
  }else
    rw <- dim - dim.null ## end if1 (dim.null)
  
  if ( dim.null >= dim || rw <= 0)
      stop("sample.improper.normal:dimension null space >= no of rows")
### matrix random nos with dim= rw(rows) x n.sample(cols)
  c <- rnorm(rw * n.sample, mean=0, sd=1)
  c.mat <- matrix(c, nrow=rw)   
  if (dim.null > 0 || any(diag(wmat) == 0)) { ## if2 (dim.null)
### eigenvalues are returned in descending order
###    print("SVD for symmetric,semi-positive definite matrix")
    ind <- 1:rw
### reg <- La.eigen(wmat, symmetric=T, only.values=F)
    lambda <- (reg$w.eigen)[ind]
    R <- reg$w.eigv
    R.orth <- R[,ind]
### diagonal matrix (rw x rw) with values 1/lambda
    if(length(lambda) > 1)
      L.inv <- diag(1/lambda)
    else
      L.inv <- matrix(1/lambda,nrow=1,ncol=1)
    
    lambda.inv <- matrix(0, nrow=dim, ncol=dim)
    lambda.inv[1:rw,1:rw] <- L.inv 
    x <- (1/sqrt(theta)) * (R.orth %*% sqrt(L.inv) %*% c.mat)
    covar <- (1/theta) * (R %*% lambda.inv %*% t(R))
  }else{   ## if2 (dim.null==0)
### symmetric,positive-definite (no null space) matrix wmat
    messout("Choleski for symmetric, positive definite matrix", verbose)
    CH <- chol(wmat,pivot=T)
    rw1 <- attr(CH, "rank")
    if(rw1 != rw){
      messout("sample.improper.normal, warning:ranks (chol and code) do not agree", verbose)
    print(c(dim,rw,rw1))
    }
    CH.inv <- backsolve(CH, diag(dim))
    x <- CH.inv %*% c.mat /sqrt(theta)
    covar <- 1/theta * (CH.inv %*% t(CH.inv)) } ## end if2 length(dim.null)
  return(list(x =x, covar=covar,rank=rw))
  }
  

  

##########################################################################################
##########################################################################################
##
## FUNCTION NAME:      first.diff
##
## PLACE:    
##
## DESCRIPTION: it returns a matrix which can be used to compute the vector
##              of first differences of another vector.
##             (EV: M differentiates vector,v,as y = M %*% v;  
##                  same as using default R-function diff(v, 1,1) = y)
##
## FORMAT: M <- first.diff(d)
##
## INPUT:   d, scalar, dimension of the vector for which we want to compute the first differences
##        
## OUTPUT:   (d-1) x d matrix of the form:
## 
##         -1    1    0    0    0 ...
## M =      0   -1    1    0    0 ...
##          0    0   -1    1    0 ...
##          0    0    0   -1    1 ...
##          .........................
##
## If this matrix is applied to a d x 1 vector (x_1, x_2, x_3, x_4, ...)
## it returns the (d-1) x 1 vector of first differences: (x_2 - x_1, x_3 - x_2, x_4 - x_3 ...)
## This is matrix is NOT the discretization of the derivative operator, because it is not
## square, however it is similar to it in spirit (the difference is due to the edge effects)
## and it can be used to define a smoothness functional. If x is the vector of the values of
## a function at equally spaced points then |Mx|^2 is like a smoothness functional which
## penalizes the first derivative (| | stands for the usual Euclidean norm).
## If we define a probability density P(x) = K exp (-|Mx|^2) = K exp (- x'M'Mx) then we obtain
## an improper normal distribution, whose inverse covariance is W = M'M (' is the transpose)
##
##     
##
## WRITTEN BY: Federico Girosi 
##             fgirosi@latte.harvard.edu
##             CBRSS, Harvard University
## 
## Last modified: 06/14/2003
## 
##########################################################################################
##########################################################################################



first.diff <- function(d){
 M <- diag(-1,nrow=d-1, ncol=d)
 M[col(M) == row(M) + 1] <- 1
 return(invisible(M));
}


##########################################################################################
##########################################################################################
##
## FUNCTION NAME:      n.diff
##
## PLACE:    
##
## DESCRIPTION: it returns a matrix which can be used to compute the vector
##              of n-th differences of another vector.
## NOTE: input matrix for find.zero.eigen(wmat); wmat = t(M) %*% M
##
## FORMAT: M <- n.diff(n,d)
##
## INPUT:   d, scalar, dimension of the vector for which we want to compute the n-th differences
##          n, scalar
##        
## OUTPUT:   (d-n) x d matrix. It simply the matrix obtained by iterating the first difference
##           operator n times. It resembles the discretization of the $n$-th derivative
##           and it can be used to define smoothness functionals.
##           If x is the vector of the values of function at equally spaced points
##           then |Mx|^2 is like a smoothness functional which penalizes the n-th derivative.
##           If we define a probability density P(x) = K exp (-|Mx|^2) = K exp (- x'M'Mx) then
##           we obtain an improper normal distribution, whose inverse covariance is
##           W = M'M (' is the transpose)
##
##     
##
## WRITTEN BY: Federico Girosi 
##             fgirosi@latte.harvard.edu
##             CBRSS, Harvard University
## 
## Last modified: 06/14/2003
## 
##########################################################################################
##########################################################################################


n.diff <-  function(n,d){
  M <- first.diff(d);
  if (n > 1) {
    for (i in 1:(n-1)){M <- first.diff(nrow(M))%*%M}
  }
  return(invisible(M));
}



derivative.prior.old <- function(d,der.v,k=0){
  W <- matrix(0,nrow=d,ncol=d);
  for (i in 1:length(der.v)){
    mat <- n.diff(i,d);
    weight <- seq(1,nrow(mat))^k;
    W <- W + der.v[i]*t(mat)%*%diag(weight)%*%mat;
  }
  return(invisible(W));
}




normalize.prior <- function(W,sigma=1){
sample <- sample.improper.normal(W,1,1);
covar <- sample$cov;
ave.var <- mean(diag(covar));
W <- W*ave.var/sigma^2;
}


mean.age.profile <- function(ebase){
  ewho <- get("env.who", envir=ebase)
  whoinsampy <- get("whoinsampy", envir=ewho)
  age.vec <- get("age.vec", envir=ewho)
  n.age <- length(age.vec)
  insampy <- list.by.cntry(whoinsampy);
  ap.mean <- rep(0,n.age);
  for (i in 1:length(insampy)){
    y <- insampy[[i]];
    ap <- colMeans(as.data.frame(y),na.rm=TRUE);
    ap.mean <- ap.mean + ap;
  }
  ap.mean <- ap.mean/length(insampy);
  return(invisible(ap.mean));
}


derif <-  function(order,d){
  if (order == 0) return(diag(d));
  D <- n.diff(order,d);
  D.first <- D[1,];
  D.last <- D[nrow(D),];
  if (order%%2 >= 1) {
    D <- rbind(D.first,D);
    order <-  order-1;
  }
  k <-  order/2;
  if (k > 0){
  for (i in 1:k){
    D <- rbind(D.first,D,D.last)
  }
}
  return(invisible(D));
}



##########################################################################################
##########################################################################################
##
## FUNCTION NAME:  derivative.prior
##
## DESCRIPTION: it returns the matrix W that defines a mixed smoothness prior
##
## IMPORTED FUNCTIONS: derif
##
## FORMAT: W <- derivative.prior(d,der.v=c(0,0,1),weight=NA)
##
## INPUT:       d: scalar, dimension of the desired matrix W
##
##          der.v: list, weight of each derivative in the mixed smoothness functional,
##                 starting from the derivative of order 0 (the identity operator).
##                 Example 1: der.v = c(0,0,1) corresponds to a smoothness functional which
##                 penalizes the 2nd derivative (0*identity + 0*1st derivative + 1*2nd derivative)
##                 Example 2: der.v = c(0,1,1) corresponds to a smoothness functional
##                 which penalizes equally the 1st and 2nd derivative 
##                 (0*identity + 1*1st derivative + 1*2nd derivative)
##
##         weight: scalar or d x 1 vector. The smoothness functional is an integral and
##                 weight defines the measure in the integral. If weight is a scalar the
##                 measure dx in the integral is taken to be x^weight dx. If weight = NA
##                 it is taken as 0, so the measure is uniform. If weight is a d x 1 vector
##                 it is taken as the discretization of the measure and used "as is".
##                 
##
## OUTPUT:      W: d x d matrix, the matrix which defines the smoothness prior.
##
##     
##
## WRITTEN BY: Federico Girosi 
##             fgirosi@latte.harvard.edu
##             CBRSS, Harvard University
## 
## Last modified: 07/28/2003
## 
##########################################################################################
##########################################################################################

derivative.prior <- function(d,der.v=c(0,0,1),weight=NA){
  if (length(weight) == 0) stop("weight vector has 0 length in derivative.prior");
  W <- matrix(0,nrow=d,ncol=d);
  if (any(is.na(weight)) == TRUE) weight <-  0;
  if (length(weight)==1) {weight <- seq(1,d)^weight;
                          weight <- weight/sum(weight);
                        }
  if (length(weight) != d) stop("weight vector has wrong length");
  for (i in 1:length(der.v)){
    if (der.v[i] != 0){
      mat <- derif(i-1,d);
      W <- W + der.v[i]*t(mat)%*%diag(weight)%*%mat;
    }
  }
  return(invisible(W));
}



##########################################################################################
##########################################################################################
##
## FUNCTION NAME:  sm.func.age.time.mu
##
## DESCRIPTION: it computes the values of a mixed age/time smoothness functional
##              for all the countries in a cross-sectional time series
##
## IMPORTED FUNCTIONS: 
##
## FORMAT: sm <- sm.func.age.time.mu(y,prior,n.age)
##
## INPUT:     y:  (list) a cross-sectional time series over countries and age groups
##                 missing values not allowed (since derivatives must be taken)
##
##        prior: (list) list of 4 parameters which uniquely define a mixed age/time
##                     smoothness functional
##                     prior$a.deriv: list for derivatives with respect to age;
##
##                     prior$a.weight: weights for age part of the functional;
##
##                     prior$t.deriv: list for derivatives with respect to time
##
##                     prior$t.weight: weights for time part of the functional;
##
##              n.age: number of age groups;
## 
## OUTPUT:   sm: (list) a list of length = # of countries 
##
##     
##
## WRITTEN BY: Federico Girosi 
##             fgirosi@latte.harvard.edu
##             CBRSS, Harvard University
## 
## Last modified: 08/5/2003
## 
##########################################################################################
##########################################################################################


sm.func.age.time.mu <- function(y,prior,n.age,prior.mean.by.cntry=NA, ebase){
  ewho <- get("env.who", envir=ebase)
  age.vec <- get("age.vec", envir=ewho)
  n.age <- length(age.vec)
    
    if (!is.na(prior.mean.by.cntry)){
      m <- make.average.age.profile(y,bycntry=prior.mean.by.cntry);
      y <-  add.subtract.age.profile(y,m,subtract=TRUE)
    }
  a.deriv <- prior$a.deriv;
  t.deriv <-  prior$t.deriv;
  a.weight <- prior$a.weight;
  t.weight <- prior$t.weight;
  if (is.null(a.weight)) a.weight <-  0;
  if (is.null(t.weight)) t.weight <-  0;
  if (is.null(t.deriv)) t.deriv <-  1;
  if (is.null(a.deriv)) a.deriv <-  1;  
  W.age <- derivative.prior(n.age,a.deriv,a.weight);  
  yc <- list.by.cntry(y);
  func <- function(mu,W.age,t.deriv,t.weight){
    if (any(is.na(mu))) stop("Missing values in sm.func.age.time.mu")
    W.time <- derivative.prior(nrow(mu),t.deriv,t.weight);  
    return(sum(diag(W.age%*%t(mu)%*%W.time%*%mu)));
  }
  h <- lapply(yc,FUN=func,W.age,t.deriv,t.weight);
}


##########################################################################################
##########################################################################################
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
##########################################################################################
##########################################################################################

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



##########################################################################################
##########################################################################################
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
##     
##
## WRITTEN BY: Federico Girosi 
##             fgirosi@latte.harvard.edu
##             CBRSS, Harvard University
## 
## Last modified: 08/14/2003
## 
##########################################################################################
##########################################################################################

impute.csts <- function(y){
         s <- lapply(y,FUN=impute.ts);
return(invisible(s))
}



##########################################################################################
##########################################################################################
##
## FUNCTION NAME:  smooth.vector.faster
##
## DESCRIPTION: it smooths a vector of uniformly spaced observations according
##              to a prescribed smoothness functional, and finds the optimal
##              parameter using Generalized Cross Validation (GCV). This
##              requires as input the eigenvector/eigenvalue decomposition
##              of the matrix defining the smoothness functional and
##              therefore is quite fast.
##
## IMPORTED FUNCTIONS: none
##
## FORMAT: s <- smooth.vector.faster(y,R,e,rank,m=NA,trace.factor=1,vector.only=FALSE)
##
## INPUT:   y: (array, vector) an N x 1 vector or array of equally spaced observations.
##             Missing values are not allowed.
##
##          R: (array) the rotation matrix of the eigenvalue/eigenvector decomposition
##             of the symmetric, positive definite matrix W which represents the prior.
##             In formulas W = R L R'
##
##          e: (vector) the vector of eigenvalues of the matrix W above (the diagonal
##             elements of W)
##
##       rank: (scalar) the rank of W above (number of non zero elements of e)
##
##          m: (vector or NA) and N x 1 vector specifying the mean of the prior
##             if missing the prior is 0 mean
##    
##   trace.factor: a factor, usually > 1, which multiplies the trace in the denominator
##                 of the GCV score. If 1 this corresponds to ordinary GCV. If set
##                 to a number larger than 1 it will make GCV select a larger
##                 smoothness parameter than it would have otherwise chosen.
##                 A number larger than (say 1.5) is useful to correct the tendency
##                 of GCV to undersmooth the data.
##
##    vector.only: (logical) it controls what is returned by the function
## 
## OUTPUT:      s: (list or array) if vector.only == TRUE then we return an N x 1 array
##                 with the smoothed version of y. Otherwise this is a list with 3 elements:
##
##                 s$y.hat is the smoothed version of y
##                 s$lambda is the optimal smoothed parameter selected by GCV
##
## SEE ALSO: smooth.vector, which does not require the eigenvector/eigenvalue
##           decomposition of the prior matrix.
##     
##
## WRITTEN BY: Federico Girosi 
##             fgirosi@latte.harvard.edu
##             CBRSS, Harvard University
## 
## Last modified: 08/16/2003
## 
##########################################################################################
##########################################################################################

smooth.vector.faster <- function(y,R,e,rank,m=NA,trace.factor=1,vector.only=FALSE){
  if (!is.array(y)) y <- as.array(y);
  rn  <- rownames(y);
  n <- nrow(y);
  if (length(dim(y)) !=2) y <- array(y,dim=c(n,1),dimnames=(list(rn,NULL)))

  if (!any(is.na(m))){
    if (!is.array(m))  m <- as.array(m);
    if (nrow(m) != nrow(y)) stop("m has wrong nrow in smooth.vector.faster");
    if (length(dim(m)) !=2) m <- array(m,dim=c(nrow(m),1));
    y <- y - m;
  }

  if (any(is.na(y)))  y <- approx(1:n,y,xout=1:n,rule=2)$y;
  tR <- t(R);
  ty <- tR %*% y;  
  inter <- range(1/e[1:rank]);
  min.obj <-  optimize(f=gcv,interval=inter,e=e,ty=ty,trace.factor=trace.factor);
  lambda.opt <- min.obj$minimum;
  S <-  R %*% diag(1/(1+lambda.opt*e)) %*% tR 
  yhat <- S%*%y;
  if (!any(is.na(m))) yhat <- yhat + m;
  rownames(yhat) <- rn;
  if (vector.only==FALSE){
    return(list(y.hat=yhat,lambda=lambda.opt));
  } else {
    return(yhat)
  }
}




##########################################################################################
##########################################################################################
##
## FUNCTION NAME:  smooth.vector
##
## DESCRIPTION: it smooths a vector of uniformly spaced observations according
##              to a prescribed smoothness functional, and finds the optimal
##              parameter using Generalized Cross Validation (GCV). 
##
## IMPORTED FUNCTIONS: smooth.vector.faster, derivative.prior
##
## FORMAT: s <- smooth.vector(y,prior,m=NA,trace.factor=1,vector.only=FALSE)
##
## INPUT:   y: (array, vector) an N x 1 vector or array (missing values allowed),
##             that we wish to smooth. Missing values are allowed, although they
##             are imputed before the smoothing with a piecewise linear interpolation,
##             which is not the optimal smoother.
##
##      prior: (list) a list with two elements, specifying the smoothness prior
##
##              prior$deriv: (vector) the "mixing vector" in the mixed smoothness
##                           prior (see manual and documentation for derivative.prior)
##             prior$weight: (scalar or vector) specifies the measure in the integral
##                           of the smoothness functions (see manual and documentation
##                           for derivative.prior)
##
##          m: (vector or NA) and N x 1 vector specifying the mean of the prior
##             if missing the prior is 0 mean
##    
##   trace.factor: a factor, usually > 1, which multiplies the trace in the denominator
##                 of the GCV score. If 1 this corresponds to ordinary GCV. If set
##                 to a number larger than 1 it will make GCV select a larger
##                 smoothness parameter than it would have otherwise chosen.
##                 A number larger than (say 1.5) is useful to correct the tendency
##                 of GCV to undersmooth the data.
##
##    vector.only: (logical) it controls what is returned by the function
## 
## OUTPUT:      s: (list or array) if vector.only == TRUE then we return an N x 1 array
##                 with the smoothed version of y. Otherwise this is a list with 3 elements:
##
##                 s$y.hat is the smoothed version of y
##                 s$y.imp is the imputed version of y
##                 s$lambda is the optimal smoothed parameter selected by GCV
##
## SEE ALSO: smooth.vector.faster
##     
##
## WRITTEN BY: Federico Girosi 
##             fgirosi@latte.harvard.edu
##             CBRSS, Harvard University
## 
## Last modified: 08/16/2003
## 
##########################################################################################
##########################################################################################

smooth.vector <- function(y,prior,m=NA,trace.factor=1,vector.only=FALSE){
  prior.deriv <- prior$deriv;
  prior.weight <- prior$weight;

  ###
  ### The following is not strictly necessary: it simply makes sure that
  ### that the argument y to smooth.vector.faster is an array of dimension
  ### dim(y) = c(n,1). This is because smooth.vector.faster was getting
  ### confused between vectors of length(n), arrays of dimension dim(y) = n and arrays
  ### of dimensions dim(y) = c(n,1). This should not be necessary anymore, since
  ### the problem has been fixed in smooth.vector.faster, and it could be eliminated
  ### if it slows things down
  ###
  
  if (!is.array(y)) y <- as.array(y);
  rn  <- rownames(y);
  n <- nrow(y);
  if (length(dim(y)) !=2) {
    y <- array(y,dim=c(n,1),dimnames=(list(rn,NULL)))
  }

  ###
  ### Here we define the quantities relating to the prior
  ### Notice that we use qr to figure out the rank. This is unnecessary
  ### since we have the eigenvalues already and we could just count
  ### the number of 0s. We do this is so that we can delegate to qr
  ### the problem of deciding what is 0 and what is not. Alternatively
  ### we could use the function find.zero.eigen: we need to find out
  ### which is faster.
  
  W <- derivative.prior(n,prior.deriv,prior.weight);
  #w <- La.eigen(W,symmetric=TRUE);
  w <- list(vectors=NULL, values=NULL)
  R <- w$vectors;
  e <- w$values;
  rank <- qr(W)$rank;

  ###
  ### This is the call which actually peform all the calculations using GCV
  ###
  
  yhat <- smooth.vector.faster(y,R,e,rank,m,trace.factor,vector.only)
  return(yhat);  
}



##########################################################################################
##########################################################################################
##
## FUNCTION NAME: gcv
##
## DESCRIPTION: it computes the Generalized Cross-Validation (GCV) score
##              for a given value of lambda. To be used in a fast search
##              for the optimal smoothness parameter.
##
## IMPORTED FUNCTIONS: none
##
## FORMAT: s <- gcv(lambda,e,ty,trace.factor=1)
##
## INPUT:  lambda: (vector) a vector of values  for the smoothness parameter.
##
##              e: (vector) the eigenvalues of the symmetric, positive definite
##                  matrix W which represents the prior. In formulas W = R diag(e) R'
##
##             ty: (array) the rotated data vector. If y is the data vector,
##                 ty = R' y, where R is given above.
##    
##   trace.factor: a factor, usually > 1, which multiplies the trace in the denominator
##                 of the GCV score. If 1 this corresponds to ordinary GCV. If set
##                 to a number larger than 1 it will make GCV select a larger
##                 smoothness parameter than it would have otherwise chosen.
##                 A number larger than (say 1.5) is useful to correct the tendency
##                 of GCV to undersmooth the data.
##
## 
## OUTPUT:      s: (vector) the GCV scores of the smoothness parameters stored
##                 in lambda
##
##
## SEE ALSO: smooth.vector.faster
##     
##
## WRITTEN BY: Federico Girosi 
##             fgirosi@latte.harvard.edu
##             CBRSS, Harvard University
## 
## Last modified: 08/16/2003
## 
##########################################################################################
##########################################################################################


gcv <- function(lambda,e,ty,trace.factor=1){
  func <- function(l,e,ty,trace.factor){
  N <- length(e);
  f1 <- sapply(e,FUN=function(x,lambda){(x/(1+lambda*x))^2},l);
  f2 <- l^2*crossprod(f1,ty^2);
  f3 <- sapply(e,FUN=function(x,lambda){1/(1+lambda*x)},l);
  f4 <- (N-trace.factor*sum(f3))^2;
  f5 <- N*f2/f4;
  return(f5)}
  g <- sapply(lambda,FUN=func,e,ty,trace.factor);
  return(g);
}



##########################################################################################
##########################################################################################
##
## FUNCTION NAME: smooth.csts.time
##
## DESCRIPTION: it smooths a cross-sectional time series over time, using
##              Generalized Cross-Validation (GCV) to select an optimal
##              smoothness parameter.
##
## IMPORTED FUNCTIONS: smooth.vector
##
## FORMAT: s <- smooth.csts.time(y,time.prior=list(deriv=c(0,0,1),weight=0),trace.factor=1)
##
## INPUT:         y: (list) a cross-sectional time series
##
##
##       time.prior: (list) a list with two elements, specifying the smoothness
##                   functional over time.
##
##         time.prior$deriv: (vector) the "mixing vector" in the mixed smoothness
##                           prior (see manual and documentation for derivative.prior)
##         time.prior$weight: (scalar or vector) specifies the measure in the integral
##                           of the smoothness functions (see manual and documentation
##                           for derivative.prior)
##
##    
##   trace.factor: a factor, usually > 1, which multiplies the trace in the denominator
##                 of the GCV score. If 1 this corresponds to ordinary GCV. If set
##                 to a number larger than 1 it will make GCV select a larger
##                 smoothness parameter than it would have otherwise chosen.
##                 A number larger than (say 1.5) is useful to correct the tendency
##                 of GCV to undersmooth the data.
##
## 
## OUTPUT:      s: (list) the smoothed version of the time series y
##     
##
## WRITTEN BY: Federico Girosi 
##             fgirosi@latte.harvard.edu
##             CBRSS, Harvard University
## 
## Last modified: 08/16/2003
## 
##########################################################################################
##########################################################################################


smooth.csts.time <- function(y,time.prior=list(deriv=c(0,0,1),weight=0),trace.factor=1, verbose=T){
  ebase <- try(get("env.base", envir=parent.frame()),silent=T)
  if(class(ebase)!="try-error")
    verbose <- get("verbose", envir=ebase)
  messout("Smoothing over time", verbose);      
################# BEGIN LOCALLY DEFINED FUNCTION #############################
###
### func is the function which smooths one cross-section over time 
###
  
  func <- function(x,time.prior,trace.factor){
    s <- smooth.vector(x,time.prior,NA,trace.factor);
    if (!is.array(s$y.hat)) stop("s is not array in smooth.csts.time");
    rownames(s$y.hat) <- rownames(x);
    return(s)
  }
################# BEGIN LOCALLY DEFINED FUNCTION #############################
  

  y <- lapply(y,func,time.prior,trace.factor);
  return(invisible(y))
}



##########################################################################################
##########################################################################################
##
## FUNCTION NAME: smooth.csts.ages
##
## DESCRIPTION: it smooths a cross-sectional time series over time, using
##              Generalized Cross-Validation (GCV) to select an optimal
##              smoothness parameter.
##
## IMPORTED FUNCTIONS: smooth.vector
##
## FORMAT: s <- smooth.csts.ages(y,age.prior=list(deriv=c(0,0,1),weight=0),trace.factor=1)
##
## INPUT:         y: (list) a cross-sectional time series
##
##
##        age.prior: (list) a list with two elements, specifying the smoothness
##                   functional over time.
##
##          age.prior$deriv: (vector) the "mixing vector" in the mixed smoothness
##                           prior (see manual and documentation for derivative.prior)
##          age.prior$weight: (scalar or vector) specifies the measure in the integral
##                           of the smoothness functions (see manual and documentation
##                           for derivative.prior)
##
##    
##   trace.factor: a factor, usually > 1, which multiplies the trace in the denominator
##                 of the GCV score. If 1 this corresponds to ordinary GCV. If set
##                 to a number larger than 1 it will make GCV select a larger
##                 smoothness parameter than it would have otherwise chosen.
##                 A number larger than (say 1.5) is useful to correct the tendency
##                 of GCV to undersmooth the data.
##
## 
## OUTPUT:      s: (list) the smoothed version of the time series y
##
##     
##
## WRITTEN BY: Federico Girosi 
##             fgirosi@latte.harvard.edu
##             CBRSS, Harvard University
## 
## Last modified: 08/16/2003
## 
##########################################################################################
##########################################################################################


smooth.csts.ages <- function(y,age.prior=list(deriv=c(0,0,1),weight=0),m=NA,
                             trace.factor=1, verbose=T){
  ebase <- try(get("env.base", envir=parent.frame()),silent=T)
  if(class(ebase)!="try-error")
    verbose <- get("verbose", envir=ebase)
messout("Smoothing over age", verbose);      
################# BEGIN LOCALLY DEFINED FUNCTION #############################
###
### func is the function which smooths one country-specific matrix over age
###
  
  func <- function(y,a.prior,m,trace.factor){
###
### Here y is a country specific logmortality matrix
###

    prior.deriv <- a.prior$deriv;
    prior.weight <- a.prior$weight;
    n <- ncol(y);
    W <- derivative.prior(n,prior.deriv,prior.weight);
    # w <- La.eigen(W,symmetric=TRUE);
  w <- list(vectors=NULL, values=NULL)
    R <- w$vectors;
    e <- w$values;
    rank <- qr(W)$rank;
    y.smooth <-  apply(y,MARGIN=1,FUN=function(x,R,e,rank,m,trace.factor){smooth.vector.faster(x,R,e,rank,m,trace.factor)$y.hat},R,e,rank,m,trace.factor)
    y.smooth <- t(y.smooth);
    colnames(y.smooth) <- colnames(y);
    return(y.smooth)
  }

################# END LOCALLY DEFINED FUNCTION #############################
  
###
### first we impute missing values
###

y <- impute.csts(y);

###
### then we convert to a list of country matrices
###

y <- list.by.cntry(y);

###
### now we smooth
### We need to do things differently depending on whether m is a list of
### country specific mean age profiles or not. 
###

    if (is.list(m)){

      ###
      ### we need to pass to each smoother a different vector
      ###

      ind <- as.list(1:length(y));
      names(ind) <- names(y);
      y <- lapply(ind,FUN=function(n,y,age.prior,m,trace.factor){yc <- y[[n]];
                                                    mc <- m[[n]];
                                                    ys <- func(yc,age.prior,mc,trace.factor);
                                                    return(ys)},y,age.prior,m,trace.factor)
    } else {
      
      ###
      ### Here m is the same for all cross-sections
      ###
      
      y <-  lapply(y,FUN=func,age.prior,m,trace.factor);      
    }

###
### and we convert to a cross-sectional list with countries and ages
###

y <- list.by.csid(y);
return(invisible(y));
}



## ######################################################################################
## ######################################################################################
##
## FUNCTION NAME: smooth.csts.age.time
##
## DESCRIPTION: it smooths a cross-sectional time series over time, using
##              Generalized Cross-Validation (GCV) to select an optimal
##              smoothness parameter.
##
## IMPORTED FUNCTIONS: smooth.vector.faster,digitpull,derivative.prior,
##                     find.zero.eigen,impute.csts,list.by.cntry.long,
##                     unlist.by.cntry.long,conv.char
##
## GLOBALS:   who.digit.first,who.cntry.digits,who.age.digits,who.year.digits 
##
##
## FORMAT: s <- smooth.csts.age.time(y,
##                  prior=list(a.deriv=c(0,0,1),a.weight=0,t.deriv=c(0,0,1),t.weight=0),
##                  trace.factor=1)
##
## INPUT:         y: (list) a cross-sectional time series
##
##
##            prior: (list) a list with 4 elements, specifying the smoothness
##                   functional over time and ages.
##
##               prior$a.deriv: (vector) the "mixing vector" in the mixed smoothness
##                              prior over ages (see manual and documentation for
##                              derivative.prior)
##               prior$a.weight: (scalar or vector) specifies the measure in the integral
##                               over ages of the smoothness functions (see manual and
##                               documentation for derivative.prior)
##               prior$t.deriv: (vector) the "mixing vector" in the mixed smoothness
##                              prior over time (see manual and documentation for
##                              derivative.prior)
##               prior$t.weight: (scalar or vector) specifies the measure in the integral
##                               over time of the smoothness functions (see manual and
##                               documentation for derivative.prior)
##
##    
##     trace.factor: a factor, usually > 1, which multiplies the trace in the denominator
##                   of the GCV score. If 1 this corresponds to ordinary GCV. If set
##                   to a number larger than 1 it will make GCV select a larger
##                   smoothness parameter than it would have otherwise chosen.
##                   A number larger than (say 1.5) is useful to correct the tendency
##                   of GCV to undersmooth the data.
##
## 
## OUTPUT:      s: (list) the smoothed version of the time series y
##
##     
##
## WRITTEN BY: Federico Girosi 
##             fgirosi@latte.harvard.edu
##             CBRSS, Harvard University
## 
## Last modified: 08/16/2003
## 
## ######################################################################################
## ######################################################################################


smooth.csts.age.time <- function(y,
                                 prior=list(a.deriv=c(0,0,1),a.weight=0,t.deriv=c(0,0,1),t.weight=0),
                                 prior.mean.by.cntry=NA,
                                 trace.factor=1, ebase){
  ewho <- get("env.who", envir=ebase);
  who.digit.first  <- get("who.digit.first", envir=ewho)
  who.cntry.digits <- get("who.cntry.digits", envir=ewho)
   who.age.digits  <- get("who.age.digits", envir=ewho)
   who.year.digits <- get("who.year.digits", envir=ewho) 
  verbose <- get("verbose", envir=ewho)
  messout("Smoothing over age/time",verbose);      
    digit.cntry.begin <- who.digit.first + 1 
    digit.cntry.end   <- who.digit.first + who.cntry.digits
    digit.age.begin   <- digit.cntry.end + 1
    digit.age.end     <- digit.cntry.end + who.age.digits
    digit.year.begin  <- digit.age.end   + 1 
    digit.year.end    <- digit.age.end   + who.year.digits 

################# BEGIN LOCALLY DEFINED FUNCTION #############################
###
### func is the function which smooths one long vector with GCV
###
  env.bae <- ebase
  func <- function(y,prior,trace.factor,n.age, ebase = env.bae){
  ewho <- get("env.who", envir=ebase);
  who.digit.first  <- get("who.digit.first", envir=ewho)
  who.cntry.digits <- get("who.cntry.digits", envir=ewho)
   who.age.digits  <- get("who.age.digits", envir=ewho)
   who.year.digits <- get("who.year.digits", envir=ewho)
  age.vec <- get("age.vec", envir=ewho)

    digit.cntry.begin <- who.digit.first + 1 
    digit.cntry.end   <- who.digit.first + who.cntry.digits
    digit.age.begin   <- digit.cntry.end + 1
    digit.age.end     <- digit.cntry.end + who.age.digits
    digit.year.begin  <- digit.age.end   + 1 
    digit.year.end    <- digit.age.end   + who.year.digits 

###
### Here y is a country specific logmortality long vector
###

###
### We use the fact that time series of different age groups have same length
###
    
    nt <- nrow(y)/n.age;
    t <- conv.char(digitpull(as.numeric(names(y)[1:nt]),digit.year.begin,digit.year.end))

###
### Auxiliary matrices needed to the prior over ages and time
###
    
    W.age <- derivative.prior(n.age,prior$a.deriv,prior$a.weight)
    W.time <- derivative.prior(nt,prior$t.deriv,prior$t.weight)

    rownames(W.age) <- 1:nrow(W.age);
    colnames(W.age) <- 1:nrow(W.age);
    rownames(W.time) <-  1:nrow(W.time);
    colnames(W.time) <-  1:nrow(W.time);

###
### W is the matrix which defines the prior over age and time
###
    
    W <- kronecker(W.age, W.time, make.dimnames=TRUE);

###
### I compute the quantities needed for the GCV smoother
###
    
    aux <- find.zero.eigen(W,only.values=FALSE);
    R <- aux$w.eigv;
    e <- aux$w.eigen;
    rank <- aux$w.rank

###
### Now we use the GCV smoother. Notice that the prior has zero mean because
### average age profiles are constant over time.
###
    
    y.smooth <-  smooth.vector.faster(y,R,e,rank,m=NA,trace.factor=trace.factor);
    return(y.smooth)
  }

################# END LOCALLY DEFINED FUNCTION #############################

###
### We read the age groups from the names of y 
###

    age.vec <- unique.default(digitpull(as.numeric(names(y)),digit.age.begin,digit.age.end));
    ages <- conv.char(age.vec);
    n.age <- length(age.vec);

###
### first we impute missing values
###

    y <- impute.csts(y);

###
###  we subtract the mean age profile (which will be 0 if prior.mean.by.cntry is NA)
###

    if (!is.na(prior.mean.by.cntry)){
      m <- make.average.age.profile(y,bycntry=prior.mean.by.cntry);
      y <-  add.subtract.age.profile(y,m,subtract=TRUE)
    }
    
    
###
### This will create a list with one long vector for each country
### The vector is the concatenation of the time series of all the age groups
###

    y <- list.by.cntry.long(y);

###
### now we smooth
###

    y <-  lapply(y,FUN=func,prior,trace.factor,n.age);
    lambda <- lapply(y,function(x){x$lambda});
    y <- lapply(y,function(x){x$y.hat});


###
### and we convert to the standard format of cross-sectional list with countries and ages
###

    y <- unlist.by.cntry.long(y);
    
###
###  we add the mean age profile (which will be 0 if prior.mean.by.cntry is NA)
###

    if (!is.na(prior.mean.by.cntry)){
      y <-  add.subtract.age.profile(y,m,subtract=FALSE)
    }
    
    return(list(y.hat=y,lambda=lambda));
}



## ######################################################################################
## ######################################################################################
##
## FUNCTION NAME: smooth.csts.gcv
##
## DESCRIPTION: it smooths a cross-sectional time series with a mixed age/time prior,
##              using Generalized Cross-Validation (GCV) to select an optimal
##              smoothness parameter.
##
## IMPORTED FUNCTIONS: smooth.vector.faster,digitpull,derivative.prior,
##                     find.zero.eigen,impute.csts,list.by.cntry.long,
##                     unlist.by.cntry.long,conv.char
##
## GLOBALS:   who.digit.first,who.cntry.digits,who.age.digits,who.year.digits 
##
##
## FORMAT: lst <- smooth.csts.gcv(y,
##                                std=NA,
##                                prior=list(a.deriv=c(0,0,1),a.weight=0,
##                                           t.deriv=c(0,0,1),t.weight=0),
##                                prior.mean.by.cntry=NA,
##                                pool.lambda=TRUE,
##                                trace.factor=1,
##                                outdir="OUTPUT/",
##                                compact=TRUE)
##
## INPUT:         y: (list) a cross-sectional time series over countries and ages
##
##              std: (list) a cross-sectional time series, with the same structure
##                   as y, of standard deviations, to be used as observation weights
##                   in the likelihood. When the prior is not zero mean the mean of
##                   the prior is weighted according to the average of the std
##                   over time, and therefore the correct result will be obtained
##                   only if the std are constant over time (which is usually the case).
##                   If missing no weighting occurs.
##
##            prior: (list) a list with 4 elements, specifying the smoothness
##                   functional over time and age.
##
##               prior$a.deriv: (vector) the "mixing vector" in the mixed smoothness
##                              prior over ages (see manual and documentation for
##                              derivative.prior)
##               prior$a.weight: (scalar or vector) specifies the measure in the integral
##                               over ages of the smoothness functions (see manual and
##                               documentation for derivative.prior)
##               prior$t.deriv: (vector) the "mixing vector" in the mixed smoothness
##                              prior over time (see manual and documentation for
##                              derivative.prior)
##               prior$t.weight: (scalar or vector) specifies the measure in the integral
##                               over time of the smoothness functions (see manual and
##                               documentation for derivative.prior)
##
##   prior.mean.by.cntry: a parameter used to choose the mean of the prior.
##                        If NA then the prior is zero mean.
##                        If FALSE then the prior has the same mean  for all countries.
##                        If TRUE then the prior has country specific means.
##
##           pool.lambda: (logical) if TRUE the regularization parameter is pooled
##                        within each country. This is meaningful only when the prior
##                        is over age only or time only, because when we smooth over
##                        age and time there is only one regularization parameter to
##                        begin with. For example if the prior is over
##                        age only and pool.lambda is TRUE then the age profiles for
##                        different years share the same lambda. Similarly, if we
##                        smooth over time only and pool.lambda is TRUE then
##                        the regularization parameter is the same in all the age groups
##                        (probably not a good idea).
##
##     trace.factor: a factor, usually > 1, which multiplies the trace in the denominator
##                   of the GCV score. If 1 this corresponds to ordinary GCV. If set
##                   to a number larger than 1 it will make GCV select a larger
##                   smoothness parameter than it would have otherwise chosen.
##                   A number larger than 1 (say 1.5) is useful to correct the tendency
##                   of GCV to undersmooth the data.
##
##           outdir: (string) pathname of the directory where we want to save the graphs
##                   if missing no graphs are produced
##
##          compact: (logical) if FALSE we save one ps file for each country,
##                   otherwise we compact all the files in one multi-page ps file.
## 
## OUTPUT:      lst: (list) an object of type "cross-sectional forecast" with the
##                   following non-empty fields:
##                      disease=whodisease (standard)
##                      gender=whogender (standard)
##                      yrest=yrest (standard)
##                      age.vec=age.vec (standard)
##                      cntry.vec=cntry.vec (standard)
##                      yhatin=y.smooth: the smoothed version of y
##                      insampy=y
##
##                    This object can be plotted with the command:
##                        plot.forecast(lst,outputdir,insample.only=TRUE,modelstr,annot)
##
##                    If outdir is not missing we produces graphs which show both the
##                    data (in dotted green) and the smoothed time series
##                    (in continuous red) in the usual format of forecast plots.
##                    Both the age profiles of the data and its smoothed version are shown.
##                    If compact=FALSE each country has its own graph, whose name follows
##                    standard conventions and the suffix "SMOOTH". For example
##                    oting3_SMOOTH_Italy.ps refers to other infectious diseases in Italian
##                    females.
##                    If compact=TRUE the country specific files are compacted in one
##                    multi-page ps file, named as before but without the country names
##                    (for example oting3_SMOOTH.ps). The legend of the file will tell
##                    which kind of smoothing has been applied (age, time or age and time).
##                   
##
## WRITTEN BY: Federico Girosi 
##             fgirosi@latte.harvard.edu
##             CBRSS, Harvard University
## 
## Last modified: 08/16/2003
## 
## ######################################################################################
## ######################################################################################


smooth.csts.gcv <- function(y,std=NA,prior=list(a.deriv=c(0,0,1),a.weight=0,
                                                t.deriv=c(0,0,1),t.weight=0),
                            prior.mean.by.cntry=NA,pool.lambda=TRUE,trace.factor=1,
                            outdir="OUTPUT/",compact=TRUE, ebase){
  ewho <- get("env.who", envir=ebase)
  who.digit.first  <- get("who.digit.first", envir=ewho)
  who.cntry.digits <- get("who.cntry.digits", envir=ewho)
  who.age.digits   <- get("who.age.digits", envir=ewho)
  who.year.digits  <- get("who.year.digits", envir=ewho)
### Structure of dataset cstsid: 
  digit.cntry.begin <- who.digit.first + 1 
  digit.cntry.end   <- who.digit.first + who.cntry.digits
  digit.age.begin   <- digit.cntry.end + 1
  digit.age.end     <- digit.cntry.end + who.age.digits
  digit.year.begin  <- digit.age.end   + 1 
  digit.year.end    <- digit.age.end   + who.year.digits 

###
### Some useful quantities
###

  age.vec <- unique.default(digitpull(as.numeric(names(y)),digit.age.begin,digit.age.end));
  cntry.vec <- unique.default(digitpull(as.numeric(names(y)),digit.cntry.begin,digit.cntry.end));

###
### The insample period is the same for all cross-sections. Therefore I can get
### yrest, the last year of the insample period, from the 1st cross-section.
### yrest is the same as the global whoyrest, I am just minimizing the number
### of globals I use.
###
  insampy <-  y;
  yrest <- max(digitpull(as.numeric(rownames(y[[1]])),digit.year.begin,digit.year.end))


###
### If requested we weight y with the inverse of the standard deviations
###
    if (!any(is.na(std))){
      if (length(std) != length(y)) stop("smooth.csts.age.time: length of std and y differ");
      ind <- 1:length(y);
      y <- lapply(ind,FUN=function(n,y,std){y[[n]]/std[[n]]},y,std);
      names(y) <-  names(insampy);
    }

###
### First we need to figure out what kind of prior we are dealing with
### In principle we could use smooth.csts.age.time for all types, but it
### it would much slower whenever we smooth over age only, or time only
###
prior.type <- "age and time";
if (length(prior$t.deriv)==1) prior.type <- "age";
if (length(prior$a.deriv)==1) prior.type <- "time";
if (prior.type == "age and time") pool.lambda <- TRUE;
  
###
### If the prior type is "age", we need to figure out whether the prior
##  is zero mean or not, and if it is zero mean whether the mean
### varies by country
###
### if prior.mean.by.cntry is NA then the prior is zero mean
### if prior.mean.by.cntry is FALSE then the prior has global mean
### (the same for all countries)
### if prior.mean.by.cntry is TRUE then the prior has country specific means

if (prior.type!="age") prior.mean.by.cntry <-  NA;
  
if (is.na(prior.mean.by.cntry)) {
  m <- NA;
} else {
  m <-  make.average.age.profile(insampy,bycntry=prior.mean.by.cntry)
}

###
### if we weight observations with standard deviations we need to weight the mean
### age profile too. Evene if the mean average profile is the same for all countries
### its weighted version will not be the same for all countries. So m must
### be converted to a list which runs over countries.
###

  if (!is.na(prior.mean.by.cntry)){
### I make std of the same type as insampy so I can use list.by.cntry
    ind <- 1:length(std);
    names(ind) <- names(insampy);      
    std <- lapply(ind,FUN=function(n,std,y){std[[n]] <- as.matrix(std[[n]]);
                                            rownames(std[[n]]) <- rownames(insampy[[n]]);
                                            rownames(std[[n]]) <- rownames(insampy[[n]]);
                                            colnames(std[[n]]) <- "std";
                                            return(std[[n]])},std,insampy);
    std.age.profile <- make.average.age.profile(std,bycntry=TRUE);
    if (is.list(m)){
    ind <- 1:length(m);
    names(ind) <- names(m);      
    m <-  lapply(ind,FUN=function(n,m,s){m[[n]] / s[[n]]},m,std.age.profile);
    } else {
    ind <- 1:length(std.age.profile);
    names(ind) <- names(std.age.profile);
    m <-  lapply(ind,FUN=function(n,m,s){m / s[[n]]},m,std.age.profile);    
    }
  }
  
  
###
### Now we smooth according to the type of prior (age,time or age and time)
###

    
if (prior.type=="age" && pool.lambda==FALSE){
  a.prior <- list(deriv=prior$a.deriv,weight=prior$a.weight);
  y.smooth <- smooth.csts.ages(y,a.prior,m,trace.factor);
}

if (prior.type=="time" && pool.lambda==FALSE){
    t.prior <- list(deriv=prior$t.deriv,weight=prior$t.weight);
    y.smooth <- smooth.csts.time(y,t.prior,trace.factor);
    y.smooth <- lapply(y.smooth,FUN=function(x){x$y.hat});
}
  
lambda <- NA;
if (pool.lambda == TRUE){
    s <- smooth.csts.age.time(y,prior,prior.mean.by.cntry,trace.factor);
    y.smooth <- s$y.hat;
    lambda <- s$lambda;
}

###
### we remove the effect of weighting by standard deviations
### 

    if (!any(is.na(std))){
      ind <- 1:length(y.smooth);
      y.smooth <- lapply(ind,FUN=function(n,y,std){y[[n]]*std[[n]]},y.smooth,std)
      names(y.smooth) <-  names(insampy);
    }

###
### Now we build an object which is suitable for plotting
###

lst <- list(yhatin=y.smooth,insampy=insampy,lambda=lambda,outsapmy = NULL)
### plot.forecast does not exist in the package
###  if (!is.na(outdir)){
###    plot.forecast(lst,outputdir=outdir,insample.only=TRUE,
###                  modelstr="SMOOTH",annot=paste("Smoothing over",prior.type))

###
### If we saved to ps files and the option compact is TRUE we condense all the printed files
### in one, multipage ps file, with 8 ps files in each page and delete the single files
###
    whodisease <- get("whodisease", envir=ewho)
    whogender <- get("whogender", envir=ewho)
    
    if (compact==TRUE){
      outps <- paste(outdir,whodisease,"g",whogender,"_SMOOTH*",".ps",sep="");
      outps.compact <- paste(outdir,whodisease,"g",whogender,"_SMOOTH",sep="");
      system(paste("/bin/rm ", outps.compact,".ps",sep=""));
      mpage.str <- paste("mpage -c1",outps,"> ",outps.compact);
      system(mpage.str);
      system(paste("/bin/rm",outps));
      system(paste("/bin/mv",outps.compact,paste(outps.compact,".ps",sep="")));
    }

return(invisible(lst))
}


## ######################################################################################
## ######################################################################################
##
## FUNCTION NAME: smooth.csts.iterative
##
## DESCRIPTION: it smooths a cross-sectional time series with a prior specified by 
##              several smoothness functionals by iteratively smoothing using one
##              prior at the time and using GCV to select the optimal smoothing
##              parameter. This attempts to solve the problem that GCV with multiple
##              smoothing parameters is hard to do efficiently.
##
## IMPORTED FUNCTIONS: smooth.csts.gcv
##
## GLOBALS: none
##
##
## FORMAT: lst <- smooth.csts.iterative(y,
##                                      prior.list,
##                                      n.iter=1,
##                                      std=NA,
##                                      prior.mean.by.cntry=NA,
##                                      pool.lambda=FALSE,
##                                      trace.factor=1)
##
## INPUT:         y: (list) a cross-sectional time series over countries and ages
##
##       prior.list: (list) a list of priors in standard format (as in smooth.csts.gcv)
##
##           n.iter: (scalar) how many times we want to go through the list of the
##                   priors and smooth accordingly
##
##       std,prior.mean.by.cntry,pool.lambda,trace.factor: see documentation for
##                                                         smooth.csts.gcv    
## 
## OUTPUT:      lst: (list) an object of type "cross-sectional forecast" with the
##                   following non-empty fields:
##                      yrest=yrest (standard)
##                      age.vec=age.vec (standard)
##                      cntry.vec=cntry.vec (standard)
##                      yhatin=y.smooth: the smoothed version of y
##                      insampy=y
##
##                    This object can be plotted with the command:
##                        plot.forecast(lst,outputdir,insample.only=TRUE,modelstr,annot)
##
##
## WRITTEN BY: Federico Girosi 
##             fgirosi@latte.harvard.edu
##             CBRSS, Harvard University
## 
## Last modified: 09/13/2003
## 
## ######################################################################################
## ######################################################################################

smooth.csts.iterative <- function(y,prior.list,n.iter,std=NA,prior.mean.by.cntry=NA,
                                  pool.lambda=FALSE,trace.factor=1, ebase){
  
### ewho <- get("env.who", envir=ebase)

  for (i in 1:n.iter){
    for (j in 1:length(prior.list)){
      print("***");
      prior <-  prior.list[[j]];
      lst <- smooth.csts.gcv(y,std=std,prior=prior,prior.mean.by.cntry=prior.mean.by.cntry,
                             trace.factor=trace.factor,pool.lambda=pool.lambda);
      y <- lst$yhatin;
    }
  }
return(lst)
}


## ######################################################################################
## ######################################################################################

prior.matrix.age.time <- function(prior,n.age=NULL,n.time, ebase){
    ewho <- get("env.who", envir=ebase)
    if(length(n.age) <= 0)
    n.age <- length(get("age.vec", envir=ewho))
         
    W.age <- derivative.prior(n.age,prior$a.deriv,prior$a.weight)
    W.time <- derivative.prior(n.time,prior$t.deriv,prior$t.weight)

###
### I give names so I can check that the matrix W is built correctly
###

    age.vec <- 1:n.age;
    time.vec <- 1:n.time;    
    rownames(W.age) <- age.vec;
    colnames(W.age) <- age.vec;
    rownames(W.time) <-  time.vec;
    colnames(W.time) <-  time.vec;

###
### W is the matrix which defines the prior over age and time
###
    
    W <- kronecker(W.age, W.time, make.dimnames=TRUE);
    return(invisible(W))
}

## ######################################################################################
## ######################################################################################


add.subtract.age.profile <- function(y,m,subtract=TRUE){

  y <-  list.by.cntry(y);

  ind <- as.list(1:length(y));
  names(ind) <- names(y);
  y <-  lapply(ind,FUN=function(n,y,m,subtract){
    yc <- y[[n]];
    mc <-  m;
    if (is.list(m)) mc <- m[[n]];
    aux <-  -1;
    if (subtract==FALSE) aux <-  1;
    mc <- aux*mc;
    t(apply(yc,1,FUN=function(x,mc){x+mc},mc))},y,m,subtract)
  y <-  list.by.csid(y);
  return(invisible(y));
}


#############################################################################################
## the following function is there for compatibility with old versions


make.W <-  function(n){
  s <- matrix(0,nrow=n,ncol=n);
  for (i in 2:(n-1)){s[i,i+1] <-  1;
                     s[i,i-1] <-  1};
  s[1,2] <-  1;
  s[n,n-1] <- 1;
  W <- diag(apply(s,1,sum))-s;
  s <- 2*s/sum(diag(W));
  W <- 2*W/sum(diag(W));
  return(s,W);
}


