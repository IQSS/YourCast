### this function returns the standard deviation of the predictions of the MAP algorithm
### it has the same format as whocov: a list with one element for each csid, and each element
### is a time series which spans both insample and out of sample.
### it takes as input ecxc, the environment of cxc.model, but it really only needs whocov and
### S.list, a list of matrices which is computed in cxc.model. 
###
### AUTHOR Federico Girosi
###        Rand, CA
###        girosi@rand.org

std.MAP <- function(S.list){
  ebase <- get("env.base", envir=parent.frame())
  env.base <- ebase
  ewho <- get("env.who", envir=ebase)
  whocov <- get("whocov",ewho)
  who.age.digits <- get("who.age.digits", envir=ewho)
  who.cntry.digits <- get("who.cntry.digits", envir=ewho)
  who.digit.first <- get("who.digit.first", envir=ewho)
  digit.cntry.begin <- who.digit.first + 1
  digit.cntry.end <- who.cntry.digits + who.digit.first
  digit.age.begin <- who.digit.first + who.cntry.digits + 1;
  digit.age.end <- who.digit.first + who.cntry.digits + who.age.digits;
  
  std.list.by.cntry <- S.list

  cntries <- digitpull(as.numeric(names(whocov)),digit.cntry.begin,digit.cntry.end)
  cntrylist <- unique.default(digitpull(as.numeric(names(whocov)),digit.cntry.begin,digit.cntry.end))
###  print(cntrylist)
  for (c in cntrylist){
    zlist <- whocov[cntries == c]
    alist <- digitpull(as.numeric(names(zlist)),digit.age.begin,digit.age.end)
    if (!identical(order(alist),1:length(alist)))
      stop("mismatch in function which computes variance of estimates")
    Z <- diag.mat(zlist)
    S <- S.list[[as.character(c)]]
    std.list.by.cntry[[as.character(c)]] <- NULL
###  tmp  <- try(sqrt(diag(Z%*%S%*%t(Z))), silent=T)
    tmp <-try (diag(Z%*%S%*%t(Z)), silent=TRUE)
    if(class(tmp)!= "try-error"){
    tmp <- sqrt(abs(tmp))
    std.list.by.cntry[[as.character(c)]] <- tmp
  }
  }

### now convert to same format as whocov

  std.list <- whocov
  for (i in 1:length(std.list)){
    cntry <- digitpull(names(whocov)[i],digit.cntry.begin,digit.cntry.end)
    s <- std.list.by.cntry[[as.character(cntry)]]
    rnames <- rownames(whocov[[i]])
    std.list[[i]] <- s[rnames]
  }
  return(std.list)
}

