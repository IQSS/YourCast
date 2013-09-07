### DESCRIPTION: it takes the formula, ff, and split it into three elements
###              the left, "~", and rigt sides.  Splits the right-side into
###              components according to the separation sign="+"
###              It uses findregx(ch) to obtain the nude variable for every element  
###              in the right-hand-side of formula
###
### INPUT : the formula ff <- log(dth) ~ lag(gdp^3, -7) + log(tobacco)^0.55 + time
###         a boolean ops, if ops = F then it returns "gdp*income", otherwise
###         returns "gdp" "income" as atomic elements. 
###         
### OUTPUT: The right and left hand side atomic components
###         gdp, tobacco, dth
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
parseff <- function(ff, ops=FALSE)
{
  fc  <- as.character(ff)
  rs <- fc[3]
  ls  <- fc[2]
  
  if(!identical(trim.blanks(fc[1]), "~") )
    stop("Bad parsing in formula")
  
  vec <- strsplit(rs,"\\+")[[1]]
  vec <- sapply(vec, FUN="trim.blanks")
  vec <- lapply(vec, FUN="findregex")
  dth <- findregex(ls);
  if(ops){
  vec <- lapply(vec, FUN="split.ops")
  dth <- lapply(dth, FUN="split.ops")
  
  }
    
  lst <- list(rightside = vec, leftside=dth)
  return(lst)
}

### DESCRIPTION: It takes a char or string to find the arguments 
###              It returns its atomic components
###              ch <- "lag(log(gdp)^4, -5)"
###              find.regex(ch) returns "gdp"
###              ch <- sqrt(lag(gdp^0.5, -10)
###              returns "gdp"
###              ch <- lag(log(gdp/time)^4, -5)
###              returns "gdp/time"
###
### Elena Villalon
### November 2005

 findregex <- function(ch)
  {
    
    if(is.na(ch))
      return(ch)
    
    ch <- trim.blanks(ch)
    ffc <- paste("1", "~", ch)
    ff <- as.formula(ffc)
    return(all.vars(ff)) ###internal R-function 
    

  }

    
## ************************************************************************
##
## FUNCTION NAME:      trim.blanks
##
##
## DESCRIPTION: It returns a string of characters with blanks eliminated  
##              at the beginning or end of the word only.
##
## FORMAT:  trim.blanks(word), with word a character string
##
## OUTPUT: word with the blanks at the beginning and end eliminated
##         For example word <- "   my house in Concord   "
##         trim.blanks(word) <- "my house in Concord". 
##
## WRITTEN BY: Elena Villalon 
##             evillalon@iq.harvard.edu
##             IQS, Harvard University
##
## LAST MODIFIED: 28/09/2005
## 
## ************************************************************************
  
  trim.blanks <- function(x) {
### at the beginning of string"^" gets anny number (+) of white spaces
    if(x=="" || is.na(x))
      return(x)
    
     x <- sub('^ +', '', x)
### at the ending of string"$" gets anny number (+) of white spaces
     x <- sub(' +$', '', x)
     return(x)
    }

### helper function check.ops(c, cols, opvec)
### for any element, c, of a formula and for corresponding columns of
### given matrix of covariates and a vector of operations. It checks if
### c contains the operation in opvec (basically, multiplication and division), 
### and after anayzing its componets, if
### all components are included as elemenst of cols.
### For example gdp * hc *fat, will find elements gdp,hc, fat and check
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
split.ops <- function(c, opvec=c("\\*", "/")){
  vec <-  c
  for(op in opvec)
    vec <- unlist(sapply(vec, strsplit,op))
  if(length(vec) <= 0)
    return(c)
  
  vec <- sapply(vec,FUN="trim.blanks")
  names(vec) <-  NULL
  return(vec)
}

  

                
  
