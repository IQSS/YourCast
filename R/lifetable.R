############################################################################
### Function to return life table from yourcast object with mortality rates
### Konstantin Kashin
### August 2013
############################################################################

lifetable <- function(obj,ax=0.5,a0=0.07,nx=NULL,dv.log=NA){
	# get information about geographies
	# count number of "g" in index.code
 	index.code.list <- unlist(strsplit(obj$aux$index.code,""))
	N.g <- sum(grepl("g",index.code.list))
  
  	# vector of unique geo areas (csids)
  	csid.vec <- substring(names(obj$yhat),1,N.g)
  	names(csid.vec) <- names(obj$yhat)
	csid.unique <- unique(csid.vec)
	
  	# match dv.log arg and check if dependent variable is logged
  	dv.log <- as.character(dv.log)
  	dv.log.arg <- match.arg(dv.log, c(NA,"TRUE","FALSE"))
  	if(is.na(dv.log.arg)){
  		if(length(obj$call$formula)<=1){
  			stop("Unable to verify if dependent variable is logged because formula is missing. Please set dv.log argument to either 'TRUE' or 'FALSE'.")
  		} else{
  		dv.log <- any(grepl("log",obj$call$formula[[2]]))
  		}
  	} else{
  		dv.log <- as.logical(dv.log.arg)
  		}
  	
  	# iterate over each geography
  	geolist <- lapply(csid.unique,function(csid){
		y.mat <- sapply(obj$yhat[grep(csid,csid.vec)],function(x) x[,1]) 
		yhat.mat <- sapply(obj$yhat[grep(csid,csid.vec)],function(x) x[,2])
		colnames(y.mat) <- colnames(yhat.mat) <-as.character(as.integer(gsub(csid,"",colnames(y.mat))))
		out <- list(obs=lifetableCalc(mat=y.mat,ax=ax,a0=a0,nx=nx,dv.log=dv.log), pred=lifetableCalc(mat=yhat.mat,ax=ax,a0=a0,nx=nx,dv.log=dv.log))
		out
		})
	names(geolist) <- csid.unique
	if(length(geolist)==1){return(geolist[[1]])} else{
		return(geolist)
	}
}
