################################################
### Function to create arrays from yourcast object
### Konstantin Kashin
### August 2013
################################################

array.yourcast <- function(x,unlog=FALSE){
	# get information about geographies
	# count number of "g" in index.code
 	index.code.list <- unlist(strsplit(x$aux$index.code,""))
	N.g <- sum(grepl("g",index.code.list))
  
  	# vector of unique geo areas (csids)
  	csid.vec <- substring(names(x$yhat),1,N.g)
  	names(csid.vec) <- names(x$yhat)
	csid.unique <- unique(csid.vec)
	sample.begin <- x$aux$sample.frame[1]
	sample.end <- x$aux$sample.frame[2]
	forecast.begin <- x$aux$sample.frame[3]
	forecast.end <- x$aux$sample.frame[4]
		
	geolist <- lapply(csid.unique,function(csid){
		y.mat <- sapply(x$yhat[grep(csid,csid.vec)],function(x) x[,1]) 
		yhat.mat <- sapply(x$yhat[grep(csid,csid.vec)],function(x) x[,2])
		colnames(y.mat) <- colnames(yhat.mat) <- as.character(as.integer(gsub(csid,"",colnames(y.mat))))
		if(unlog){
			y.mat <- exp(y.mat)
			yhat.mat <- exp(yhat.mat)
		}
		array.out <- array(data=NA,dim=c(dim(y.mat),3),dimnames=c(dimnames(y.mat),list(c("y","yhat","comb"))))
		array.out[,,"y"] <- y.mat
		array.out[,,"yhat"] <- yhat.mat
		array.out[1:which(rownames(y.mat)==as.character(sample.end)),,"comb"] <- y.mat[1:which(rownames(y.mat)==as.character(sample.end)),]
		array.out[which(rownames(y.mat)==as.character(forecast.begin)):nrow(y.mat),,"comb"] <- yhat.mat[which(rownames(y.mat)==as.character(forecast.begin)):nrow(y.mat),]
		array.out
	})
	names(geolist) <- csid.unique
	if(length(geolist)==1){return(geolist[[1]])} else{
		return(geolist)
	}
}
