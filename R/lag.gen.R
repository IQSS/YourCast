lag.gen <- function(geolist.i,data,split.index,N.g,lag,years.insamp,
                    years.pred,years.lag,years.tot,sample.frame,response,
                    covars,index.add,vars.nolag) {
  index.geo.i <- grep(paste("^",geolist.i,sep=""),names(data))
  geo.i <- data[index.geo.i]
  csids.i <- names(geo.i)

  # Determine age groups in geo area i
  # Parse index.code to find out how to extract age groups
  N.a <- length(split.index[split.index=="a"])
  ages <- as.numeric(sapply(csids.i,substr,start=N.g+1,
                                stop=N.g+N.a))
  # Find out where lagged data for each age group is stored,
  # unless group is younger than lag time
  ages.lag <- ages-lag
  ages.lag[ages.lag<min(ages)] <- NA
  agelag.loc <- c()
  for(a in 1:length(ages)) {
    if(is.na(ages.lag[a])){agelag.loc <- NA}
    if(!is.na(ages.lag[a])){
      find.age <- which(ages==ages.lag[a])
      if(length(find.age)>0){agelag.loc[a] <- find.age}
      if(length(find.age)==0){age.close <- which.min(abs(ages-ages.lag[a]))
                              agelag.loc[a] <- age.close
                              warning(paste("Could not find age",
                                            ages.lag[a],
                                            "data for age",ages[a],
                                            "cohort in covariates to be lagged. Used data from age",
                                            ages[age.close],"cohort instead."))
                            }
    }
  }

  # Isolate the response for each cross section in insample years
  # and place NAs in prediction years
  response.i <- lapply(geo.i,
                       function(x,response){
                         # Try to include response data from
                         # years.pred in case this is a validation
                         # exercise
                         year.max <- max(as.numeric(rownames(x)))
                         if(identical(year.max,sample.frame[2])) {
                         out <- x[as.character(years.insamp),
                                  paste(response)]
                         out <- c(out,rep(NA,length(years.pred))) }
                         else if(year.max>=sample.frame[4]) {
                           out <- x[as.character(years.tot),
                                  paste(response)]}
                         else {out <- x[as.character(years.insamp),
                                  paste(response)]
                         out <- c(out,rep(NA,length(years.pred)))}
                         out <- as.matrix(out)
                         colnames(out) <- response
                         rownames(out) <- years.tot
                         return(out)},
                       response=response)

  # Now need to extract covariates from relevant lag years
  # (depends on lag time)
  covars.i <- lapply(geo.i,function(x,index.add){
                # Find out which columns match a covariate in the formula
                grab <- unlist(sapply(covars,function(y){which(colnames(x)==y)}))
                out <- x[as.character(years.lag),grab]
                year.lag <- as.numeric(rownames(out))
                out <- cbind(out,years.lag)
                rownames(out) <- years.tot
                colnames(out) <- c(covars,"year.lag")
                return(out)},index.add=index.add)

  # Now take covariates from each cross section from lag period years
  # and put them in spot necessary to match with relevant response data
  covars.i.lag <- covars.i[agelag.loc]

  # Deal with any covariates that will not be lagged
  covars.nolag.i <- NULL
  if(!is.null(vars.nolag)){
    covars.nolag.i <- lapply(geo.i,function(x,vars.nolag){
      grab <- unlist(sapply(vars.nolag,function(y){which(colnames(x)==y)}))
      year.max <- max(as.numeric(rownames(x)))
      if(year.max>=sample.frame[4]){
        out <- x[as.character(years.tot),grab]}
      if(year.max<sample.frame[4]){
        years.have <- sample.frame[1]:year.max
        nyears.miss <- length(years.tot)-length(years.have)
        out <- x[as.character(years.have),grab]
        na.add <- matrix(NA,nrow=nyears.miss,ncol=length(grab))
        out <- cbind(out,na.add)}
       out <- as.matrix(out)
      colnames(out) <- vars.nolag
      rownames(out) <- years.tot
      return(out)},
                             vars.nolag=vars.nolag)
  }

  # Now merge response data and lagged covariates
  if(!is.null(covars.nolag.i)){data.geo <- mapply(cbind,response.i,covars.i.lag,
                   covars.nolag.i,SIMPLIFY=FALSE)}
  if(is.null(covars.nolag.i)){data.geo <- mapply(cbind,response.i,covars.i.lag,
                   SIMPLIFY=FALSE)}
  # If 'index' variable included in formula, add an index to every
  # cross section
  if(index.add) {data.geo <- lapply(X=data.geo,FUN=cbind,
                                    index=1:length(years.tot))}
  
  return(data.geo)
}
