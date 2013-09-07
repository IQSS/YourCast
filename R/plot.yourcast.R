########################################################
### YourCast S3 plot method
### Author: Konstantin Kashin
### August 2013
###
###
### Function based upon previous version of YourCast S3
### plotting method written by Jon Bischof
########################################################

##### START OF PLOT FUNCTION #####
plot.yourcast <- function(x, plots=c("age","time"),
						 title=NULL, subtitle=NULL, age.opts=list(), time.opts=list(),
                         threedim.opts=list(),totcount.opts=list(), ex0.opts=list(),
                         print="device",
                         print.args = list(),
                         dv.log=NA,...){
  
  	# get plot types and do checks
  	types <- match.arg(plots,c("age","time","totcount","ex0","threedim"),several.ok=TRUE)
  	if(length(types)>3){
  		stop("You cannot select more than 3 plot types.")
  	}
  	if("totcount" %in% types & "ex0" %in% types){
  		stop("You cannot plot both total counts and life expectancy at birth.")
  	}
  	
  	# set default options for plotting
    age.options <- list(xlab = NULL, ylab = NULL, insamp.obs = TRUE, insamp.predict = TRUE, age.select=NULL,time.select=NULL,unlog=FALSE)
    time.options <- list(xlab = NULL, ylab = NULL, insamp.obs = FALSE, insamp.predict = TRUE, age.select=NULL,time.select=NULL,unlog=FALSE)   
    totcount.options <- list(xlab = NULL, ylab = NULL, insamp.obs = TRUE, insamp.predict = TRUE)    
    ex0.options <- list(xlab = NULL, ylab = NULL, insamp.obs = TRUE, insamp.predict = TRUE, ax=0.5, a0=0.07,nx=NULL)
    threedim.options <- list(xlab=NULL,ylab=NULL,zlab=NULL,insamp.predict=TRUE,
    				 unlog=FALSE,args.wireframe=list(), screen=list(z=-40, x=-60, y=0))
    print.options <- list(height=NULL,
                         width=NULL,filename=NULL,dpath=getwd())
                         
    # get names of option lists
    age.names <- names(age.options)
    time.names <- names(time.options)
    totcount.names <- names(totcount.options)
    ex0.names <- names(ex0.options)
    threedim.names <- names(threedim.options)
    print.names <- names(print.options)
    
    # update default options with user specified options
    age.options[(nam.age.opts <- names(age.opts))] <- age.opts
    time.options[(nam.time.opts <- names(time.opts))] <- time.opts
    totcount.options[(nam.totcount.opts <- names(totcount.opts))] <- totcount.opts
    ex0.options[(nam.ex0.opts <- names(ex0.opts))] <- ex0.opts
    threedim.options[(nam.threedim.opts <- names(threedim.opts))] <- threedim.opts
    print.options[(nam.print.args <- names(print.args))] <- print.args		
    
    # return warnings if option names do not match
    if (length(noNames <- nam.age.opts[!nam.age.opts %in% age.names]))
        warning("Unknown names in age plot options: ", paste(noNames, collapse = ", "), call.=FALSE)
    if (length(noNames <- nam.time.opts[!nam.time.opts %in% time.names])) 
        warning("Unknown names in time plot options: ", paste(noNames, collapse = ", "), call.=FALSE) 
    if (length(noNames <- nam.totcount.opts[!nam.totcount.opts %in% totcount.names])) 
        warning("Unknown names in total counts plot options: ", paste(noNames, collapse = ", "), call.=FALSE)
    if (length(noNames <- nam.ex0.opts[!nam.ex0.opts %in% ex0.names])) 
        warning("Unknown names in life expectancy plot options: ", paste(noNames, collapse = ", "), call.=FALSE)     
    if (length(noNames <- nam.threedim.opts[!nam.threedim.opts %in% threedim.names])) 
        warning("Unknown names in 3D plot options: ", paste(noNames, collapse = ", "), call.=FALSE)    
   if (length(noNamesPrint <- nam.print.args[! nam.print.args %in% print.names])) 
        warning("Unknown names in plot print options: ", paste(noNamesPrint, collapse = ", "), call.=FALSE)    		 
    				 
  	# loading aux files from 'yourcast' object
  	G.names <- NULL
  	if(!is.null(x$aux$G.names) && length(x$aux) > 0)
  	  {G.names <- x$aux$G.names}
  	
  	# match dv.log arg and check if dependent variable is logged
  	dv.log <- as.character(dv.log)
  	dv.log.arg <- match.arg(dv.log, c(NA,"TRUE","FALSE"))
  	if(is.na(dv.log.arg)){
  		if(length(x$call$formula)<=1){
  			stop("Unable to verify if dependent variable is logged because formula is missing. Please set dv.log argument to either 'TRUE' or 'FALSE'.")
  		} else{
  		dv.log <- any(grepl("log",x$call$formula[[2]]))
  		}
  	} else{
  		dv.log <- as.logical(dv.log.arg)
  		}
  	
  	# set unlog options for the graphs
  	# unless both unlog==TRUE and the DV is logged, set unlog to FALSE
  	if(!(age.options$unlog==TRUE & dv.log)){
  		age.options$unlog <- FALSE
  	}
  	if(!(time.options$unlog==TRUE & dv.log)){
  		time.options$unlog <- FALSE
  	}
  	# for ex0 and totcount plots, set 'unlog' to TRUE if DV is logged
  	# since quantities must be unlogged for these plots
  	if(dv.log){
  		ex0.options$unlog <- totcount.options$unlog <- TRUE
  	} else{
  		ex0.options$unlog <- totcount.options$unlog <- FALSE
  	}
  	
  	# count number of "g", "a", "t", and "-" in index.code
 	index.code.list <- unlist(strsplit(x$aux$index.code,""))
	N.g <- sum(grepl("g",index.code.list))
	N.a <- sum(grepl("a",index.code.list))
	N.t <- sum(grepl("t",index.code.list))
	N.dash <- sum(grepl("-",index.code.list))
  
  	# vector of unique geo areas (csids)
  	csid.vec <- substring(names(x$yhat),1,N.g)
  	names(csid.vec) <- names(x$yhat)
	csid.unique <- unique(csid.vec)

	# make a geolist per each csid and assign to environment as geoGGGG where GGGG is code
	# geolist is a list, with first element matrix y (observed values) and second element yhat (predicted values)
	# rows are times, columns are ages

	# geonames is a vector of geolists by name
 	geonames <- c()
 	
 	for(i in 1:length(csid.unique)){
	csid.i <- csid.unique[i]
	y.mat <- sapply(x$yhat[grep(csid.i,csid.vec)],function(x) x[,1]) 
	yhat.mat <- sapply(x$yhat[grep(csid.i,csid.vec)],function(x) x[,2])
	colnames(y.mat) <- colnames(yhat.mat) <- as.character(as.integer(gsub(csid.i,"",colnames(y.mat))))
	# create name
	geolist <- list(y=y.mat, yhat=yhat.mat)
	# assign to environment as geoGGGG and store name
	assign(paste("geo",csid.unique[i],sep=""),geolist)
	geonames[i] <- paste("geo",csid.unique[i],sep="")
	rm(csid.i,y.mat,yhat.mat,geolist)
	}
	
  # If printing to pdf and filenames provided, make sure vector of
  # strings the correct length
  filename <- print.options$filename
  if(!is.null(filename)){
  	if(length(filename)==length(geonames)){
  		filename.vec <- filename
  		} else{
  		filename.vec <- rep(filename,length.out=length(geonames))	
  		}
  	}
   
  # set plot size dimensions (if NULL)
  if(is.null(print.options$width)){
  	width <- switch(length(types),"1"=6,"2"=10,"3"=12)}
  if(is.null(print.options$height)){height <- 5}
  
  # extract sample frame
  sf <- x$aux$sample.frame
  
  # iterate over geonames and plot
  for(i in 1:length(geonames)) {	
      # Use either user provided file name or default
      if(!is.null(filename)){
      	if(length(filename)==length(geonames)){
      		file <- paste(print.options$dpath,"/",filename.vec[i],
                             ".pdf",sep="")	
      	} else{
      		file <- paste(print.options$dpath,"/",filename.vec[i],csid.unique[i],
                             ".pdf",sep="")		
      		}
      	} else{
      		file <- paste(print.options$dpath,"/","plot",csid.unique[i],
                             ".pdf",sep="")
      	}
             
      # set country name (geography name)
      cntryname <- NULL
              
      # if there is a G.names object, pull the csid names from it
      # set the label for ith csid to cntryname
      if(!is.null(G.names)) {
        cntryname <- G.names[grep(csid.unique[i],
                                  G.names[,1]),2]}
     
     # define main = title for entire device
     # format for main is "title, cntryname" if both are provided; else either "title", "ctryname" or NULL
     main <- paste(if(!is.null(title)){title},
                   if(!is.null(title) & !is.null(cntryname)){", "},
                   if(!is.null(cntryname)){cntryname},sep="")
                             
	if(print=="device"){dev.new(height=5,width=width)}
	if(print=="pdf"){pdf(file=file,height=5,width=width)}
	
	### Call plots here
	plot.list <- lapply(types,function(x){ switch(x,
		age=ageplot(get(geonames[i]),age.opts=age.options, sample.frame=sf),
		time=timeplot(get(geonames[i]),time.opts=time.options, sample.frame=sf),
		totcount=totalplot(get(geonames[i]),tot.opts=totcount.options, sample.frame=sf),
		ex0=ex0plot(get(geonames[i]),ex0.opts=ex0.options, sample.frame=sf),
		threedim=threedimplot(get(geonames[i]),threedim.opts=threedim.options, sample.frame=sf))})
	
	### Pass arguments to grid
	grid.args <- c(plot.list,list(nrow=1,ncol=length(types),main=main,sub=subtitle))
	do.call(grid.arrange,grid.args)
	if(print=="pdf"){dev.off()}
	if(print=="device"){if(i != length(geonames)) {
    user.option(i,csid.unique,cntryname,G.names)}}
  } # end of loop over geonames
} # end of plot fxn

############################################
# user prompt function
user.option <- function(i,csid.unique,cntryname,G.names) {
       ANSWER <-
         readline(paste("Plot for '",
                        ifelse(is.null(cntryname),
                               paste("geo unit",csid.unique[i+1]),
                               paste(G.names[grep(csid.unique[i+1],
                                                  G.names[,1]),2])),
                        "' (plot ",i+1," of ",length(csid.unique),
                        ") next. Type 'y' to continue or 'n' to quit: ",
                        sep=""))
       if (substr(ANSWER, 1, 1) == "n")
         stop("Function terminated by user.")
     }

####################################################
##### Define global variables so as to appease codetools for variables passed to aes() in ggplot2
if(getRversion() >= "2.15.1")  utils::globalVariables(c("y", "yhat","age","time","ex.pred","ex.obs"))

####################################################

### Age plot: time on x axis, dep var on y axis, groupings are ages
ageplot <- function(el,age.opts, sample.frame){
	if(is.null(age.opts$xlab)) {xlab <- "Time"}
    if(is.null(age.opts$ylab)) {ylab <- "Data and Forecasts"}
	
	times <- as.integer(rownames(el$y)) 
	ages <- as.integer(colnames(el$y))
	time.min <- min(times)
	time.max <- max(times)
	age.min <- min(ages)
	age.max <- max(ages)
	sample.end <- sample.frame[2]
	forecast.begin <- sample.frame[3]
	
	# melt dataframes
	y.melt <- melt(el$y)
	colnames(y.melt) <- c("time", "age", "y")		
	yhat.melt <- melt(el$yhat)
	colnames(yhat.melt) <- c("time", "age", "yhat")
	
	# if insamp.predict==FALSE, remove the in-sample predictions from the dataframe
	if(!age.opts$insamp.predict){
		yhat.melt <- yhat.melt[yhat.melt$time>sample.end,]
	}
	
	# if age.select option is TRUE:
	if(!is.null(age.opts$age.select)){
		if(!is.numeric(age.opts$age.select)){
			warning("'age.select' must be a numeric vector. Ignoring option.", call.=FALSE)
			age.select <- ages
		} else{
			if(length(age.opts$age.select)>1){
				age.select <- age.opts$age.select	
			} else{
				age.select <- unique(c(seq(age.min,age.max,by=as.integer(age.opts$age.select)),age.max))
			} 
		}
		y.melt <- y.melt[y.melt$age %in% age.select,]
		yhat.melt <- yhat.melt[yhat.melt$age %in% age.select,]
	}
	# if time.select option is TRUE:
	if(!is.null(age.opts$time.select)){
		if(!is.numeric(age.opts$time.select)){
			warning("'time.select' must be a numeric vector. Ignoring option.", call.=FALSE)
			time.select <- times
		} else{
			if(length(age.opts$time.select)>1){
				time.select <- age.opts$time.select	
			} else{
				time.select <- unique(c(seq(time.min,time.max,by=as.integer(age.opts$time.select)),time.max))
			} 
		}
		y.melt <- y.melt[y.melt$time %in% time.select,]
		yhat.melt <- yhat.melt[yhat.melt$time %in% time.select,]
	}	
	
	# if 'unlog'=TRUE, take exp()
	if(age.opts$unlog){
		y.melt$y <- exp(y.melt$y)
		yhat.melt$yhat <- exp(yhat.melt$yhat)
	}
	
	#plot
	plot.age <- ggplot(yhat.melt, aes(x=time, y=yhat, color=age, group=age)) + geom_line() + theme_bw()  + scale_x_continuous(xlab) + scale_y_continuous(ylab)
	plot.age <- plot.age + scale_color_gradientn("Age",colours=rainbow(7)) + theme(legend.margin=unit(-0.02,"npc"),legend.text=element_text(size=8))
	# if insamp.obs==TRUE, plot in-sample observed values:
	if(age.opts$insamp.obs){
	plot.age <- plot.age + geom_path(data=y.melt,aes(x=time,y=y,color=age,group=age), linetype="dashed",na.rm=TRUE)
		}
	plot.age
}

####################################################

### time plot: age on x axis, dep var on y axis, groupings are times	
timeplot <- function(el,time.opts,sample.frame){
    if(is.null(time.opts$xlab)) {xlab <- "Age"}
    if(is.null(time.opts$ylab)) {ylab <- "Forecasts"}
	
	times <- as.integer(rownames(el$y)) 
	ages <- as.integer(colnames(el$y))
	time.min <- min(times)
	time.max <- max(times)
	age.min <- min(ages)
	age.max <- max(ages)
	sample.end <- sample.frame[2]
	
	# melt dataframes
	y.melt <- melt(el$y)
	colnames(y.melt) <- c("time", "age", "y")		
	yhat.melt <- melt(el$yhat)
	colnames(yhat.melt) <- c("time", "age", "yhat")
	
	# if insamp.predict==FALSE, remove the in-sample predictions from the dataframe
	if(!time.opts$insamp.predict){
		yhat.melt <- yhat.melt[yhat.melt$time>sample.end,]
	}
	
	# if age.select option is TRUE:
	if(!is.null(time.opts$age.select)){
		if(!is.numeric(time.opts$age.select)){
			warning("'age.select' must be a numeric vector. Ignoring option.", call.=FALSE)
			age.select <- ages
		} else{
			if(length(time.opts$age.select)>1){
				age.select <- time.opts$age.select	
			} else{
				age.select <- unique(c(seq(age.min,age.max,by=as.integer(time.opts$age.select)),age.max))
			} 
		}
		y.melt <- y.melt[y.melt$age %in% age.select,]
		yhat.melt <- yhat.melt[yhat.melt$age %in% age.select,]
	}
	# if time.select option is TRUE:
	if(!is.null(time.opts$time.select)){
		if(!is.numeric(time.opts$time.select)){
			warning("'time.select' must be a numeric vector. Ignoring option.", call.=FALSE)
			time.select <- times
		} else{
			if(length(time.opts$time.select)>1){
				time.select <- time.opts$time.select	
			} else{
				time.select <- unique(c(seq(time.min,time.max,by=as.integer(time.opts$time.select)),time.max))
			} 
		}
		y.melt <- y.melt[y.melt$time %in% time.select,]
		yhat.melt <- yhat.melt[yhat.melt$time %in% time.select,]
	}	
	
	# if 'unlog'=TRUE, take exp()
	if(time.opts$unlog){
		y.melt$y <- exp(y.melt$y)
		yhat.melt$yhat <- exp(yhat.melt$yhat)
	}
		
	plot.time <- ggplot(yhat.melt, aes(x=age, y=yhat, color=time, group=time)) + geom_line() + theme_bw() + scale_x_continuous(xlab) + scale_y_continuous(ylab)
	plot.time <- plot.time + scale_color_gradientn("Time",colours=rainbow(7)) + theme(legend.justification=c(0,1),legend.position=c(0.05,1),legend.direction="horizontal", legend.text=element_text(angle=45), legend.title.align=1, legend.background = element_rect(fill="transparent"))
	if(time.opts$insamp.obs){
	plot.time <- plot.time + geom_line(data=na.omit(y.melt),aes(x=age,y=y,color=time,group=time), linetype="dashed") 
	}	
	plot.time
}

####################################################

### total plot (only makes sense for counts)
totalplot <- function(el,tot.opts,sample.frame){
	times <- as.integer(rownames(el$y)) 
	sample.end <- sample.frame[2]
	
	# if 'unlog'=TRUE, take exp()
	if(tot.opts$unlog){
		y <- exp(el$y)
		yhat <- exp(el$yhat)
	}
	
	if(!is.null(tot.opts$age.select)){
		y <- y[which(times %in% as.character(tot.opts$age.select)),]
		yhat <- yhat[which(times %in% as.character(tot.opts$age.select)),]
	}

	totals.df <- data.frame(time=times,y=rowSums(y),yhat=rowSums(yhat))
	totals.df.pred <- totals.df[totals.df$time>sample.end,]
	
	plot.total <- ggplot(data=switch(as.character(tot.opts$insamp.predict),"TRUE"=totals.df,"FALSE"=totals.df.pred), aes(x=time, y=yhat)) +
    	geom_line() + theme_bw() + scale_x_continuous("Time") +
    	scale_y_continuous("Data and Forecasts (Totals)")
    	if(tot.opts$insamp.obs){
    		plot.total <- plot.total + geom_point(data=na.omit(totals.df),aes(x=time,y=y))
    	}
	plot.total
}

####################################################

### life expectancy plots
ex0plot <- function(el,ex0.opts,sample.frame){
	times <- as.integer(rownames(el$y)) 
	sample.end <- sample.frame[2]
	
	# calculate life expectancies at birth
	ex0.obs <- lifetableCalc(el$y,ax=ex0.opts$ax,a0=ex0.opts$a0,nx=ex0.opts$nx, dv.log=ex0.opts$unlog)[,"0","ex"]
	ex0.pred <- lifetableCalc(el$yhat,ax=ex0.opts$ax,a0=ex0.opts$a0, nx=ex0.opts$nx, dv.log=ex0.opts$unlog)[,"0","ex"]
	
	# plot life expectancy
	ex.df <- data.frame(time=as.numeric(times),ex.obs=ex0.obs,ex.pred=ex0.pred)
	ex.df.pred <- ex.df[ex.df$time>sample.end,]
	
	plot.ex0 <- ggplot(data=switch(as.character(ex0.opts$insamp.predict),"TRUE"=ex.df,"FALSE"=ex.df.pred), aes(x=time, y=ex.pred)) +
    	geom_line() + theme_bw() + scale_x_continuous("Time") +
    	scale_y_continuous("Life Expectancy at Birth")
    	if(ex0.opts$insamp.obs){
    		plot.ex0 <- plot.ex0 + geom_point(data=na.omit(ex.df),aes(x=time,y=ex.obs))
    	} 
    	plot.ex0
}

####################################################
### 3D plotter
threedimplot <- function(el,threedim.opts,sample.frame) {
  # Set up labels 
  if(is.null(threedim.opts$xlab)) {threedim.opts$xlab <- "Time"}
  if(is.null(threedim.opts$ylab)) {threedim.opts$ylab <- "Age"}
  if(is.null(threedim.opts$zlab)) {threedim.opts$zlab <- "Forecasts"}
  
  ages <- colnames(el$yhat)
  times <- rownames(el$yhat)
  sample.end <- sample.frame[2]
  
  # melt dataframes
  y.melt <- melt(el$y)
  yhat.melt <- melt(el$yhat)
  colnames(y.melt) <- colnames(yhat.melt) <-  c("x", "y", "z")	
	
  if(threedim.opts$insamp.predict){
   data <- yhat.melt
  }

  if(!threedim.opts$insamp.predict) {
  	data <- rbind(y.melt[y.melt$x<=sample.end,],yhat.melt[y.melt$x>sample.end,])
  }	
  
  # if 'unlog'=TRUE, take exp()
	if(threedim.opts$unlog){
		data$z <- exp(data$z)
	}
	
  # create lattice object
  args.wireframe <- threedim.opts$args.wireframe
  args.wireframe$x <- z ~ x * y
  args.wireframe$data <- data
  args.wireframe$xlab <- threedim.opts$xlab
  args.wireframe$ylab <- threedim.opts$ylab
  args.wireframe$zlab <- threedim.opts$zlab
  args.wireframe$screen <- threedim.opts$screen
  if(is.null(args.wireframe$shade)){args.wireframe$shade <- TRUE}
  if(is.null(args.wireframe$scales)){args.wireframe$scales <-
                                       list(arrows=F)}
  out <- do.call("wireframe",args=args.wireframe)
  out
}
