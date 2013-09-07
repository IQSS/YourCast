### function to convert mortality rates to conditional probabilities of death
mx.to.qx <- function(mx,ax,nx)
  return((nx*mx)/(1+(nx-ax)*mx))

### Helper function that calculates life tables from a matrix of mortality rates
# input is matrix of either y or yhat
# output is lifetable array
lifetableCalc <- function(mat,ax,a0,nx,dv.log=TRUE){
	# extract times and ages
	times <- as.integer(rownames(mat))
	ages <- as.integer(colnames(mat))
	N.times <- length(times)
	N.ages <- length(ages)
	
	if(is.matrix(ax)){
		if(!(ncol(ax)==N.ages) | !(nrow(ax)==N.times) ){
			stop("ax must either be a scalar, a vector of length equal to number of age intervals, or a matrix with the number of rows equal to the number of times and the number of columns equal to the number of age intervals.", call.=FALSE)
		} else{
			ax.mat <- ax
			}
	} else if(is.vector(ax) & length(ax)>1){
		if(!(length(ax)==N.ages)){
			stop("ax must either be a scalar, a vector of length equal to number of age intervals, or a matrix with the number of rows equal to the number of times and the number of columns equal to the number of age intervals.", call.=FALSE)
		} else{
			ax.mat <- matrix(ax,nrow=N.times,ncol=N.ages,byrow=TRUE)
		}
	} else{
		ax.mat <- matrix(data=ax,nrow=N.times,ncol=N.ages)
		# fix a0
		ax.mat[,1] <- rep(a0,N.times)
	}
	
	# set age interval width (nx)
	if(is.null(nx)){
		nx.vec <- c(diff(ages),NA)
	} else{
		if(length(nx)>1){
			if(length(nx)!=N.ages){
				stop("If manually set, nx must either be a scalar or a vector of length equal to number of age intervals.",call.=FALSE)
			} else{
				nx.vec <- nx
			}
		 } else{
		 	nx.vec <- rep(nx,N.ages)
		 }
		}
	
	# create T x A version matrix version of age interval width (nx)
	nx.mat <- matrix(nx.vec,nrow=N.times,ncol=N.ages,byrow=TRUE)
	
	# define mx (mortality rates)
	if(dv.log){
		mx <- exp(mat)
	} else{
	    mx <- mat	
	}
	
	
	# calculate conditional probabilities of death
	qx <- sapply(1:N.ages,function(i) mx.to.qx(mx=mx[,i],ax=ax.mat[,i],nx=nx.vec[i]))
	
	# need to set open-ended conditional probability of death to 1 and ax to 1/mx
	qx[,N.ages] <- 1
	ax.mat[,N.ages] <- 1/mx[,N.ages]
	# px is the probability of surviving from age x to x+n
	px <- 1-qx

	# lx is nuber of age surviror (where there are 100k at the beggining)	
	# lxpn is proportion of survivors at end of age interval x
	lx <- t(apply(px,MARGIN=1,function(x) cumprod(c(1, x))))*100000
	lxpn <- lx[,-1]
	
	# dx is the number of deaths in the age interval
	dx <- t(apply(lx,MARGIN=1,function(x) -diff(x)))
	
	# Lx is person-years lived by population in interval [x,x+n)
	nx.mat.noNA <- nx.mat
	nx.mat.noNA[is.na(nx.mat.noNA)] <- 0
	Lx <- lxpn*nx.mat.noNA + dx * ax.mat
	
	# Tx is person-years remaining for individuals at start of age interval x
	Tx <- t(apply(Lx,MARGIN=1,function(x) rev(cumsum(rev(x)))))
	
	# ex is remaining life expectancy at start of age interval x 
	ex <- Tx/lx[,1:N.ages]

	# create array
	vars <- c("x","nx","ax","qx","px","lx","dx","Lx","Tx","ex")
	array.out <- array(data=NA,dim=c(length(times),N.ages, length(vars)), dimnames=list(times,ages,vars))
	array.out[,,"x"] <- matrix(ages,nrow=length(times),ncol=N.ages, byrow=TRUE)
	array.out[,,"nx"] <- nx.mat
	array.out[,,"ax"] <- round(ax.mat, 4)
	array.out[,,"qx"] <- round(qx, 4)
	array.out[,,"px"] <- round(px, 4)
	array.out[,,"lx"] <- round(lx[,1:N.ages], 4)
    array.out[,,"dx"] <- round(dx, 4)
    array.out[,,"Lx"] <- round(Lx, 4)
    array.out[,,"Tx"] <- round(Tx,4) 
    array.out[,,"ex"] <- round(ex, 4)
	array.out
	}
