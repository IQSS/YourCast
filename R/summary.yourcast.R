# summary.yourcast function

summary.yourcast <- function(object, ...) {

if(!is.null(object$call)) {
#retrieving sample.frame
sample.frame <- object$aux$sample.frame

# establish how many "g", "a", "t", and "-" used in index.code
indexcheck <- unlist(strsplit(object$aux$index.code,split=""))
g <- length(grep("g",indexcheck))
a <- length(grep("a",indexcheck))
t <- length(grep("t",indexcheck))
dash <- length(grep("-",indexcheck))

# looking for groups of geo areas: create a list of the unique
# geo areas
geolist <- vector("list",length(names(object$yhat)))
  geolist <- sapply(names(object$yhat),substr,1,g)
  uniquecsid <- unique(geolist)

out <- vector("list")
out$sample.frame <- sample.frame
out$params <- object$params
out$model <- object$call$model
out$formula <- object$call$formula
out$numcs <- length(object$yhat)
out$cntry.codes <- uniquecsid

G.names <- object$aux$G.names
if(!is.null(G.names)) {
for(i in 1:length(uniquecsid)) {
out$cntry.names[i] <- G.names[grep(uniquecsid[i],
                          G.names[,1]),2]}}

if(object$call$model == "OLS") {
out$coef <- object$coeff
# Applying colnames to coefficients
for(i in 1:length(out$coef))
  {colnames(out$coef[[i]]) <- substr(names(out$coef)[i],g+1,g+a)}
}

if(object$call$model == "LC") {
out$coef <- object$coeff[[2]]
names(out$coef) <- uniquecsid
for(i in length(out$coef)) {
out$coef[[i]] <- as.matrix(out$coef[[i]])
rownames(out$coef[[i]]) <- sapply(names(object$yhat)
                      [grep(paste("^",uniquecsid[i],
                                      sep=""),
                         names(object$yhat))],substr,g+1,g+a)
colnames(out$coef[[i]]) <- "Lee-Carter Coefficients"
}
}

if(object$call$model == "map") {
out$coef <- object$coeff
}
}

else {out <- vector("list")
    out$params <- object$summary}

class(out) <- "summary.yourcast"

return(out)
}

