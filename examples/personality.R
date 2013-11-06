dmat <- function(lw,up)
    {
        sub <- subset(bfi,age>lw & age < up)
        X <- sub[,1:25]
        complete <- do.call('c',alply(X,1,function(v) !any(is.na(v))))
        XC <- as.matrix(X[complete,])
        as.matrix(dist(t(XC)))/sqrt(nrow(XC))
    }
