eye <- function(n,sparse=F)
  {
    if (sparse)
      {
        sparseMatrix(i=1:n,j=1:n,x=rep(1,n))
      }
    else
      {
        diag(rep(1,n))
      }
  }

a.n <- as.numeric

is.square <- function(m) dim(m)[1] == dim(m)[2]

center.conf <- function(XL)
            {
                llply(XL,function(X) t(scale(t(X),scale=FALSE)))
            }

Penalty.M <- function(t , periodic=FALSE)
  {
    tmp = matrix(0,t,t)
    tmp[seq(1,t-1),seq(2,t)] = diag(t-1)
    L= -1*diag(t) + tmp
    Delta = L %*% L
    Delta = Delta[-c(t-1,t),,drop=F]
    
    if (periodic)
      {
        Delta[,1] <- Delta[,1] + 1
        Delta[,t] <- Delta[,t] + 1
      }
    
    M <- t(Delta) %*% Delta
    M
  }

penalty.per.curve.list <- function(DL,XL,params)
  {
    Dist <- laply(DL, function(D) D)
    X <- laply(XL, function(X) X, .drop = FALSE)
    
    pen <- function(Dist,x,params) params$l* O(Dist, x , params)
    
    penalties <- aaply( X , 3 , function(x) pen(Dist,array(x,c(params$T,params$D,1)),params))
    
    penalties <- data.frame(index=1:dim(X)[3],penalties=penalties)
    penalties
  }

C.L <- function(DL , XL , params)
  {
    Dist <- laply(DL, function(D) D)
    Dist <- array(Dist,c(params$T,params$N,params$N))
    X <- laply(XL, function(X) X, .drop = FALSE)

    if (params$T >= 3){
      L(Dist , X , params) +  params$l * O(Dist , X , params)
    } else {
      L(Dist, X , params)
    }      
  }

L.L <- function(DL, XL, params)
  {
    N <- params$N
    T <- params$T

    Dist <- laply(DL, function(D) D)
    Dist <- array(Dist,c(T,N,N))
    X <- laply(XL, function(X) X, .drop = FALSE)
    
    D.X <- Dist.matrix(X, params) # TxNxN
    
    stopifnot(dim(D.X) == c(T,N,N))# "Ds dimension error")
    
    sum((D.X - Dist)^2)
  }

O.L <- function(DL , XL , params)
  {
    T <- params$T
    D <- params$D
    N <- params$N

    Dist <- laply(DL, function(D) D)
    X <- laply(XL, function(X) X, .drop = FALSE)
    
    #stopifnot( dim(M.DN) == c( T*D , T*D ) )# "M.large dimension error")
    #stopifnot( dim(X) == c( T , D , N ) )
    
    n <- dim(X)[3] # if called from f, this is =1 instead of N

    if (n == 1)
      {
        M <- params$M.D
        stopifnot( dim(M) == c( T*D , T*D ) )# "M.large dimension error")
      }
    else
      {
        M <- params$M.DN
        stopifnot( dim(M) == c( T*D*N , T*D*N ) )# "M.large dimension error")
      }

    Y <- array( X , c( T*D*n,1 ) )

    #dim of a.m(Y[,i]): T*Dx1

    pen <- t(Y) %*% M %*% Y #length N

    #stopifnot( length(pen.per.curve) == n  )# "error in calculating penalty of one curve" )
    pen
  }


L <- function(Dist , X , params)
  {
    N <- params$N
    T <- params$T
    
    D.X <- Dist.matrix( X , params) # TxNxN
    
    stopifnot( dim(D.X) == c(T,N,N) )# "Ds dimension error")
    
    sum( (D.X - Dist)^2 )
  }

O <- function(Dist , X , params)
  {
    T <- params$T
    D <- params$D
    N <- params$N
    
    #stopifnot( dim(M.DN) == c( T*D , T*D ) )# "M.large dimension error")
    #stopifnot( dim(X) == c( T , D , N ) )

    n <- dim(X)[3] # if called from f, this is =1 instead of N

    if (n == 1)
      {
        M <- params$M.D
        stopifnot( dim(M) == c( T*D , T*D ) )# "M.large dimension error")
      }
    else
      {
        M <- params$M.DN
        stopifnot( dim(M) == c( T*D*N , T*D*N ) )# "M.large dimension error")
      }

    Y <- array( X , c( T*D*n,1 ) )

    #dim of a.m(Y[,i]): T*Dx1

    pen <- t(Y) %*% M %*% Y #length N

    #stopifnot( length(pen.per.curve) == n  )# "error in calculating penalty of one curve" )
    pen
  }

Dist.matrix <- function(X,params)
  {
    T <- params$T
    D <- params$D
    N <- params$N
    stopifnot( dim(X) == c(T , D , N) )
    array(aaply(X,1,function(Xt) as.matrix(dist(t(matrix(Xt,D,N))))),c(T,N,N))
  }

transform.res <- function(res){
  Dist <- laply(res$DL, function(D) D)
  X <- laply(res$XL, function(X) X, .drop = FALSE)    
  D.X <- Dist.matrix(X, res$params)
  list(D = Dist, X = X, D.X = D.X)
}

x2list <- function(x,params){
  T <- params$T
  D <- params$D
  N <- params$N
  X <- array(x,c(T,D,N))
  XL <- alply(X, 1, function(x) matrix(x,D,N))
  XL
}

XL2num <- function(XL){
  a.n(laply(XL,function(x) as.numeric(x)))
}

Cost <- function(x,DL,params){
  XL <- x2list(x,params)
  C.L(DL, XL, params)
}

check.local.minimum <- function(res){
  x <- XL2num(res$XL)
  n <- length(x)
  g <- grad(Cost, x, DL = res$DL, params = res$params)
  
  if (all.equal(rep(0,n),g)){
    cat("The result is a local minimum of the cost function.\n")
  } else {
    cat("The automatic test for local minimum failed. \n You are now in a browser to check the actual value of the numerical gradient, g.\n")
    browser()
  }
}

Dist.List <- function(XL,params){
  X <- laply(XL, function(X) X, .drop = FALSE)  
  D.X <- Dist.matrix(X, params)
  DL.X <- alply(D.X,1,function(x) x)
  DL.X
}
  
