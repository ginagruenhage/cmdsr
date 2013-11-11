#' @title
#' cmds computes a smooth embedding in k dimensions
#'
#' @description
#' For a given set of distance matrices, \code{cmds} finds a smooth embedding in \code{k} dimensions.
#'
#' @details
#' The algorithm is based on multi-dimensional scaling (MDS). It solves multiple MDS problems simultaneously and connects the solutions smoothly via a smoothing penalty. The result is a list of coordinates of a set of \code{N} points, such that at each timestep, the euclidean distances between the points are as close as possible to the distances in the distance matrix. The smoothness can be regulated and custom weights for the MDS cost function can be used.
#' \code{cmds} provides several MDS variants via the parameter \code{W}. \code{W} can be set to
#' \itemize{
#' \item \code{NULL} The algorithm uses an unweighted MDS cost function.
#' \item \code{kamada-kawai} The weights are set to \eqn{w_{ij} = 1/(d_{ij}^2)}.
#' \item A list of custom weight matrices.
#' }
#'
#' @param DL A list of \code{T} distance matrices of the same size, where each list item holds the distance matrix of size \code{NxN} for one timestep, where \code{N} is the size of the underlying dataset. The distance matrices should be positive, symmetric and have zero diagonal.
#' @param k An integer defining the embedding dimension. Defaults to one.
#' @param l The regularization parameter. Should be a positive real number.
#' @param W An optional list of weight matrices. Should be of the same size as \code{D}. It can be used to implement variants of MDS.
#' @param v verbose. If set to TRUE the function outputs information about the convergence during runtime.
#' @param per periodic. If set to TRUE the penalty will be adjusted to enforce periodic embeddings.
#' @param M An optional custom penalty matrix of size \code{TxT}.
#' @param init The intialization method. Defaults to \code{average}, meaning that the algorithm is initialized with constant curves based on the average distance matrix. An alternative method is \code{random}.
#' @param eps The accepted deviation from the previous iteration for convergence checking.
#' @export
#' @return res A list with the following elements:
#'   \item{XL}{A list of Nxk matrices, whose rows contain the coordinates of the points for a given timestep.}
#'   \item{DL}{The input list of distance matrices.}
#'   \item{XL.init}{The initial configuration for the algorithm.}
#'   \item{params}{A list of parameters used by the algorithm.}
#'   \item{con}{A list of convergence characteristics}
#' @examples 
#' res <- cmds(QuadCurves)
#' res <- cmds(ExpandingTriangle, k = 2, v = TRUE)
#' res <- cmds(QuadCurves, l = 10)

cmds <- function(DL, k = 1, l = 0, W = "NULL", v = FALSE, per = FALSE, M = "NULL", init = "average", eps = 1e-2) {

  ## input checks on D
  if (class(DL) != "list") stop("'D' is not a list.")

  llply(DL, function(d) {
    if (!is.square(d)) stop("All elements of 'D' must be square matrices.")
    if (!isSymmetric(d)) stop("All elements of 'D' must be symmetric matrices.")
    if (any(diag(d) != 0)) stop("All elements of 'D' must have zero diagonal.")
  })

  N <- dim(DL[[1]])[1]
  T <- length(DL)

  llply(DL, function(d) {
    tmp <- aaply(d, 1, function(v) (all(v == 0) | sum(is.na(v))==(N-1) ) )
    tmp[is.na(tmp)] <- FALSE
    if (any(tmp)) warning("Some rows of 'D' are equal to zero or NA.")
    if (any(is.infinite(d))) stop("'D' contains infinite values.")
    if (any(d[!is.na(d)] < 0)) stop("'D' contains negative values.")
  })
  
  llply(DL, function(d) {
    if (dim(d)[1] != N) stop("All elements of 'D' must have the same dimensions.")
  })

  ## weight checking
  if (class(W) == "character") {
    if (W == "NULL") {
      WL <- llply(replicate(T,list(matrix(1,N,N))), function(W){
        diag(W) <- 0;
        W
      })
    }
    else if (W == "kamada-kawai") {
      WL <- llply(DL,function(d) { 
        W <- 1/(d*d);
        diag(W) <- 0;
        W
      })
    }
  } else {
    if (class(W) == "list" & length(W) == T) {
      tmp <- aaply(seq_len(T), 1, function(i) dim(DL[[i]]) == dim(W[[i]]))
      if (!all(tmp)) stop("Your custom weights should be of the same dimensions as D.")
    } else {
      stop("W must be either 'classic', 'kamada-kawai' or a custom list of weights.")
    }
  }
  
  ## check D for NA's and set the corresponding weights to zero
  WL <- alply(seq_len(T), 1, function(i) {
    if (any(is.na(DL[[i]]))) {
      warning("There are some missing values in 'D'.")
      WL[[i]][is.na(DL[[i]])] <- 0
    }
    WL[[i]]
  })

  ## row-normalize weight matrices
  WL <- llply(WL,function(W){
    sweep(W,1,rowSums(W),FUN="/")})

  ## check k
  k <- as.integer(k)
  if (k > (N-1) | k < 1) stop("'k' must be in 1,2, ..., N-1.")

  ## check l
  if (l < 0) stop("l must be positive.")

  ## check M
  if (T >=3 ){
  if (class(M) == "character") {
    if (M == "NULL") {
      M <- Penalty.M(T, periodic = per)
    } 
  } else {
    if (all(dim(as.matrix(M)) != c(T,T))) stop("The custom penalty matrix must be of size TxT.")
  }
} else { M = 1 }
  ## check init
  if (class(init) == "list") {
    tmp <- aaply(seq_len(T), 1, function(i) dim(init[[i]]) == c(N,k))
    if(!all(tmp)) stop("Your custom initialization should be a list with matrices of size Nxk.")
  } else {
    if (class(init) == "character") {
      if (init != "average" & init != "random") stop("'init' must be 'average' or 'random' or a custom list.")
      if (init == "random" & l == 0) warning("We recommend to use at least a little bit of smoothing with random initialization.")
    } else
      stop("'init' must be 'average' or 'random' or a custom list.")
  }
  
  ## set params
  max.iter = 50 ## number of outer loop iterations, multiple of 5
  
  params <- list(N = N, D = k , T = T, l = l, weighted = TRUE, WL = WL, Regress = solve(l*M + eye(T)), eps = eps, M = M, M.D = diag(k) %x% M , M.DN = diag(k*N) %x% M)
  
  ## function to compute embedding
  do.cmds <- function(DL, XL, params, max.iter){
    con <- list()
    con$trace <- data.frame(iter = c(seq(0,5),seq(10,max.iter,by=5),max.iter + 1),C=NA)
    
    for (i in seq(1,max.iter)) {
      if (i == 1) {
        e <- C.L(DL, XL, params)/N
        if (v == TRUE) cat("Initialization: Total cost C: ",e,"\n")
        con$trace[1,]$C <- e
      }      
      OuterLoop(DL, XL, params,1)
      if (i %in% c(seq(1,5),seq(10, max.iter, by=5))) {
        e <- C.L(DL, XL, params)/N
        if (v == TRUE) cat("Iteration: ", i, ", total cost C: ",e,"\n")
        con$trace[con$trace$iter == i,]$C <- e
      }
    }  
    OuterLoop(DL, XL, params,1)
    e0 <- C.L(DL, XL, params)/N
    cat("Total cost C: ",e0,"\n")
    con$trace[con$trace$iter == max.iter + 1,]$C <- e0
    list(DL = DL, XL = XL, params = params, con = con, e = e, e0 = e0)
  }

  ## Initialize
  add <- function(x) Reduce("+", x)

  if (class(init) == "character"){
    if (init == "random"){
      ## initialization is done later      
    } else if (init == "average") {
      MeanD <- add(DL)/T
      w <- matrix(1,N,N)
      w[is.na(MeanD)] <- 0
      MeanD[is.na(MeanD)] <- 0
      X <- as.matrix(smacofSym(MeanD,ndim=k)$conf,N,k)
      X <- raply(T,X,.drop=FALSE)
      XL <- alply(X,1,function(x) t(x))
      attributes(XL) <- NULL
      ## store a hard copy of the initial configuration
      XL.init <- llply(XL,function(X) X+0)

      #plot.cmds(list(DL = DL,XL = XL, params = params))
    }
  } else {
    XL <- init
    ## store a hard copy of the initial configuration
    XL.init <- llply(XL,function(X) X+0)
    
    plot.cmds(list(DL = DL,XL = XL, params = params))
  }
  
  ## compute embedding
  if (class(init) == "character"){
    if (init == "random"){
      RL <- list()
      for (t in seq(1,40)){
        XL <- alply(seq_len(T), 1, function(i) x <- matrix(rnorm(N*k,0,max(DL[[i]])),k,N))
        XL.init <- llply(XL,function(X) X+0)
        plot.cmds(list(DL = DL,XL = XL, params = params))
        RL[[t]] <- do.cmds(DL, XL.init, params, max.iter)
      }
      res <- RL[[which.min(laply(RL,function(res) res$con$trace[16,2]))]]
      res$XL.init <- XL.init
    } else {
      res <- do.cmds(DL, XL, params, max.iter)
      res$XL.init <- XL.init
    }
  } else {
    res <- do.cmds(DL, XL, params, max.iter)
    res$XL.init <- XL.init
  }
  
  #plot.cmds(res)  

  ## convergence checking
  delta <- (res$e - res$e0)
  rel.diff <- delta / res$e0

  res$con <- c(res$con,list(eps = eps, delta = delta, rel.diff = rel.diff))
  
  if ((delta >= 0) & delta < eps) {
    cat("The algorithm converged. (delta = ", delta,")\n")
  } else if (delta < 0 & abs(delta) < 1e-3) {
    cat("The algorithm converged with small increase. (delta = ", delta,", relative difference = ",rel.diff,")\n")
  } else {
    cat("The algorithm did not converge. (delta = ", delta, ", relative difference = ",rel.diff,")\n")
  }

  ## center configuration
  res$XL <- llply(res$XL,function(X) t(scale(t(X),scale=FALSE)))
  
  res
}

  
  
  
