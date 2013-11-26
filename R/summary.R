#' View a summary of cmds results
#'
#' \code{summary} prints a summary of cmds results and statistics.
#'
#' After running \code{cmds}, the summary function prints information such as whether the algorithm converged, the embedding dimension, the distortion and total error of the embedding and the runtime. It also prints the distortion per timestep. 
#'
#' @param res The results from running \code{cmds}.
#'
#' @export

summary.cmds <- function(res) {

  DL <- res$DL
  XL <- res$XL
  params <- res$params
  con <- res$con

  ans <- list()
  
  ans$Error <- a.n(C.L(DL,XL,params) / params$sum.D)
  ans$Distortion <- a.n(L.L(DL,XL,params) / params$sum.D)
  ans$Penalty <- a.n(ans$Error-ans$Distortion)

  DL.X <- Dist.List(XL,params)

  ans$Distortion_per_timestep <- aaply(1:params$T,1, function(i){
    sum((DL.X[[i]]-DL[[i]])^2) / sum(DL[[i]]^2) })
  
  con$trace$rate <- aaply(seq(1,dim(con$trace)[1]), 1, function(i){
    if (i==1) Inf
    else (con$trace[i-1,]$C - con$trace[i,]$C)
  })

  sub <- subset(con$trace, rate < 1e-6)
  
  ans$Convergence <- cat("The algorithm converged after", sub$iter[1], "iterations.\n")

  ans
}
