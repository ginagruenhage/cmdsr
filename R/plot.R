#' Plot cmds results
#'
#' This function gives a basic visualization of cmds results.
#'
#' @param res The results from running \code{cmds}.
#' @param embedding Defaults to TRUE. If nothing else is specified, only the embedding results of cmds are plotted. 
#' @param animation If set to TRUE, an animated display for \code{k=2} is created.
#' @param delay This is an optional parameter for the animation, setting the delay between two animation frames.
#' @param convergence If set to TRUE, a trace of cost values is plotted.
#' @param shepard If set to TRUE, a shapard plot is created.
#' @export
#' @examples
#' plot.cmds(cmds(QuadCurves))
#' plot.cmds(cmds(QuadCurves), convergence = TRUE)
#' plot.cmds(cmds(QuadCurves), shepard = TRUE)
#' plot.cmds(cmds(ExpandingTriangle, k = 2))

plot.cmds <- function(res, embedding = TRUE,  animation = FALSE, delay = .1, convergence = FALSE, shepard = FALSE) {
  DL <- res$DL
  XL <- res$XL
  params <- res$params
  DL.X <- Dist.List(XL,params)
  D <- params$D
  T <- params$T
  
  if (embedding){
    if (D == 1) {
      par(mfrow=c(1,1))
      lim <- max(laply(DL,max))
      plot(1:T, 1:T, type="n", ylim = c(-lim, lim), xlab="alpha", ylab="Dimension 1")
      if (T > 1){
        apply(sapply(center.conf(XL),function(v) v[1,]),1,function(v) lines(v))
      } else {
        points(rep(1,params$N),center.conf(XL)[[1]])
      }
      
    } else {
      if (animation == FALSE){
        par(mfrow = c(ceiling(T/3),3))
        lim <- max(laply(DL,max))
        XLC <- center.conf(XL)
        for (i in 1:T){
          par(pty = "s")
          plot(0, 0, type="n", ylim = c(-lim, lim), xlim = c(-lim, lim), xlab = "Dimension 1", ylab = "Dimension 2")          
          points(XLC[[i]][1,], XLC[[i]][2,], pch = 19)
        } 
      } else {
        show.anim <- function(delay=1)
          {
            lim <- max(laply(DL,max))
            XLC <- center.conf(XL)
            for (ind in 1:T)
                {
                  par(mfrow = c(1,1))
                  plot(0, 0, type="n", ylim = c(-lim, lim), xlim = c(-lim, lim), xlab = "Dimension 1", ylab = "Dimension 2")
                  points(t(XLC[[ind]]), pch = 19)
                  ani.pause(delay)
                }
          }
        show.anim(delay)
      }
    }
  }
  
  if (convergence){
    plot.new()
    par(mfrow=c(1,2))
    
    plot(res$con$trace$iter, res$con$trace$C, xlab = "Iteration", ylab = "Cost")
    plot(res$con$trace$iter[-1], res$con$trace$C[-1], xlab = "Iteration", ylab = "Cost")
  }

  if (shepard){
    plot.new()
    par(mfrow = c(ceiling(T/3),3))
    lim <- max(laply(DL,max))
    
    for (i in seq(1,T)){
      par(pty = "s")
      plot(0, 0, type="n", ylim = c(0,lim), xlim = c(0,lim), xlab = "Original Distances", ylab = "Embedding Distances")      
      abline(0,1)
      points(a.n(DL[[i]]), a.n(DL.X[[i]]))
    }
  }
    
  
}
  
