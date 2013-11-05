#' Continuous embeddings using cmdsr
#'
#' The package cmdsr provides functions to compute and analyze continuous embeddings. Continuous embeddings are an extension of Multi-Dimensional Scaling (MDS). They are the solution of multiple MDS problems that are solved simultaneously together with a smoothing penalty. This ensures that the individual MDS solutions are similar local minima, so that they can be connected smoothly. More information about each function can be found in the function documentation.
#'
#' Computing continuous embeddings
#'
#' The function \code{\link{cmds}} computes the continuous embedding. Its only obligatory parameter is the list of distance matrices \code{DL}. The most commonly used optional parameters are \code{k}, the dimension of the MDS embedding and \code{l}, the smoothing parameter.
#'
#' Plotting continuous embeddings
#'
#' The function \code{\link{plot.cmds}} can be used to plot the results from \code{\link{cmds}}. \code{plot.cmds(res)} plots the continuous embedding. It automatically detects the embedding dimension. Setting the optional parameter \code{shepard} to \code{TRUE} yields a Shepard plot, a common diagnostic tool for MDS. A plot of the convergence trace is computed when calling \code{plot.cmds} with \code{convergence = TRUE}. Two plots are created. The second is similar to the first, but skips the first iteration. This results in better scaling of the y-axis.
#'
#' Summarizing continuous embeddings
#'
#' The function \code{\link{summary.cmds}} gives a summary of the \code{cmds} results. 
#'
#' @useDynLib cmdsr
#' @name cmdsr
#' @docType package
NULL

#' Constant Triangle
#' 
#' A list of distance matrices representing a constant equilateral triangle at five timesteps.
#' 
#' @docType data
#' @keywords datasets
#' @format A list of five 3x3 matrices.
#' @name ConstantTriangle
NULL

#' Expanding Circle
#' 
#' A list of distance matrices representing ten points arranged on a circle with equal spacing. The circle has linearly increasing radius. There is a bit gaussian noise on the position of the points. 
#' 
#' @docType data
#' @keywords datasets
#' @format A list of ten 10x10 matrices.
#' @name ExpandingCircle
NULL

#' Expanding Triangle
#' 
#' A list of distance matrices representing three points on the corner of an equilateral triangle. The side lenghts of the triangle are linearly increasing.
#' 
#' @docType data
#' @keywords datasets
#' @format A list of five 3x3 matrices.
#' @name ExpandingTriangle
NULL

#' Lines
#' 
#' A list of distance matrices representing fifty lines with linearly increasing slopes.
#' 
#' @docType data
#' @keywords datasets
#' @format A list of ten 50x50 matrices.
#' @name Lines
NULL

#' QuadCurves
#' 
#' A list of distance matrices representing five quadratic curves with linearly increasing slopes.
#' 
#' @docType data
#' @keywords datasets
#' @format A list of eleven 5x5 matrices.
#' @name QuadCurves
NULL
