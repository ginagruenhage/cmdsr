\name{cmds}
\alias{cmds}
\title{cmds computes a smooth embedding in k dimensions}
\usage{
  cmds(DL, k = 1, l = 0, W = "NULL", v = FALSE,
    per = FALSE, M = "NULL", init = "average", eps = 0.01)
}
\arguments{
  \item{DL}{A list of \code{T} distance matrices of the
  same size, where each list item holds the distance matrix
  of size \code{NxN} for one timestep, where \code{N} is
  the size of the underlying dataset. The distance matrices
  should be positive, symmetric and have zero diagonal.}

  \item{k}{An integer defining the embedding dimension.
  Defaults to one.}

  \item{l}{The regularization parameter. Should be a
  positive real number.}

  \item{W}{An optional list of weight matrices. Should be
  of the same size as \code{D}. It can be used to implement
  variants of MDS.}

  \item{v}{verbose. If set to TRUE the function outputs
  information about the convergence during runtime.}

  \item{per}{periodic. If set to TRUE the penalty will be
  adjusted to enforce periodic embeddings.}

  \item{M}{An optional custom penalty matrix of size
  \code{TxT}.}

  \item{init}{The intialization method. Defaults to
  \code{average}, meaning that the algorithm is initialized
  with constant curves based on the average distance
  matrix. An alternative method is \code{random}.}

  \item{eps}{The accepted deviation from the previous
  iteration for convergence checking.}
}
\value{
  res A list with the following elements: \item{XL}{A list
  of Nxk matrices, whose rows contain the coordinates of
  the points for a given timestep.} \item{DL}{The input
  list of distance matrices.} \item{XL.init}{The initial
  configuration for the algorithm.} \item{params}{A list of
  parameters used by the algorithm.} \item{con}{A list of
  convergence characteristics}
}
\description{
  For a given set of distance matrices, \code{cmds} finds a
  smooth embedding in \code{k} dimensions.
}
\details{
  The algorithm is based on multi-dimensional scaling
  (MDS). It solves multiple MDS problems simultaneously and
  connects the solutions smoothly via a smoothing penalty.
  The result is a list of coordinates of a set of \code{N}
  points, such that at each timestep, the euclidean
  distances between the points are as close as possible to
  the distances in the distance matrix. The smoothness can
  be regulated and custom weights for the MDS cost function
  can be used. \code{cmds} provides several MDS variants
  via the parameter \code{W}. \code{W} can be set to
  \itemize{ \item \code{NULL} The algorithm uses an
  unweighted MDS cost function. \item \code{kamada-kawai}
  The weights are set to \eqn{w_{ij} = 1/(d_{ij}^2)}. \item
  A list of custom weight matrices. }
}
\examples{
res <- cmds(QuadCurves)
res <- cmds(ExpandingTriangle, k = 2, v = TRUE)
res <- cmds(QuadCurves, l = 10)
}

