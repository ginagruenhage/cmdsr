\name{plot.cmds}
\alias{plot.cmds}
\title{Plot cmds results}
\usage{
  plot.cmds(res, embedding = TRUE, animation = FALSE,
    delay = 0.1, convergence = FALSE, shepard = FALSE)
}
\arguments{
  \item{res}{The results from running \code{cmds}.}

  \item{embedding}{Defaults to TRUE. If nothing else is
  specified, only the embedding results of cmds are
  plotted.}

  \item{animation}{If set to TRUE, an animated display for
  \code{k=2} is created.}

  \item{delay}{This is an optional parameter for the
  animation, setting the delay between two animation
  frames.}

  \item{convergence}{If set to TRUE, a trace of cost values
  is plotted.}

  \item{shepard}{If set to TRUE, a shapard plot is
  created.}
}
\description{
  This function gives a basic visualization of cmds
  results.
}
\examples{
plot.cmds(cmds(QuadCurves))
plot.cmds(cmds(QuadCurves), convergence = TRUE)
plot.cmds(cmds(QuadCurves), shepard = TRUE)
plot.cmds(cmds(ExpandingTriangle, k = 2))
}

