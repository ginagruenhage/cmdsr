\name{summary.cmds}
\alias{summary.cmds}
\title{View a summary of cmds results}
\usage{
  summary.cmds(res)
}
\arguments{
  \item{res}{The results from running \code{cmds}.}
}
\description{
  \code{summary} prints a summary of cmds results and
  statistics.
}
\details{
  After running \code{cmds}, the summary function prints
  information such as whether the algorithm converged, the
  embedding dimension, the distortion and total error of
  the embedding, the runtime, the weights that were used by
  the algorithm. It also prints the value of optional
  parameters.
}

