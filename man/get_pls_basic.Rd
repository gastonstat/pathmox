\name{get_pls_basic}
\alias{get_pls_basic}
\title{Basic PLS-PM algorithm}
\usage{
  get_pls_basic(DT, IDM, blocks, modes, scheme, scaled,
    tol, iter)
}
\arguments{
  \item{DT}{data table}

  \item{IDM}{inner design matrix}

  \item{blocks}{blocks of manifest variables}

  \item{modes}{vector of measurement modes}

  \item{scheme}{inner weighting scheme}

  \item{scaled}{logical indicating whether to scale the
  data}

  \item{tol}{convergence threshold}

  \item{iter}{maximum number of iterations}
}
\description{
  Internal function. \code{get_pls_basic} is called by
  \code{pathmox}, \code{techmox}, \code{fix.pathmox},
  \code{fix.techmox}, \code{treemox.pls}, and
  \code{treemox.boot}. This function uses internal
  functions of the super cool package \code{plspm}
}
\keyword{internal}

