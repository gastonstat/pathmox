\name{get_pls_basic}
\alias{get_pls_basic}
\title{Basic PLS-PM algorithm}
\usage{
  get_pls_basic(DT, path_matrix, blocks, specs)
}
\arguments{
  \item{DT}{data table}

  \item{path_matrix}{inner design matrix}

  \item{blocks}{blocks of manifest variables}

  \item{specs}{list with pls algorithm specifications}
}
\description{
  Internal function. \code{get_pls_basic} is called by
  \code{pathmox}, \code{techmox}, \code{fix.pathmox},
  \code{fix.techmox}, \code{treemox.pls}, and
  \code{treemox.boot}.
}
\keyword{internal}

