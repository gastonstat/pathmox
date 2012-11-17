\name{get_xexeloa}
\alias{get_xexeloa}
\title{Core algorithm for pathmox and techmox}
\usage{
  get_xexeloa(pls, DT, EXEV, type.exev, elemnod, nv, size,
    mox)
}
\arguments{
  \item{pls}{object of class plspm}

  \item{DT}{data table}

  \item{EXEV}{data frame with segmentation variables}

  \item{type.exev}{vector with types of categorical
  variables (ordinal, nominal)}

  \item{elemnod}{element in node nv}

  \item{nv}{number of node}

  \item{size}{minimum size of elements inside a node}

  \item{mox}{string indicating type of algorithm "pathmox"
  or "techmox"}
}
\description{
  Internal function. \code{get_xexeloa} is called by
  \code{pathmox}
}
\keyword{internal}

