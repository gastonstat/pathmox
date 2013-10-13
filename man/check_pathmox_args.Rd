\name{check_pathmox_args}
\alias{check_pathmox_args}
\title{Check arguments for \code{pathmox} and \code{pathmox.fit}}
\usage{
  check_pathmox_args(pls, EXEV, X, signif, size, deep,
    tree)
}
\arguments{
  \item{pls}{An object of class \code{"plspm"} returned by
  \code{\link{plspm}}.}

  \item{EXEV}{A data frame of factors contaning the
  segmentation variables.}

  \item{X}{Optional dataset (matrix or data frame) used
  when argument \code{dataset=NULL} inside \code{pls}.}

  \item{signif}{A numeric value indicating the significance
  threshold of the F-statistic. Must be a decimal number
  between 0 and 1.}

  \item{size}{numeric value indicating the minimum size of
  elements inside a node.}

  \item{deep}{integer indicating the depth level of the
  tree. Must be an integer greater than 1.}

  \item{tree}{whether the tree should be displayed
  (\code{TRUE} by default).}
}
\value{
  list of validated arguments
}
\description{
  Check arguments for \code{pathmox} and \code{pathmox.fit}
}
\details{
  Internal function. \code{check_pathmox_args} is called by
  \code{pathmox}.
}
\keyword{internal}

