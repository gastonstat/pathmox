\name{techmox}
\alias{fix.techmox}
\alias{techmox}
\title{TECHMOX Algorithm: Segmentation Trees in PLS Path Modeling}
\usage{
  techmox(pls, EXEV, X = NULL, signif = 0.05, size = 0.1,
    deep = 2, tree = TRUE)
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

  \item{size}{A numeric value indicating the minimum size
  of elements inside a node.}

  \item{deep}{An integer indicating the depth level of the
  tree. Must be an integer greater than 1.}

  \item{tree}{A logical value indicating if the tree should
  be displayed (\code{TRUE} by default).}
}
\value{
  An object of class \code{"treemox"}. Basically a list
  with the following results:

  \item{MOX}{Data frame containing the results of the
  segmentation tree.}

  \item{FT}{Data frame containing the results of the
  F-tests for each node partition.}

  \item{candidates}{List of data frames containing the
  candidate splits of each node partition.}

  \item{list.nodes}{List of elements for each node.}
}
\description{
  The function \code{techmox} calculates a binary
  segmentation tree following the TECHMOX algorithm.  In
  contrast, \code{fix.techmox} obtains a supervised TECHMOX
  tree in the sense of allowing the user to interactively
  fix the partitions along the construction process of the
  tree. \cr
}
\details{
  The argument \code{EXEV} must be a data frame containing
  segmentation variables as factors (see
  \code{\link{factor}}). The number of rows in \code{EXEV}
  must be the same as the number of rows in the data used
  in \code{pls}. \cr

  The argument \code{size} can be defined as a decimal
  value (i.e. proportion of elements inside a node), or as
  an integer (i.e. number of elements inside a node). \cr

  When the object \code{pls} does not contain a data matrix
  (i.e. \code{pls$data=NULL}), the user must provide the
  data matrix or data frame in \code{X}.
}
\examples{
\dontrun{
 ## example of PLS-PM in customer satisfaction analysis
 ## model with seven LVs and reflective indicators
 data(csimobile)

 # select manifest variables
 data_mobile = csimobile[,8:33]

 # define path matrix (inner model)
 IMAG = c(0, 0, 0, 0, 0, 0, 0)
 EXPE = c(1, 0, 0, 0, 0, 0, 0)
 QUAL = c(0, 1, 0, 0, 0, 0, 0)
 VAL = c(0, 1, 1, 0, 0, 0, 0)
 SAT = c(1, 1, 1, 1, 0, 0, 0)
 COM = c(0, 0, 0, 0, 1, 0, 0)
 LOY = c(1, 0, 0, 0, 1, 1, 0)
 mob_path = rbind(IMAG, EXPE, QUAL, VAL, SAT, COM, LOY)

 # blocks of indicators (outer model)
 mob_blocks = list(1:5, 6:9, 10:15, 16:18, 19:21, 22:24, 25:26)
 mob_modes = rep("A", 7)

 # apply plspm
 mob_pls = plspm(data_mobile, mob_path, mob_blocks, modes = mob_modes,
                 scheme = "factor", scaled = FALSE)

 # re-ordering those segmentation variables with ordinal scale
 # (Age and Education)
 csimobile$Education = factor(csimobile$Education,
     levels=c("basic","highschool","university"),
     ordered=TRUE)

 # select the segmentation variables
 seg_vars = csimobile[,1:7]

 # Techmox Analysis
 mob_techmox = techmox(mob_pls, seg_vars, signif=.10, size=.10, deep=2)
 }
}
\author{
  Gaston Sanchez
}
\references{
  Sanchez, G. (2009) \emph{PATHMOX Approach: Segmentation
  Trees in Partial Least Squares Path Modeling.} Doctoral
  Dissertation.

  \url{http://www.gastonsanchez.com/Pathmox_Approach_Thesis_Gaston_Sanchez.pdf}
}
\seealso{
  \code{\link{pathmox}}, \code{\link{plot.treemox}}
  \code{\link{treemox.pls}}.
}

