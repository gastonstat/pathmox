\name{plspmox}
\alias{plspmox}
\title{PLS Path Model of a node from a PATHMOX or TECHMOX tree}
\usage{
  plspmox(pls, treemox, X = NULL, node, boot.val = FALSE,
    br = NULL, dataset = FALSE)
}
\arguments{
  \item{pls}{An object of class \code{"plspm"}}

  \item{treemox}{An object of class \code{"treemox"}}

  \item{X}{Optional argument for data table}

  \item{node}{An integer value indicating the number of the
  node. Must be an integer larger than 1.}

  \item{boot.val}{A logical value indicating whether
  bootstrap validation is performed (\code{FALSE} by
  default).}

  \item{br}{An integer indicating the number bootstrap
  resamples. Used only when \code{boot.val=TRUE}.}

  \item{dataset}{A logical value indicating whether the
  data matrix should be included in the list of results
  (\code{FALSE} by default).}
}
\value{
  An object of class \code{"plspm"}
}
\description{
  Calculates a PLS Path Model on a selected node from a
  pathmox or techmox tree
}
\details{
  Performs a PLS-PM analysis with the elements contained in
  \code{node}, by calling the function \code{\link{plspm}}.
  The rest of the parameters to perform the PLS-PM analysis
  (i.e. \code{inner}, \code{outer}, \code{modes},
  \code{scheme}, \code{scaled}, \code{tol}, \code{iter})
  are inherited from the object in argument \code{pls}. \cr
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

 # define inner model
 IMAG = c(0, 0, 0, 0, 0, 0, 0)
 EXPE = c(1, 0, 0, 0, 0, 0, 0)
 QUAL = c(0, 1, 0, 0, 0, 0, 0)
 VAL = c(0, 1, 1, 0, 0, 0, 0)
 SAT = c(1, 1, 1, 1, 0, 0, 0)
 COM = c(0, 0, 0, 0, 1, 0, 0)
 LOY = c(1, 0, 0, 0, 1, 1, 0)
 mob_inner = rbind(IMAG, EXPE, QUAL, VAL, SAT, COM, LOY)

 # outer model
 mob_outer = list(1:5, 6:9, 10:15, 16:18, 19:21, 22:24, 25:26)
 mob_modes = rep("A", 7)

 # apply plspm
 mob_pls = plspm(data_mobile, mob_inner, mob_outer, mob_modes, scheme="factor", scaled=FALSE)

 # re-ordering those segmentation variables with ordinal scale (Age and Education)
 csimobile$Education = factor(csimobile$Education,
     levels=c("basic","highschool","university"),
     ordered=TRUE)

 # select the segmentation variables
 seg_vars = csimobile[,1:7]

 # Pathmox Analysis
 mob_pathmox = pathmox(mob_pls, seg_vars, signif=.10, size=.10, deep=2)

 # get PLS-PM of nodes 2 and 3
 node2 = plspmox(mob_pls, mob_pathmox, node=2)
 node3 = plspmox(mob_pls, mob_pathmox, node=3)
 }
}
\author{
  Gaston Sanchez
}
\seealso{
  See Also as \code{\link{plspm}}, \code{\link{pathmox}}
}

