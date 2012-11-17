\name{treemox.pls}
\alias{treemox.pls}
\title{PLS-PM results of terminal nodes from a PATHMOX or TECHMOX tree}
\usage{
  treemox.pls(pls, treemox, X = NULL)
}
\arguments{
  \item{pls}{An object of class \code{"plspm"} returned by
  \code{\link{plspm}}.}

  \item{treemox}{An object of class \code{"treemox"}
  returned by \code{\link{pathmox}} or
  \code{\link{techmox}}.}

  \item{X}{Optional dataset (matrix or data frame) used
  when argument \code{dataset=NULL} inside \code{pls}.}
}
\value{
  An object of class \code{"treemox.pls"}. Basically a list
  with the following results:

  \item{weights}{Matrix of outer weights for each terminal
  node}

  \item{loadings}{Matrix of loadings for each terminal
  node}

  \item{paths}{Matrix of path coefficients for each
  terminal node}

  \item{r2}{Matrix of r-squared coefficients for each
  terminal node}
}
\description{
  Calculates basic PLS-PM results for the terminal nodes of
  PATHMOX and TECHMOX trees
}
\details{
  The argument \code{pls} must be the same used for
  calculating the \code{treemox} object.  When the object
  \code{pls} does not contain a data matrix (i.e.
  \code{pls$data=NULL}), the user must provide the data
  matrix or data frame in \code{X}.
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

 # applying function treemox.pls
 mob_nodes = treemox.pls(mob_pls, mob_pathmox)

 # comparative barplots
 plot(mob_nodes)
 }
}
\author{
  Gaston Sanchez
}
\seealso{
  \code{\link{pathmox}}, \code{\link{techmox}},
  \code{\link{plot.treemox}}.
}

