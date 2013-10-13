\name{treemox.boot}
\alias{treemox.boot}
\title{Bootstrapping validation for PATHMOX or TECHMOX trees}
\usage{
  treemox.boot(pls, treemox, X = NULL, br = 100)
}
\arguments{
  \item{pls}{An object of class \code{"plspm"} returned by
  \code{\link{plspm}}.}

  \item{treemox}{An object of class \code{"treemox"}
  returned by either \code{\link{pathmox}} or
  \code{\link{techmox}}.}

  \item{X}{Optional dataset (matrix or data frame) used
  when argument \code{dataset=NULL} inside \code{pls}.}

  \item{br}{An integer indicating the number bootstrap
  resamples (\code{br=100} by default).}
}
\value{
  An object of class \code{"bootnodes"}. Basically a list
  with the following results:

  \item{PC}{Matrix of original path coefficients for the
  root node and the terminal nodes.}

  \item{PMB}{Matrix of bootstrap path coefficients (mean
  value) for the root node and the terminal nodes.}

  \item{PSB}{Matrix of bootstrap standard errors of path
  coefficients for the root node and the terminal nodes.}

  \item{PP05}{Matrix of 0.05 bootstrap percentile of path
  coefficients for the root node and the terminal nodes.}

  \item{PP95}{Matrix of 0.95 bootstrap percentile of path
  coefficients for the root node and the terminal nodes.}
}
\description{
  Performs bootstrapping validation on path coefficients of
  terminal nodes from a PATHMOX or TECHMOX tree
}
\details{
  The default number of re-samples is 100. However,
  \code{br} can be specified in a range from 50 to 500. \cr
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

 # Pathmox Analysis
 mob_pathmox = pathmox(mob_pls, seg_vars, signif=.10, size=.10, deep=2)

 # applying function treemox.pls
 mob_nodes_boot = treemox.boot(mob_pls, mob_pathmox)

 # plot of results for path coefficient number 12
 plot(mob_nodes_boot, pc=12)
 }
}
\author{
  Gaston Sanchez
}
\seealso{
  \code{\link{pathmox}}, \code{\link{techmox}},
  \code{\link{treemox.pls}}.
}

