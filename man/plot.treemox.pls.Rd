\name{plot.treemox.pls}
\alias{plot.treemox.pls}
\title{Comparative plot between nodes from a PATHMOX or TECHMOX tree}
\usage{
  \method{plot}{treemox.pls} (x, comp.by = "nodes",
    nodes.names = NULL, ordered = TRUE, decreasing = FALSE,
    color = NULL, show.box = TRUE, border = NA,
    cex.names = 0.75, cex.axis = 0.75, short.labs = TRUE,
    short.min = NULL, ...)
}
\arguments{
  \item{x}{An object of class \code{"treemox.pls"} returned
  by \code{\link{treemox.pls}}.}

  \item{comp.by}{One of "nodes" or "latents". This argument
  indicates the type of barplots comparison.}

  \item{nodes.names}{Optional vector of names for the
  terminal nodes (must be a vector of length equal to the
  number of terminal nodes).}

  \item{ordered}{A logical value indicating whether the
  barplots are shown in increasing ordered.}

  \item{decreasing}{A logical value indicating if the sort
  order should be increasing or decreasing.}

  \item{color}{Optional vector of colors for the bars. When
  \code{color=NULL} rainbow colors are used.}

  \item{show.box}{A logical value indicating whether a box
  is drawn around each barplot.}

  \item{border}{The color to be used for the border of the
  bars. Use \code{border=NA} to omit borders.}

  \item{cex.names}{Expansion factor for axis names (bar
  labels).}

  \item{cex.axis}{Expansion factor for numeric axis
  labels.}

  \item{short.labs}{Logical value indicating if the labels
  of the barplots should be abbreviated (\code{TRUE} by
  default).}

  \item{short.min}{Integer number indicating the minimum
  length of the abbreviations for the labels. Only used
  when \code{short.labs=TRUE}.}

  \item{\dots}{Arguments to be passed to/from other
  methods.}
}
\description{
  Plot method for objects of class \code{"treemox.pls"}.
  Barplots of path coefficients of terminal nodes with
  respect to those of the global (root) model
}
\details{
  This function aims to visualize the comparison between
  path coefficients of the terminal nodes against the path
  coefficients of the global model in the root node. \cr
  When \code{comp.by="nodes"} a graphic window is displayed
  for each endogenous latent variable of the PLS model, and
  barplots of nodes are shown. \cr When
  \code{comp.by="latents"} a graphic window is displayed
  for each endogenous relationship of the PLS model, and
  barplots of independent latent variables are shown.
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

 # re-ordering those segmentation variables with ordinal scale (Age and Education)
 csimobile$Education = factor(csimobile$Education,
     levels=c("basic","highschool","university"),
     ordered=TRUE)

 # select the segmentation variables
 seg_vars = csimobile[,1:7]

 # Pathmox Analysis
 mob_pathmox = pathmox(mob_pls, seg_vars, signif=.10, size=.10, deep=2)

 # applying function treemox.pls
 mob_nodes <- treemox.pls(mob_pls, mob_pathmox)

 # default plot
 # comparative barplots of endogenous latent variables between nodes
 plot(mob_nodes, comp.by="nodes")

 # comparative barplots of nodes between latent regressors
 plot(mob_nodes, comp.by="latents", decreasing=TRUE)
 }
}
\seealso{
  \code{\link{treemox.pls}}, \code{\link{pathmox}},
  \code{\link{techmox}}
}

