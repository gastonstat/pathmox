\name{plot.treemox}
\alias{plot.treemox}
\title{Plot PATHMOX and TECHMOX trees}
\usage{
  \method{plot}{treemox} (x, root.col = "#eeeeee",
    root.bor = "#cccccc", root.txt = "#757575",
    root.cex = 0.8, root.lwd = 3, root.shadow = "gray40",
    node.col = "#feb769", node.bor = "#FE9929",
    node.txt = "#555555", node.cex = 0.7, node.lwd = 3,
    node.shadow = "gray30", leaf.col = "#93c4e5",
    leaf.bor = "#5a99c5", leaf.txt = "#555555",
    leaf.cex = 0.7, leaf.lwd = 3, leaf.shadow = "gray30",
    shadow.size = 0, arr.lwd = 3, lcol = "#ddddddbb",
    arr.col = "gray95", seg.cex = 0.7, seg.col = "#2cb0a7",
    cat.cex = 0.8, cat.col = "#555555", show.pval = TRUE,
    pval.col = "#2cb0a7", main = NULL, cex.main = 1,
    col.main = "gray50", ...)
}
\arguments{
  \item{x}{An object of class \code{"treemox"} returned by
  \code{\link{pathmox}} or \code{\link{techmox}}.}

  \item{root.col}{Fill color of root node.}

  \item{root.bor}{Border color of root node.}

  \item{root.txt}{Text color of root node.}

  \item{root.cex}{magnification to be used for text in root
  node.}

  \item{root.lwd}{Line width of border in the root node.}

  \item{root.shadow}{Color of shadow of root node.}

  \item{node.col}{Fill color of child nodes.}

  \item{node.bor}{Border color of child nodes.}

  \item{node.txt}{Text color of child nodes.}

  \item{node.cex}{magnification to be used for text in
  child nodes.}

  \item{node.lwd}{Line width of border in child nodes.}

  \item{node.shadow}{Color of shadow of child nodes.}

  \item{leaf.col}{Fill color of leaf nodes.}

  \item{leaf.bor}{Border color of leaf nodes.}

  \item{leaf.txt}{Text color of leaf nodes.}

  \item{leaf.cex}{magnification to be used for text in leaf
  nodes.}

  \item{leaf.lwd}{Line width of border in leaf nodes.}

  \item{leaf.shadow}{Color of shadow of leaf nodes.}

  \item{shadow.size}{Relative size of shadows.}

  \item{arr.lwd}{Line width of the tree branches.}

  \item{lcol}{color of lines}

  \item{arr.col}{color of arrows}

  \item{seg.cex}{A numerical value indicating the
  magnification to be used for plotting text.}

  \item{seg.col}{The color to be used for the labels of the
  segmentation variables.}

  \item{cat.cex}{magnification to be used for the
  categories}

  \item{cat.col}{The color to be used for the labels of the
  categories}

  \item{show.pval}{Logical value indicating whether the
  p-values should be plotted.}

  \item{pval.col}{The color to be used for the labels of
  the p-values.}

  \item{main}{A main title for the plot.}

  \item{cex.main}{The magnification to be used for the main
  title.}

  \item{col.main}{Color to be used for the main title}

  \item{\dots}{Further arguments are ignored.}
}
\description{
  The function \code{plot.treemox} allows to display binary
  trees of PATHMOX and TECHMOX analyses. If
  \code{shadow.size=0}, no shadows are drawn.
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
 mob_outer = list(1:5, 6:9, 10:15, 16:18, 19:21, 22:24, 25:26)
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

 # plot pathmox tree
 plot(mob_pathmox, root.col="lightblue", node.col="turquoise", leaf.col="skyblue3",
      shadow.size=0, seg.col="blue2", pval.col="magenta")
 }
}

