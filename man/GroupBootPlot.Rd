\name{GroupBootPlot}
\alias{GroupBootPlot}
\title{
Create a plot of bootstrapped confidence intervals for all edges of a Joint Graphical Lasso model.
}
\description{
This function plots output from bootstrapped networks computed with \link{GroupNetworkBoot}.
}
\usage{
GroupBootPlot(BootOut, GroupNames, edges.x, edges.y, 
labels = TRUE, transparency = 0.15, point.size = 1.5, line.size = 1, scales = "fixed", 
legend.position = "none", GroupNamesCheck = FALSE)
}
\arguments{
  \item{BootOut}{The output from \link{GroupNetworkBoot}.
}
  \item{GroupNames}{
A vector of optional group names that will be printed as facet labels in plot. By default, names of the networks are taken. If specified, GroupNames should match the alphabetical order of names of network groups. If unsure, you can check the matching of names by setting \code{GroupNamesCheck = TRUE}.
}
\item{edges.x}{
  If only a subset of edge combinations is of interest for the plot, this subset can be specified by setting \code{edges.x} and \code{edges.y}. Specifically, node names can be specified as vectors for \code{edges.x} and \code{edges.y} and all unique combinations of \code{edges.x} and \code{edges.y} will be plotted. For example, \code{edges.x = c("a", "b")} and \code{edges.y = "c"} will plot edges a-c and b-c but not a-b.
}
\item{edges.y}{
  See \code{edges.x}.
}
  \item{labels}{
  Logical, should edge labels be included in plots. Default is \code{labels = TRUE}.
}
\item{transparency}{
  Set \pkg{ggplot2} alpha channel (transparency) for confidence interval ribbon in plot.
}
\item{point.size}{
  Set point size. 
}
\item{line.size}{
  Set line size. 
}
\item{scales}{
  Set \pkg{ggplot2} facet scales. Default is \code{scale = "fixed"}. See ?facet_grid in ggplot2 for details.
}
\item{legend.position}{
  Define legend position to indicate colour for sample and bootstrap means. See ?theme in ggplot2.
}
\item{GroupNamesCheck}{
  Option to print match of indicated GroupNames to console. Only prints if GroupNames is specified. See \code{GroupNames} for details.
}
}

\details{
The code for the Joint Graphical Lasso procedure was adapted from the R package \pkg{JGL}. Some of the code for the cross-validation procedure was adapted from package \pkg{parcor}. Some of the code was inspired by package \pkg{qgraph}.
GroupBootPlot automatically calls \link{BootTable} to format \link{GroupNetworkBoot} output, so see \link{BootTable} for completely independent plotting. 
}
\value{
The output of GroupBootPlot returns a plot based on \pkg{ggplot2} with the bootstrapped confidence intervals of edges across groups.
}

\references{
Danaher, P., Wang, P., and Witten, D. M. (2014). The joint graphical lasso for inverse covariance estimation across multiple classes. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 76(2), 373-397. http://doi.org/10.1111/rssb.12033
}
\author{
Nils Kappelmann <n.kappelmann@gmail.com>, Giulio Costantini
}
\seealso{
\pkg{JGL}, \pkg{qgraph}, \pkg{parcor}
}
\keyword{graphs}
\keyword{multivariate}