\name{GroupNetworkBoot}
\alias{GroupNetworkBoot}
\title{
Compute bootstrap networks for a Joint Graphical Lasso model on data collected on observations from different groups.
}
\description{
This bootstrapping function resamples from intial dataframes to compute bootstrapping intervals for edges estimated using EstimateGroupNetwork.
}
\usage{
                     
GroupNetworkBoot(data_list, groupNetwork, nboots = 100, bootSeed, ...)                  
                     
}
\arguments{
 \item{data_list}{A list of dataframes, one by group. Each dataframe must be structured in the same way (the same variables for each group). This needs to be the same input as was used for the original Joint Graphical Lasso network estimated with \code{EstimateGroupNetwork}.
}
  \item{groupNetwork}{
The to-be-bootstrapped network estimated with the \code{EstimateGroupNetwork} function. Importantly, the initial Joint Graphical Lasso needs to be estimated with \code{simplifyOutput = FALSE}.
}
  \item{nboots}{
The number of bootstraps to-be-conducted.
}
  \item{bootSeed}{
An optional random seed for ensuring replicability of the results.
}
  \item{...}{
  All further arguments need to be specified as done for the initial computation of the \code{EstimateGroupNetwork} function. Here all arguments apply and have the default values of function \code{EstimateGroupNetwork}, with the exceptions being the arguments \code{inputType = "list.of.dataframes"}, \code{simplifyOutput = FALSE}, and  \code{labels}, as node labels are taken directly from the original network. These arguments are set by default.
}
}

\details{
Some of the code for the cross-validation procedure was adapted from package \pkg{parcor}. Some of the code was inspired by package \pkg{qgraph}.
}
\value{
The output of GroupNetworkBoot returns a list with the following elements:
\item{data}{The original list of dataframes supplied to the function}
\item{sample}{A list including the original output from \code{EstimateGroupNetwork}}
\item{boot}{A list of matrices, each including a bootstrapped network}
}

\references{
Danaher, P (2013). JGL: Performs the Joint Graphical Lasso for sparse inverse covariance estimation on multiple classes. R package version 2.3. https://CRAN.R-project.org/package=JGL
Danaher, P., Wang, P., and Witten, D. M. (2014). The joint graphical lasso for inverse covariance estimation across multiple classes. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 76(2), 373-397. http://doi.org/10.1111/rssb.12033
}
\author{
Nils Kappelmann <n.kappelmann@gmail.com>, Giulio Costantini
}
\seealso{
\pkg{JGL}, \pkg{qgraph}
}
\examples{

\dontrun{

## Load packages:
library("psych")
library("EstimateGroupNetwork")

# In this example, the BFI network of males and females are compared for the subset of
# Agreeableness items
# Load BFI data
data(bfi)

## The bfi data is subset to Agreeableness items only for the first 500 individuals to decrease
# computational time
bfi <- bfi[, c(paste("A", 1:5, sep = ""), "gender")]

# remove observations with missing values on items or gender
bfi <- na.omit(bfi)

# Create list split by gender
bfi_list <- list(males = bfi[bfi$gender == 1, 1:5],
                 females = bfi[bfi$gender == 2, 1:5])

# Estimate JGL:
bfi_net <- EstimateGroupNetwork(bfi_list, inputType = "list.of.dataframes", simplifyOutput = FALSE)

# Bootstrap network 10 times (this will take a few minutes)
boot_bfi_net <- GroupNetworkBoot(data_list = bfi_list, groupNetwork = bfi_net, 
                                 nboots = 10, bootSeed = 1234, ncores = 1)

# use BootTable to obtain a table with information for each boostrapped edge
BootTable(boot_bfi_net)

## Use GroupBootPlot to obtain plots as a list with each group plot as one element
GroupBootPlot(boot_bfi_net) 

## Get plot for a subset of edges (here: all edges including A1). Also check Groupnames
GroupBootPlot(boot_bfi_net, edges.x = "A1", edges.y = c("A2", "A3", "A4", "A5"), 
              GroupNames = c("Females", "Males"), GroupNamesCheck = TRUE, legend.position = "top")

}
}
\keyword{graphs}
\keyword{multivariate}