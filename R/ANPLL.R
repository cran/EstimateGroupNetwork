
# negative predictive log likelihood for one pairs of S and O
# notice that the sums of the elementwise products of S*O is the same as the trace of the matrix produce S%*%O, but it requires much less computational time
npll <- function(S, O) sum(S*O) - log(det(O))

# average negative predictive log likelihood function, as in Guo et al. (2011), Joint estimation of multiple graphical models. Biometrika, 98(1): 1-15. doi:10.1093/biomet/asq060
anpll <- function(S, O)   
{
  nvar <- ncol(S[[1]])
  nclasses <- length(S)
  if(length(S) != length(O))
    stop("S and O must be of the same length")
  if (!all(sapply(S, is.square.matrix))) 
    stop("each element in S must be a square matrix")
  if (!all(sapply(O, is.square.matrix))) 
    stop("all elements in O and S must be square matrices of the same size")
  if (!all(c(sapply(O, dim), sapply(S, dim)) == nvar))
    stop("all elements in O and S must be square matrices of the same size")
  # return the sum of the negative predictive log likelihood for all pairs in S and O
  pred_neg_log_lik <- sapply(1:nclasses, function(cls) npll(S[[cls]], O[[cls]]))
  sum(pred_neg_log_lik)
}

is.square.matrix <- function(m)
{
  is.matrix(m) & dim(m)[1] == dim(m)[2]
}

