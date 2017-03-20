EstimateGroupNetwork <- function(X,
                                 inputType = c("dataframe", "list.of.dataframes", "list.of.covariance.matrices"), n, covfun = covNoBessel, groupID, labels,
                                 method = c("InformationCriterion", "crossvalidation"), strategy = c("sequential", "simultaneous"),
                                 nlambda1 = 100, lambda1.min.ratio = .01, logseql1 = TRUE,
                                 nlambda2 = 100, lambda2.min.ratio = .01, logseql2 = TRUE,
                                 k = 10, seed,
                                 criterion = c("ebic", "bic", "aic"), count.unique = FALSE, gamma = .5, dec = 5, optimize = TRUE, optmethod = "CG",
                                 penalty = c("fused", "group"), weights = c("equal", "sample.size"), penalize.diagonal = FALSE, maxiter = 500, rho = 1, truncate = 1e-5,
                                 ncores = 1, simplifyOutput = TRUE)

{
  # Check input data #
  ####################
  # This section of the code checks data in input and convert different types of data to
  # a list of covariance matrices S and a vector of sample sizes n
  criterion <- match.arg(criterion)
  penalty <- match.arg(penalty)
  method <- match.arg(method)
  strategy <- match.arg(strategy)
  l2max <- 1

  if(!is.list(X)) stop('data argument must be a list of dataframes or a list of variance-covariance matrices')

  # Identify the type of input
  if(missing(inputType))
  {
    if(is.data.frame(X))
    {
      inputType <- "dataframe"
    } else if (is.list(X)) {
      dm <- sapply(X, dim)
      if(all(dm[1,] == dm[2,]))
      {
        inputType <- "list.of.covariance.matrices"
      } else if(all(sapply(X, is.data.frame))) {
        inputType <- "list.of.dataframes"
      } else {
        stop("Cannot determine the type of input X. Please specify 'X' and 'inputType' consistently. see  ?EstimateGroupNetwork")
      }
    }
  } else {
    inputType <- match.arg(inputType)
  }

  if(inputType == "list.of.covariance.matrices" & method == "crossvalidation")
  {
    stop ("To use crossvalidation, the original dataset is necessary. If you only have covariance matrices,
          select another 'method' for determining tuning parameters")
  }
  # Convert any kind of input to a list of covariance matrices S
  switch(inputType,
         "dataframe" = {
           if(missing(groupID)) stop("if 'X' is a single dataframe, the 'groupID' factor must be specified. See ?EstimateGroupNetwork")
           if(is.character(groupID))
           {
             if(!groupID %in% names(X))
             {
               stop("'groupID' must be a string corresponding to a column name in 'X' or a number corresponding to a column in 'X'. See ?EstimateGroupNetwork")
             } else {
               grp <- as.factor(X[,groupID])
               X <- split(X[, !names(X) %in% groupID], grp)
               groups <- levels(grp)
             }
           } else if(is.numeric(groupID)) {
             if(ncol(X) < groupID)
             {
               stop("'groupID' must be a string corresponding to a column name in 'X' or a number corresponding to a column in 'X'. See ?EstimateGroupNetwork")
             } else {
               grp <- as.factor(X[,groupID])
               X <- split(X[, -groupID], grp)
               groups <- levels(grp)
             }
           } else {
             stop("'groupID' must be a string corresponding to a column name in 'X' or a number corresponding to a column in 'X'. See ?EstimateGroupNetwork")
           }
           if(!missing(n)) warning("argument 'n' is ignored if the input type is a dataframe. The sample sizes are extracted from the input data.")

           S <- lapply(X, covfun)
           S <- lapply(S, cov2cor)
           n <- sapply(X, nrow)
         },

         "list.of.dataframes" = {
           if(!missing(n)) warning("argument 'n' is ignored if the input type is a list of dataframes. The sample sizes are extracted from the input data.")
           n <- sapply(X, nrow)
           ncl <- sapply(X, ncol)
           allequal <- length(unique(ncl)) == 1
           if(!allequal) stop("The input data are a list of dataframes with a different number of variables per group. Please check the input and try again. See ?EstimateGroupNetwork")

           S <- lapply(X, covfun)
           S <- lapply(S, cov2cor)
           n <- sapply(X, nrow)
           groups <- names(X)
         },

         "list.of.covariance.matrices" = {
           if(missing(n)) stop("if 'X' encodes a covariance/correlation matrix, sample sizes 'n' must be specified. See ?EstimateGroupNetwork")
           ncl <- sapply(X, ncol)
           allequal <- length(unique(ncl)) == 1
           if(!allequal) stop("The input data are a list of dataframes with a different number of variables per group. Please check the input and try again. See ?EstimateGroupNetwork")
           S <- X
           S <- lapply(S, cov2cor)
           groups <- names(X)
         })

  # number of variables
  p <- unique(sapply(S, ncol))
  if(length(p) > 1) stop("all datasets must include the same number of variables")
  K <- length(S) # number of groups
  # set weights
  wei <- weights <- match.arg(weights)
  switch(weights,
         equal = {weights = rep(1,K)},
         sample.size = {weights = n/sum(n)})
  # set argument list for myJGL
  tol <- 10^(-dec)
  myJGLarglist <-list("penalty" = penalty,
                      "weights" = weights,
                      "penalize.diagonal" = penalize.diagonal,
                      "maxiter" = maxiter,
                      "tol" = tol,
                      "rho" = rho,
                      "truncate" = truncate)

  # assign names to variables if none exist
  if(missing(labels))
  {
    if(!is.null(colnames(S[[1]])))
    {
      labels <- colnames(S[[1]])
    } else
    {
      labels <- paste0("V", seq(p))
      for(i in seq(length(S))) colnames(S[[i]]) <- rownames(S[[i]]) <- labels
    }
  } else if(length(labels) != p)
  {
    warning("The length of argument 'labels' must match that of the data in 'X'. Argument 'labels' was ignored")
    labels <- paste0("V", seq(p))
    for(i in seq(length(S))) colnames(S[[i]]) <- rownames(S[[i]]) <- labels
  }

  #####################################################
  # Perform a search of the values of lambda 1 and 2  #
  #####################################################

  ## method 1: k-fold crossvalidation
  #####################################
  if(method == "crossvalidation")
  {
    # define folds for k-fold cross validation
    if(!missing(seed)) set.seed(seed)
    all.folds <- lapply(n, function(x) split(sample(seq(x)), rep(seq(k), length = x)))

    if(strategy == "simultaneous")
    {
      # Candidate l1 and l2 values
      l1cand <- l1sequence(S, nlambda1, lambda1.min.ratio, logseql1)
      l2cand <- l2sequence(S, n, l1cand, nlambda2, lambda2.min.ratio, logseql2, l2max, ncores, myJGLarglist)

      # put all candidate lambda1 and lambda2 pairs in a matrix
      lambdas <- matrix(ncol = 2, nrow = nlambda1 * nlambda2)
      lambdas[,1] <- sort(rep(l1cand, nlambda2))
      lambdas[,2] <- l2cand

      # use function crossval to compute a matrix with the average negative predictive log likelihood function (cost function) in each split
      anpllmat <- sapply(seq(k), function(i) crossval(i,  X = X, all.folds = all.folds, K = K, ncores = ncores, lambdas = lambdas, covfun = covfun,
                                                      weights, penalty = penalty, penalize.diagonal=penalize.diagonal, maxiter=maxiter, tol=tol, rho=rho, truncate=truncate))
      # average across folds
      anpllvec <- apply(anpllmat, 1, sum)

      # the pair of lambdas values with the minimum cost function in the validation set
      la <- lambdas[which.min(anpllvec), ]
      lambda1 <- unlist(la[1])
      lambda2 <- unlist(la[2])

    } else if(strategy == "sequential")
    {
      ###################
      # select l1 first #
      ###################
      # Candidate l1 values
      l1cand <- l1sequence(S, nlambda1, lambda1.min.ratio, logseql1)
      # Select a very low value of l2 as a candidate
      l2cand <- l2sequence(S, n, l1cand[1], nlambda2, lambda2.min.ratio, logseql2, l2max, ncores, myJGLarglist)
      l2cand <- l2cand[1,]

      # put all candidate lambda1 and lambda2 pairs in a matrix
      lambdas <- matrix(ncol = 2, nrow = nlambda1)
      lambdas[,1] <- l1cand
      lambdas[,2] <- l2cand


      # use function crossval to compute a matrix with the average negative predictive log likelihood function (cost function) in each split
      anpllmat <- sapply(seq(k), function(i) crossval(i,  X = X, all.folds = all.folds, K = K, ncores = ncores, lambdas = lambdas, covfun = covfun,
                                                      weights, penalty = penalty, penalize.diagonal=penalize.diagonal, maxiter=maxiter, tol=tol, rho=rho, truncate=truncate))
      # average across folds
      anpllvec <- apply(anpllmat, 1, sum)

      # select the lambda1 associated with the minimum cost for that minimum value of l2
      lambda1 <- lambdas[which.min(anpllvec), 1]

      ##################
      # select l2 then #
      ##################
      l2cand <- l2sequence(S, n, lambda1, nlambda2, lambda2.min.ratio, logseql2, l2max, ncores, myJGLarglist)

      # put all candidate lambda1 and lambda2 pairs in a matrix
      lambdas <- matrix(ncol = 2, nrow = nlambda2)
      lambdas[,1] <- lambda1
      lambdas[,2] <- l2cand

      # use function crossval to compute a matrix with the average negative predictive log likelihood function (cost function) in each split
      anpllmat <- sapply(seq(k), function(i) crossval(i,  X = X, all.folds = all.folds, K = K, ncores = ncores, lambdas = lambdas, covfun = covfun,
                                                      weights, penalty = penalty, penalize.diagonal=penalize.diagonal, maxiter=maxiter, tol=tol, rho=rho, truncate=truncate))
      # average across folds
      anpllvec <- apply(anpllmat, 1, sum)

      # select the lambda1 associated with the minimum cost for that minimum value of l2
      lambda2 <- lambdas[which.min(anpllvec), 2]
    }
  }

  ## method 2: use an information criterion
  ##########################################
  if(method == "InformationCriterion")
  {
    if(strategy == "sequential")
    {
      ###################
      # select l1 first #
      ###################
      # Candidate l1 values
      l1cand <- l1sequence(S, nlambda1, lambda1.min.ratio, logseql1)
      # Select a very low value of l2 as a candidate
      l2cand <- l2sequence(S, n, l1cand[1], nlambda2, lambda2.min.ratio, logseql2, l2max, ncores, myJGLarglist)
      l2cand <- l2cand[1,]

      # select the value of l1 with minimal value of an information criterion among the candidates, considering a minimal value of l2
      fun1 <- function(lambda1, lambda2, S, n, myJGLarglist, criterion, count.unique, gamma, dec)
      {
        theta <- myJGL(S = S, n = n, lambda1 = lambda1, lambda2 = lambda2,
                       penalty = myJGLarglist$penalty, weights = myJGLarglist$weights,
                       penalize.diagonal = myJGLarglist$penalize.diagonal,
                       maxiter = myJGLarglist$maxiter, tol = myJGLarglist$tol, rho = myJGLarglist$rho,
                       truncate = myJGLarglist$truncate)
        theta <- theta$concentrationMatrix
        InformationCriterion(theta = theta, S = S, n = n, criterion = criterion, count.unique = count.unique, gamma = gamma, dec = dec)
      }
      if(ncores > 1)
      {
        cl <- makeCluster(ncores)
        clusterExport(cl, list("fun1", "S", "n", "l2cand", "myJGLarglist", "criterion", "count.unique", "gamma", "dec"), envir = environment())
        ICl1 <- parSapply(cl = cl, l1cand, fun1, lambda2 = l2cand, S = S, n = n, myJGLarglist = myJGLarglist, criterion = criterion, count.unique = count.unique, gamma = gamma, dec = dec)
        stopCluster(cl)
      } else {
        ICl1 <- sapply(l1cand, fun1, lambda2 = l2cand, S = S, n = n, myJGLarglist = myJGLarglist, criterion = criterion, count.unique = count.unique, gamma = gamma, dec = dec)
      }

      lambda1 <- l1cand[which.min(ICl1)]

      if(optimize)
      {
        l1candOp <- c(l1cand[1]*(1-lambda1.min.ratio ), l1cand, l1cand[nlambda1]*(1+lambda1.min.ratio))
        l1min <- l1candOp[which.min(ICl1)]
        l1max <- l1candOp[which.min(ICl1) + 2]
        # see if AIC can be further improved in the neighborhood of the chosen value
        opt <- optimize(f = fun1, interval = c(l1min, l1max), lambda2 = l2cand, S = S, n = n, myJGLarglist = myJGLarglist, criterion = criterion, count.unique = count.unique, gamma = gamma, dec = dec)

        lambda1 <- ifelse(which.min(c(min(ICl1), opt$objective)) == 1,
                          lambda1,
                          opt$minimum)
      }

      ##################
      # Select lambda2 #
      ##################
      l2cand <- l2sequence(S, n, lambda1, nlambda2, lambda2.min.ratio, logseql2, l2max, ncores, myJGLarglist)

      # select the value of l2 with minimal value of an information criterion among the candidates, considering a minimal value of l1
      fun2 <- function(lambda2, lambda1, S, n, myJGLarglist, criterion, count.unique, gamma, dec)
      {
        theta <- myJGL(S = S, n = n, lambda1 = lambda1, lambda2 = lambda2,
                       penalty = myJGLarglist$penalty, weights = myJGLarglist$weights,
                       penalize.diagonal = myJGLarglist$penalize.diagonal,
                       maxiter = myJGLarglist$maxiter, tol = myJGLarglist$tol, rho = myJGLarglist$rho,
                       truncate = myJGLarglist$truncate)
        theta <- theta$concentrationMatrix
        InformationCriterion(theta = theta, S = S, n = n, criterion = criterion, count.unique = count.unique, gamma = gamma, dec = dec)
      }

      if(ncores > 1)
      {
        cl <- makeCluster(ncores)
        clusterExport(cl, list("fun2", "S", "n", "l2cand", "lambda1", "myJGLarglist", "criterion", "count.unique", "gamma", "dec"), envir = environment())
        ICl2 <- parSapply(cl = cl, l2cand, fun2, lambda1 = lambda1, S = S, n = n, myJGLarglist = myJGLarglist, criterion = criterion, count.unique = count.unique, gamma = gamma, dec = dec)
        stopCluster(cl)
      } else {
        ICl2 <- sapply(l2cand, fun2, lambda1 = lambda1, S = S, n = n, myJGLarglist = myJGLarglist, criterion = criterion, count.unique = count.unique, gamma = gamma, dec = dec)
      }

      lambda2 <- l2cand[which.min(ICl2)]

      if(optimize)
      {
        l2candOp <- c(l2cand[1]*(1-lambda2.min.ratio ), l2cand, l2cand[nlambda1]*(1+lambda2.min.ratio))
        l2min <- l2candOp[which.min(ICl2)]
        maxl2 <- l2candOp[which.min(ICl2) + 2]
        # see if AIC can be further improved in the neighborhood of the chosen value
        opt <- optimize(f = fun2, interval = c(l2min, maxl2), lambda1 = lambda1, S = S, n = n, myJGLarglist = myJGLarglist, criterion = criterion, count.unique = count.unique, gamma = gamma, dec = dec)

        lambda2 <- ifelse(which.min(c(min(ICl2), opt$objective)) == 1,
                          lambda2,
                          opt$minimum)

    }
    } else if(strategy == "simultaneous")
    {
      # Candidate l1 and l2 values
      l1cand <- l1sequence(S, nlambda1, lambda1.min.ratio, logseql1)
      l2cand <- l2sequence(S, n, l1cand, nlambda2, lambda2.min.ratio, logseql2, l2max, ncores, myJGLarglist)

      # put all candidate lambda1 and lambda2 pairs in a matrix
      lambdas <- matrix(ncol = 2, nrow = nlambda1 * nlambda2)
      lambdas[,1] <- sort(rep(l1cand, nlambda2))
      lambdas[,2] <- l2cand

      # among each pair of l1 and l2 values, search for the one that results in the minimum value of the information criterion
      fun3 <- function(lambdas, S, n, myJGLarglist, criterion, count.unique, gamma, dec)
      {
        theta <- myJGL(S = S, n = n, lambda1 = lambdas[1], lambda2 = lambdas[2],
                       penalty = myJGLarglist$penalty, weights = myJGLarglist$weights,
                       penalize.diagonal = myJGLarglist$penalize.diagonal,
                       maxiter = myJGLarglist$maxiter, tol = myJGLarglist$tol, rho = myJGLarglist$rho,
                       truncate = myJGLarglist$truncate)
        theta <- theta$concentrationMatrix
        InformationCriterion(theta = theta, S = S, n = n, criterion = criterion, count.unique = count.unique, gamma = gamma, dec = dec)
      }

      if(ncores > 1)
      {
        cl <- makeCluster(ncores)
        clusterExport(cl, list("fun3", "S", "n", "lambdas", "myJGLarglist", "criterion", "count.unique", "gamma", "dec"), envir = environment())
        IC <- parApply(cl = cl, lambdas, 1, fun3, S = S, n = n, myJGLarglist = myJGLarglist, criterion = criterion, count.unique = count.unique, gamma = gamma, dec = dec)
        stopCluster(cl)
      } else
      {
        IC <- apply(lambdas, 1, fun3, S = S, n = n, myJGLarglis = myJGLarglist, criterion = criterion, count.unique = count.unique, gamma = gamma, dec = dec)
      }

      lambda1 <- lambdas[which.min(IC), 1]
      lambda2 <- lambdas[which.min(IC), 2]

      if(optimize)
      {
        opt <- optim(c(lambda1, lambda2), fun3, S = S, n = n, myJGLarglis = myJGLarglist, criterion = criterion, count.unique = count.unique, gamma = gamma, dec = dec,
                     method = optmethod)

        if(opt$value < min(IC))
        {
          lambda1 <- opt$par[1]
          lambda2 <- opt$par[2]
        }

        if(lambda1 < 0) lambda1 <- 0
        if(lambda2 < 0) lambda2 <- 0
      }
    }
  }

  # Put the final output together
  theta <- myJGL(S = S, n = n, lambda1 = lambda1, lambda2 = lambda2,
                 penalty = myJGLarglist$penalty, weights = myJGLarglist$weights,
                 penalize.diagonal = myJGLarglist$penalize.diagonal,
                 maxiter = myJGLarglist$maxiter, tol = myJGLarglist$tol, rho = myJGLarglist$rho,
                 truncate = myJGLarglist$truncate)

  for(i in seq(length(theta$network)))
    colnames(theta$network[[i]]) <- rownames(theta$network[[i]]) <-
    colnames(theta$concentrationMatrix[[i]]) <- rownames(theta$concentrationMatrix[[i]]) <-labels
  out <- lapply(theta$network, as.matrix)
  names(out) <- groups

  if(!simplifyOutput)
  {
    out <- list()
    out$network <- lapply(theta$network, as.matrix)
    names(out$network) <- groups
    out$concentrationMatrix <- lapply(theta$concentrationMatrix, as.matrix)
    names(out$concentrationMatrix) <- groups
    out$correlationMatrix <- lapply(S, as.matrix)
    names(out$correlationMatrix) <- groups
    out$InformationCriteria <- c("AIC" = InformationCriterion(theta = theta$concentrationMatrix, S = S, n = n, criterion = "aic", count.unique = count.unique, gamma = gamma, dec = dec),
                                 "BIC" = InformationCriterion(theta = theta$concentrationMatrix, S = S, n = n, criterion = "bic", count.unique = count.unique, gamma = gamma, dec = dec),
                                 "eBIC" = InformationCriterion(theta = theta$concentrationMatrix, S = S, n = n, criterion = "ebic", count.unique = count.unique, gamma = gamma, dec = dec),
                                 "Count Unique" = count.unique,
                                 "gamma for eBIC" = gamma)
    out$TuningParameters <- c("lambda1" = lambda1,
                              "lambda2" = lambda2)
    out$n <- n
    out$miscellaneous <- c("method" = method,
                           "strategy" = strategy,
                           "k for crossvalidation" = k,
                           "seed for crossvalidation" = ifelse(missing("seed"), "none", seed),
                           "information criterion used" = ifelse(criterion == "InformationCriterion", criterion, "none"),
                           "nlambda1" = nlambda1,
                           "lambda1.min.ratio" = lambda1.min.ratio,
                           "logseql1" = logseql1,
                           "nlambda2" = nlambda2,
                           "lambda2.min.ratio" = lambda2.min.ratio,
                           "logseql2" = logseql2,
                           "optimize" = optimize,
                           "optmethod" = optmethod,
                           "type of penalty" = ifelse(penalty == "fused", "Fused Graphical Lasso", "Group Graphical Lasso"),
                           "weights" = wei,
                           "lambda 1 penalized diagonal?" = penalize.diagonal,
                           "maxiter" = maxiter,
                           "rho" = rho,
                           "truncate" = truncate)



  }

  out
}


InformationCriterion <- function(theta, S, n, criterion = c("aic", "bic", "ebic"), count.unique = TRUE, gamma = .5, dec = 5)
{
  # theta = output of myJGL
  # n = vector of sample sizes corresponding to each element of jgl (same order)
  # S = list of covariance matrices corresponding to each element of jgl (same order)

  criterion <- match.arg(criterion)

  p <- ncol(theta[[1]])
  K <- length(theta)

  # log likelihoods
  Lik <- numeric(K)
  for (i in seq(K)){
    Lik[i] <- logGaus(S[[i]], theta[[i]], n[i])
  }
  L <- sum(Lik)

  # number of parameters
  k <- npar(theta, count.unique, dec)

  switch(criterion,
         aic = {
           out = 2 * k -2 * L
         },
         bic = {
           out = k * log(sum(n)) -2 * L
         },
         ebic = {
           out = k * log(sum(n)) -2 * L + 4 * k * gamma * log(p)
         })
  out
}


npar <- function(theta, count.unique = TRUE, dec = 5)
{
  if(count.unique == FALSE) {
    E <- sum(sapply(theta, function(x) x[lower.tri(x)] != 0))
  } else if (count.unique == TRUE)
  {
    p <- ncol(theta[[1]])
    K <- length(theta)

    # Number of unique pars:
    arr <- array(unlist(theta), dim = c(p, p, K))
    arr <- round(arr, dec)
    Nuni <- apply(arr, 1:2, function(x) {
      un <- unique(x)
      length(un[un!=0])
    })

    # Number of pars:
    E <- sum(Nuni[upper.tri(Nuni)])
  }

  E
}

logGaus <- function(S, theta, n)
{
  theta <- as.matrix(theta)
  S <- as.matrix(S)
  trThS = sum(theta * S)
  return(n/2 * (log(det(theta)) - trThS)  )
}

l1sequence <- function(S, nlambda1, lambda1.min.ratio, logseql1)
{
  # find l1_theormax, the value of l1 that would cause all edges in at least one network to be missing
  maxl1 <- function(S) max(max(S - diag(nrow(S))), -min(S - diag(nrow(S))))
  l1_theormax  <- min(sapply(S, maxl1))

  # define the sequence of l1 values
  l1min <- l1_theormax * lambda1.min.ratio
  if(logseql1 == TRUE)
  {
    l1cand <- exp(seq(log(l1min), log(l1_theormax), length = nlambda1))
  } else
  {
    l1cand <- seq(from = 0, to = l1_theormax, length.out = nlambda1)
  }
  l1cand
}


l2sequence <- function(S, n, l1cand, nlambda2, lambda2.min.ratio, logseql2, l2max, ncores, myJGLarglist)
{
  nlambda1 <- length(l1cand)
  # For each candidate l1, identify the lowest value of lambda2 that makes all the elements of the precision matix
  #  equal among different classes (max difference allowed = tol), using the stats::optimize function

  fun4 <- function(lambda2, lambda1, S, n, myJGLarglist)
  {
    #cat(paste(lambda1, lambda2, "\n"))
    tol <- myJGLarglist$tol
    maxdiff <- c()
    thetas <- myJGL(S = S, n = n, lambda2 = lambda2, lambda1 = lambda1, penalty = myJGLarglist$penalty,
                    weights = myJGLarglist$weights, penalize.diagonal = myJGLarglist$penalize.diagonal,
                    maxiter = myJGLarglist$maxiter, tol = myJGLarglist$tol, rho = myJGLarglist$rho,
                    truncate = myJGLarglist$truncate)
    thetas <- thetas$concentrationMatrix

    for(i in seq(2, length(thetas)))
    {
      for(j in seq(i-1))
      {
        maxdiff <- c(maxdiff, max(abs(thetas[[i]] - thetas[[j]])))
      }
    }
    maxdiff <- max(maxdiff)

    ifelse(maxdiff > tol, l2max+1, lambda2)
  }

  if(ncores > 1 & nlambda1 > 1)
  {
    cl <- makeCluster(ncores)
    clusterExport(cl, list("fun4", "S", "n", "l2max", "l1cand", "myJGLarglist"), envir = environment())
    l2_theormax <- parSapply(cl = cl, seq(nlambda1), function(i) optimize(f = fun4, S = S, n = n, lambda1 = l1cand[i], myJGLarglist = myJGLarglist, interval = c(0, l2max))$minimum)
    stopCluster(cl)
  } else {
    l2_theormax <-  sapply(seq(nlambda1), function(i) optimize(f = fun4, lambda1 = l1cand[i], S= S, n = n, myJGLarglist = myJGLarglist, interval = c(0, l2max))$minimum)
  }

  l2min <- l2_theormax * lambda2.min.ratio

  if(logseql2 == TRUE)
  {
    logl2min <- log(l2min)
    logl2max <- log(l2_theormax)
    l2cand <- sapply(seq_along(l2min), function(i) exp(seq(logl2min[i], logl2max[i], length = nlambda2)))
  } else
  {
    l2cand <- sapply(seq_along(l2min), function(i) seq(l2min[i], l2_theormax[i], length = nlambda2))
  }
  l2cand
}

# this way of computing the covariance matrix was taken from function JGL::JGL
covNoBessel <- function(x, ...)
{
  n <- nrow(x)
  (cov(x, ...)*(n-1))/n
}



crossval <- function(i, X, all.folds, K, ncores, lambdas, covfun, ...)
{
  # some of the code for the k-fold cross validation has been taken from the R package parcor, function adalasso

  # X = input data as a list of dataframes
  # for a given fold, function crossval computes Sigmas, the covaraince matrices of the training set for each class,
  # and Omegas, a list that includes, for each pair of lambdas, the estimated concentration matrices for each class
  # finally, using function anpll, for each pair of lambdas it computes the average negative predictive log likelihood function
  # in each training sample
  # i is the fold, an integer between 1 and k
  # ... -> arguments to be passed to myJGL
  # The function is parallelized for Windows
  omit <- lapply(all.folds, function(q) q[[i]]) # the indices of the observations in the ith fold in each class
  dat_train <- lapply(1:K, function(k) X[[k]][-omit[[k]], , drop = FALSE]) # for each class, the data after omitting the ith fold = training set
  dat_test <- lapply(1:K, function(k) X[[k]][omit[[k]], , drop = FALSE]) # for each class, the data after keeping only the ith fold = validation set
  SigmaTrain <- lapply(dat_train, covfun) # for each class, the covariance matrices in the validation set
  SigmaTest <- lapply(dat_test, covfun) # for each class, the covariance matrices in the validation set

  # for each pair of l1 and l2, compute the concentration matrix of the training set and use anpll to
  # validate it against the sample covariance matrix of the validation set
  n <- sapply(dat_train, nrow)
  if(ncores > 1)
  {
    cl <- makeCluster(ncores)
    clusterExport(cl, list("SigmaTrain", "n", "lambdas", "myJGL"), envir=environment())
    Omegas <- parLapply(cl = cl, seq(nrow(lambdas)), function(la) myJGL(S = SigmaTrain, n = n, lambda1 = lambdas[la, 1], lambda2 = lambdas[la, 2], ...)$concentrationMatrix)
    stopCluster(cl)

  } else  {
    Omegas <- lapply(seq(nrow(lambdas)), function(la) myJGL(S = SigmaTrain, n = n, lambda1 = lambdas[la, 1], lambda2 = lambdas[la, 2], ...)$concentrationMatrix)
  }

  fits <- sapply(Omegas, function(O) anpll(S = SigmaTest, O = O))
  fits
}


