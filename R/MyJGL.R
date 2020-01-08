myJGL <-
  function(S, n, weights, penalty = c("fused", "group"), lambda1, lambda2, warm, penalize.diagonal=FALSE, maxiter=500, tol=1e-5, rho=1, truncate=1e-5)
  {
    p <- ncol(S[[1]])
    K <- length(S)

    connected <- rep(TRUE,p)
    lam1 <- lambda1
    lam2 <- lambda2
    ## examine criteria:  (value 0 where S allows theta=0 to satisfy KKT; value 1 where theta must be connected)
    if(penalty=="fused")
    {
      if(K==2)  #use bi-conditional screening rule to identify block structure exactly
      {
        crit1 = list()
        for(k in 1:K) { crit1[[k]] =  abs(S[[k]])*weights[k] > lam1 + lam2 }
        S.sum = matrix(0,sum(connected),sum(connected))
        for(k in 1:K) {S.sum = S.sum + weights[k]*S[[k]]}
        S.sum = abs(S.sum)
        crit2 = S.sum > 2*lam1
      }

      if(K>2)  #use sufficient screening rule to identify larger-grained block structure
      {
        crit1 = list()
        for(k in 1:K) { crit1[[k]] =  abs(S[[k]])*weights[k] > lam1 }
        crit2 = matrix(0,sum(connected),sum(connected))
      }

      # are both criteria met?
      critboth = crit2
      for(k in 1:K) {critboth = critboth + crit1[[k]]}
      critboth = (critboth!=0)
      diag(critboth) = 1
    }

    if(penalty=="group")
    {
      ## examine criteria:  (value 0 where S allows theta=0 to satisfy KKT; value 1 where theta must be connected)
      tempsum = matrix(0,sum(connected),sum(connected))
      for(k in 1:K) {tempsum = tempsum + (pmax(weights[k]*abs(S[[k]]) - lam1,0))^2 }
      critboth = tempsum > lam2^2
      diag(critboth) = 1
    }

    ## now identify block structure using igraph:
    g1 <- graph.adjacency(critboth)
    cout = clusters(g1)
    blocklist = list()
    # identify unconnected elements, and get blocks:
    unconnected = c()

    # adapt cout$membership to start with index 1:
    if(min(cout$membership)==0){cout$membership=cout$membership+1}
    for(i in 1:(cout$no))
    {
      if(sum(cout$membership==i)==1) { unconnected <- c(unconnected,which(cout$membership==i)) }
      if(sum(cout$membership==i)>1) { blocklist[[length(blocklist)+1]] <- which(cout$membership==i) }
    }

    # final set of connected nodes
    connected[unconnected] = FALSE
    # connected indices of connected nodes:  0 for unconnected nodes, and 1:length(connected) for the rest.
    # maps features 1:p to their order in the connected features
    connected.index = rep(0,p)
    connected.index[connected] = 1:sum(connected)
    # regular indices of connected nodes: map connected nodes onto 1:p indexing:

    # redefine unconnected as !connected (up until now it's been extra nodes caught as unconnected)
    unconnected=!connected

    ## define theta on all connected:   (so theta is really theta.connected).
    theta = list()
    for(k in 1:K)
    {
      theta[[k]] = matrix(0,sum(connected),sum(connected))
      if(sum(connected)>0)
      {
        dimnames(theta[[k]])[[1]]=dimnames(theta[[k]])[[2]]=dimnames(S[[k]])[[2]][connected]
      }
    }

    ## get solution on unconnected nodes
    Su <- lapply(S, function(x) unname(diag(x)[unconnected]))
    # penalty vectors:
    # note: for admm.iters.unconnected, we use the penalize.diagonal argument before calling the function.  for admm.iters, we use it IN the function.
    if(length(lambda1)==1) { lam1.unconnected = lambda1 }
    if(length(lambda1)>1) { lam1.unconnected = diag(lambda1)[unconnected] }
    if(length(lambda2)==1) { lam2.unconnected = lambda2 }
    if(length(lambda2)>1) { lam2.unconnected = diag(lambda2)[unconnected] }
    # if penalize.diagonal==FALSE, then set the appropriate penalty vectors to zero:
    if(!penalize.diagonal)
    {
      lam1.unconnected = lam1.unconnected * 0
      if(penalty=="group") {lam2.unconnected = lam2.unconnected * 0}
    }

    # get the unconnected portion of theta:
    if(sum(unconnected)>0)
    {
      theta.unconnected = myadmm.iters.unconnected(S = Su, lambda1=lam1.unconnected,lambda2=lam2.unconnected,penalty=penalty,rho=rho,weights=weights,maxiter=maxiter,tol=tol)$Z
      for(k in 1:K) { names(theta.unconnected[[k]])=dimnames(S[[k]])[[2]][!connected] }
    }
    if(sum(unconnected)==0) {theta.unconnected = NULL}

    ## now run JGL on each block of the connected nodes to fill in theta:
    if(length(blocklist)>0){
      for(i in 1:length(blocklist)){
        # the variables in the block
        bl <- blocklist[[i]]
        Sbl = list()
        # get the data on only those variables
        for(k in 1:K)
        {
          Sbl[[k]] = S[[k]][bl, bl]
        }
        # penalty matrices:
        if(length(lambda1)==1) { lam1.bl = lambda1 }
        if(length(lambda1)>1) { lam1.bl = lambda1[bl,bl] }
        if(length(lambda2)==1) { lam2.bl = lambda2 }
        if(length(lambda2)>1) { lam2.bl = lambda2[bl,bl] }
        # initialize lambdas:
        lam1.bl = penalty.as.matrix(lam1.bl,dim(Sbl[[1]])[2],penalize.diagonal=penalize.diagonal)
        if(penalty=="fused") {lam2.bl = penalty.as.matrix(lam2.bl,dim(Sbl[[1]])[2],penalize.diagonal=TRUE)}
        if(penalty=="group") {lam2.bl = penalty.as.matrix(lam2.bl,dim(Sbl[[1]])[2],penalize.diagonal=penalize.diagonal)}

        # implement warm start if desired
        if(missing(warm)) {warm.bl = NULL
        } else {
          warm.bl = list()
          for(k in 1:K) { warm.bl[[k]] = warm[[k]][bl,bl] }
        }
        # run JGL on the block:
        Thetabl = myadmm.iters(Sbl, lam1.bl,lam2.bl,penalty=penalty,rho=rho,weights=weights,penalize.diagonal=TRUE,maxiter=maxiter,tol=tol,warm=warm.bl)
        # update Theta with Thetabl's results:
        for(k in 1:K) {theta[[k]][connected.index[bl],connected.index[bl]] = Thetabl$Z[[k]]}
      }}

     # round very small theta entries down to zero:
    if(dim(theta[[1]])[1]>0)
    {
      for(k in 1:K)
      {
        rounddown = abs(theta[[k]])<truncate; diag(rounddown)=FALSE
        theta[[k]]=theta[[k]]*(1-rounddown)
      }}

    # return output: theta on connected nodes, diagonal theta on unconnected nodes, and the identities of the connected nodes
    whole.theta = list()
      for(k in 1:K)
      {
        whole.theta[[k]] = matrix(0,p,p)
        diag(whole.theta[[k]])[unconnected] = theta.unconnected[[k]]
        whole.theta[[k]][connected,connected] = theta[[k]]
        dimnames(whole.theta[[k]])[[1]] = dimnames(whole.theta[[k]])[[2]] = dimnames(S[[k]])[[2]]
      }

    theta <- whole.theta

    list("network" = lapply(theta, wi2net), "concentrationMatrix" = theta)
  }



myadmm.iters.unconnected = function(S, lambda1,lambda2,penalty="fused",rho=1,rho.increment=1,weights,maxiter = 1000,tol=1e-5,warm=NULL)
{
  K = length(S)
  p = length(S[[1]])
  n = weights

  # initialize theta:
  theta = list()
  for(k in 1:K){theta[[k]] = 1/S[[k]]}
  # initialize Z:
  Z = list(); for(k in 1:K){Z[[k]] = rep(0,p)}
  # initialize W:
  W = list();	for(k in 1:K){W[[k]] = rep(0,p)}

  # initialize lambdas:
  lam1 = lambda1
  lam2 = lambda2

  # iterations:
  iter=0
  diff_value = 10
  while((iter==0) || (iter<maxiter && diff_value > tol))
  {
    # update theta to minimize -logdet(theta) + <S,theta> + rho/2||theta - Z + W ||^2_2:
    theta.prev = theta
    for(k in 1:K)
    {
      B = n[k]*S[[k]] - rho*(Z[[k]] - W[[k]])
      theta[[k]] = 1/(2*rho) * ( -B + sqrt(B^2 + 4*rho*n[k]) )
      #		B = S[[k]] - rho/n[k]*(Z[[k]] - W[[k]])
      #		theta[[k]] = n[k]/(2*rho) * ( -B + sqrt(B^2 + 4*rho/n[k]) )  # written as in our paper
    }

    # update Z:
    # define A matrices:
    A = list()
    for(k in 1:K){ A[[k]] = theta[[k]] + W[[k]] }
    if(penalty=="fused")
    {
      # use flsa to minimize rho/2 ||Z-A||_F^2 + P(Z):
      if(K==2){Z = flsa2(A,rho,lam1,lam2,penalize.diagonal=TRUE)}
      if(K>2){Z = flsa.general(A,rho,lam1,lam2,penalize.diagonal=TRUE)}
    }
    if(penalty=="group")
    {
      # minimize rho/2 ||Z-A||_F^2 + P(Z):
      Z = dsgl(A,rho,lam1,lam2,penalize.diagonal=TRUE)
    }

    # update the dual variable W:
    for(k in 1:K){W[[k]] = W[[k]] + (theta[[k]]-Z[[k]])}

    # bookkeeping:
    iter = iter+1
    diff_value = 0
    for(k in 1:K) {diff_value = diff_value + sum(abs(theta[[k]] - theta.prev[[k]])) / sum(abs(theta.prev[[k]]))}
    # increment rho by a constant factor:
    rho = rho * rho.increment
  }
  diff = 0; for(k in 1:K){diff = diff + sum(abs(theta[[k]]-Z[[k]]))}
  out = list(theta=theta,Z=Z,diff=diff,iters=iter)
  return(out)
}


### ADMM for FGL:
### ADMM for FGL:
myadmm.iters = function(S, lambda1,lambda2,penalty="fused",rho=1,rho.increment=1,weights,penalize.diagonal,maxiter = 1000,tol=1e-5,warm=NULL)
{
  K = length(S)
  p = ncol(S[[1]])
  n=weights
  
  warned = FALSE
  
  # initialize theta:
  theta = list()
  for(k in 1:K){theta[[k]] = diag(1/diag(S[[k]]))}
  # initialize Z:
  Z = list(); for(k in 1:K){Z[[k]]=matrix(0,p,p)}
  # initialize W:
  W = list();	for(k in 1:K) {W[[k]] = matrix(0,p,p) }

  # initialize lambdas:  (shouldn't need to do this if the function is called from the main wrapper function, JGL)
  lam1 = penalty.as.matrix(lambda1,p,penalize.diagonal=penalize.diagonal)
  if(penalty=="fused") {lam2 = penalty.as.matrix(lambda2,p,penalize.diagonal=TRUE)}
  if(penalty=="group") {lam2 = penalty.as.matrix(lambda2,p,penalize.diagonal=penalize.diagonal)}
  # iterations:
  iter=0
  diff_value = 10

  while((iter==0) || (iter<maxiter && diff_value > tol))
  {
    # reporting
    #	if(iter%%10==0)
    # if(FALSE)
    # {
    #   print(paste("iter=",iter))
    #   if(penalty=="fused")
    #   {
    #     print(paste("crit=",crit(theta,S,n=rep(1,K),lam1,lam2,penalize.diagonal=penalize.diagonal)))
    #     print(paste("crit=",crit(Z,S,n=rep(1,K),lam1,lam2,penalize.diagonal=penalize.diagonal)))
    #   }
    #   if(penalty=="group"){print(paste("crit=",gcrit(theta,S,n=rep(1,K),lam1,lam2,penalize.diagonal=penalize.diagonal)))}
    # }

    # update theta:
    theta.prev = theta

    for(k in 1:K){
      edecomp = eigen(S[[k]] - rho*Z[[k]]/n[k] + rho*W[[k]]/n[k])
      D = edecomp$values
      V = edecomp$vectors
      D2 = n[k]/(2*rho) * ( -D + sqrt(D^2 + 4*rho/n[k]) )
      theta[[k]] = V %*% diag(D2) %*% t(V)
    }

    # update Z:
    # define A matrices:
    A = list()
    for(k in 1:K){ A[[k]] = theta[[k]] + W[[k]] }
    
    # added by GC to handle complex eigenvalues - thanks to Jonas Haslbeck for finding this error
    if(any(sapply(A, is.complex))) 
    {
      A <- lapply(A, Re)
      
      if(!warned)
      {
      warning(paste0("complex eigenvalues were returned during the ADMM optimization for the following values of lambda:",
                     "
                     lambda1 = ", 
                     lambda1[1,2],
                     "
                     lambda2 = ",
                     lambda2[1,2],
                     ".
                     Set simplifyOutput = FALSE to know the selected values of lambda (see TuningParameters). If these values of lambda were selected, results might be unreliable. Solutions could include: changing the parameter covfun (es. covfun = cor) or changing the lambda parameters (arguments nlambda1, lambda1.min.ratio, nlambda2, lambda2.min.ratio, logseql1, logseql2) to a different value"))
      warned <- TRUE
        }
      
    }
    ## end GC mod
    
    if(penalty=="fused")
    {
      # use flsa to minimize rho/2 ||Z-A||_F^2 + P(Z):
      if(K==2){Z = flsa2(A,rho,lam1,lam2,penalize.diagonal=TRUE)}
      if(K>2){Z = flsa.general(A,rho,lam1,lam2,penalize.diagonal=TRUE)}  # the option to not penalize the diagonal is exercised when we initialize the lambda matrices
    }

    if(penalty=="group")
    {
      #  minimize rho/2 ||Z-A||_F^2 + P(Z):
      Z = dsgl(A,rho,lam1,lam2,penalize.diagonal=TRUE)
    }

    # update the dual variable W:
    for(k in 1:K){W[[k]] = W[[k]] + (theta[[k]]-Z[[k]])}

    # bookkeeping:
    iter = iter+1
    diff_value = 0
    for(k in 1:K) {diff_value = diff_value + sum(abs(theta[[k]] - theta.prev[[k]])) / sum(abs(theta.prev[[k]]))}
    # increment rho by a constant factor:
    rho = rho*rho.increment
  }

  
  
  
  diff = 0; for(k in 1:K){diff = diff + sum(abs(theta[[k]]-Z[[k]]))}
  out = list(theta=theta,Z=Z,diff=diff,iters=iter)
  
  return(out)
}


flsa2 <-
  function(A,L,lam1,lam2,penalize.diagonal)  #A is a list of 2 matrices from which we apply an L2 penalty to departures
  {
    # 1st apply fused penalty:
    S1 = abs(A[[1]]-A[[2]])<=2*lam2/L
    X1 = (A[[1]]+A[[2]])/2
    Y1 = X1

    S2 = (A[[1]] > A[[2]]+2*lam2/L)
    X2 = A[[1]] - lam2/L
    Y2 = A[[2]] + lam2/L

    S3 = (A[[2]] > A[[1]]+2*lam2/L)
    X3 = A[[1]] + lam2/L
    Y3 = A[[2]] - lam2/L

    X = soft(a = S1*X1 + S2*X2 + S3*X3, lam = lam1/L, penalize.diagonal=penalize.diagonal)
    Y = soft(a = S1*Y1 + S2*Y2 + S3*Y3, lam = lam1/L, penalize.diagonal=penalize.diagonal)

    return(list(X,Y))
  }




flsa.general <-
  function(A,L,lam1,lam2,penalize.diagonal)
  {
    trueA = A
    if(is.matrix(A[[1]])) {p=dim(A[[1]])[1]}
    if(is.vector(A[[1]])) {p=length(A[[1]])}
    K = length(A)
    # results matrices:
    X = list()
    #for(k in 1:K) {X[[k]] = matrix(NA,p,p)}
    for(k in 1:K) {X[[k]] = A[[1]]*NA}
    if(is.matrix(A[[1]])) {fusions = array(FALSE,dim=c(K,K,p,p))}
    if(is.vector(A[[1]])) {fusions = array(FALSE,dim=c(K,K,p,1))}

    # get starting newc: list of matrices.  newc[[k]][i,j] gives the (2*ordermats[[k]]-K-1) adjustment for lam2 at that element.  Basically, how many rank higher vs how many rank lower?
    newc = list()
    for(k in 1:K)
    {
      others = setdiff(1:K,k)
      others.smaller.k = 1:(k-1)
      newc[[k]] = A[[1]]*0
      for(o in others) {newc[[k]] = newc[[k]] + (A[[o]]-A[[k]]< -1e-4) - (A[[o]]-A[[k]]>1e-4)}
    }

    ######### start the loop here:
    for(iter in 1:(K-1))
    {

      # create order matrices:
      ordermats = list()
      for(k in 1:K)
      {
        others = setdiff(1:K,k)
        others.smaller.k = 1:(k-1)
        ordermats[[k]] = A[[1]]*0
        for(o in others) {ordermats[[k]] = ordermats[[k]] + (A[[k]]-A[[o]]>1e-4)}
        # to deal with ties, also add a unit to ordermat[[k]] if a class with a lower k has a tie at element i,j:
        if(k>1)
        {
          for(o in others.smaller.k) {ordermats[[k]] = ordermats[[k]] + (abs(A[[o]]-A[[k]])<1e-4)}
        }
        ordermats[[k]] = ordermats[[k]] + 1
      }

      # create beta.g matrices, holding the solution to Holger's "unconstrained problem"
      #  (prending we're not constraining the order of the solution to match the order of the A matrices)
      betas.g = list()
      for(k in 1:K)
      {
        betas.g[[k]] = A[[k]] - lam2/L*newc[[k]]
      }

      # identify and fuse all elements for which the betas.g are out of order:
      new.ordermats = list()
      for(k in 1:K)
      {
        others = setdiff(1:K,k)
        others.smaller.k = 1:(k-1)
        new.ordermats[[k]] = A[[1]]*0
        for(o in others) {new.ordermats[[k]] = new.ordermats[[k]] + (betas.g[[k]]-betas.g[[o]]>1e-4)}
        # to deal with ties, also add a unit to ordermat[[k]] if a class with a lower k has a tie at element i,j:
        if(k>1)
        {
          for(o in others.smaller.k) {new.ordermats[[k]] = new.ordermats[[k]] + (abs(betas.g[[o]]-betas.g[[k]])<1e-4)}
        }
        new.ordermats[[k]] = new.ordermats[[k]] + 1
      }

      # identify neighboring fusions:  "fusions": K x K x p x p array: K x K matrices, T/F for fusions
      for(k in 1:K){
        for(kp in 1:K){
          #given k,kp, declare a fusion when their ordermats entries are adjacent, and their new.ordermats entries have the opposite direction:
          fusions[k,kp,,] = fusions[k,kp,,]+
            ((ordermats[[k]]-1==ordermats[[kp]])&(new.ordermats[[k]]<new.ordermats[[kp]]))+
            ((ordermats[[k]]+1==ordermats[[kp]])&(new.ordermats[[k]]>new.ordermats[[kp]]))+
            (abs(A[[k]]-A[[kp]])<1e-4)
          #(existing fusions, neighboring fusions, and ties)
          fusions = (fusions>0)*1
        }}


      # now we've noted fusions between all entries with adjacent ordermats entries and reversed new.ordermats entries
      # next: extend fusions to non-adjecent entries: if a-b and b-c, then connect a-c:
      for(k in 1:K){
        for(kp in 1:K){
          others = setdiff(1:K,c(k,kp))
          for(o in others)
          {
            #identify elements in o which are fused with the same element in both k and kp, then add them to the list of k-kp fusions:
            bothfused = fusions[k,o,,] & fusions[kp,o,,]
            fusions[k,kp,,] = fusions[k,kp,,] | bothfused
          }
        }}

      # now recalculate A with the new fused entries:
      # to recalculate A, for each non-zero entry, identify the classes k with which it must be fused, and get their average:
      for(k in 1:K)
      {
        others = setdiff(1:K,k)
        #fusemean and denom: the mean value of all the trueA to be fused, and the number of values to be fused:
        fusemean = trueA[[k]]
        denom = A[[1]]*0+1
        for(o in others)
        {
          fusemean = fusemean+fusions[k,o,,]*trueA[[o]]  #add the values of the elements which must be fused to fusemean
          denom = denom+fusions[k,o,,]
        }
        # now redefine A[[k]]: unchanged from trueA if there's no fusion, and the mean of the fused elements when there is fusion:
        A[[k]] = fusemean/denom
      }

      #newc: list of matrices.  newc[[k]][i,j] gives the (2*ordermats[[k]]-K-1) adjustment for lam2 at that element.  Basically, how many rank higher vs how many rank lower?
      newc = list()
      for(k in 1:K)
      {
        others = setdiff(1:K,k)
        others.smaller.k = 1:(k-1)
        newc[[k]] = A[[1]]*0
        for(o in others) {newc[[k]] = newc[[k]] + (A[[o]]-A[[k]]< -1e-4) - (A[[o]]-A[[k]]>1e-4)}
      }

    } #end loop here

    # final version of betas.g:
    for(k in 1:K)
    {
      betas.g[[k]] = A[[k]] - lam2/L*newc[[k]]
    }
    # now soft-threshold the solution matices:
    for(k in 1:K)
    {
      X[[k]] = soft(betas.g[[k]],lam=lam1/L,penalize.diagonal=penalize.diagonal)
    }
    return(X)
  }


penalty.as.matrix <-
  function(lambda,p,penalize.diagonal)
  {
    # for matrix penalties:  check dim and symmetry:
    if(is.matrix(lambda))
    {
      if(sum(lambda!= t(lambda))>0) {stop("error: penalty matrix is not symmetric")}
      if(sum(abs(dim(lambda)-p))!=0 ) {stop("error: penalty matrix has wrong dimension")}
    }
    # for scalar penalties: convert to matrix form:
    if(length(lambda)==1) {lambda=matrix(lambda,p,p)}
    # apply the penalize.diagonal argument:
    if(!penalize.diagonal) {diag(lambda)=0}
    return(lambda)
  }


soft <-
  function(a,lam,penalize.diagonal){ # if last argument is FALSE, soft-threshold a matrix but don't penalize the diagonal
    out <- sign(a)*pmax(0, abs(a)-lam)
    if(!penalize.diagonal) diag(out) <- diag(a)
    return(out)
  }



dsgl <-
  function(A,L,lam1,lam2,penalize.diagonal)
  {
    lam1 = lam1*1/L
    lam2 = lam2*1/L

    if(is.matrix(A[[1]])) {p=dim(A[[1]])[1]}
    if(is.vector(A[[1]])) {p=length(A[[1]])}
    K=length(A)
    softA = A
    for(k in 1:K) {softA[[k]] = soft(A[[k]],lam1,penalize.diagonal=penalize.diagonal) }   #if penalize.diagonal=FALSE was used in ggl(), then this will not penalize the diagonal.
    normsoftA = A[[1]]*0
    for(k in 1:K) {normsoftA = normsoftA + (softA[[k]])^2}

    normsoftA = sqrt(normsoftA)

    notshrunk = (normsoftA>lam2)*1
    # reset 0 elements of normsoftA to 1 so we don't get NAs later.
    normsoftA = normsoftA + (1-notshrunk)

    out = A
    for(k in 1:K)
    {
      out[[k]] = softA[[k]]*(1-lam2/normsoftA)
      out[[k]] = out[[k]]*notshrunk
    }
    return(out)
  }

