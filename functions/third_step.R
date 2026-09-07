# ------------------------------------------------------------------------------
# Third step of the three-step procedure.
#
# ThirdStep() estimates the initial and transition probabilities of the hidden
# Markov model as functions of covariates, by weighted maximum likelihood, under
# one of three estimators:
#
#   3S      the weights are the posterior probabilities of the second step
#   3S-IMP  the second and third steps are iterated with the first step held
#           fixed, so that the weights are progressively updated
#           (Bartolucci, Montanari and Pandolfi, 2015)
#   BCH     the weights are corrected through the classification error matrix
#           (Bolck, Croon and Hagenaars, 2004; Di Mari, Oberski and Vermunt, 2016)
#
# The covariates entering each of the two models are chosen by name. Any
# recoding of the covariates is left to the caller: the data are used as given.
#
# est_multilogit() and prob_multilogit(), at the end of this file, are adapted
# from the LMest package.
#
# Requires: MASS
# ------------------------------------------------------------------------------


## Arguments
##
##   data       data frame with one row per observation, in the same order as
##              the rows of the responses analysed in the first step
##   Second     output of SecondStep()
##   cov.ini    names of the columns entering the model for the initial
##              probabilities
##   cov.trans  names of the columns entering the model for the transition
##              probabilities
##   method     estimator to be used
##   intercept  TRUE to prepend a column of ones to both design matrices
##   id, time   names of the identifier and of the occasion columns
##   tol        convergence tolerance of 3S-IMP, on the largest absolute change
##              of the regression coefficients between two iterations
##   maxit      maximum number of iterations of 3S-IMP
##   fort       passed to est_multilogit()
##   verbose    TRUE to report the estimator, the information it uses and, for
##              3S-IMP, the iterations
##
## Value
##
##   ini, trans   output of est_multilogit() for the two models
##   Be, Ga       estimates, standard errors and p-values, in long format
##   P.ini        fitted initial probabilities, one row per subject observed at
##                the first occasion
##   P.trans      fitted transition probabilities, k rows per pair of
##                consecutive occasions, in the order of Pairs
##   Pairs        row indices of the pairs of consecutive occasions
##   iter, conv   iterations performed and convergence flag of 3S-IMP

ThirdStep <- function(data, Second, cov.ini, cov.trans,
                      method    = c("3S", "3S-IMP", "BCH"),
                      intercept = TRUE,
                      id      = "id",
                      time    = "time",
                      tol     = 1e-6,
                      maxit   = 25,
                      fort    = FALSE,
                      verbose = TRUE) {
  
  method <- match.arg(method)
  
  Pp <- Second$Pp
  k  <- ncol(Pp)
  
  stopifnot(nrow(data) == nrow(Pp),
            all(c(id, time, cov.ini, cov.trans) %in% names(data)))
  
  
  ## ---------------------------------------------------------------------------
  ## 1. Observations entering each of the two models
  ## ---------------------------------------------------------------------------
  
  Id <- data[[id]]; Tm <- data[[time]]
  
  ## the initial probabilities use the subjects observed at the first occasion
  ind.ini <- which(Tm == min(Tm))
  
  ## the transition probabilities use the pairs of consecutive occasions only
  Pairs <- do.call(rbind, lapply(split(seq_len(nrow(data)), Id), function(ii) {
    ii <- ii[order(Tm[ii])]
    j  <- which(diff(Tm[ii]) == 1)
    if (length(j) == 0) return(NULL)
    cbind(prev = ii[j], curr = ii[j+1])
  }))
  if (is.null(Pairs)) stop("No pair of consecutive occasions is available.")
  nP <- nrow(Pairs)
  
  
  ## ---------------------------------------------------------------------------
  ## 2. Design arrays
  ## ---------------------------------------------------------------------------
  
  X.ini.raw <- as.matrix(data[, cov.ini,   drop = FALSE])
  X.trs.raw <- as.matrix(data[, cov.trans, drop = FALSE])
  if (intercept) {
    X.ini.raw <- cbind("(Intercept)" = 1, X.ini.raw)
    X.trs.raw <- cbind("(Intercept)" = 1, X.trs.raw)
    cov.ini   <- c("(Intercept)", cov.ini)
    cov.trans <- c("(Intercept)", cov.trans)
  }
  J.ini <- ncol(X.ini.raw); J.trs <- ncol(X.trs.raw)
  
  ## class 1 is the reference category of the initial probabilities
  Ref   <- diag(k)[, -1, drop = FALSE]
  X.ini <- array(NA, c(k, J.ini*(k-1), length(ind.ini)))
  for (i in seq_along(ind.ini))
    X.ini[, , i] <- Ref %x% t(X.ini.raw[ind.ini[i], ])
  
  ## persistence in the current state is the reference category of the
  ## transition probabilities; the covariates are taken at the arrival occasion
  X.trans <- array(NA, c(k, k*(k-1)*J.trs, nP*k))
  h <- 0
  for (p in 1:nP) for (u in 1:k) {
    h   <- h + 1
    Sel <- matrix(0, 1, k); Sel[u] <- 1
    X.trans[, , h] <- (Sel %x% diag(k)[, -u, drop = FALSE]) %x%
      t(X.trs.raw[Pairs[p, "curr"], ])
  }
  
  
  if (verbose) {
    Lab <- c("3S"     = "weights given by the posterior probabilities",
             "3S-IMP" = "second and third steps iterated, first step held fixed",
             "BCH"    = "weights corrected through the classification error matrix")
    cat("\nThird step, ", method, ": ", Lab[[method]], "\n", sep = "")
    cat("  ", nrow(data), " observations | ", length(ind.ini),
        " subjects at the first occasion | ", nP,
        " pairs of consecutive occasions\n", sep = "")
    cat("  ", J.ini, " covariates on the initial and ", J.trs,
        " on the transition probabilities",
        if (intercept) " (intercept included)" else "", "\n", sep = "")
  }
  
  
  ## ---------------------------------------------------------------------------
  ## 3. Weights and estimation
  ## ---------------------------------------------------------------------------
  
  ## builds the two response matrices from a matrix of observation weights
  MakeY <- function(Wg) {
    Yt <- matrix(NA, nP*k, k)
    for (p in 1:nP) for (u in 1:k)
      Yt[(p-1)*k + u, ] <- Wg[Pairs[p, "prev"], u] * Wg[Pairs[p, "curr"], ]
    list(ini = Wg[ind.ini, , drop = FALSE], trans = Yt)
  }
  
  Fit <- function(Y) list(
    ini   = est_multilogit(Y = Y$ini,   Xdis = X.ini,   fort = fort),
    trans = est_multilogit(Y = Y$trans, Xdis = X.trans, fort = fort)
  )
  
  Iter <- NA_integer_; Conv <- NA
  
  if (method == "3S") {
    
    Out <- Fit(MakeY(Pp))
    
  } else if (method == "BCH") {
    
    Out <- Fit(MakeY(Second$Wbch))
    
  } else {
    
    Lik <- Second$Lik
    
    ## lambda(u | x) on arbitrary rows, given the coefficients Be
    LamIni <- function(rows, Be) {
      Eta <- cbind(0, X.ini.raw[rows, , drop = FALSE] %*% Be)
      E   <- exp(Eta - apply(Eta, 1, max))
      E / rowSums(E)
    }
    
    ## model-implied marginals, propagated along the consecutive occasions; the
    ## propagation restarts at the beginning of every consecutive run, because
    ## participation is intermittent and the covariates are not observed within
    ## the gaps
    PropMarg <- function(Be, PI) {
      Lam <- matrix(NA, nrow(data), k)
      for (p in 1:nP) {
        Pr <- Pairs[p, "prev"]; Cu <- Pairs[p, "curr"]
        if (is.na(Lam[Pr, 1])) Lam[Pr, ] <- LamIni(Pr, Be)
        Lam[Cu, ] <- as.vector(Lam[Pr, , drop = FALSE] %*% PI[((p-1)*k+1):(p*k), ])
      }
      Miss <- which(is.na(Lam[, 1]))
      if (length(Miss) > 0) Lam[Miss, ] <- LamIni(Miss, Be)
      Lam[ind.ini, ] <- LamIni(ind.ini, Be)
      Lam
    }
    
    ## initialised at the 3S estimates
    Out    <- Fit(MakeY(Pp))
    Be.cur <- matrix(Out$ini$be, J.ini)
    PI.cur <- Out$trans$P
    Par    <- c(Out$ini$be, Out$trans$be)
    Conv   <- FALSE
    
    for (Iter in 1:maxit) {
      
      ## --- the second step, repeated with the current structural parameters ---
      Lam <- PropMarg(Be.cur, PI.cur)
      Wg  <- Lam * Lik; Wg <- Wg / rowSums(Wg)          # q_iu^(t)
      
      ## --- the third step, on the updated weights ---
      Yt <- matrix(NA, nP*k, k)
      for (p in 1:nP) {
        M <- matrix(NA, k, k)
        for (u in 1:k)
          M[u, ] <- Wg[Pairs[p, "prev"], u] * PI.cur[(p-1)*k + u, ] *
            Lik[Pairs[p, "curr"], ]
        Yt[((p-1)*k+1):(p*k), ] <- M / sum(M)           # sum_{u,v} q_iuv^(t) = 1
      }
      
      Out    <- Fit(list(ini = Wg[ind.ini, , drop = FALSE], trans = Yt))
      Be.cur <- matrix(Out$ini$be, J.ini)
      PI.cur <- Out$trans$P
      
      ## the weights change from one iteration to the next, so the two
      ## log-likelihoods are not comparable across iterations: convergence is
      ## assessed on the regression coefficients instead
      Par.new <- c(Out$ini$be, Out$trans$be)
      Delta   <- max(abs(Par.new - Par))
      Par     <- Par.new
      
      if (verbose)
        cat(sprintf("  iteration %2d: max |change| = %.3e\n", Iter, Delta))
      if (Delta < tol) { Conv <- TRUE; break }
    }
    
    if (!Conv)
      warning("3S-IMP did not converge in ", maxit, " iterations; ",
              "the largest change was ", format(Delta), ".")
  }
  
  
  if (verbose) {
    if (method == "3S-IMP")
      cat("  ", if (Conv) "converged" else "NOT converged", " after ", Iter,
          " iterations\n", sep = "")
    cat("  log-likelihood: ", format(Out$ini$lk, digits = 8), " (initial), ",
        format(Out$trans$lk, digits = 8), " (transition)\n", sep = "")
  }
  
  
  ## ---------------------------------------------------------------------------
  ## 4. Estimates, standard errors and p-values
  ## ---------------------------------------------------------------------------
  
  Pv <- function(e, s) 2 * (1 - pnorm(abs(e / s)))
  
  Be   <- matrix(Out$ini$be, J.ini, k-1)
  seBe <- matrix(sqrt(diag(MASS::ginv(Out$ini$Fi))), J.ini, k-1)
  Tab.Be <- data.frame(
    covariate = rep(cov.ini, k-1),
    class     = rep(2:k, each = J.ini),
    estimate  = as.vector(Be), se = as.vector(seBe),
    p.value   = as.vector(Pv(Be, seBe)), row.names = NULL)
  
  Ga   <- array(Out$trans$be, c(J.trs, k-1, k))
  seGa <- array(sqrt(diag(MASS::ginv(Out$trans$Fi))), c(J.trs, k-1, k))
  Tab.Ga <- do.call(rbind, lapply(1:k, function(u) {
    To <- setdiff(1:k, u)
    do.call(rbind, lapply(seq_along(To), function(v) data.frame(
      covariate = cov.trans, from = u, to = To[v],
      estimate  = Ga[, v, u], se = seGa[, v, u],
      p.value   = Pv(Ga[, v, u], seGa[, v, u]), row.names = NULL)))
  }))
  
  list(ini = Out$ini, trans = Out$trans, Be = Tab.Be, Ga = Tab.Ga,
       P.ini = Out$ini$P, P.trans = Out$trans$P,
       Pairs = Pairs, ind.ini = ind.ini,
       method = method, iter = Iter, conv = Conv)
}


## ---------------------------------------------------------------------------
## Multinomial logit model, adapted from the LMest package
## ---------------------------------------------------------------------------

# Multinomial logit model fitted by Fisher scoring.
#
# The rows of Y need not be indicator vectors: in the third step of the procedure they hold
# the posterior probabilities of the latent states, so the fit is a weighted maximum
# likelihood one.
#
#   Y      n x k matrix of (possibly fractional) responses
#   Xdis   design array, k x ncov x ndis; a two-dimensional Xdis is expanded internally
#   label  maps each row of Y to a slice of Xdis, for repeated design matrices
#   be     optional starting values for the regression coefficients
#   ex     TRUE to stop after computing score and information, without updating be
#   fort   TRUE to use the compiled Fortran routines of LMest, FALSE (the default) for
#          the pure R code, which makes this file self-contained
#
# Returns the estimated coefficients, the fitted probabilities, the score and the
# information matrix.
est_multilogit <- function(Y,Xdis,label=1:n,be=NULL,Pdis=NULL,dis=FALSE,
                           fort=FALSE,ex=FALSE,tol=10^-8){
  
  # preliminaries
  n = dim(Y)[1]
  k = dim(Y)[2]
  nbe = dim(Xdis)[2]
  ndis = max(label)
  
  # correct covariance matrix
  if(length(dim(Xdis))==2) Xdis = aperm(array(Xdis,c(ndis,nbe,1)),c(3,2,1))
  if(dim(Xdis)[1]<k){
    Xdis0 = Xdis
    Xdis = array(0,c(k,nbe,ndis))
    Xdis[2:k,,] = Xdis0
  }
  
  # starting values
  if(is.null(be)){
    be = rep(0,nbe)
    Pdis = NULL}
  if(is.null(Pdis)){
    Pdis = prob_multilogit(Xdis,be,label,fort)$Pdis
  }
  Ydis = matrix(0,ndis,k)
  if(fort==F){
    for(i in 1:ndis){
      li = (label==i)
      if(sum(li)==1) Ydis[i,] = as.numeric(Y[li, ]) else Ydis[i,] = colSums(Y[li,])
    }
  } else {
    out = .Fortran("sum_Y",Ydis=Ydis,Y=Y,label=as.integer(label),ndis=as.integer(ndis),ns=as.integer(n),k=as.integer(k))
    Ydis = out$Ydis
  }
  ny = rowSums(Ydis)
  lk = sum(Ydis*log(Pdis))
  if(dis) print(c(0,lk))
  # iterate until convergence
  it = 0; lko = lk
  while(((lk-lko)/max(abs(lko),1)>tol & it<100) | it==0){
    it = it+1; lko = lk
    if(fort==F){
      sc = 0; Fi = 0
      for(i in 1:ndis){
        pdis = Pdis[i,]
        sc = sc+t(Xdis[,,i])%*%(Ydis[i,]-ny[i]*pdis)
        Fi = Fi+ny[i]*t(Xdis[,,i])%*%(diag(pdis)-pdis%o%pdis)%*%Xdis[,,i]
      }
    } else {
      out = .Fortran("nr_multilogit",Xdis=Xdis,be=be,Pdis=Pdis,Ydis=Ydis,ny=ny,k=as.integer(k),ndis=as.integer(ndis),ncov=as.integer(nbe),
                     sc=rep(0,nbe),Fi=matrix(0,nbe,nbe))
      sc = out$sc; Fi = out$Fi
    }
    
    if(ex==FALSE){
      dbe = as.vector(ginv(Fi)%*%sc)
      mdbe = max(abs(dbe))
      if(mdbe>0.5) dbe = dbe/mdbe*0.5
      be0 = be
      flag = TRUE
      while(flag){
        be = be0+dbe
        Eta = matrix(0,ndis,k)
        for(i in 1:ndis){
          if(nbe==1) Eta[i,] = Xdis[,,i]*c(be)
          else Eta[i,] = Xdis[,,i]%*%be
        }
        if(max(abs(Eta))>100){
          dbe = dbe/2
          flag = TRUE
        } else {
          flag = FALSE
        }
      }
    }
    
    # compute again probabilities
    out = prob_multilogit(Xdis,be,label,fort)
    P = out$P; Pdis = out$Pdis
    lk = sum(Ydis*log(Pdis))
    if(dis) print(c(it,lk,lk-lko))
  }
  out = list(be=be, P=P, Pdis=Pdis, sc=sc, Fi=Fi, lk=lk)
}


# Fitted probabilities implied by the coefficients `be`, returned both per distinct design
# matrix (Pdis) and per observation (P). With der = TRUE the derivatives with respect to be
# are also returned.
prob_multilogit <- function(Xdis, be, label, fort = FALSE, der = FALSE) {
  
  k = dim(Xdis)[1]
  ndis = max(label)
  n = length(label)
  ncov = length(be)
  Pdis = matrix(0,ndis,k); P = matrix(0,n,k)
  if(fort == FALSE){
    for(i in 1:ndis){
      if(ncov==1) pdis = exp(Xdis[,,i]*be) 
      else pdis = exp(Xdis[,,i]%*%be) 
      pdis = pdis/sum(pdis)
      Pdis[i,] = pdis
      mul = sum(label==i)
      P[label==i,] = rep(1,mul)%o%pdis
    }
  } else {
    out = .Fortran("prob_multilogif",Xdis=Xdis,be=be,label=as.integer(label),Pdis=Pdis,P=P,k=as.integer(k),ndis=as.integer(ndis),
                   ns=as.integer(n),ncov=as.integer(ncov))
    Pdis = out$Pdis; P = out$P
  }
  if(der){
    lbe = length(be)	
    dPdis = array(0,c(ndis,k,lbe)); dP = array(0,c(n,k,lbe))
    for(i in 1:ndis){
      pdis = Pdis[i,]
      Op = diag(pdis)-pdis%o%pdis
      dPdis[i,,] = Op%*%Xdis[,,i]
      ind = which(label==i)
      for(j in ind) dP[j,,] = dPdis[i,,]
    }
  } else {
    dPdis = NULL; dP = NULL
  }
  out = list(Pdis=Pdis,P=P,dPdis=dPdis,dP=dP)
  return(out)
}