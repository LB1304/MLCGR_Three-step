# ------------------------------------------------------------------------------
# Second step of the three-step procedure.
#
# SecondStep() takes the model fitted in the first step and returns the
# quantities that describe the classification of the observations, together with
# everything the third step needs under any of its three estimators.
#
# No package is required.
# ------------------------------------------------------------------------------


## Arguments
##
##   Model   the fitted model returned by FirstStep() in its Model component,
##           or any object providing Pp (posterior probabilities, one row per
##           observation) and piv (class weights)
##
## Value
##
##   Pp        posterior probabilities p(U = u | y_h), one row per observation
##   piv       class weights of the measurement model
##   Lik       normalised likelihood of the responses, proportional to
##             p(y_h | U = u); used by the 3S-IMP estimator
##   MAP       state assigned to each observation by the maximum-a-posteriori
##             rule
##   W         modal assignment in indicator form, one row per observation
##   Entropy   normalised Shannon entropy of each observation, between 0 and 1
##   Summary   summary statistics of Entropy, as reported in the paper
##   D         classification error matrix, D[u, w] = p(W = w | U = u)
##   Wbch      weights of the BCH estimator, the rows of W times D^(-1)

SecondStep <- function(Model) {

  Pp  <- as.matrix(Model$Pp)
  piv <- as.vector(Model$piv)
  k   <- ncol(Pp)

  stopifnot(length(piv) == k, all(piv > 0))


  ## ---------------------------------------------------------------------------
  ## 1. Likelihood of the responses
  ## ---------------------------------------------------------------------------

  ## p(U = u | y_h) is proportional to p(y_h | U = u) * piv_u, so dividing by the
  ## class weights and rescaling returns the responses' likelihood up to a
  ## constant that cancels in the third step
  Lik <- Pp / rep(piv, each = nrow(Pp))
  Lik <- Lik / rowSums(Lik)


  ## ---------------------------------------------------------------------------
  ## 2. Modal assignment and its uncertainty
  ## ---------------------------------------------------------------------------

  MAP <- max.col(Pp)
  W   <- matrix(0, nrow(Pp), k)
  W[cbind(seq_len(nrow(Pp)), MAP)] <- 1

  ## normalised Shannon entropy, with the convention 0 * log(0) = 0
  Ent     <- Pp * log(Pp); Ent[Pp == 0] <- 0
  Entropy <- -rowSums(Ent) / log(k)

  Summary <- c(Minimum = min(Entropy), Q1 = unname(quantile(Entropy, 0.25)),
               Median = median(Entropy), Mean = mean(Entropy),
               Q3 = unname(quantile(Entropy, 0.75)), Maximum = max(Entropy))


  ## ---------------------------------------------------------------------------
  ## 3. Classification error
  ## ---------------------------------------------------------------------------

  ## D[u, w] is the share of the expected number of observations in state u that
  ## the modal rule assigns to state w
  D <- (t(Pp) %*% W) / colSums(Pp)
  dimnames(D) <- list(true = 1:k, assigned = 1:k)

  if (abs(det(D)) < 1e-10)
    stop("The classification error matrix is nearly singular: the BCH estimator ",
         "cannot be used. Its determinant is ", format(det(D)), ".")

  Wbch <- W %*% solve(D)

  list(Pp = Pp, piv = piv, Lik = Lik, MAP = MAP, W = W,
       Entropy = Entropy, Summary = Summary, D = D, Wbch = Wbch)
}
