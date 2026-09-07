# ------------------------------------------------------------------------------
# Simulation of multivariate longitudinal ordinal data from a multidimensional
# latent class graded response model whose latent process is a first-order
# Markov chain with covariates.
#
# The items are arranged in blocks. Each block is assigned a vector of support
# points, one per latent class; blocks sharing the same vector measure the same
# latent trait, so the number of latent traits and the partition of the blocks
# are determined by the Support argument and are returned with the data. They
# can then be compared with the structure recovered by GradedModelSearch().
#
# Participation is intermittent, subjects may leave the study permanently, and
# the blocks are not administered simultaneously.
# ------------------------------------------------------------------------------


## Arguments
##
##   Support     matrix with one row per block of items and one column per
##               latent class, holding the support points; rows that coincide
##               identify blocks measuring the same latent trait
##   block.size  number of items in each block; one entry per row of Support
##   n, TT       sample size and number of time occasions
##   n.cat       number of response categories, common to all items
##   eta         discrimination parameter, either common to all items or one
##               per item
##   beta        difficulty parameters, of length n.cat - 1, common to all items
##   n.cov.ini   number of time-fixed binary covariates entering the initial
##               probabilities; each is drawn from a Bernoulli(0.50)
##   n.cov.trans number of time-varying binary covariates entering the
##               transition probabilities; each is drawn from a Bernoulli(0.35)
##   Be          coefficients of the initial probabilities: a matrix with
##               1 + n.cov.ini rows, the first being the intercept, and k-1
##               columns, for classes 2, ..., k against class 1. NULL uses a
##               default.
##   Ga          coefficients of the transition probabilities: an array of
##               dimension (1 + n.cov.trans) x (k-1) x k, where Ga[, , u]
##               refers to the k-1 states other than u, in the order given by
##               setdiff(1:k, u), against persistence in u. NULL uses a
##               default.
##   p.resp      response probability at the first occasion
##   decay       factor by which the response probability declines per occasion
##   p.drop      probability of permanently leaving the study, per occasion
##   p.block     probability that a block is administered, given a response
##   seed        passed to set.seed(); NULL leaves the generator untouched
##
## Value
##
##   Y.Red           id, time and the item responses
##   X.Red           id, time and the covariates, named X.fixed.1, ... and
##                   X.timevar.1, ...
##   Full.Multi      block structure, to be passed to GradedModelSearch()
##   True.Multi      trait structure used to generate the data
##   True.Partition  blocks grouped by latent trait
##   D               number of latent traits
##   Be, Ga          coefficients actually used
##   State           latent state underlying each row of Y.Red

SimulateData <- function(Support,
                         block.size,
                         n       = 500,
                         TT      = 10,
                         n.cat   = 4,
                         eta     = 1,
                         beta    = c(-1, 0, 1),
                         n.cov.ini   = 1,
                         n.cov.trans = 1,
                         Be      = NULL,
                         Ga      = NULL,
                         p.resp  = 0.75,
                         decay   = 0.90,
                         p.drop  = 0.12,
                         p.block = 0.70,
                         seed    = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  Support <- as.matrix(Support)
  B <- nrow(Support)
  k <- ncol(Support)
  r <- sum(block.size)
  
  n.ini   <- n.cov.ini            # time-fixed covariates
  n.trans <- n.cov.trans          # time-varying covariates
  
  stopifnot(length(block.size) == B, B >= 1, k >= 2, TT >= 2, n.cat >= 2,
            length(beta) == n.cat - 1, length(eta) %in% c(1, r),
            n.ini >= 1, n.trans >= 1)
  
  Eta   <- if (length(eta) == 1) rep(eta, r) else eta
  Block <- rep(1:B, block.size)
  
  
  ## -------------------------------------------------------------------------
  ## 1. Latent traits implied by Support
  ## -------------------------------------------------------------------------
  
  Key   <- apply(Support, 1, paste, collapse = "_")
  Trait <- as.integer(factor(Key, levels = unique(Key)))
  D     <- max(Trait)
  Xi    <- Support[!duplicated(Key), , drop = FALSE]        # D x k
  Item.Trait <- Trait[Block]
  
  Pad <- function(rows) {
    W <- max(sapply(rows, length))
    matrix(unlist(lapply(rows, function(z) c(z, rep(0, W - length(z))))),
           nrow = length(rows), byrow = TRUE)
  }
  Full.Multi     <- Pad(split(1:r, Block))                  # one row per block
  True.Multi     <- Pad(split(1:r, Item.Trait))             # one row per trait
  True.Partition <- split(1:B, Trait)
  
  
  ## -------------------------------------------------------------------------
  ## 2. Default coefficients of the latent process
  ## -------------------------------------------------------------------------
  
  ## Persistence decreases with the distance between states
  ## The covariate pushes towards the higher states
  ## every covariate is given the same effect by default
  if (is.null(Be))
    Be <- rbind(rep(0, k-1),
                matrix(seq(0.4, 0.8, length.out = k-1), n.ini, k-1, byrow = TRUE))
  if (is.null(Ga)) {
    Ga <- array(0, c(1 + n.trans, k-1, k))
    for (u in 1:k) {
      Oth <- setdiff(1:k, u)
      Ga[1, , u] <- -1.0 - 1.5 * (abs(Oth - u) - 1)
      Ga[-1, , u] <- matrix(0.6 * sign(Oth - u), n.trans, k-1, byrow = TRUE)
    }
  }
  Be <- as.matrix(Be)
  rownames(Be) <- c("intercept", sprintf("x%d", 1:n.ini))
  stopifnot(dim(Be) == c(1 + n.ini, k-1), dim(Ga) == c(1 + n.trans, k-1, k))
  
  
  ## -------------------------------------------------------------------------
  ## 3. Response probabilities of the graded response model
  ## -------------------------------------------------------------------------
  
  ## p(Y_j >= y | U = xi) = plogis(eta_j * xi - beta_y), y = 1, ..., n.cat - 1
  Prob <- array(NA, c(r, k, n.cat))
  for (j in 1:r) for (u in 1:k) {
    Ge <- plogis(Eta[j] * Xi[Item.Trait[j], u] - beta)
    Prob[j, u, ] <- diff(-c(1, Ge, 0))
  }
  
  
  ## -------------------------------------------------------------------------
  ## 4. Covariates and latent process
  ## -------------------------------------------------------------------------
  
  x.ini   <- matrix(rbinom(n * n.ini, 1, 0.50), n, n.ini)
  x.trans <- array(rbinom(n * TT * n.trans, 1, 0.35), c(n, TT, n.trans))
  
  Softmax <- function(z) { z <- z - max(z); exp(z) / sum(exp(z)) }
  
  State <- matrix(NA_integer_, n, TT)
  for (i in 1:n) {
    Lin <- c(0, Be[1, ] + drop(x.ini[i, ] %*% Be[-1, , drop = FALSE]))
    State[i, 1] <- sample(k, 1, prob = Softmax(Lin))
    for (t in 2:TT) {
      u   <- State[i, t-1]
      Oth <- setdiff(1:k, u)
      Gu  <- matrix(Ga[-1, , u], n.trans, k-1)
      Lin <- numeric(k)
      Lin[Oth] <- Ga[1, , u] + drop(x.trans[i, t, ] %*% Gu)
      State[i, t] <- sample(k, 1, prob = Softmax(Lin))
    }
  }
  
  
  ## -------------------------------------------------------------------------
  ## 5. Participation
  ## -------------------------------------------------------------------------
  
  Keep <- vector("list", n)
  for (i in 1:n) {
    tt <- integer(0)
    for (t in 1:TT) {
      if (runif(1) < p.resp * decay^(t-1)) tt <- c(tt, t)
      if (t < TT && runif(1) < p.drop) break          # permanent drop-out
    }
    Keep[[i]] <- tt
  }
  
  Rows <- data.frame(id   = rep(sprintf("P%04d", 1:n), sapply(Keep, length)),
                     time = unlist(Keep), stringsAsFactors = FALSE)
  
  ## Blocks administered at each retained occasion
  ## Occasions at which no block is administered are discarded
  Adm  <- matrix(runif(nrow(Rows)*B) < p.block, nrow(Rows), B)
  Ok   <- rowSums(Adm) > 0
  Rows <- Rows[Ok, , drop = FALSE]
  Adm  <- Adm[Ok, , drop = FALSE]
  N    <- nrow(Rows)
  
  
  ## -------------------------------------------------------------------------
  ## 6. Item responses
  ## -------------------------------------------------------------------------
  
  Idx <- match(Rows$id, sprintf("P%04d", 1:n))
  Y   <- matrix(NA_integer_, N, r,
                dimnames = list(NULL, sprintf("item%02d", 1:r)))
  for (h in 1:N) {
    u <- State[Idx[h], Rows$time[h]]
    for (j in which(Adm[h, Block]))
      Y[h, j] <- sample(0:(n.cat-1), 1, prob = Prob[j, u, ])
  }
  
  
  ## -------------------------------------------------------------------------
  ## 7. Output
  ## -------------------------------------------------------------------------
  
  list(Y.Red = cbind(Rows, as.data.frame(Y)),
       X.Red = setNames(data.frame(Rows,
                                   x.ini[Idx, , drop = FALSE],
                                   vapply(1:n.trans,
                                          function(m) x.trans[cbind(Idx, Rows$time, m)],
                                          numeric(N)),
                                   stringsAsFactors = FALSE),
                        c("id", "time", sprintf("X.fixed.%d", 1:n.ini),
                          sprintf("X.timevar.%d", 1:n.trans))),
       Full.Multi     = Full.Multi,
       True.Multi     = True.Multi,
       True.Partition = True.Partition,
       D              = D,
       Be             = Be,
       Ga             = Ga,
       State          = State[cbind(Idx, Rows$time)])
}