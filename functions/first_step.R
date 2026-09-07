# ------------------------------------------------------------------------------
# First step of the three-step procedure.
#
# FirstStep() either estimates a multidimensional latent class graded response
# model for a given trait structure, or selects the structure first through the
# stepwise aggregation of the blocks of items and then estimates the selected
# model. The remaining functions in this file implement the search and are not
# meant to be called directly.
#
# Requires: MultiLCIRT, future.apply, progress, RcppAlgos
# ------------------------------------------------------------------------------


## Arguments
##
##   S           matrix of item responses, one row per observation, coded from
##               0; NA for the items not administered
##   k           number of latent classes
##   multi       trait structure: one row per latent trait when search = FALSE,
##               one row per block of items when search = TRUE, padded with
##               zeros. SimulateData() returns both, as True.Multi and
##               Full.Multi.
##   search      TRUE to select the number of latent traits before estimating
##   criterion   information criterion driving the selection; AIC, BIC and ICL
##               are computed and reported in any case
##   tol1, tol2  tolerances of the search and of the final estimate
##   N.Init      starting values per model
##   parallel    TRUE to distribute the starting values over workers
##   quiet       TRUE to suppress the output printed by est_multi_poly
##
## Value
##
##   Model       the fitted model, with standard errors
##   multi       trait structure of the fitted model, renumbered so that the
##               items measuring the same trait are adjacent
##   Ord         permutation applied to the columns of S
##   Search      output of the search, or NULL when search = FALSE
##   Criteria    all the models evaluated during the search, or NULL

FirstStep <- function(S, k, multi, search = FALSE, criterion = c("AIC", "BIC", "ICL"),
                      tol1 = 1e-6, tol2 = 1e-8, N.Init = 20,
                      parallel = FALSE, workers = 1, quiet = TRUE) {

  criterion <- match.arg(criterion)


  ## ---------------------------------------------------------------------------
  ## 1. Trait structure
  ## ---------------------------------------------------------------------------

  if (search) {
    Res <- GradedModelSearch(S = S, k = k, Full.Multi = multi,
                             tol1 = tol1, tol2 = tol2, N.Init = N.Init,
                             parallel = parallel, workers = workers,
                             criterion = criterion, quiet = quiet)

    ## the search stops when no merge improves on the current model, so at the
    ## last step the best row is the current model, whose fit is the best of
    ## the previous step
    Last  <- length(Res)
    Which <- which.min(Res[[Last]]$Res[[criterion]])
    Sel   <- if (Which == 1) Res[[Last-1]]$Mod[[1]]$multi else Res[[Last]]$Mod[[Which]]$multi
    Crit  <- CollectCriteria(Res)
  } else {
    Res <- NULL; Sel <- multi; Crit <- NULL
  }


  ## ---------------------------------------------------------------------------
  ## 2. Items reordered so that those measuring the same trait are adjacent
  ## ---------------------------------------------------------------------------

  Ord   <- as.vector(t(Sel)); Ord <- Ord[Ord != 0]
  S.Ord <- S[, Ord]

  Multi <- t(Sel)
  Mask  <- Multi != 0
  Multi[Mask] <- 1:sum(Mask)
  Multi <- t(Multi)


  ## ---------------------------------------------------------------------------
  ## 3. Final estimate
  ## ---------------------------------------------------------------------------

  Model <- ModEst(S = S.Ord, k = k, multi = Multi, tol1 = tol2, tol2 = tol2,
                  N.Init = N.Init, parallel = parallel, workers = workers,
                  output = TRUE, out_se = TRUE, quiet = quiet)

  stopifnot(nrow(Model$Pp) == nrow(S))

  list(Model = Model, multi = Multi, Ord = Ord, Search = Res, Criteria = Crit)
}


# Stepwise selection of the number of latent traits.
#
# Starts from the model in which each block of items defines its own latent trait and, at
# each iteration, estimates every model obtained by merging a pair of traits. The merge
# attaining the best value of `criterion` is retained; the procedure stops when no merge
# improves on the current model. AIC, BIC and ICL are always computed and reported for
# every model, whichever criterion drives the selection.
#
# Returns a list with one element per iteration, each holding the results table (`Res`) and
# the fitted models (`Mod`); the criterion used is stored as an attribute.
GradedModelSearch <- function(S, k, Full.Multi, tol1 = 1e-6, tol2 = 1e-8, N.Init, parallel, workers,
                              criterion = c("AIC", "BIC", "ICL"), quiet = TRUE) {
  criterion <- match.arg(criterion)

  # the pool of workers is opened once for the whole search and restored on
  # exit; ModEst is then called with set.plan = FALSE
  Old.Plan <- future::plan()
  on.exit(future::plan(Old.Plan), add = TRUE)
  if (parallel) plan(multisession, workers = workers) else plan(sequential)
  Step <- 1
  Old.Multi <- Full.Multi
  Old.Model <- paste(paste0("(", 1:nrow(Old.Multi), ")"), collapse = " - ")
  cat(paste("Step", Step, ": Estimation of the model with", nrow(Old.Multi), "dimensions."))
  Old.Est <- ModEst(S = S, k = k, multi = Old.Multi, tol1 = tol1, tol2 = tol2, 
                    N.Init = N.Init, parallel = parallel, workers = workers,
                    quiet = quiet, set.plan = FALSE)
  cat("\n\n")
  df2print <- data.frame(Model = Old.Model, LogLik = Old.Est$lk, NumPar = Old.Est$np, AIC = Old.Est$aic, BIC = Old.Est$bic,
                         Ent = Old.Est$ent, ICL = GetICL(Old.Est))
  print(FormatRes(df2print), row.names = F)
  cat("\n")
  cat(FitNote(attr(Old.Est, "NumOK"), attr(Old.Est, "NumFail"), attr(Old.Est, "NumInvErr")))
  cat("\n")
  cat(paste0("   [Selection criterion: ", criterion, "]\n\n"))
  
  All.Results <- vector(mode = "list", length = nrow(Old.Multi))
  # slot 1 holds the initial model; Mod is a list, as in every other slot, so that
  # All.Results[[s]]$Mod[[1]] has the same meaning at every step
  All.Results[[1]] <- list(Res = NULL,
                           Mod = list(append(Old.Est, list(multi = Old.Multi))))
  GoOn <- TRUE
  while (GoOn) {
    Step <- Step + 1
    cat(paste("Step", Step, ": Estimation of the models with", nrow(Old.Multi)-1, "dimensions."))
    cat("\n"); cat("\n")
    flush.console()
    Results <- data.frame(matrix(NA, ncol = 7, nrow = 1 + choose(nrow(Old.Multi), 2)))
    colnames(Results) <- c("Model", "LogLik", "NumPar", "AIC", "BIC", "Ent", "ICL")
    
    Results$Model[1] <- Old.Model
    Results$LogLik[1] <- Old.Est$lk; Results$NumPar[1] <- Old.Est$np
    Results$AIC[1] <- Old.Est$aic; Results$BIC[1] <- Old.Est$bic
    Results$Ent[1] <- Old.Est$ent; Results$ICL[1] <- GetICL(Old.Est)
    
    Diag <- c(0, 0, 0)
    D <- nrow(Old.Multi)
    Comb <- RcppAlgos::comboGeneral(v = D, m = 2)
    List.Models <- vector(mode = "list", length = nrow(Comb))
    pb <- progress_bar$new(
      format = "   Estimation: [:bar]:percent; Model: :current/:total; Execution time: :elapsedfull; Remaining time: :eta.", 
      total = nrow(Comb), clear = FALSE
    )
    pb$tick(len = 0) 
    Sys.sleep(0.5)
    pb$tick(len = 0)
    for (i in 1:nrow(Comb)) {
      New.Multi <- GetNewMulti(Old.Multi, Ind2Bind = Comb[i, ])
      Aux.Model <- strsplit(Old.Model, split = " - ")[[1]]
      New.Model <- paste0("(", paste(gsub("[()]", "", Aux.Model[Comb[i, ]]), collapse = ","), ")")
      New.Model <- paste(c(New.Model, Aux.Model[-Comb[i, ]]), collapse = " - ")
      New.Est <- ModEst(S = S, k = k, multi = New.Multi, tol1 = tol1, tol2 = tol2, 
                        N.Init = N.Init, parallel = parallel, workers = workers,
                        quiet = quiet, set.plan = FALSE)
      pb$tick()
      Diag <- Diag + c(attr(New.Est, "NumOK"), attr(New.Est, "NumFail"), attr(New.Est, "NumInvErr"))
      
      Results$Model[i+1] <- New.Model
      Results$LogLik[i+1] <- New.Est$lk; Results$NumPar[i+1] <- New.Est$np
      Results$AIC[i+1] <- New.Est$aic; Results$BIC[i+1] <- New.Est$bic
      Results$Ent[i+1] <- New.Est$ent; Results$ICL[i+1] <- GetICL(New.Est)
      List.Models[[i]] <- append(New.Est, list(multi = New.Multi))
      names(List.Models)[i] <- New.Model
    }
    
    Crit <- Results[[criterion]]
    First <- Crit[1]
    if (is.na(First)) stop("The current model has no valid value of ", criterion, ".")
    if (all(is.na(Crit[-1]))) stop("No candidate model produced a valid value of ", criterion, ".")
    Best <- min(Crit, na.rm = TRUE)
    
    if (Best < First) {
      # A model with fewer traits is better: iterate
      Ind.Best <- which.min(Crit)
      Old.Multi <- GetNewMulti(Old.Multi, Ind2Bind = Comb[Ind.Best-1, ])
      Old.Model <- Results$Model[Ind.Best]
      Old.Est <- List.Models[[Ind.Best-1]]
      
      Ord <- order(Crit, decreasing = F)
      Results <- Results[Ord, ]
      # a NULL placeholder keeps Mod[[j]] aligned with row j of Res: the current model,
      # which has no fitted candidate of its own, takes the NULL slot
      List.Models <- c(list(NULL), List.Models)[Ord]
      
      cat("\n\n")
      print(FormatRes(Results), row.names = F)
      cat("\n")
      cat(FitNote(Diag[1], Diag[2], Diag[3]))
      cat("\n\n")
      
      GoOn <- TRUE
    } else if (Best == First) {
      # The current model is the best: stop
      Ord <- order(Crit, decreasing = F)
      Results <- Results[Ord, ]
      # a NULL placeholder keeps Mod[[j]] aligned with row j of Res: the current model,
      # which has no fitted candidate of its own, takes the NULL slot
      List.Models <- c(list(NULL), List.Models)[Ord]
      cat("\n\n")
      print(FormatRes(Results), row.names = F)
      cat("\n")
      cat(FitNote(Diag[1], Diag[2], Diag[3]))
      cat("\n\n")
      GoOn <- FALSE
    }
    All.Results[[length(All.Results) - D + 2]] <- list(Res = Results, Mod = List.Models)
  }
  All.Results[sapply(All.Results, is.null)] <- NULL
  attr(All.Results, "criterion") <- criterion
  return(All.Results)
}


# ICL in the ICL-BIC form of Biernacki, Celeux and Govaert (2000):
#   ICL = -2 * logLik + np * log(n) + 2 * EN = BIC + 2 * ent
# where EN is the classification entropy. est_multi_poly already returns this quantity in
# $ent, so no additional computation or refitting is needed. Lower values indicate a better
# model, consistently with $aic and $bic.
GetICL <- function(Est) {
  if (is.null(Est$ent)) return(NA_real_)
  Est$bic + 2 * Est$ent
}


# Rounds the criteria for display only; the values stored in the results are kept at full
# precision.
FormatRes <- function(Res) {
  for (cl in c("LogLik", "AIC", "BIC", "Ent", "ICL")) {
    if (!is.null(Res[[cl]])) Res[[cl]] <- round(Res[[cl]], 1)
  }
  Res
}


# One-line summary of the fits suppressed by quiet = TRUE. NumInvErr counts the
# "inversion error in P" messages emitted by the M step of the EM algorithm, not by the
# standard error computation, which is not performed when out_se = FALSE. A non-zero value
# signals numerical difficulties and should not be ignored.
FitNote <- function(NumOK, NumFail, NumInvErr) {
  paste0("   [Starting values: ", NumOK, " converged, ", NumFail, " failed]")
}


# Number of difficulty parameters counted in excess by est_multi_poly.
#
# The package assumes l = max(S) + 1 categories for every item and counts (l - 1) free
# thresholds each, whereas items with c_j < l categories only have (c_j - 1) identified;
# the remaining ones diverge and are not free parameters. The difference, sum_j (l - c_j),
# is constant across the number of traits and of classes, so it does not affect model
# comparisons, but it must be subtracted from the absolute values reported.
#
# Also warns about items with unobserved intermediate categories, a different and more
# serious problem.
CountOvercount <- function(S) {
  l  <- max(S, na.rm = TRUE) + 1
  Cj <- apply(S, 2, function(x) max(x, na.rm = TRUE)) + 1
  Gap <- apply(S, 2, function(x) {
    u <- sort(unique(na.omit(x)))
    length(u) < (max(x, na.rm = TRUE) + 1)
  })
  if (any(Gap)) warning("Items with unobserved intermediate categories: ",
                        paste(which(Gap), collapse = ", "))
  sum(l - Cj)
}


# Recomputes np, AIC, BIC and ICL after subtracting the overcount returned by
# CountOvercount. `n` is the number of observations used by the package in the BIC, that is
# sum(yv); with yv = 1 this is the number of rows of S after removing units with no
# responses.
AdjustCriteria <- function(Res, Over, n) {
  Res$NumParAdj <- Res$NumPar - Over
  Res$AICadj <- -2 * Res$LogLik + 2 * Res$NumParAdj
  Res$BICadj <- -2 * Res$LogLik + Res$NumParAdj * log(n)
  Res$ICLadj <- Res$BICadj + 2 * Res$Ent
  Res
}


# Collects, in a single data frame, every model evaluated at each iteration, with the three
# criteria side by side and the differences from the best model of that iteration. Not
# needed by the search itself, but convenient for the sensitivity analysis.
CollectCriteria <- function(All.Results) {
  Out <- do.call(rbind, lapply(seq_along(All.Results), function(s) {
    Res <- All.Results[[s]]$Res
    if (is.null(Res)) return(NULL)
    data.frame(Iter = s - 1, Res, row.names = NULL)
  }))
  if (is.null(Out)) return(NULL)
  for (cr in c("AIC", "BIC", "ICL")) {
    Out[[paste0("d", cr)]] <- ave(Out[[cr]], Out$Iter, FUN = function(x) x - min(x, na.rm = TRUE))
  }
  Out
}


# Builds the trait structure obtained by merging the two rows of `Old.Multi` indexed by
# `Ind2Bind`. The merged row is placed first and the remaining rows follow, padded with
# zeros to a common number of columns.
GetNewMulti <- function(Old.Multi, Ind2Bind) {
  NewRow <- sort(as.vector(Old.Multi[Ind2Bind, ]))
  NewRow <- NewRow[NewRow > 0]
  # drop = FALSE keeps OldRow a matrix also when a single trait is left unmerged,
  # so that ncol() below is always well defined
  OldRow <- Old.Multi[-Ind2Bind, , drop = FALSE]
  
  NumR <- nrow(Old.Multi) - 1
  NumC <- max(length(NewRow), ncol(OldRow))
  New.Multi <- matrix(0, nrow = NumR, ncol = NumC)
  
  New.Multi[1, 1:length(NewRow)] <- NewRow
  # when all traits are merged into one there is no remaining row to copy
  if (NumR > 1) New.Multi[2:NumR, 1:ncol(OldRow)] <- OldRow
  
  return(New.Multi)
}


# Fits the MLCGR model for a given trait structure from N.Init starting values (the first
# deterministic, the others random) and returns the fit attaining the highest
# log-likelihood, refined with a tighter tolerance.
#
# With quiet = TRUE the output that est_multi_poly prints unconditionally at every call is
# captured and discarded; only the counts of converged starts, failed starts and inversion
# messages are kept, and returned as attributes of the result.
ModEst <- function(S, k, multi, tol1, tol2, N.Init, output = FALSE, out_se = FALSE,
                   parallel = FALSE, workers = 1, quiet = TRUE, set.plan = TRUE) {
  Init_vec <- c(0, rep(1, N.Init-1))
  
  # set.plan = FALSE when the caller has already opened the pool of workers, so
  # that it is not torn down and rebuilt at every model of the search
  if (set.plan) {
    Old.Plan <- future::plan()
    on.exit(future::plan(Old.Plan), add = TRUE)
    if (parallel) plan(multisession, workers = workers) else plan(sequential)
  }
  
  Out <- future_lapply(1:length(Init_vec), function(j) {
    Init <- Init_vec[j]
    
    if (quiet) {
      Log <- utils::capture.output(
        Fit <- try(est_multi_poly(S = S, k = k, start = Init, tol = tol1, 
                                  link = 1, disc = 1, difl = 0, multi = multi), silent = TRUE)
      )
    } else {
      Log <- character(0)
      Fit <- try(est_multi_poly(S = S, k = k, start = Init, tol = tol1, 
                                link = 1, disc = 1, difl = 0, multi = multi), silent = TRUE)
    }
    list(Fit = Fit, NumInvErr = sum(grepl("inversion error", Log, fixed = TRUE)))
  }, future.seed = TRUE)
  
  Est <- lapply(Out, "[[", "Fit")
  NumInvErr <- sum(unlist(lapply(Out, "[[", "NumInvErr")))
  IsErr <- unlist(lapply(Est, function(x) any(class(x) == "try-error")))
  NumFail <- sum(IsErr)
  Est <- Est[!IsErr]
  if (length(Est) == 0) {
    Msg <- unique(sapply(Out[IsErr], function(x) attr(x$Fit, "condition")$message))
    stop("No starting value produced a valid fit. ",
         "Errors returned by est_multi_poly:\n  ", paste(Msg, collapse = "\n  "))
  }
  
  Best <- Est[[which.max(unlist(lapply(Est, "[[", "lk")))]]
  if (quiet) {
    Log <- utils::capture.output(
      Best <- est_multi_poly(S = S, k = k, start = 2, tol = tol2, output = output, out_se = out_se, 
                             link = 1, disc = 1, difl = 0, multi = multi, 
                             piv = Best$piv, Phi = Best$Phi, gac = Best$gac, De = Best$De)
    )
    NumInvErr <- NumInvErr + sum(grepl("inversion error", Log, fixed = TRUE))
  } else {
    Best <- est_multi_poly(S = S, k = k, start = 2, tol = tol2, output = output, out_se = out_se, 
                           link = 1, disc = 1, difl = 0, multi = multi, 
                           piv = Best$piv, Phi = Best$Phi, gac = Best$gac, De = Best$De)
  }
  
  attr(Best, "NumOK") <- length(Est)
  attr(Best, "NumFail") <- NumFail
  attr(Best, "NumInvErr") <- NumInvErr
  
  return(Best)
}
