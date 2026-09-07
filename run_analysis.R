# ------------------------------------------------------------------------------
# Worked example of the three-step procedure, on a simulated data set.
#
# Generates data from the model, selects the number of latent traits and
# estimates the measurement model (first step), classifies the observations
# (second step), and estimates the initial and transition probabilities under
# the three estimators 3S, 3S-IMP and BCH (third step).
#
# Because the data are simulated, the number of latent traits and the partition
# of the blocks of items are known, and the structure recovered by the search
# can be compared with them.
#
# Run from the root of the repository.
# ------------------------------------------------------------------------------

library(MultiLCIRT)
library(future.apply)
library(progress)
library(RcppAlgos)
library(MASS)

source("functions/simulate_data.R")
source("functions/first_step.R")
source("functions/second_step.R")
source("functions/third_step.R")

k       <- 3      # latent classes
N.Init  <- 20     # starting values per model
Workers <- 20     # parallel workers


## -----------------------------------------------------------------------------
## 1. Simulated data
## -----------------------------------------------------------------------------

## Four blocks of five items, with four ordered response categories. 
## Blocks 1 and 2 share the same vector of support points, 
## so the data are generated with three latent traits and 
## the true partition of the blocks is {B1, B2}{B3}{B4}.

Support <- rbind(c(-1.5,  0.2, 1.5),
                 c(-1.5,  0.2, 1.5),
                 c(-1.5,  1.0, 0.5),
                 c(-1.5, -1.0, 2.5))

Sim <- SimulateData(Support     = Support,
                    block.size  = rep(5, 4),
                    n           = 500,
                    TT          = 10,
                    n.cat       = 4,
                    n.cov.ini   = 1,
                    n.cov.trans = 1,
                    seed        = 2801)

S    <- Sim$Y.Red[, -(1:2)]
Data <- Sim$X.Red

cat("subjects:", length(unique(Sim$Y.Red$id)),
    " rows:", nrow(Sim$Y.Red),
    " items:", ncol(S),
    " true latent traits:", Sim$D, "\n\n")


## -----------------------------------------------------------------------------
## 2. First step: selection of the number of latent traits (AIC)
##                and estimation of the measurement model
## -----------------------------------------------------------------------------

First <- FirstStep(S = S, k = k, multi = Sim$Full.Multi, search = TRUE,
                   criterion = "AIC", N.Init = N.Init,
                   parallel = TRUE, workers = Workers)

cat("\nLatent traits: true =", Sim$D, " selected =", nrow(First$multi), "\n")
cat("True partition of the blocks:",
    paste(sapply(Sim$True.Partition,
                 function(z) paste0("{", paste(z, collapse = ","), "}")),
          collapse = " "), "\n")
cat("Selected structure:\n"); print(First$multi)

print(First$Criteria)


## -----------------------------------------------------------------------------
## 3. Second step: classification and its uncertainty
## -----------------------------------------------------------------------------

Second <- SecondStep(First$Model)

cat("\nNormalised entropy of the classification:\n")
print(round(Second$Summary, 3))
cat("\nClassification error matrix:\n")
print(round(Second$D, 3))


## -----------------------------------------------------------------------------
## 4. Third step: initial and transition probabilities
## -----------------------------------------------------------------------------

Cov.Ini   <- grep("^X\\.fixed",   names(Data), value = TRUE)
Cov.Trans <- grep("^X\\.timevar", names(Data), value = TRUE)

Est <- list(
  "3S"     = ThirdStep(Data, Second, Cov.Ini, Cov.Trans, method = "3S"),
  "3S-IMP" = ThirdStep(Data, Second, Cov.Ini, Cov.Trans, method = "3S-IMP",
                       tol = 1e-6, maxit = 25, verbose = FALSE),
  "BCH"    = ThirdStep(Data, Second, Cov.Ini, Cov.Trans, method = "BCH")
)


## -----------------------------------------------------------------------------
## 5. Estimated regression coefficients
## -----------------------------------------------------------------------------

Cmp.Be <- cbind(Est[["3S"]]$Be[, c("covariate", "class")],
                sapply(Est, function(x) x$Be$estimate))
Cmp.Ga <- cbind(Est[["3S"]]$Ga[, c("covariate", "from", "to")],
                sapply(Est, function(x) x$Ga$estimate))

cat("\nInitial probabilities, estimated coefficients:\n")
print(Cmp.Be, row.names = FALSE, digits = 3)
cat("\nTransition probabilities, estimated coefficients:\n")
print(Cmp.Ga, row.names = FALSE, digits = 3)


## -----------------------------------------------------------------------------
## 6. Estimated initial and transition probabilities,
##    averaged over the observed covariate patterns
## -----------------------------------------------------------------------------

From <- rep(1:k, times = nrow(Est[["3S"]]$Pairs))

for (m in names(Est)) {
  cat("\n", m, " average initial probabilities:\n", sep = "")
  print(round(colMeans(Est[[m]]$P.ini), 3))
  cat(m, " average transition probabilities:\n", sep = "")
  print(round(t(sapply(1:k, function(u)
    colMeans(Est[[m]]$P.trans[From == u, , drop = FALSE]))), 3))
}


## -----------------------------------------------------------------------------
## 7. Saving
## -----------------------------------------------------------------------------

save(Sim, First, Second, Est, Cmp.Be, Cmp.Ga, file = "Results.RData")
