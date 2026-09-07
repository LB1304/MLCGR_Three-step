# Three-step estimation of a multidimensional hidden Markov model for ordinal data

R code accompanying the paper:

> Brusa L., Pennoni F., Bartolucci F., Cheli M., Maggi L.
> *Multidimensional hidden Markov model for ordinal response variables: a three-step
> approach to map disease progression in myasthenia gravis.*

The code implements the three-step estimation procedure described in the paper: a
multidimensional latent class graded response model for the measurement component,
the classification of the observations into latent states, and the weighted maximum
likelihood estimation of the initial and transition probabilities of a hidden Markov
model with covariates.

---

## Repository structure

```
functions/
  simulate_data.R    SimulateData()  data generation
  first_step.R       FirstStep()     selection of the latent traits and measurement model
  second_step.R      SecondStep()    classification and classification error
  third_step.R       ThirdStep()     initial and transition probabilities

run_analysis.R       worked example, from the data to the three estimators
```

Each file in `functions/` is self-contained: it holds the exported function and the
routines it relies on.

---

## The functions

### `SimulateData()`

Generates multivariate longitudinal ordinal data from the model. The items are
arranged in blocks, and each block is given a vector of support points, one per
latent class: **blocks sharing the same vector measure the same latent trait**, so
the number of latent traits and the partition of the blocks are determined by the
`Support` argument and returned with the data. Participation is intermittent,
subjects may leave the study permanently, and the blocks are not administered
simultaneously.

The number of latent classes is taken from `ncol(Support)`, and the number of
covariates entering each of the two structural models from `n.cov.ini` and
`n.cov.trans`.

### `FirstStep()`

With `search = FALSE`, estimates the model for the trait structure passed in
`multi`. With `search = TRUE`, selects the structure first: starting from one latent
trait per block, at each iteration it estimates every model obtained by merging a
pair of traits and retains the merge that improves the chosen information criterion,
stopping when none does. AIC, BIC and ICL are computed and reported for every model
evaluated, whichever criterion drives the selection.

Each model is fitted from several starting values, which can be distributed over
workers; the items are then reordered so that those measuring the same trait are
adjacent, and the selected model is re-estimated with standard errors.

### `SecondStep()`

Takes the fitted model and returns the posterior probabilities, the modal
assignment, the normalised Shannon entropy of each observation with its summary
statistics, the classification error matrix, and the weights of the BCH estimator.
It is called once, and its output serves all three estimators of the third step.

### `ThirdStep()`

Estimates the two multinomial logit models by weighted maximum likelihood, under one
of three estimators:

| `method` | weights |
|---|---|
| `"3S"` | the posterior probabilities of the second step |
| `"3S-IMP"` | the second and third steps iterated with the first held fixed (Bartolucci, Montanari and Pandolfi, 2015) |
| `"BCH"` | corrected through the classification error matrix (Bolck, Croon and Hagenaars, 2004; Di Mari, Oberski and Vermunt, 2016) |

The covariates entering each of the two models are chosen by name through `cov.ini`
and `cov.trans`, and may differ. The data are used as given: any recoding, merging of
categories or treatment of the missing values is left to the caller.

The initial probabilities use the subjects observed at the first occasion; the
transition probabilities use the pairs of consecutive occasions only. Under
`3S-IMP` the iterations stop when the largest absolute change of the regression
coefficients falls below `tol`, or after `maxit` iterations. Convergence is assessed
on the coefficients rather than on the log-likelihood because the weights change from
one iteration to the next, so the two log-likelihoods are not comparable across
iterations.

---

## Requirements

R (≥ 4.0) and the packages `MultiLCIRT`, `MASS`, `future.apply`, `progress` and
`RcppAlgos`.

```r
install.packages(c("MultiLCIRT", "MASS", "future.apply", "progress", "RcppAlgos"))
```

---

## Running the example

From the root of the repository:

```r
source("run_analysis.R")
```

The example simulates four blocks of five items with four ordered response
categories, three latent classes and one covariate per structural model. Blocks 1
and 2 share the same support points, so the data are generated with three latent
traits and the true partition is `{B1, B2}{B3}{B4}`. Because the truth is known, the
structure recovered by the search can be compared with it.

The script then runs the three steps, calls the third step under each of the three
estimators, and reports the estimated coefficients side by side together with the
fitted initial and transition probabilities. Results are written to `Results.RData`.

---

## Data

The application reported in the paper uses data from the MyRealWorld MG study, which
are not publicly available: they were provided by Vitaccess Ltd under licence, and
the conditions under which they may be obtained are stated in the *Availability of
data and materials* section of the paper. This repository therefore ships a data
generator rather than the data, so that the whole procedure can be run and inspected
without them.

---

## License

`ThirdStep()` includes routines adapted from the **LMest** package, released under
the GPL. This repository is therefore distributed under the GPL-3 licence.