## ---- echo = FALSE-------------------------------------------------------
set.seed(42)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options(digits = 2)

## ----echo=TRUE, results="hide", warning=FALSE----------------------------
library(estimatr)

# Example dataset to be used throughout built using fabricatr and randomizr
library(fabricatr)
library(randomizr)
dat <- fabricate(
  N = 100,                        # sample size
  x = runif(N, 0, 1),             # pre-treatment covariate
  y0 = rnorm(N, mean = x),        # control potential outcome
  y1 = y0 + 0.35,                 # treatment potential outcome
  z = complete_ra(N),             # complete random assignment to treatment
  y = ifelse(z, y1, y0),          # observed outcome

  # We will also consider clustered data
  clust = sample(rep(letters[1:20], each = 5)),
  z_clust = cluster_ra(clust),
  y_clust = ifelse(z_clust, y1, y0)
)

head(dat)

## ----echo=FALSE----------------------------------------------------------
knitr::kable(head(dat))

## ---- lm_robust----------------------------------------------------------
res <- lm_robust(y ~ z + x, data = dat)
summary(res)

## ---- results="hide"-----------------------------------------------------
tidy(res)

## ---- echo=FALSE---------------------------------------------------------
knitr::kable(tidy(res))

## ---- echo=TRUE, results="hide"------------------------------------------
res_cl <- lm_robust(
  y_clust ~ z_clust + x,
  data = dat,
  clusters = clust
)
tidy(res_cl)

## ----echo=FALSE----------------------------------------------------------
knitr::kable(tidy(res_cl))

## ----echo=TRUE, results="hide"-------------------------------------------
res_stata <- lm_robust(
  y_clust ~ z_clust + x,
  data = dat,
  clusters = clust,
  se_type = "stata"
)
tidy(res_stata)

## ----echo=FALSE----------------------------------------------------------
knitr::kable(tidy(res_stata))

## ----margins, echo=TRUE--------------------------------------------------
library(margins)

res_int <- lm_robust(y ~ x * z, data = dat)
mar_int <- margins(res_int, vce = "delta")
summary(mar_int)

## ----texreg, echo=TRUE, eval=FALSE---------------------------------------
#  library(texreg)
#  
#  tex_int <- extract(res_int)
#  texreg(tex_int, file = "ex.tex")

## ----echo=TRUE, results="hide"-------------------------------------------
res_lin <- lm_lin(
  y ~ z,
  covariates = ~ x,
  data = dat
)
tidy(res_lin)

## ----echo=FALSE----------------------------------------------------------
knitr::kable(tidy(res_lin))

## ----echo=TRUE, results="hide"-------------------------------------------
# Simple version
res_dim <- difference_in_means(
  y ~ z,
  data = dat
)
tidy(res_dim)

## ----echo=FALSE----------------------------------------------------------
knitr::kable(tidy(res_dim))

## ----echo=TRUE, results="hide"-------------------------------------------
# Clustered version
res_dim_cl <- difference_in_means(
  y_clust ~ z_clust,
  data = dat,
  clusters = clust
)

## ----echo=FALSE----------------------------------------------------------
knitr::kable(tidy(res_dim_cl))

## ----mp------------------------------------------------------------------
data(sleep)
res_mps <- difference_in_means(extra ~ group, data = sleep, blocks = ID)
res_mps$design

## ---- results="hide"-----------------------------------------------------
# Complete random assignment declaration
crs_decl <- declare_ra(
  N = nrow(dat),
  prob = 0.5,
  simple = FALSE
)

ht_comp <- horvitz_thompson(
  y ~ z,
  data = dat,
  declaration = crs_decl
)
tidy(ht_comp)

## ---- echo=FALSE---------------------------------------------------------
knitr::kable(tidy(ht_comp))

## ---- results = "hide"---------------------------------------------------
# Clustered random assignment declaration
crs_clust_decl <- declare_ra(
  N = nrow(dat),
  clusters = dat$clust,
  prob = 0.5,
  simple = FALSE
)

ht_clust <- horvitz_thompson(
  y_clust ~ z_clust,
  data = dat,
  declaration = crs_clust_decl
)
tidy(ht_clust)

## ----echo=FALSE----------------------------------------------------------
knitr::kable(tidy(ht_clust))

## ---- results="hide"-----------------------------------------------------
# arbitrary permutation matrix
possible_treats <- cbind(
  c(1, 1, 0, 1, 0, 0, 0, 1, 1, 0),
  c(0, 1, 1, 0, 1, 1, 0, 1, 0, 1),
  c(1, 0, 1, 1, 1, 1, 1, 0, 0, 0)
)
arb_pr_mat <- permutations_to_condition_pr_mat(possible_treats)

# Simulating a column to be realized treatment
dat <- data.frame(
  z = possible_treats[, sample(ncol(possible_treats), size = 1)],
  y = rnorm(nrow(possible_treats))
)

ht_arb <- horvitz_thompson(
  y ~ z,
  data = dat,
  condition_pr_mat = arb_pr_mat
)
tidy(ht_arb)

## ----echo=FALSE----------------------------------------------------------
knitr::kable(tidy(ht_arb))

