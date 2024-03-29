#' Ordinary Least Squares with Robust Standard Errors
#'
#' @description This formula fits a linear model, provides a variety of
#' options for robust standard errors, and conducts coefficient tests
#'
#' @param formula an object of class formula, as in \code{\link{lm}}
#' @param data A \code{data.frame}
#' @param weights the bare (unquoted) names of the weights variable in the
#' supplied data.
#' @param subset An optional bare (unquoted) expression specifying a subset
#' of observations to be used.
#' @param clusters An optional bare (unquoted) name of the variable that
#' corresponds to the clusters in the data.
#' @param fixed_effects An optional right-sided formula containing the fixed
#' effects that will be projected out of the data, such as \code{~ blockID}. Do not
#' pass multiple-fixed effects with intersecting groups. Speed gains are greatest for
#' variables with large numbers of groups and when using "HC1" or "stata" standard errors.
#' See 'Details'.
#' @param se_type The sort of standard error sought. If \code{clusters} is
#' not specified the options are "HC0", "HC1" (or "stata", the equivalent),
#'  "HC2" (default), "HC3", or
#' "classical". If \code{clusters} is specified the options are "CR0", "CR2" (default), or "stata". Can also specify "none", which may speed up estimation of the coefficients.
#' @param ci logical. Whether to compute and return p-values and confidence
#' intervals, TRUE by default.
#' @param alpha The significance level, 0.05 by default.
#' @param return_vcov logical. Whether to return the variance-covariance
#' matrix for later usage, TRUE by default.
#' @param try_cholesky logical. Whether to try using a Cholesky
#' decomposition to solve least squares instead of a QR decomposition,
#' FALSE by default. Using a Cholesky decomposition may result in speed gains, but should only
#' be used if users are sure their model is full-rank (i.e., there is no
#' perfect multi-collinearity)
#'
#' @details
#'
#' This function performs linear regression and provides a variety of standard
#' errors. It takes a formula and data much in the same was as \code{\link{lm}}
#' does, and all auxiliary variables, such as clusters and weights, can be
#' passed either as quoted names of columns, as bare column names, or
#' as a self-contained vector. Examples of usage can be seen below and in the
#' \href{https://declaredesign.org/r/estimatr/articles/getting-started.html}{Getting Started vignette}.
#'
#' The mathematical notes in
#' \href{https://declaredesign.org/r/estimatr/articles/mathematical-notes.html}{this vignette}
#' specify the exact estimators used by this function.
#' The default variance estimators have been chosen largely in accordance with the
#' procedures in
#' \href{https://github.com/acoppock/Green-Lab-SOP/blob/master/Green_Lab_SOP.pdf}{this manual}.
#' The default for the case
#' without clusters is the HC2 estimator and the default with clusters is the
#' analogous CR2 estimator. Users can easily replicate Stata standard errors in
#' the clustered or non-clustered case by setting \code{`se_type` = "stata"}.
#'
#' The function estimates the coefficients and standard errors in C++, using
#' the \code{RcppEigen} package. By default, we estimate the coefficients
#' using Column-Pivoting QR decomposition from the Eigen C++ library, although
#' users could get faster solutions by setting \code{`try_cholesky` = TRUE} to
#' use a Cholesky decomposition instead. This will likely result in quicker
#' solutions, but the algorithm does not reliably detect when there are linear
#' dependencies in the model and may fail silently if they exist.
#'
#' If \code{`fixed_effects`} are specified, both the outcome and design matrix
#' are centered using the method of alternating projections (Halperin 1962; Gaure 2013). Specifying
#' fixed effects in this way will result in large speed gains with standard error
#' estimators that do not need to invert the matrix of fixed effects. This means using
#' "classical", "HC0", "HC1", "CR0", or "stata" standard errors will be faster than other
#' standard error estimators. Be wary when specifying fixed effects that may result
#' in perfect fits for some observations or if there are intersecting groups across
#' multiple fixed effect variables (e.g. if you specify both "year" and "country" fixed effects
#' with an unbalanced panel where one year you only have data for one country).
#'
#' As with \code{`lm()`}, multivariate regression (multiple outcomes) will only admit
#' observations into the estimation that have no missingness on any outcome.
#'
#' @return An object of class \code{"lm_robust"}.
#'
#' The post-estimation commands functions \code{summary} and \code{\link{tidy}}
#' return results in a \code{data.frame}. To get useful data out of the return,
#' you can use these data frames, you can use the resulting list directly, or
#' you can use the generic accessor functions \code{coef}, \code{vcov},
#' \code{confint}, and \code{predict}. Marginal effects and uncertainty about
#' them can be gotten by passing this object to
#' \code{\link[margins]{margins}} from the \pkg{margins},
#' or to \code{emmeans} in the \pkg{emmeans} package.
#'
#' Users who want to print the results in TeX of HTML can use the
#' \code{\link[texreg]{extract}} function and the \pkg{texreg} package.
#'
#' If users specify a multivariate linear regression model (multiple outcomes),
#' then some of the below components will be of higher dimension to accommodate
#' the additional models.
#'
#' An object of class \code{"lm_robust"} is a list containing at least the
#' following components:
#'   \item{coefficients}{the estimated coefficients}
#'   \item{std.error}{the estimated standard errors}
#'   \item{statistic}{the t-statistic}
#'   \item{df}{the estimated degrees of freedom}
#'   \item{p.value}{the p-values from a two-sided t-test using \code{coefficients}, \code{std.error}, and \code{df}}
#'   \item{conf.low}{the lower bound of the \code{1 - alpha} percent confidence interval}
#'   \item{conf.high}{the upper bound of the \code{1 - alpha} percent confidence interval}
#'   \item{term}{a character vector of coefficient names}
#'   \item{alpha}{the significance level specified by the user}
#'   \item{se_type}{the standard error type specified by the user}
#'   \item{res_var}{the residual variance}
#'   \item{N}{the number of observations used}
#'   \item{k}{the number of columns in the design matrix (includes linearly dependent columns!)}
#'   \item{rank}{the rank of the fitted model}
#'   \item{vcov}{the fitted variance covariance matrix}
#'   \item{r.squared}{The \eqn{R^2},
#'   \deqn{R^2 = 1 - Sum(e[i]^2) / Sum((y[i] - y^*)^2),} where \eqn{y^*}
#'   is the mean of \eqn{y[i]} if there is an intercept and zero otherwise,
#'   and \eqn{e[i]} is the ith residual.}
#'   \item{adj.r.squared}{The \eqn{R^2} but penalized for having more parameters, \code{rank}}
#'   \item{fstatistic}{a vector with the value of the F-statistic with the numerator and denominator degrees of freedom}
#'   \item{weighted}{whether or not weights were applied}
#'   \item{call}{the original function call}
#'   \item{fitted.values}{the matrix of predicted means}
#' We also return \code{terms} and \code{contrasts}, used by \code{predict}. If \code{fixed_effects} are specified, then we return \code{proj_fstatistic}, \code{proj_r.squared}, and \code{proj_adj.r.squared}, which are model fit statistics that are computed on the projected model (after demeaning the fixed effects).
#'
#' @references
#' Abadie, Alberto, Susan Athey, Guido W Imbens, and Jeffrey Wooldridge. 2017. "A Class of Unbiased Estimators of the Average Treatment Effect in Randomized Experiments." arXiv Pre-Print. \url{https://arxiv.org/abs/1710.02926v2}.
#'
#' Bell, Robert M, and Daniel F McCaffrey. 2002. "Bias Reduction in Standard Errors for Linear Regression with Multi-Stage Samples." Survey Methodology 28 (2): 169-82.
#'
#' Gaure, Simon. 2013. "OLS with multiple high dimensional category variables." Computational Statistics & Data Analysis 66: 8-1. \doi{10.1016/j.csda.2013.03.024}
#'
#' Halperin, I. 1962. "The product of projection operators." Acta Scientiarum Mathematicarum (Szeged) 23(1-2): 96-99.
#'
#' MacKinnon, James, and Halbert White. 1985. "Some Heteroskedasticity-Consistent Covariance Matrix Estimators with Improved Finite Sample Properties." Journal of Econometrics 29 (3): 305-25. \doi{10.1016/0304-4076(85)90158-7}.
#'
#' Pustejovsky, James E, and Elizabeth Tipton. 2016. "Small Sample Methods for Cluster-Robust Variance Estimation and Hypothesis Testing in Fixed Effects Models." Journal of Business & Economic Statistics. Taylor & Francis. \doi{10.1080/07350015.2016.1247004}.
#'
#' Samii, Cyrus, and Peter M Aronow. 2012. "On Equivalencies Between Design-Based and Regression-Based Variance Estimators for Randomized Experiments." Statistics and Probability Letters 82 (2). \doi{10.1016/j.spl.2011.10.024}.
#'
#' @examples
#' set.seed(15)
#' library(fabricatr)
#' dat <- fabricate(
#'   N = 40,
#'   y = rpois(N, lambda = 4),
#'   x = rnorm(N),
#'   z = rbinom(N, 1, prob = 0.4)
#' )
#'
#' # Default variance estimator is HC2 robust standard errors
#' lmro <- lm_robust(y ~ x + z, data = dat)
#'
#' # Can tidy() the data in to a data.frame
#' tidy(lmro)
#' # Can use summary() to get more statistics
#' summary(lmro)
#' # Can also get coefficients three ways
#' lmro$coefficients
#' coef(lmro)
#' tidy(lmro)$estimate
#' # Can also get confidence intervals from object or with new 1 - `alpha`
#' lmro$conf.low
#' confint(lmro, level = 0.8)
#'
#' # Can recover classical standard errors
#' lmclassic <- lm_robust(y ~ x + z, data = dat, se_type = "classical")
#' tidy(lmclassic)
#'
#' # Can easily match Stata's robust standard errors
#' lmstata <- lm_robust(y ~ x + z, data = dat, se_type = "stata")
#' tidy(lmstata)
#'
#' # Easy to specify clusters for cluster-robust inference
#' dat$clusterID <- sample(1:10, size = 40, replace = TRUE)
#'
#' lmclust <- lm_robust(y ~ x + z, data = dat, clusters = clusterID)
#' tidy(lmclust)
#'
#' # Can also match Stata's clustered standard errors
#' lm_robust(
#'   y ~ x + z,
#'   data = dat,
#'   clusters = clusterID,
#'   se_type = "stata"
#' )
#'
#' # Works just as LM does with functions in the formula
#' dat$blockID <- rep(c("A", "B", "C", "D"), each = 10)
#'
#' lm_robust(y ~ x + z + factor(blockID), data = dat)
#'
#' # Weights are also easily specified
#' dat$w <- runif(40)
#'
#' lm_robust(
#'   y ~ x + z,
#'   data = dat,
#'   weights = w,
#'   clusters = clusterID
#' )
#'
#' # Subsetting works just as in `lm()`
#' lm_robust(y ~ x, data = dat, subset = z == 1)
#'
#' # One can also choose to set the significance level for different CIs
#' lm_robust(y ~ x + z, data = dat, alpha = 0.1)
#'
#' # We can also specify fixed effects
#' # Speed gains with fixed effects are greatest with "stata" or "HC1" std.errors
#' tidy(lm_robust(y ~ z, data = dat, fixed_effects = ~ blockID, se_type = "HC1"))
#'
#' \dontrun{
#'   # Can also use 'margins' or 'emmeans' package if you have them installed
#'   # to get marginal effects
#'   library(margins)
#'   lmrout <- lm_robust(y ~ x + z, data = dat)
#'   summary(margins(lmrout))
#'
#'   # Can output results using 'texreg'
#'   library(texreg)
#'   texreg(lmrout)
#'
#'   # Using emmeans to obtain covariate-adjusted means
#'   library(emmeans)
#'   fiber.rlm <- lm_robust(strength ~ diameter + machine, data = fiber)
#'   emmeans(fiber.rlm, "machine")
#' }
#'
#' @export
lm_robust <- function(formula,
                      data,
                      weights,
                      subset,
                      clusters,
                      fixed_effects,
                      se_type = NULL,
                      ci = TRUE,
                      alpha = .05,
                      return_vcov = TRUE,
                      try_cholesky = FALSE) {
  datargs <- enquos(
    formula = formula,
    weights = weights,
    subset = subset,
    cluster = clusters,
    fixed_effects = fixed_effects
  )
  data <- enquo(data)
  model_data <- clean_model_data(data = data, datargs)

  fes <- is.character(model_data[["fixed_effects"]])
  if (fes) {
    yoriginal <- model_data[["outcome"]]
    Xoriginal <- model_data[["design_matrix"]]
    model_data <- demean_fes(model_data)
    attr(model_data$fixed_effects, "fe_rank") <- sum(model_data[["fe_levels"]]) + 1
  } else {
    Xoriginal <- NULL
    yoriginal <- NULL
  }

  return_list <-
    lm_robust_fit(
      y = model_data$outcome,
      X = model_data$design_matrix,
      yoriginal = yoriginal,
      Xoriginal = Xoriginal,
      weights = model_data$weights,
      cluster = model_data$cluster,
      fixed_effects = model_data$fixed_effects,
      ci = ci,
      se_type = se_type,
      alpha = alpha,
      return_vcov = return_vcov,
      try_cholesky = try_cholesky,
      has_int = attr(model_data$terms, "intercept"),
      iv_stage = list(0)
    )

  return_list <- lm_return(
    return_list,
    model_data = model_data,
    formula = formula
  )

  return_list[["call"]] <- match.call()

  return(return_list)
}
