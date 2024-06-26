% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimatr_horvitz_thompson.R
\name{horvitz_thompson}
\alias{horvitz_thompson}
\title{Horvitz-Thompson estimator for two-armed trials}
\usage{
horvitz_thompson(
  formula,
  data,
  blocks,
  clusters,
  simple = NULL,
  condition_prs,
  condition_pr_mat = NULL,
  ra_declaration = NULL,
  subset,
  condition1 = NULL,
  condition2 = NULL,
  se_type = c("youngs", "constant", "none"),
  ci = TRUE,
  alpha = 0.05,
  return_condition_pr_mat = FALSE
)
}
\arguments{
\item{formula}{an object of class formula, as in \code{\link{lm}}, such as
\code{Y ~ Z} with only one variable on the right-hand side, the treatment.}

\item{data}{A data.frame.}

\item{blocks}{An optional bare (unquoted) name of the block variable. Use
for blocked designs only. See details.}

\item{clusters}{An optional bare (unquoted) name of the variable that
corresponds to the clusters in the data; used for cluster randomized
designs. For blocked designs, clusters must be within blocks.}

\item{simple}{logical, optional. Whether the randomization is simple
(TRUE) or complete (FALSE). This is ignored if \code{blocks} are specified,
as all blocked designs use complete randomization, or either
\code{ra_declaration} or \code{condition_pr_mat} are passed. Otherwise, it
defaults to \code{TRUE}.}

\item{condition_prs}{An optional bare (unquoted) name of the variable with
the condition 2 (treatment) probabilities. See details. May also use a single
number for the condition 2 probability if it is constant.}

\item{condition_pr_mat}{An optional 2n * 2n matrix of marginal and joint
probabilities of all units in condition1 and condition2. See details.}

\item{ra_declaration}{An object of class \code{"ra_declaration"}, from
the \code{\link[randomizr]{declare_ra}} function in the \pkg{randomizr}
package. This is the third way that one can specify a design for this
estimator. Cannot be used along with any of \code{condition_prs},
\code{blocks}, \code{clusters}, or \code{condition_pr_mat}. See details.}

\item{subset}{An optional bare (unquoted) expression specifying a subset of
observations to be used.}

\item{condition1}{value in the treatment vector of the condition
to be the control. Effects are
estimated with \code{condition1} as the control and \code{condition2} as the
treatment. If unspecified, \code{condition1} is the "first" condition and
\code{condition2} is the "second" according to levels if the treatment is a
factor or according to a sortif it is a numeric or character variable (i.e
if unspecified and the treatment is 0s and 1s, \code{condition1} will by
default be 0 and \code{condition2} will be 1). See the examples for more.}

\item{condition2}{value in the treatment vector of the condition to be the
treatment. See \code{condition1}.}

\item{se_type}{can be one of \code{c("youngs", "constant", "none")} and corresponds
the estimator of the standard errors. Default estimator uses Young's
inequality (and is conservative) while the other uses a constant treatment
effects assumption and only works for simple randomized designs at the
moment. If "none" then standard errors will not be computed which may speed up run time if only the point estimate is required.}

\item{ci}{logical. Whether to compute and return p-values and
confidence intervals, TRUE by default.}

\item{alpha}{The significance level, 0.05 by default.}

\item{return_condition_pr_mat}{logical. Whether to return the condition
probability matrix. Returns NULL if the design is simple randomization,
FALSE by default.}
}
\value{
Returns an object of class \code{"horvitz_thompson"}.

The post-estimation commands functions \code{summary} and \code{\link{tidy}}
return results in a \code{data.frame}. To get useful data out of the return,
you can use these data frames, you can use the resulting list directly, or
you can use the generic accessor functions \code{coef} and
\code{confint}.

An object of class \code{"horvitz_thompson"} is a list containing at
least the following components:

  \item{coefficients}{the estimated difference in totals}
  \item{std.error}{the estimated standard error}
  \item{statistic}{the z-statistic}
  \item{df}{the estimated degrees of freedom}
  \item{p.value}{the p-value from a two-sided z-test using \code{coefficients} and \code{std.error}}
  \item{conf.low}{the lower bound of the \code{1 - alpha} percent confidence interval}
  \item{conf.high}{the upper bound of the \code{1 - alpha} percent confidence interval}
  \item{term}{a character vector of coefficient names}
  \item{alpha}{the significance level specified by the user}
  \item{nobs}{the number of observations used}
  \item{outcome}{the name of the outcome variable}
  \item{condition_pr_mat}{the condition probability matrix if \code{return_condition_pr_mat} is TRUE}
}
\description{
Horvitz-Thompson estimators that are unbiased for designs in
which the randomization scheme is known
}
\details{
This function implements the Horvitz-Thompson estimator for
treatment effects for two-armed trials. This estimator is useful for estimating unbiased
treatment effects given any randomization scheme as long as the
randomization scheme is known.

In short, the Horvitz-Thompson estimator essentially reweights each unit
by the probability of it being in its observed condition. Pivotal to the
estimation of treatment effects using this estimator are the marginal
condition probabilities (i.e., the probability that any one unit is in
a particular treatment condition). Pivotal to estimating the variance
whenever the design is more complicated than simple randomization are the
joint condition probabilities (i.e., the probabilities that any two units
have a particular set of treatment conditions, either the same or
different). The estimator we provide here considers the case with two
treatment conditions.

Users interested in more details can see the
\href{https://declaredesign.org/r/estimatr/articles/mathematical-notes.html}{mathematical notes}
for more information and references, or see the references below.

There are three distinct ways that users can specify the design to the
function. The preferred way is to use
the \code{\link[randomizr]{declare_ra}} function in the \pkg{randomizr}
package. This function takes several arguments, including blocks, clusters,
treatment probabilities, whether randomization is simple or not, and more.
Passing the outcome of that function, an object of class
\code{"ra_declaration"} to the \code{ra_declaration} argument in this function,
will lead to a call of the \code{\link{declaration_to_condition_pr_mat}}
function which generates the condition probability matrix needed to
estimate treatment effects and standard errors. We provide many examples
below of how this could be done.

The second way is to pass the names of vectors in your \code{data} to
\code{condition_prs}, \code{blocks}, and \code{clusters}. You can further
specify whether the randomization was simple or complete using the \code{simple}
argument. Note that if \code{blocks} are specified the randomization is
always treated as complete. From these vectors, the function learns how to
build the condition probability matrix that is used in estimation.

In the case
where \code{condition_prs} is specified, this function assumes those
probabilities are the marginal probability that each unit is in condition2
and then uses the other arguments (\code{blocks}, \code{clusters},
\code{simple}) to learn the rest of the design. If users do not pass
\code{condition_prs}, this function learns the probability of being in
condition2 from the data. That is, none of these arguments are specified,
we assume that there was a simple randomization where the probability
of each unit being in condition2 was the average of all units in condition2.
Similarly, we learn the block-level probability of treatment within
\code{blocks} by looking at the mean number of units in condition2 if
\code{condition_prs} is not specified.

The third way is to pass a \code{condition_pr_mat} directly. One can
see more about this object in the documentation for
\code{\link{declaration_to_condition_pr_mat}} and
\code{\link{permutations_to_condition_pr_mat}}. Essentially, this 2n * 2n
matrix allows users to specify marginal and joint marginal probabilities
of units being in conditions 1 and 2 of arbitrary complexity. Users should
only use this option if they are certain they know what they are doing.
}
\examples{

# Set seed
set.seed(42)

# Simulate data
n <- 10
dat <- data.frame(y = rnorm(n))

library(randomizr)

#----------
# 1. Simple random assignment
#----------
dat$p <- 0.5
dat$z <- rbinom(n, size = 1, prob = dat$p)

# If you only pass condition_prs, we assume simple random sampling
horvitz_thompson(y ~ z, data = dat, condition_prs = p)
# Assume constant effects instead
horvitz_thompson(y ~ z, data = dat, condition_prs = p, se_type = "constant")

# Also can use randomizr to pass a declaration
srs_declaration <- declare_ra(N = nrow(dat), prob = 0.5, simple = TRUE)
horvitz_thompson(y ~ z, data = dat, ra_declaration = srs_declaration)

#----------
# 2. Complete random assignment
#----------

dat$z <- sample(rep(0:1, each = n/2))
# Can use a declaration
crs_declaration <- declare_ra(N = nrow(dat), prob = 0.5, simple = FALSE)
horvitz_thompson(y ~ z, data = dat, ra_declaration = crs_declaration)
# Can precompute condition_pr_mat and pass it
# (faster for multiple runs with same condition probability matrix)
crs_pr_mat <- declaration_to_condition_pr_mat(crs_declaration)
horvitz_thompson(y ~ z, data = dat, condition_pr_mat = crs_pr_mat)

#----------
# 3. Clustered treatment, complete random assigment
#-----------
# Simulating data
dat$cl <- rep(1:4, times = c(2, 2, 3, 3))
dat$prob <- 0.5
clust_crs_decl <- declare_ra(N = nrow(dat), clusters = dat$cl, prob = 0.5)
dat$z <- conduct_ra(clust_crs_decl)
# Easiest to specify using declaration
ht_cl <- horvitz_thompson(y ~ z, data = dat, ra_declaration = clust_crs_decl)
# Also can pass the condition probability and the clusters
ht_cl_manual <- horvitz_thompson(
  y ~ z,
  data = dat,
  clusters = cl,
  condition_prs = prob,
  simple = FALSE
)
ht_cl
ht_cl_manual

# Blocked estimators specified similarly

#----------
# More complicated assignment
#----------

# arbitrary permutation matrix
possible_treats <- cbind(
  c(1, 1, 0, 1, 0, 0, 0, 1, 1, 0),
  c(0, 1, 1, 0, 1, 1, 0, 1, 0, 1),
  c(1, 0, 1, 1, 1, 1, 1, 0, 0, 0)
)
arb_pr_mat <- permutations_to_condition_pr_mat(possible_treats)
# Simulating a column to be realized treatment
dat$z <- possible_treats[, sample(ncol(possible_treats), size = 1)]
horvitz_thompson(y ~ z, data = dat, condition_pr_mat = arb_pr_mat)

}
\references{
Aronow, Peter M, and Joel A Middleton. 2013. "A Class of Unbiased Estimators of the Average Treatment Effect in Randomized Experiments." Journal of Causal Inference 1 (1): 135-54. \doi{10.1515/jci-2012-0009}.

Aronow, Peter M, and Cyrus Samii. 2017. "Estimating Average Causal Effects Under Interference Between Units." Annals of Applied Statistics, forthcoming. \url{https://arxiv.org/abs/1305.6156v3}.

Middleton, Joel A, and Peter M Aronow. 2015. "Unbiased Estimation of the Average Treatment Effect in Cluster-Randomized Experiments." Statistics, Politics and Policy 6 (1-2): 39-75. \doi{10.1515/spp-2013-0002}.
}
\seealso{
\code{\link[randomizr]{declare_ra}}
}
