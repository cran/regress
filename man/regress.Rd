\name{regress}
\alias{regress}

\title{Fit a Gaussian Linear Model with Linear Covariance Structure}

\description{Fits Gaussian linear models in which the covariance
  structure can be expressed as a linear combination of known matrices.
  For example, block effects models and spatial models that include a
  nugget effect.  Fits model by maximising the residual log likelihood,
  also known as the REML log likelihood or restricted log likelihood.
  Uses a Newton-Raphson algorithm to maximise the residual log
  likelihood.}

\synopsis{regress(formula, Vformula, identity=TRUE, start=NULL, fraction=1, pos, verbose=0, gamVals=NULL, maxcyc=50, tol=1e-4, data, print.level)}

\usage{regress(formula, Vformula, identity=TRUE, start=NULL, fraction=1,
pos, verbose=0, gamVals=NULL, maxcyc=50, tol=1e-4, data, print.level)}

\arguments{
  \item{formula}{ a symbolic description of the model to be fitted. The
    details of model specification are the same as for \code{lm}}
  \item{Vformula}{Specifies the matrices to include in the covariance
    structure.  Each term is either a symmetric matrix, or a factor.
    Independent Gaussian random effects are included by passing the
    corresponding block factor.}
  \item{identity}{Logical variable, includes the identity as the
    final matrix of the covariance structure.  Default is TRUE}
  \item{start}{Specify the variance components at which the
    Newton-Raphson algorithm starts.  Default value is
    \code{rep(var(y),k)}.}
  \item{fraction}{The proportion of each step to take.  Default value is
    1. Useful to prevent taking huge steps in the first few iterations.}
  \item{pos}{logical vector of length k, where k is the number of
    matrices in the covariance structure.  Indicates which variance
    components are positive (TRUE) and which are real (FALSE).
    Important for multivariate problems.}
  \item{verbose}{Controls level of time output, takes values 0, 1 or
    2, Default is 0, level 1 gives parameter estimates and value of log
    likelihood at each stage.}
  \item{gamVals}{When k=2, the marginial log likelihood based on the
    residual configuration statistic (see Tunnicliffe Wilson(1989)), is
    evaluated first at \code{(1-gam) V1 + gam V2} for each value of
    \code{gam} in \code{gamVals}, a set of values from the unit
    interval.  Subsequently the Newton-Raphson algorithm is started at
    variance components corresponding the the value of \code{gam} that
    has the highest marginal log likelihood.  This is overridden if
    \code{start} is specified.}
  \item{maxcyc}{Maximum number of cycles allowed.  Default value is
    50. A warning is output to the screen if this is reached before
    convergence.}
  \item{tol}{Convergence criteria.  If the change in residual log
    likelihood for one cycle is less than \code{tol} the algorithm
    finishes.  Default value is \code{1e-4}.}
  \item{data}{an optional data frame containing the variables in the model.
    By default the variables are taken from 'environment(formula)',
    typically the environment from which 'regress' is called.}
  \item{print.level}{Deprecated}
}

\details{As the code is running it outputs the variance components, and
  the residual log likelihood at each iteration.

  To avoid confusion over terminology.  I define variance components to
  be the multipliers of the matrices and variance parameters to the
  parameter space over which the Newton-Raphson algorithm is run.  I can
  force a component to be positive be defining the corresponding variance
  parameter on the log scale.
  
  All output to the screen is for variance components (i.e. the
  multiples of the matrices).  Values for \code{start} are on the
  variance component scale.  Use \code{pos} to force certain variance
  components to be positive.

  NOTE: The final stage of the algorithm converts the estimates of the
  variance components and the Fisher Information to the usual linear
  scale, i.e. as if \code{pos} were a vector of zeroes.

  NOTE: No \code{predict} functionality is provided with regress due to
  some ambiguity.  Are we predicting conditional on the observed data.
  Are we predicting observations from the fitted model itself?  It is
  all normal anyway so it is straightforward, see our paper on regress.
  }

\value{
  \item{trace}{Matrix with one row for each iteration of algorithm.
    Each row contains the residual log likelihood, marginal log
    likelihood, variance parameters and increments.}
  \item{llik}{Value of the marginal log likelihood at the point
  of convergence.}
  \item{cycle}{Number of cycles to convergence.}
  \item{rdf}{Residual degrees of freedom.}
  \item{beta}{Estimate of the linear effects.}
  \item{beta.cov}{Estimate of the covariance structure for terms in
    beta.}
  \item{beta.se}{Standard errors for terms in beta.}
  \item{sigma}{Variance component estimates, interpretation does not depend on
    value of \code{pos}}
  \item{sigma.cov}{Covariance matrix for the variance component
    estimates based on the Fisher Information at the point of
    convergence.}
  \item{W}{Inverse of covariance matrix at point of convergence.}
  \item{Q}{$I - X^T (X^T W X)^-1 X^T W$ at point of convergence.}
  \item{fitted}{$X beta$, the fitted values.}
  \item{predicted}{If \code{identity=TRUE}, decompose y into the part
  associated with the identity and that assosicated with the rest of the
  variance structure, this second part is the predicted values.  If
  $Sigma = V1 + V2$ at point of convergence then y = V1 W y + V2 W y is
  the decomposition.}
  \item{pos}{Indicator for the scale for each variance parameter.}
  \item{Vnames}{Names associated with each variance component, used in
    \code{print.regress}.}
  \item{formula}{Copy of formula}
  \item{Vformula}{Updated version of Vformula to include identity if
    necessary}
  \item{model}{Response, covariates and matrices/factors to be used for
    model fitting}
}

\references{
  G. Tunnicliffe Wilson, "On the use of marginal likelihood in time
  series model estimation."  JRSS B, Vol 51, No 1, 1989, 15-27.
  D. Clifford and P. McCullagh, "The regress Library" R News.
}

\author{David Clifford, Peter McCullagh.}

\examples{
  ## Example of Random Effects model from Venables and Ripley, page 205
  library("nlme")
  ##library("regress")
  data(Oats)
  names(Oats) <- c("B","V","N","Y")
  Oats$N <- as.factor(Oats$N)

  ## Using regress
  oats.reg <- regress(Y~N+V,~B+I(B:V),identity=TRUE,verbose=1,data=Oats)
  summary(oats.reg)

  ## Using lme
  oats.lme <- lme(Y~N+V,random=~1|B/V,data=Oats,method="REML")
  summary(oats.lme)
}

\keyword{regression}
\keyword{models}
\keyword{multivariate}

