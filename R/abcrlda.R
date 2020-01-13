
#' Asymptotically Bias-Corrected Regularized Linear Discriminant Analysis
#' for Cost-Sensitive Binary Classification
#' @description Constructs Asymptotically Bias-Corrected Regularized Linear Discriminant Analysis.
#' @param x Input matrix or data.frame of dimension \code{nobs x nvars}; each row is an feature vector.
#' @param y A numeric vector or factor of class labels. Factor should have either two levels or be
#'   a vector with two distinct values.
#'   If \code{y} is presented as a vector, it will be coerced into a factor.
#'   Length of \code{y} has to correspond to number of samples in \code{x}.
#' @param gamma Regularization parameter \eqn{\gamma}{gamma} in the ABC-RLDA discriminant function given by:
#'   \deqn{W_{ABC}^{RLDA} = \gamma (x-\frac{\bar{x}_0 +
#'   \bar{x}_1}{2})^T H (\bar{x}_0 - \bar{x}_1)
#'   - log(\frac{C_{01}}{C_{10}}) + \hat{\omega}_{opt}}{W_ABCRLDA = gamma (x - (x0 + x1)/2) H (x0 - x1) + log(C_01/C_10) + omega_opt}
#'   \deqn{H = (I_p + \gamma \hat{\Sigma})^{-1}}{H = (I_p + gamma Sigma_hat)^-1}
#'   Formulas and derivations for parameters used in above equation can be found in the article under reference section.
#' @param cost Parameter that controls the overall misclassification costs.
#'  This is a vector of length 1 or 2 where the first value is \eqn{C_{10}}{C_10} (represents the cost of assigning label 1 when the true label is 0)
#'  and the second value, if provided, is \eqn{C_{01}}{C_01} (represents the cost of assigning label 0 when the true label is 1).
#'  The default setting is c(0.5, 0.5), so both classes have equal misclassification costs
#'
#'  If a single value is provided, it should be normalized to lie between 0 and 1 (but not including 0 or 1).
#'  This value will be assigned to \eqn{C_{10}}{C_10} while
#'  \eqn{C_{01}}{C_01} will be equal to \eqn{(1 - C_{10})}{1 - C_10}.
##  In a vector of length 1, values bigger than 0.5 prioretizes correct classification of 0 class while values less than 0.5 prioretizes 1 class.
#' @param bias_correction Takes in a boolean value.
#'   If \code{bias_correction} is TRUE, then asymptotic bias correction will be performed.
#'   Otherwise, (if \code{bias_correction} is FALSE) asymptotic bias correction will not be performed and
#'   the ABCRLDA is the classical RLDA.
#'   The default is TRUE.
#' @return An object of class "abcrlda" is returned which can be used for class prediction (see predict()).
#'   \item{a}{Coefficient vector of a discriminant hyperplane: W(\strong{x}) = \strong{a}' \strong{x} + m.}
#'   \item{m}{Intercept of discriminant hyperplane: W(\strong{x}) = \strong{a}'\strong{x} + m.}
#'   \item{cost}{Vector of cost values that are used to construct ABC-RLDA.}
#'   \item{ncost}{Normilized cost such that \eqn{C_{10}}{C_10} + \eqn{C_{01}}{C_01} = 1.}
#'   \item{gamma}{Regularization parameter value used in ABC_RLDA discriminant function.}
#'   \item{lev}{Levels corresponding to the labels in y.}
#' @section Reference:
#'   A. Zollanvari, M. Abdirash, A. Dadlani and B. Abibullaev,
#'   "Asymptotically Bias-Corrected Regularized Linear Discriminant Analysis for Cost-Sensitive
#'   Binary Classification," in IEEE Signal Processing Letters, vol. 26, no. 9, pp. 1300-1304,
#'   Sept. 2019. doi: 10.1109/LSP.2019.2918485
#'   URL: \url{http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8720003&isnumber=8770167}
#' @export
#' @family functions in the package
#' @example inst/examples/example_abcrlda.R
abcrlda <- function(x, y, gamma=1, cost=c(0.5, 0.5), bias_correction=TRUE){

  ## check requirements
  if (is.null(dim(x)))
    stop("'x' is not a matrix or data.frame")

  x <- as.matrix(x)
  if (any(!is.finite(x)))
    stop("Infinite, NA or NaN values in 'x'")

  p <- ncol(x)  # number of dimensions
  n <- nrow(x)  # number of samples

  if (n != length(y))
    stop("nrow(x) and length(grouping) are different")

  if (!is.factor(y))
    y <- as.factor(y)

  lev <- levels(y)
  k <- nlevels(y)

  if (k != 2)
    stop("number of groups != 2, this is binary classifier")

  x0 <- x[y == lev[1], ]
  x1 <- x[y == lev[2], ]
  n0 <- nrow(x0) # number of samples in x0
  n1 <- nrow(x1) # number of samples in x1
  S0 <- stats::cov(x0)
  S1 <- stats::cov(x1)
  S <- ( (n0 - 1) * S0 + (n1 - 1) * S1) / (n0 + n1 - 2)
  Hinv <- (diag(ncol(x)) + gamma * S)
  H <- solve(Hinv)
  # ----------------------------------------------------
  if (length(cost) == 1){
    if (cost >= 1 | cost <= 0)
      stop("While providing single valued vector
           cost should be between 0 and 1 (not including)")
    cost <- c(cost, 1 - cost)
  }
  if (length(cost) != 2)
    stop("cost vector should be of length 1 or 2, this is binary classifier")

  m0 <- colMeans(x0)
  m1 <- colMeans(x1)
  mdif <- m0 - m1
  msum <- m1 + m0
  a <- H %*% mdif
  mt <- t(a) %*% msum
  m_rlda <- -0.5 * mt - log( cost[2] / cost[1]) / gamma
  # ------- omega optimal calculation -------------
  trace_h <- sum(diag(H))  # sum of diagonal elements in H
  deltahat <- (p / (n0 + n1 - 2) - trace_h / (n0 + n1 - 2)) /
              (gamma * (1 - p / (n0 + n1 - 2) + trace_h / (n0 + n1 - 2)))
  G0 <- 0.5 * t(m0 - m1) %*% H %*% mdif - log( cost[2] / cost[1]) / gamma
  G1 <- 0.5 * t(m1 - m0) %*% H %*% mdif - log( cost[2] / cost[1]) / gamma
  Ghat0 <- G0 - ( (n0 + n1 - 2) / n0) * deltahat
  Ghat1 <- G1 + ( (n0 + n1 - 2) / n1) * deltahat
  D <- t(a) %*% S %*% a
  # D <- t(mdif) %*% H %*% S %*% H %*% mdif  # no difference, more explicit
  Dhat <- D * (1 + gamma * deltahat) ^ 2
  omegaopt <- gamma * (Dhat * log( cost[2] / cost[1]) / (Ghat1 - Ghat0) -
              0.5 * (Ghat0 + Ghat1))
  if (bias_correction)
    m_abcrlda <- as.numeric(m_rlda + omegaopt / gamma)
  else
    m_abcrlda <- as.numeric(m_rlda)
  return(structure(list(a = a,
                        m = m_abcrlda,
                        cost = cost,
                        ncost = c(cost[1] / sum(cost), cost[2] / sum(cost)),
                        gamma = gamma,
                        D = D,
                        G0 = G0,
                        G1 = G1,
                        Ghat0 = Ghat0,
                        Ghat1 = Ghat1,
                        Dhat = Dhat,
                        omegaopt = omegaopt,
                        bias_correction = bias_correction,
                        lev = lev), class = "abcrlda"))

}


#' Class Prediction for abcrlda objects
#' @description Classifies observations based on a given abcrlda object.
#' @param object An object of class "abcrlda".
#' @param newx Matrix of new values for x at which predictions are to be made.
# @param out_type Determines a type of output. Two type of input could be provided.
#   If "class" value is provided this will return factor with levels corresponding to lev stored in object.
#   If "raw" value is provided this will return numeric vector with values obtained from discriminant function.
#' @param ... Argument used by generic function predict(object, x, ...).
#'
#' @return
#'  Returns factor vector with predictions (i.e., assigned labels) for each observation. Factor levels are inhereted from the object variable.
#' @export
#' @family functions in the package
#'
#' @example inst/examples/example_abcrlda.R
#' @inheritSection abcrlda Reference
predict.abcrlda <- function(object, newx, ...){
  ## check requirements
  if (class(object) != "abcrlda")
    stop("object has to be of type abcrlda")

  if (!is.vector(newx) && !is.matrix(newx) &&
      !is.data.frame(newx) && !is.data.frame(newx))
    stop("'x' has to be a vector, matrix or data.frame")

  newx <- as.matrix(newx, drop = FALSE)
  pred <- as.numeric(newx %*% object$a + object$m <= 0) + 1
  cl <- object$lev[pred]
  cl <- factor(cl, levels=object$lev)
  # cl <- as.factor(cl)
  # levels(cl) <- object$lev
  return(cl)
}


#' Risk Calculate
#' @description Estimates risk and error by applying a constructed classifier (an object of class abcrlda) to a given set of observations.
#' @param object An object of class "abcrlda".
#' @param x_true Matrix of values for x for which true class labels are known.
#' @param y_true A numeric vector or factor of true class labels. Factor should have either two levels or be a vector with two distinct values.
#'   If \code{y_true} is presented as a vector, it will be coerced into a factor.
#'   Length of \code{y_true} has to correspond to number of samples in \code{x_test}.
#'
#' @return A list of parameters where
#'   \item{actual_err0}{Error rate for class 0.}
#'   \item{actual_err1}{Error rate for class 1.}
#'   \item{actual_errTotal}{Error rate overall.}
#'   \item{actual_normrisk}{Risk value normilized to be between 0 and 1.}
#'   \item{actual_risk}{Risk value without normilization.}
#' @export
#' @family functions in the package
#'
#' @example inst/examples/example_risk_calculate.R
risk_calculate <- function(object, x_true, y_true){
  if (!is.factor(y_true))
    y_true <- as.factor(y_true)

  if (!all(levels(y_true) == object$lev))
    stop("object and y_true have different levels")

  test0 <- x_true[y_true == object$lev[1], ]
  test1 <- x_true[y_true == object$lev[2], ]
  res0 <- as.numeric(stats::predict(object, test0)) - 1
  res1 <- as.numeric(stats::predict(object, test1)) - 1
  nerr0 <- sum(res0)
  nerr1 <- sum(!res1)
  err0 <- nerr0 / length(res0)
  err1 <- nerr1 / length(res1)
  return(structure(list(actual_err0 = err0,
                        actual_err1 = err1,
                        actual_errTotal = (nerr0 + nerr1) / length(y_true),
                        actual_normrisk = err0 * object$ncost[1] + err1 * object$ncost[2],
                        actual_risk = err0 * object$cost[1] + err1 * object$cost[2])))
}
