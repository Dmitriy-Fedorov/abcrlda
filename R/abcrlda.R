
#' Asymptotically Bias-Corrected Regularized Linear Discriminant Analysis for Cost-Sensitive Binary Classification
#' @description Performs symptotically Bias-Corrected Regularized Linear Discriminant Analysis
#'
#' @param x Matrix or data.frame of observations.
#' @param grouping Grouping variable. A vector of numeric values 0 and 1 is recommended. Length has to correspond to nrow(x).
#' @param gamma regularization parameter
#' @param cost_10 parameter that controls prioretization of classes.
#' It's value should be between 0 and 1 (0 < cost_10 < 1)
#' Values bigger than 0.5 prioretizes correct classification of 0 class while values less than 0.5 prioretizes 1 class
#'
#' @return An object of class "rrlda" is returned which can be used for class prediction (see predict())
#' \describe{
#'   \item{a}{Slope of a discriminant hyperplane. W(x) = a'x + m}
#'   \item{m}{Bias term. W(x) = a'x + m}
#'   \item{cost_10}{asd}
#'   \item{gamma}{asd}
#'   \item{Ghat0, Ghat1}{asd}
#'   \item{Dhat}{asd}
#'   \item{omegaopt}{asd}
#'   \item{lev}{Levels. Corresponds to the groups.}
#' }
#' @export
#' @examples
#' data(iris)
#' traindata = iris[ which(iris[,ncol(iris)]=='virginica' |
#'                         iris[,ncol(iris)]=="versicolor"), 1:4]
#' trainlabel = factor(iris[ which(iris[,ncol(iris)]=='virginica' |
#'                                 iris[,ncol(iris)]=="versicolor"), 5])
#' rr <- abcrlda(traindata, trainlabel, gamma = 0.5, cost_10 = 0.75)
#' predict(rr, traindata)
abcrlda <- function(x, grouping, gamma, cost_10 = 0.5){  # cost_01 = 1 -  cost_10

  ## check requirements
  if (is.null(dim(x)))
    stop("'x' is not a matrix or data.frame")

  x <- as.matrix(x)
  if (any(!is.finite(x)))
    stop("Infinite, NA or NaN values in 'x'")



  p <- ncol(x)  # number of dimensions
  n <- nrow(x)  # number of samples

  if (n != length(grouping))
    stop("nrow(x) and length(grouping) are different")

  if(!is.factor(grouping))
    grouping <- as.factor(grouping)

  lev <- levels(grouping)
  k <- nlevels(grouping)

  if (k != 2)
    stop("number of groups != 2, this is binary classifier")

  x0 <- x[grouping == lev[1],]
  x1 <- x[grouping == lev[2],]

  n0 <- nrow(x0) # number of samples in x0
  n1 <- nrow(x1) # number of samples in x1

  S0 <- cov(x0)
  S1 <- cov(x1)
  S <- ((n0-1)*S0+(n1-1)*S1)/(n0+n1-2)
  Hinv <- (diag(ncol(x)) + gamma*S)
  H <- solve(Hinv)
  # ------ -- - - -- - - - - -
  m0 <- colMeans(x0)
  m1 <- colMeans(x1)
  mdif <- m0 - m1
  msum <- m1 + m0
  a <- H %*% mdif
  mt <- t(a) %*% msum
  m_rlda <- -0.5*mt - log((1-cost_10)/cost_10)/gamma
  # ------- omega optimal calculation -------------
  traceH <- sum(diag(H))  # sum of diagonal elements in H
  deltahat <- (p/(n0+n1-2)-traceH/(n0+n1-2)) / (gamma*(1 - p/(n0+n1-2) + traceH/(n0+n1-2)))
  G0 <- 0.5*t(m0-m1) %*% H %*% mdif - log((1-cost_10)/cost_10)/gamma
  G1 <- 0.5*t(m1-m0) %*% H %*% mdif - log((1-cost_10)/cost_10)/gamma
  Ghat0 <- G0 - ((n0+n1-2)/n0) * deltahat
  Ghat1 <- G1 + ((n0+n1-2)/n1) * deltahat
  D <- t(a) %*% S %*% a
  # D <- t(mdif) %*% H %*% S %*% H %*% mdif  # no difference, more explicit
  Dhat <- D*(1+gamma*deltahat)^2
  omegaopt <- gamma * (Dhat*log((1-cost_10)/cost_10)/(Ghat1-Ghat0) - 0.5*(Ghat0+Ghat1))
  m_abcrlda <- as.numeric(m_rlda + omegaopt/gamma)
  return(structure(list(a = a,
                        m = m_abcrlda,
                        cost_10 = cost_10,
                        gamma = gamma,
                        Ghat0 = Ghat0,
                        Ghat1 = Ghat1,
                        Dhat = Dhat,
                        omegaopt = omegaopt,
                        lev=lev), class="abcrlda"))

}


#' Class Prediction for abcrlda objects
#' @description Computes class predictions for new data based on a given abcrlda object
#' @param object An object of class "abcrlda".
#' @param x New data for which the classes are to predict
#' @param ... Argument used by generic function predict(object, x, ...).
#'
#' @return
#' @param class Class prediction for each observation.
#' @param raw Raw values.
#' @export
#'
#' @examples
#' data(iris)
#' traindata = iris[ which(iris[,ncol(iris)]=='virginica' |
#'                         iris[,ncol(iris)]=="versicolor"), 1:4]
#' trainlabel = factor(iris[ which(iris[,ncol(iris)]=='virginica' |
#'                                 iris[,ncol(iris)]=="versicolor"), 5])
#' rr <- abcrlda(traindata, trainlabel, gamma = 0.5, cost_10 = 0.75)
#' predict(rr, traindata)
predict.abcrlda <- function(object, x, ...){
  ## check requirements
  if(class(object) != "abcrlda")
    stop("object has to be of type abcrlda")

  if(!is.vector(x) && !is.matrix(x) && !is.data.frame(x) && !is.data.frame(x)) #  && !is.data.frame(x)
    stop("'x' has to be a vector, matrix or data.frame")

  x <- as.matrix(x)
  pred <- as.numeric(x %*% object$a + object$m <= 0)
  cl <- object$lev[pred + 1]

  return(list(class=cl, raw=pred))  # object$cost

}



