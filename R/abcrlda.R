
#' Title
#'
#' @param x
#' @param grouping
#' @param gamma
#' @param cost_10
#'
#' @return
#' @export
#'
#' @examples
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


#' Title
#'
#' @param object
#' @param x
#'
#' @return
#' @export
#'
#' @examples
predict.abcrlda <- function(object, x){
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



