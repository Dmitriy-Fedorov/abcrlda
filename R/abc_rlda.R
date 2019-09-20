
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

  x0 <- x[grouping == 0,]
  x1 <- x[grouping == 1,]
  p <- ncol(x0)  # number of dimensions
  n0 <- nrow(x0) # number of samples in x0
  n1 <- nrow(x1) # number of samples in x1
  X0 <- as.matrix(x0)
  X1 <- as.matrix(x1)
  S0 <- cov(X0)
  S1 <- cov(X1)
  S <- ((n0-1)*S0+(n1-1)*S1)/(n0+n1-2)
  Hinv <- (diag(ncol(x)) + gamma*S)
  H <- solve(Hinv)
  # ------ -- - - -- - - - - -
  m0 <- colMeans(X0)
  m1 <- colMeans(X1)
  mdif <- m0 - m1
  msum <- m1 + m0
  at <- H %*% mdif
  mt <- t(at) %*% msum
  m <- -0.5*mt - log((1-cost_10)/cost_10)/gamma
  # ------- omega optimal calculation -------------
  traceH <- sum(diag(H))  # sum of diagonal elements in H
  deltahat <- (p/(n0+n1-2)-traceH/(n0+n1-2)) / (gamma*(1 - p/(n0+n1-2) + traceH/(n0+n1-2))) #+
  G0 <- 0.5*t(m0-m1) %*% H %*% mdif - log((1-cost_10)/cost_10)/gamma  # in paper here is gamma
  G1 <- 0.5*t(m1-m0) %*% H %*% mdif - log((1-cost_10)/cost_10)/gamma  # in paper here is gamma
  Ghat0 <- G0 - ((n0+n1-2)/n0) * deltahat
  Ghat1 <- G1 + ((n0+n1-2)/n1) * deltahat
  D <- t(at) %*% S %*% at
  # D <- t(mdif) %*% H %*% S %*% H %*% mdif  # no difference, more explicit
  Dhat <- D*(1+gamma*deltahat)^2
  omegaopt <- gamma * (Dhat*log((1-cost_10)/cost_10)/(Ghat1-Ghat0) - 0.5*(Ghat0+Ghat1))
  m <- as.numeric(m + omegaopt/gamma)
  return(structure(list(a = at,
                        m = m,
                        cost_10 = cost_10,
                        gamma = gamma,
                        Ghat0 = Ghat0,
                        Ghat1 = Ghat1,
                        Dhat = Dhat,
                        omegaopt = omegaopt), class="abcrlda"))

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
predict_abcrlda <- function(object, x){
  ## check requirements
  if(class(object) != "abcrlda")
    stop("object has to be of type abcrlda")

  if(!is.vector(x) && !is.matrix(x) && !is.data.frame(x) && !is.data.frame(x)) #  && !is.data.frame(x)
    stop("'x' has to be a vector, matrix or data.frame")

  return(as.numeric(x %*% object$a + object$m <= 0))  # object$cost

}



