# ----------------- Parameter estimation -----------------------
grid_search <- function(range_gamma, range_C_01, range_C_10){

}

cross_validation <- function(gamma, C_01, C_10){

}

equation29 <- function(gamma, C_01, C_10){

}

train <- function(){

}

# ----------------------------------------------------------------

#' Title
#'
#' @param x0
#' @param x1
#' @param gamma
#' @param cost_10
#'
#' @return
#' @export
#'
#' @examples
abcrlda <- function(x0, x1, gamma, cost_10, kappa = 1){  # cost_01 = 1 -  cost_10
  p = ncol(x0)  # number of dimensions
  n0 = nrow(x0) # number of samples in x0
  n1 = nrow(x1) # number of samples in x1
  X0 = t(as.matrix(x0))
  X1 = t(as.matrix(x1))
  S0 = cov(t(X0))
  S1 = cov(t(X1))
  S = ((n0-1)*S0+(n1-1)*S1)/(n0+n1-2)
  Hinv = (diag(dim(X1)[1])+gamma*S)
  H = solve(Hinv)
  # ------ -- - - -- - - - - -
  m1 = rowMeans(X0)
  m2 = rowMeans(X1)
  mdif = m1 - m2
  msum = m2 + m1
  at = kappa * H %*% mdif   # according to paper we should multiply by kappa?
  mt = t(at) %*% msum  # should be scalar
  m = -0.5*mt - log((1-cost_10)/cost_10)/kappa
  # ------- omega optimal calculation -------------
  traceH = sum(diag(H))  # sum of diagonal elements in H
  deltahat = (p/(n0+n1-2)-traceH/(n0+n1-2)) / (gamma*(1 - p/(n0+n1-2) + traceH/(n0+n1-2)))
  G0 = 0.5*t(m1-m2) %*% H %*% mdif - log((1-cost_10)/cost_10)/kappa  # in paper here is gamma
  G1 = 0.5*t(m2-m1) %*% H %*% mdif - log((1-cost_10)/cost_10)/kappa  # in paper here is gamma
  Ghat0 = G0 - ((n0+n1-2)/n0) * deltahat
  Ghat1 = G1 + ((n0+n1-2)/n1) * deltahat
  D = t(at) %*% S %*% at
  Dhat = D*(1+gamma*deltahat)^2
  omegaopt = gamma * (Dhat*log((1-cost_10)/cost_10)/(Ghat1-Ghat0) - 0.5*(Ghat0+Ghat1))
  m = as.numeric(m + omegaopt/gamma) # add the bias term
  return(structure(list(a=at, m=m, cost=cost_10), class="abcrlda"))

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

  if(!is.vector(x) && !is.matrix(x) && !is.data.frame(x)) #  && !is.data.frame(x)
    stop("'x' has to be a vector, matrix or data.frame")

  if(!is.data.frame(x))
    x = as.matrix(x)

  return(as.numeric(x %*% object$a + object$m <= object$cost))

}
