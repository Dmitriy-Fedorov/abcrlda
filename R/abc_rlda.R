#' Title
#'
#' @param X1 asd
#' @param X2 asd
#' @param H asd
#' @param S asd
#' @param c asd
#' @param gamma asd
#' @param kappa asd
#'
#' @return asd
#' @export
#'
#' @examples
abc_rlda <- function(X1,X2,H,S,c,gamma,kappa) {
	m1 = rowMeans(X1)
	m2 = rowMeans(X2)
	md = m1 - m2
	ms = m2 + m1
	at = H%*%md
	mt = t(at)%*%ms
	# Mat = t(at)%*%sig.true%*%at
	m = -0.5*mt - log((1-c)/c)/kappa
	p = dim(X1)[1]
	n1 = dim(X1)[2]
	n2 = dim(X2)[2]
	trH = sum(diag(H))
	s = 1-p/(n1+n2-2)+trH/(n1+n2-2)
	deltah = (p/(n1+n2-2)-trH/(n1+n2-2))/(gamma*s)
	G1 = 0.5*t(m1-m2)%*%H%*%md - ((n1+n2-2)/n1) * deltah - log((1-c)/c)/kappa
	G2 = 0.5*t(m2-m1)%*%H%*%md + ((n1+n2-2)/n2) * deltah - log((1-c)/c)/kappa
	D = t(at)%*%S%*%at
	D = D*(1+gamma*deltah)^2
	omegaopt = gamma*(D*log((1-c)/c)/(G2-G1)-0.5*(G1+G2))
  m = m + omegaopt/gamma #add the bias term
	av = c(m,at)
	return(av)
}


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

new_abcrlda <- function(x0, x1, gamma, cost_01, cost_10){

  X0 = t(as.matrix(x0))
  X1 = t(as.matrix(x1))
  S0 = cov(t(X0))
  S1 = cov(t(X1))

  n0 = nrow(x0)
  n1 = nrow(x1)

  S = ((n0-1)*S0+(n1-1)*S1)/(n0+n1-2)
  Hinv = (diag(dim(X1)[1])+gamma*S)
  H = solve(Hinv)
}

# ----------------------------------------------------------------
# predicts group membership for new data x and object from abcrlda
predict.abcrlda <- function(object, x, ...){
  ## check requirements
  if(class(object) != "abcrlda")
    stop("object has to be of type abcrlda")

  ## get dimension
  p = ncol(object$means)

  if(!is.vector(x) && !is.matrix(x) && !is.data.frame(x))
    stop("'x' has to be a vector, matrix or data.frame")

  n <- 0
  if(is.vector(x)){
    if(length(x) != p)
      stop("'x' has to be of length ", p)

    y <- NULL
    for(i in 1:length(x))
      y <- cbind(y, x[i])
    x <- y
    n <- 1
  }else{
    if(ncol(x) != p)
      stop("'x' has to be of dimension ", p)

    x <- as.matrix(x)
    n <- nrow(x)
  }

  probs <- t(-2*object$means %*% object$covi %*% t(x) + diag(object$means %*% object$covi %*% t(object$means)) - 2*log(object$prior[[1]]))
  cl = object$lev[apply(probs, 1, which.min)]

  ret = list(class=cl, posterior=probs)
  ret
}

#' Title
#'
#' @param X
#' @param a
#' @param m
#' @param c
#'
#' @return
#' @export
#'
#' @examples
test_lda<-function(X,a,m,c) {
  l = X %*% a
  if ((l + m) > c) {
    d = 1
  }else{
    d = 2
  }
  return(d)
}
