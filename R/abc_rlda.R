
#' Title
#'
#' @param x
#' @param grouping
#' @param range_gamma
#' @param range_C_10
#' @param method
#' @param k_fold
#'
#' @return
#' @export
#'
#' @examples
grid_search <- function(x, grouping, range_gamma, range_C_10, method="estimator",
                        k_fold=10){

  list_gamma = numeric()
  list_C_10 = numeric()
  list_estimates = numeric()

  if(method == "estimator"){
    for(gamma in range_gamma){
      for(C_10 in range_C_10){
        abcrlda_model <- abcrlda(x, grouping, gamma, C_10)
        list_gamma = c(list_gamma, gamma)
        list_C_10 = c(list_C_10, C_10)
        list_estimates = c(list_estimates, risk_estimate_20(abcrlda_model))
      }
    }
  }else if(method == "cross"){
    for(gamma in range_gamma){
      for(C_10 in range_C_10){
        list_gamma = c(list_gamma, gamma)
        list_C_10 = c(list_C_10, C_10)
        list_estimates = c(list_estimates, cross_validation(x, grouping, gamma, C_10, k_fold))
      }
    }
  }
  best_param_index <- which(list_estimates == min(list_estimates))
  return(structure(list(gamma = list_gamma[best_param_index],
                        C_10 = list_C_10[best_param_index],
                        e = list_estimates[best_param_index]
                        )))
}

#' Title
#'
#' @param x
#' @param grouping
#' @param gamma
#' @param C_10
#' @param kfolds
#'
#' @return
#' @export
#'
#' @examples
cross_validation <- function(x, grouping, gamma, C_10, kfolds=10){
  shufled_index <- sample(nrow(x))
  x <- x[shufled_index,]
  grouping <- grouping[shufled_index]
  folds <- cut(seq(1,nrow(x)),breaks=kfolds,labels=FALSE)
  e_cross <- numeric()

  # print(paste("____________", gamma, C_10))
  for(i in 1:kfolds){
    # print(paste(i, "fold"))
    #Segement your data by fold using the which() function
    testIndexes <- which(folds==i, arr.ind=TRUE)
    testData <- x[testIndexes, ]
    testLabel <- grouping[testIndexes]
    trainData <- x[-testIndexes, ]
    trainLabel <- grouping[-testIndexes]

    abcrlda_model <- abcrlda(trainData, trainLabel, gamma, C_10)
    test0 <- testData[testLabel == 0,]
    test1 <- testData[testLabel == 1,]
    res0 = predict_abcrlda(abcrlda_model, test0)
    res1 = predict_abcrlda(abcrlda_model, test1)
    nerr0 = sum(res0)
    nerr1 = sum(!res1)
    err0 = nerr0/length(res0)
    err1 = nerr1/length(res1)

    error <- err0*C_10 + err1*(1 - C_10)

    # predictions <- predict_abcrlda(abcrlda_model, testData)
    # error <- sum(predictions != testLabel) / length(testLabel)
    # print(error)
    e_cross <- c(e_cross, error)
  }
  return(mean(e_cross))
}

# ----------------------------------------------------------------

#' Title
#'
#' @param x
#' @param grouping
#' @param gamma
#' @param cost_10
#' @param kappa
#'
#' @return
#' @export
#'
#' @examples
abcrlda <- function(x, grouping, gamma, cost_10 = 0.5, kappa = 1){  # cost_01 = 1 -  cost_10

  x0 <- x[grouping == 0,]
  x1 <- x[grouping == 1,]
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
  # D = t(mdif) %*% H %*% S %*% H %*% mdif  # no difference, more explicit
  Dhat = D*(1+gamma*deltahat)^2
  omegaopt = gamma * (Dhat*log((1-cost_10)/cost_10)/(Ghat1-Ghat0) - 0.5*(Ghat0+Ghat1))
  m = as.numeric(m + omegaopt/gamma) # add the bias term, ? should it be kappa
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

  if(!is.vector(x) && !is.matrix(x) && !is.data.frame(x)) #  && !is.data.frame(x)
    stop("'x' has to be a vector, matrix or data.frame")

  if(!is.data.frame(x))
    x = as.matrix(x)

  return(as.numeric(x %*% object$a + object$m <= 0))  # object$cost

}

#' Title
#'
#' @param object
#'
#' @return
#' @export
#'
#' @examples
risk_estimate_20 <- function(object){
  ## check requirements
  if(class(object) != "abcrlda")
    stop("object has to be of type abcrlda")

  e0 = error_estimate_29(object, 0)
  e1 = error_estimate_29(object, 1)
  return(object$cost_10 * e0 + (1 - object$cost_10) * e1)
 }


#' Title
#'
#' @param object
#' @param i
#'
#' @return
#' @export
#'
#' @examples
error_estimate_29 <- function(object, i){
  ## check requirements
  if(class(object) != "abcrlda")
    stop("object has to be of type abcrlda")

  Ghat = c(object$Ghat0, object$Ghat1)
  X = (-1)^i * (Ghat[1+i] + object$omegaopt/object$gamma) / sqrt(object$Dhat)
  return(1 - pnorm(X))
}


