
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

  list_gamma <- numeric()
  list_C_10 <- numeric()
  list_estimates <- numeric()

  if(method == "estimator"){
    for(gamma in range_gamma){
      for(C_10 in range_C_10){
        abcrlda_model <- abcrlda(x, grouping, gamma, C_10)
        list_gamma <- c(list_gamma, gamma)
        list_C_10 <- c(list_C_10, C_10)
        list_estimates <- c(list_estimates, risk_estimate_20(abcrlda_model))
      }
    }
  }else if(method == "cross"){
    for(gamma in range_gamma){
      for(C_10 in range_C_10){
        list_gamma <- c(list_gamma, gamma)
        list_C_10 <- c(list_C_10, C_10)
        list_estimates <- c(list_estimates, cross_validation(x, grouping, gamma, C_10, k_fold))
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
    res0 <- predict_abcrlda(abcrlda_model, test0)
    res1 <- predict_abcrlda(abcrlda_model, test1)
    nerr0 <- sum(res0)
    nerr1 <- sum(!res1)
    err0 <- nerr0/length(res0)
    err1 <- nerr1/length(res1)

    error <- err0*C_10 + err1*(1 - C_10)

    # predictions <- predict_abcrlda(abcrlda_model, testData)
    # error <- sum(predictions != testLabel) / length(testLabel)
    # print(error)
    e_cross <- c(e_cross, error)
  }
  return(mean(e_cross))
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

  e0 <- error_estimate_29(object, 0)
  e1 <- error_estimate_29(object, 1)
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

  Ghat <- c(object$Ghat0, object$Ghat1)
  X <- (-1)^i * (Ghat[1+i] + object$omegaopt/object$gamma) / sqrt(object$Dhat)
  return(1 - pnorm(X))
}
