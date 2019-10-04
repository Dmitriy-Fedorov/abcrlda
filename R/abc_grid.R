
#' Grid Search
#' @description Performs grid search based on cross validation or error estimation formula.
#'
#' @param x Matrix or data.frame of observations.
#' @param grouping Grouping variable. A vector of numeric values 0 and 1 is recommended. Length has to correspond to nrow(x).
#' @param range_gamma vector of gamma values to check
#' @param range_C_10 vector of cost values to check
#' @param range_C_01 vector of cost values to check
#' @param method selects method to evaluete error. "estimator" and "cross"
#' @param k_fold number of fold to use with cross-validation
#'
#' @return List of best founded parameters
#' @export
#'
#' @example inst/examples/example_grid.R

grid_search <- function(x, grouping, range_gamma, range_C_10, range_C_01=NULL,
                        method="estimator", k_fold=10){

  list_gamma <- numeric()
  list_C_10 <- numeric()
  list_estimates <- numeric()

  if (method == "estimator"){
    for (gamma in range_gamma){
      for (C_10 in range_C_10){
        abcrlda_model <- abcrlda(x, grouping, gamma, C_10)
        list_gamma <- c(list_gamma, gamma)
        list_C_10 <- c(list_C_10, C_10)
        list_estimates <- c(list_estimates, risk_estimate_20(abcrlda_model))
      }
    }
  }else if (method == "cross"){
    for (gamma in range_gamma){
      for (C_10 in range_C_10){
        list_gamma <- c(list_gamma, gamma)
        list_C_10 <- c(list_C_10, C_10)
        list_estimates <- c(list_estimates,
                            cross_validation(x, grouping, gamma, C_10, k_fold))
      }
    }
  }
  best_param_index <- which(list_estimates == min(list_estimates))
  return(structure(list(gamma = list_gamma[best_param_index],
                        C_10 = list_C_10[best_param_index],
                        e = list_estimates[best_param_index]
  )))
}

#' Cross Validation
#' @inheritParams grid_search
#' @param gamma regularization parameter
#' @param C_10 parameter that controls prioretization of classes.
#' It's value should be between 0 and 1 (0 < cost_10 < 1)
#' Values bigger than 0.5 prioretizes correct classification of 0 class while values less than 0.5 prioretizes 1 class
#' @param kfolds Number of for cross validation algorithm
#' @return Returns average error of cross validation
#' @export
#' @family abcrlda binary classifier
#'
#' @example inst/examples/example_cross.R
cross_validation <- function(x, grouping, gamma=1, C_10=0.5, kfolds=10){
  shufled_index <- sample(nrow(x))
  x <- x[shufled_index, ]
  grouping <- grouping[shufled_index]
  folds <- cut(seq(1, nrow(x)), breaks = kfolds, labels = FALSE)
  e_cross <- numeric()

  # print(paste("____________", gamma, C_10))
  for (i in 1:kfolds){
    # Segement your data by fold using the which() function
    test_indexes <- which(folds == i, arr.ind = TRUE)
    test_data <- x[test_indexes, ]
    test_label <- grouping[test_indexes]
    train_data <- x[-test_indexes, ]
    train_label <- grouping[-test_indexes]

    abcrlda_model <- abcrlda(train_data, train_label, gamma, C_10)
    test0 <- test_data[test_label == 0, ]
    test1 <- test_data[test_label == 1, ]
    res0 <- stats::predict(abcrlda_model, test0)
    res1 <- stats::predict(abcrlda_model, test1)
    nerr0 <- sum(res0$raw)
    nerr1 <- sum(!res1$raw)
    err0 <- nerr0 / length(res0$raw)
    err1 <- nerr1 / length(res1$raw)

    error <- err0 * C_10 + err1 * (1 - C_10)
    # print(paste(i, error, mean(e_cross)))
    e_cross <- c(e_cross, error)
  }
  return(mean(e_cross))
}



#' Risk Estimator
#' @description Calculates weighted error based on normalized cost values
#' @inheritParams predict.abcrlda
#'
#' @return Weighted error based on "abcrlda" object
#' @export
#' @example inst/examples/example_risk.R
risk_estimate_20 <- function(object){
  ## check requirements
  if (class(object) != "abcrlda")
    stop("object has to be of type abcrlda")

  error_estimate_29 <- function(object, i){
    Ghat <- c(object$Ghat0, object$Ghat1)
    X <- (-1) ^ i * (Ghat[1 + i] + object$omegaopt / object$gamma) /
         sqrt(object$Dhat)
    return(1 - stats::pnorm(X))
  }

  e0 <- error_estimate_29(object, 0)
  e1 <- error_estimate_29(object, 1)
  return(object$cost_10 * e0 + (1 - object$cost_10) * e1)
}
