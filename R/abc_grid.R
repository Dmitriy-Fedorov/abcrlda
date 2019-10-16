
#' Grid Search
#' @description Performs grid search based on cross validation or error estimation formula.
#'
# @param x Matrix or data.frame of observations.
# @param y Grouping variable. A vector of numeric values 0 and 1 is recommended.
# Length has to correspond to nrow(x).
#' @param range_gamma vector of gamma values to check
#' @param range_cost nobs x 1 vector (values should be between 0 and 1) or
#'   nobs x 2 matrix (each row is cost pair value c(\eqn{C_{10}}{C_10}, \eqn{C_{01}}{C_01}))
#'   of cost values to check
#' @param method selects method to evaluete error. "estimator" and "cross"
#' @param nfolds number of fold to use with cross-validation. Default is 10.
#' @inheritParams abcrlda
#' @return List of best founded parameters
#' @export
#' @family functions in the package
#'
#' @example inst/examples/example_grid.R

grid_search <- function(x, y, range_gamma, range_cost,
                        method="estimator", nfolds=10){

  list_gamma <- numeric()
  list_cost <- numeric()
  list_estimates <- numeric()

  range_cost <- as.matrix(range_cost)

  if (method == "estimator"){
    for (gamma in range_gamma){
      for (i in 1:nrow(range_cost)){
        cost <- range_cost[i, ]
        abcrlda_model <- abcrlda(x, y, gamma, cost)
        list_gamma <- c(list_gamma, gamma)
        list_cost <- c(list_cost, cost)
        list_estimates <- c(list_estimates, da_risk_estimator(abcrlda_model))
      }
    }
  }else if (method == "cross"){
    for (gamma in range_gamma){
      for (i in 1:nrow(range_cost)){
        cost <- range_cost[i, ]
        list_gamma <- c(list_gamma, gamma)
        list_cost <- c(list_cost, cost)
        list_estimates <- c(list_estimates,
                            cross_validation(x, y, gamma, cost, nfolds))
      }
    }
  }
  best_param_index <- which(list_estimates == min(list_estimates))
  return(structure(list(gamma = list_gamma[best_param_index],
                        cost = list_cost[best_param_index],
                        e = list_estimates[best_param_index])))
}

#' Cross Validation
#' @inheritParams grid_search
#' @inheritParams abcrlda
# @param gamma regularization parameter
# @param cost parameter that controls prioretization of classes.
# It's value should be between 0 and 1 (0 < cost_10 < 1)
# Values bigger than 0.5 prioretizes correct classification of 0 class while values less than 0.5 prioretizes 1 class
# @param nfolds Number of for cross validation algorithm
#' @return Returns average error of cross validation
#' @export
#' @family functions in the package
#'
#' @example inst/examples/example_cross.R
cross_validation <- function(x, y, gamma=1, cost=c(0.5, 0.5), nfolds=10){

  shufled_index <- sample(nrow(x))
  x <- x[shufled_index, ]
  y <- y[shufled_index]
  folds <- cut(seq(1, nrow(x)), breaks = nfolds, labels = FALSE)
  e_cross <- numeric()

  if (!is.factor(y))
    y <- as.factor(y)
  lev <- levels(y)
  k <- nlevels(y)

  if (k != 2)
    stop("number of groups != 2, this is binary classifier")

  if (length(cost) == 1){
    if (cost >= 1 | cost <= 0)
      stop("While providing single valued vector
           cost should be between 0 and 1 (not including)")
    cost <- c(cost, 1 - cost)
  }
  if (length(cost) != 2)
    stop("cost vector should be of length 1 or 2, this is binary classifier")

  for (i in 1:nfolds){
    # Segement your data by fold using the which() function
    test_indexes <- which(folds == i, arr.ind = TRUE)
    test_data <- x[test_indexes, ]
    test_label <- y[test_indexes]
    train_data <- x[-test_indexes, ]
    train_label <- y[-test_indexes]

    abcrlda_model <- abcrlda(train_data, train_label, gamma, cost)
    # print(test_label)
    test0 <- test_data[test_label == lev[1], ]
    test1 <- test_data[test_label == lev[2], ]
    # print(test0)
    # print(test1)
    res0 <- stats::predict(abcrlda_model, test0, type = "raw") - 1
    res1 <- stats::predict(abcrlda_model, test1, type = "raw") - 1
    # print(res0)
    # print(res1)
    nerr0 <- sum(res0)
    nerr1 <- sum(!res1)
    # print(c(nerr0, nerr1))
    err0 <- nerr0 / length(res0)
    err1 <- nerr1 / length(res1)

    error <- err0 * cost[1] + err1 * cost[2]
    e_cross <- c(e_cross, error)
  }
  # print(e_cross)
  return(mean(e_cross))
}


# da_risk_estimator double asymptotic
#' Double Asymptotic Risk Estimator
#' @description Calculates weighted error based on normalized cost values
#' @inheritParams predict.abcrlda
#'
#' @return Weighted error based on "abcrlda" object
#' @export
#' @family functions in the package
#' @example inst/examples/example_risk.R
#' @inheritSection abcrlda Reference
da_risk_estimator <- function(object){
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
  return(object$cost[1] * e0 + object$cost[2] * e1)
}
