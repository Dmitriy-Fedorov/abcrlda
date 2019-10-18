
#' Grid Search
#' @description Performs grid search based on double asymptotic risk estimation or cross validation.
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
  risk_estimates <- numeric()

  range_cost <- as.matrix(range_cost)

  if (method == "estimator"){
    for (gamma in range_gamma){
      for (i in 1:nrow(range_cost)){
        cost <- range_cost[i, ]
        abcrlda_model <- abcrlda(x, y, gamma, cost)
        list_gamma <- c(list_gamma, gamma)
        list_cost <- rbind(list_cost, cost)
        risk_estimates <- c(risk_estimates, da_risk_estimator(abcrlda_model))
      }
    }
  }else if (method == "cross"){
    for (gamma in range_gamma){
      for (i in 1:nrow(range_cost)){
        cost <- range_cost[i, ]
        list_gamma <- c(list_gamma, gamma)
        list_cost <- rbind(list_cost, cost)
        risk_estimates <- c(risk_estimates,
                            cross_validation(x, y, gamma, cost, nfolds)$ecr)
      }
    }
  }
  best_param_index <- risk_estimates == min(risk_estimates)
  return(structure(list(gamma = list_gamma[best_param_index],
                        cost = list_cost[best_param_index, ],
                        risk = risk_estimates[best_param_index])))
}

#' Cross Validation for separate sampling adjasted for cost
#' @description
#' @inheritParams grid_search
#' @inheritParams abcrlda
# @param gamma regularization parameter
# @param cost parameter that controls prioretization of classes.
# It's value should be between 0 and 1 (0 < cost_10 < 1)
# Values bigger than 0.5 prioretizes correct classification of 0 class while values less than 0.5 prioretizes 1 class
# @param nfolds Number of for cross validation algorithm
#' @return Returns average risk after cross validation
#' @export
#' @family functions in the package
#'
#' @example inst/examples/example_cross.R
cross_validation <- function(x, y, gamma=1, cost=c(0.5, 0.5), nfolds=10){

  x <- as.matrix(x)
  shufled_index <- sample(nrow(x))
  x <- x[shufled_index, ]
  y <- y[shufled_index]

  if (!is.factor(y))
    y <- as.factor(y)
  lev <- levels(y)
  k <- nlevels(y)

  x0 <- x[y == lev[1], , drop = FALSE]
  x1 <- x[y == lev[2], , drop = FALSE]
  y0 <- y[y == lev[1]]
  y1 <- y[y == lev[2]]

  fold0 <- cut(seq(1, nrow(x0)), breaks = nfolds, labels = FALSE)
  fold1 <- cut(seq(1, nrow(x1)), breaks = nfolds, labels = FALSE)

  e0 <- numeric()
  e1 <- numeric()

  for (i in 1:nfolds){
    for (j in 1:nfolds){
      test_index0 <- which(fold0 == i, arr.ind = TRUE)
      test_index1 <- which(fold1 == j, arr.ind = TRUE)

      test_data0 <- x0[test_index0, , drop = FALSE]
      test_data1 <- x1[test_index1, , drop = FALSE]
      test_label0 <- y[test_index0]
      test_label1 <- y[test_index1]

      train_data <- rbind(x0[-test_index0, , drop = FALSE],
                          x1[-test_index1, , drop = FALSE])
      train_label <- factor(c(y0[-test_index0], y1[-test_index1]),
                            levels=1:k, labels=lev)

      model <- abcrlda(train_data, train_label, gamma, cost)

      res0 <- as.numeric(test_data0 %*% model$a + model$m <= 0)
      res1 <- as.numeric(test_data1 %*% model$a + model$m <= 0)

      e0 <- c(e0, sum(res0) / length(res0))
      e1 <- c(e1, sum(!res1) / length(res1))

    }
  }

  return(list(ecross = mean(e0) * cost[1] + mean(e1) * cost[2],
              e0 = mean(e0),
              e1 = mean(e1)))
}

#' Double Asymptotic Risk Estimator
#' @description Generalized consistent estimator of risk
#' @inheritParams predict.abcrlda
#'
#' @return Calculates risk based on estimated class error rates and misclassification costs
#'   \deqn{\Re = \varepsilon_0 * cost_{10} + \varepsilon_1 * cost_{01}}{R = e0 * cost_10 + e1 * cost_01)}
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

