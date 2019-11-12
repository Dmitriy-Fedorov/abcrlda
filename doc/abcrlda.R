## ----eval=F, include=F---------------------------------------------------
#  For some application cost of making mistake
#  For example, in case of some diseases it is more important to
#  decrease number of false negatives (illness is not detected) even at expence

## ------------------------------------------------------------------------
library(abcrlda)
data(iris)

## ------------------------------------------------------------------------
train_data <- iris[which(iris[, ncol(iris)] == "versicolor" |
                         iris[, ncol(iris)] == "virginica"), 1:4]
train_label <- factor(iris[which(iris[, ncol(iris)] == "versicolor" |
                                 iris[, ncol(iris)] == "virginica"), 5])

## ------------------------------------------------------------------------
model <- abcrlda(train_data, train_label)
predict(model, train_data)

## ------------------------------------------------------------------------
risk_calculate(model, train_data, train_label)

## ------------------------------------------------------------------------
model <- abcrlda(train_data, train_label, cost = c(0.75, 0.25))


test_versicolor <- iris[which(iris[, ncol(iris)] == "versicolor"), 1:4]
test_virginica <- iris[which(iris[, ncol(iris)] == "virginica"), 1:4]

r <- predict(model, test_versicolor)
sum(r != "versicolor") / nrow(test_virginica)  # resubstitution error for first class

r <- predict(model, test_virginica)
sum(r != "virginica") / nrow(test_virginica)  # resubstitution error for second class


## ------------------------------------------------------------------------
model <- abcrlda(train_data, train_label, cost = c(0.25, 0.75))


test_versicolor <- iris[which(iris[, ncol(iris)] == "versicolor"), 1:4]
test_virginica <- iris[which(iris[, ncol(iris)] == "virginica"), 1:4]

r <- predict(model, test_versicolor)
sum(r != "versicolor") / nrow(test_virginica)  # resubstitution error for first class

r <- predict(model, test_virginica)
sum(r != "virginica") / nrow(test_virginica)  # resubstitution error for second class

## ------------------------------------------------------------------------
model <- abcrlda(train_data, train_label, gamma = 0.5, cost = 0.75)
a <- predict(model, train_data)

## ------------------------------------------------------------------------
# same params but more explicit
model <- abcrlda(train_data, train_label, gamma = 0.5, cost = c(0.75, 0.25))
b <- predict(model, train_data)

## ------------------------------------------------------------------------
# same class costs ratio
model <- abcrlda(train_data, train_label, gamma = 0.5, cost = c(3, 1))
c <- predict(model, train_data)

# all this model will give the same predictions
all(a == b & a == c & b == c)

## ----eval=F, include=T---------------------------------------------------
#  cost_range <- seq(0.1, 0.9, by = 0.2)
#  gamma_range <- c(0.1, 1, 10, 100, 1000)
#  
#  gs <- grid_search(train_data, train_label,
#                    range_gamma = gamma_range,
#                    range_cost = cost_range,
#                    method = "estimator")
#  model <- abcrlda(train_data, train_label,
#                   gamma = gs$gamma, cost = gs$cost)
#  

## ----eval=F, include=T---------------------------------------------------
#  cost_range <- matrix(1:10, ncol = 2)
#  gamma_range <- c(0.1, 1, 10, 100, 1000)
#  
#  gs <- grid_search(train_data, train_label,
#                    range_gamma = gamma_range,
#                    range_cost = cost_range,
#                    method = "cross")
#  model <- abcrlda(train_data, train_label,
#                   gamma = gs$gamma, cost = gs$cost)

