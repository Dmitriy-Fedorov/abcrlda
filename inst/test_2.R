library(abcrlda)
library(caret)
library(mvtnorm)
library(Matrix)
library(matrixcalc)
library(corpcor)
# bayes_error = pnorm(2)
# N = 1000
# N_test = 10000
# p0 = c(0, 0, 0.5)  # mean x, mean y, standart deviation
# p1 = c(1, 1, 0.5)
# Sys.which("pdflatex")
# devtools::build_manual(path = "inst")
# lintr::lint_package()

# glmnet

generate_train_test <- function(N=100, N_test=10000,
                                m0=c(0, 0), m1=c(0, 1), s0=diag(0.5,nrow = 2), s1=diag(0.5,nrow = 2)){
  train0 <- rmvnorm(N, mean = m0, sigma = s0)
  train1 <- rmvnorm(N, mean = m1, sigma = s1)
  test0 <- rmvnorm(N_test, mean = m0, sigma = s0)
  test1 <- rmvnorm(N_test, mean = m1, sigma = s1)
  train <- rbind(train0, train1)
  train_label <- c(replicate(N, 0), replicate(N, 1))
  test <- rbind(test0, test1)
  test_label <- c(replicate(N_test, 0), replicate(N_test, 1))
  return(structure(list(train=train,
                        train_label=train_label,
                        test=test,
                        test_label=test_label)))
}


err <- function(model, test, grouping){
  test0 <- test[grouping == model$lev[1],]
  test1 <- test[grouping == model$lev[2],]
  res0 = as.numeric(predict(model, test0)) - 1
  res1 = as.numeric(predict(model, test1)) - 1
  nerr0 = sum(res0)
  nerr1 = sum(!res1)
  err0 = nerr0 / length(res0)
  err1 = nerr1 / length(res1)
  return(t(data.frame(actual_err0 = err0, actual_err1 = err1,
                      actual_errTotal = (nerr0 + nerr1) / length(grouping),
                      actual_normrisk = err0*model$ncost[1] + err1*model$ncost[2],
                      actual_risk = err0*model$cost[1] + err1*model$cost[2])))
}

plott <- function(model, x, grouping, start = -2, finish = 2){
  p1 <- -(start * model$a[1] / model$a[2]) - model$m / model$a[2]
  p2 <- -(finish * model$a[1] / model$a[2]) - model$m / model$a[2]
  separator <- matrix(c(start, finish, p1, p2), ncol = 2)
  index <- which(grouping == 0)
  grouping[index] <- rgb(red = 1, green = 0, blue = 0, alpha = 0.5)
  grouping[-index] <- rgb(red = 0, green = 0, blue = 1, alpha = 0.5)
  plot(x, col = grouping,
       xlim=c(-2, 2), ylim = c(-1.5,3.5), pch = 16,cex = .8)
  lines(separator, lwd=3)
  legend('topright', c("1", "0"), col=c('blue', 'red'), pch=c(16,16))
}


# -------------------------------- 2D test manual --------------------------------
gen <- generate_train_test(N=1000, m0=c(0, 0), m1=c(0, 2))
costtt <- c(0.7, 0.3)
model <- abcrlda::abcrlda(gen$train, gen$train_label, gamma = 1, cost = costtt, bias_correction = 1)
e = err(model, gen$test, gen$test_label)
rbind(e,
  t(data.frame(
    estimated_risk = abcrlda::da_risk_estimator(model),
    cross = abcrlda::cross_validation(gen$train, gen$train_label, gamma = 1, cost=costtt, nfold=10))
  ))

abcrlda::cross_validation(gen$train, gen$train_label, gamma = 1, cost=costtt, nfold=10)


confusionMatrix(reference = as.factor(gen$test_label),
                data = (predict(model, gen$test)))
asd = predict(model, gen$test)
asd
# asd2 = predict(model, gen$test, type="class")
# all(asd == as.numeric(asd2))
# class(predict(model, gen$test, type="class"))
plott(model, gen$train, gen$train_label)
plott(model, gen$test, gen$test_label)
# -------------------------------- 2D test grid_search --------------------------------
gen <- generate_train_test(m0=c(0, 0), m1=c(0, 1))
crange <- seq(0.4, 0.7, by=0.05)
grange <- (c(5,10) %o% 10^(-1:3))
gs_e <- grid_search(gen$train, gen$train_label,
                  range_gamma = grange,
                  range_cost = crange,
                  method = "estimator")

gs_c <- grid_search(gen$train, gen$train_label,
                   range_gamma = grange,
                   range_cost = crange,
                   method = "cross")


model_est = abcrlda(gen$train, gen$train_label, gamma = gs_e$gamma, cost = gs_e$cost)
model_crs = abcrlda(gen$train, gen$train_label, gamma = gs_c$gamma, cost = gs_c$cost)

er_e = err(model_est, gen$test, gen$test_label)
er_c = err(model_crs, gen$test, gen$test_label)

cbind(a = rbind(t(as.data.frame(gs_e)), er_e),
      b = rbind(t(as.data.frame(gs_c)), er_c))
# rbind(t(as.data.frame(gs_e)), er_e)
# rbind(t(as.data.frame(gs_c)), er_c)

confusionMatrix(reference = as.factor(gen$test_label),
                data = as.factor((predict(model_est, gen$test))),
                mode='everything')
confusionMatrix(reference = as.factor(gen$test_label),
                data = as.factor((predict(model_crs, gen$test))),
                mode='everything')

plott(model_est, gen$train, gen$train_label)
plott(model_crs, gen$train, gen$train_label)


# -------------------------------- multivariate test --------------------------------

set.seed(9)

gen <- generate_train_test(N=100,N_test=1000,
                           m0 = c(0,0,0,0,0,0,0,0,0,0),
                           m1 = c(1,1,1,1,1,1,1,1,1,1),
                           s0 = make.positive.definite(as.matrix(abs(
                             forceSymmetric(Matrix(rnorm(100),10)))), tol=1e-3),
                           s1 = make.positive.definite(as.matrix(abs(
                             forceSymmetric(Matrix(rnorm(100),10)))), tol=1e-3))

crange <- seq(0.3, 0.7, by=0.05)
grange <- (c(5,10) %o% 10^(-1:3))
gs_e <- grid_search(gen$train, gen$train_label,
                    range_gamma = grange,
                    range_cost = crange,
                    method = "estimator")

gs_c <- grid_search(gen$train, gen$train_label,
                    range_gamma = grange,
                    range_cost = crange,
                    method = "cross")


model_est = abcrlda(gen$train, gen$train_label, gamma = gs_e$gamma, cost = gs_e$cost)
model_crs = abcrlda(gen$train, gen$train_label, gamma = gs_c$gamma, cost = gs_c$cost)

er_e = err(model_est, gen$test, gen$test_label)
er_c = err(model_crs, gen$test, gen$test_label)

cbind(a = rbind(t(as.data.frame(gs_e)), er_e),
      b = rbind(t(as.data.frame(gs_c)), er_c))
# rbind(t(as.data.frame(gs_e)), er_e)
# rbind(t(as.data.frame(gs_c)), er_c)

confusionMatrix(reference = as.factor(gen$test_label),
                data = as.factor((predict(model_est, gen$test))))
confusionMatrix(reference = as.factor(gen$test_label),
                data = as.factor((predict(model_crs, gen$test))))


# -------------------------------- iris dataset --------------------------------

data(iris)

traindata = iris[ which(iris[,ncol(iris)]=='virginica' |
                          iris[,ncol(iris)]=="versicolor"), 1:4]
trainlabel = factor(iris[ which(iris[,ncol(iris)]=='virginica' |
                                  iris[,ncol(iris)]=="versicolor"), 5])
bias = 1

rr <- abcrlda(traindata, trainlabel, gamma=1, cost = 0.75, bias_correction = bias)
# stats::predict(rr, traindata)
e = err(rr, traindata, trainlabel)

rbind(e,
      t(data.frame(
        bias_correction = bias,
        restimate = abcrlda::da_risk_estimator(rr),
        cross = abcrlda::cross_validation(traindata, trainlabel, cost=c(0.75,0.25)))
      ))
# abcrlda::cross_validation(traindata, trainlabel, cost=c(0.75, 0.25), nfolds = 3)
# asd = predict(rr, traindata, type="raw")
# asd2 = predict(rr, traindata, type="class")
# all(asd == as.numeric(asd2))

# -------------------------------- s5 dataset --------------------------------
train = read.csv(file="A:/Dropbox/9th/R/Code/s5/SigPar_Train.csv", header=FALSE, sep=",")
test = read.csv(file="A:/Dropbox/9th/R/Code/s5/SigPar_Test.csv", header=FALSE, sep=",")

train_data <- train[, 1:ncol(train) - 1]
train_label <- train[, ncol(train)]
test_data <- test[, 1:ncol(train) - 1]
test_label <- test[, ncol(train)]

costtt <- c(0.7, 0.3)
bias = FALSE
model <- abcrlda::abcrlda(train_data, train_label, gamma = 1, cost = costtt, bias_correction = bias)
e = err(model, test_data, test_label)
rbind(e,
      t(data.frame(
        bias_correction = bias,
        estimated_risk = abcrlda::da_risk_estimator(model),
        cross = abcrlda::cross_validation(train_data, train_label, gamma = 1, cost=costtt, nfold=10))
      ))

abcrlda::cross_validation(train_data, train_label, gamma = 1, cost=costtt, nfold=10)


# -------------------------------- CRAN ---------------------------------
usethis::use_release_issue()
# -------------------------------- other --------------------------------
data(iris)
train_data <- iris[which(iris[, ncol(iris)] == "virginica" |
                           iris[, ncol(iris)] == "versicolor"), 1:4]
train_label <- factor(iris[which(iris[, ncol(iris)] == "virginica" |
                                   iris[, ncol(iris)] == "versicolor"), 5])
model <- abcrlda(train_data, train_label, gamma = 0.5, cost = 0.75)
a <- predict(model, train_data)
# same params but more explicit
model <- abcrlda(train_data, train_label, gamma = 0.5, cost = c(0.75, 0.25))
b <- predict(model, train_data)
# same class costs ratio
model <- abcrlda(train_data, train_label, gamma = 0.5, cost = c(3, 1))
c <- predict(model, train_data)


cost_range <- matrix(1:10, ncol = 2)
gamma_range <- c(0.1, 1, 10, 100, 1000)

gs <- grid_search(train_data, train_label,
                  range_gamma = gamma_range,
                  range_cost = cost_range,
                  method = "cross")
model <- abcrlda(train_data, train_label,
                 gamma = gs$gamma, cost = gs$cost)
predict(model, train_data)


data(iris)
x <- iris[which(iris[, ncol(iris)] == "virginica" |
                           iris[, ncol(iris)] == "versicolor"), 1:4]
y <- factor(iris[which(iris[, ncol(iris)] == "virginica" |
                                   iris[, ncol(iris)] == "versicolor"), 5])

cross_validation(x, y, gamma = 5, nfolds = 100)

model <- abcrlda(x, y, gamma = 5)
da_risk_estimator(model)
model$a


folds = cut(seq(1, nrow(x)), breaks = 100, labels = FALSE)


test_indexes <- which(folds == 50, arr.ind = TRUE)
test_indexes
x <- as.matrix(x)
x
test_data <- x[test_indexes, , drop = FALSE]
test_data
dim(test_data)
typeof(test_data)
predict(model, test_data)


test_label <- y[test_indexes]
train_data <- x[-test_indexes, ]
train_label <- y[-test_indexes]

abcrlda_model <- abcrlda(train_data, train_label)


a = c(1,0,0,0,1,1,0)
b <- as.factor(a)
b[b == 0]
b[b == levels(b)[2]]
ad = c(b[b == levels(b)[1]], b[b == levels(b)[2]])
factor(ad, levels=1:nlevels(b), labels=levels(b))
