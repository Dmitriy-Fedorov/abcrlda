library(abcrlda)
library(caret)
library(mvtnorm)
# bayes_error = pnorm(2)
# N = 1000
# N_test = 10000
# p0 = c(0, 0, 0.5)  # mean x, mean y, standart deviation
# p1 = c(1, 1, 0.5)

generate_train_test <- function(N=100, N_test=10000,
                                p0=c(0, 0, 0.5), p1=c(0, 1, 0.5)){

  train0 <- matrix(c(rnorm(N, mean = p0[1], sd = p0[3]),
                     rnorm(N, mean = p0[2], sd = p0[3])), ncol=2)
  train1 <- matrix(c(rnorm(N, mean = p1[1], sd = p1[3]),
                     rnorm(N, mean = p1[2], sd = p1[3])), ncol=2)

  test0 <- matrix(c(rnorm(N_test, mean = p0[1], sd = p0[3]),
                    rnorm(N_test, mean = p0[2], sd = p0[3])), ncol=2)
  test1 <- matrix(c(rnorm(N_test, mean = p1[1], sd = p1[3]),
                    rnorm(N_test, mean = p1[2], sd = p1[3])), ncol=2)
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
  test0 <- test[grouping == 0,]
  test1 <- test[grouping == 1,]
  res0 = predict(model, test0)
  res1 = predict(model, test1)
  nerr0 = sum(res0$raw)
  nerr1 = sum(!res1$raw)
  err0 = nerr0/length(res0$raw)
  err1 = nerr1/length(res1$raw)
  return(structure(list(err0 = err0, err1 = err1,
                        errTotal = (nerr0 + nerr1)/length(grouping),
                        errCost = err0*model$cost_10 + err1*(1-model$cost_10))))
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

gen <- generate_train_test(p1=c(1, 1, 0.5))

# -----------------------------------------------------------------------------------------
crange <- seq(0.05, 0.95, by=0.05)
grange <- (c(5,10) %o% 10^(-1:3))
gs0 <- grid_search(gen$train, gen$train_label,
                  range_gamma = grange,
                  range_C_10 = crange,
                  method = "estimator")

gs1 <- grid_search(gen$train, gen$train_label,
                   range_gamma = grange,
                   range_C_10 = crange,
                   method = "cross")


model_s0 = abcrlda(gen$train, gen$train_label, gamma = gs0$gamma[1], cost_10 = gs0$C_10[1])
model_s1 = abcrlda(gen$train, gen$train_label, gamma = gs1$gamma[1], cost_10 = gs1$C_10[1])

e0 = err(model_s0, gen$test, gen$test_label)
e1 = err(model_s1, gen$test, gen$test_label)


gs0
e0
gs1
e1

confusionMatrix(reference = as.factor(gen$test_label),
                data = as.factor((predict(model_s0, gen$test))$raw),
                mode='everything')
confusionMatrix(reference = as.factor(gen$test_label),
                data = as.factor((predict(model_s1, gen$test))$raw),
                mode='everything')

plott(model_s0, gen$train, gen$train_label)
plott(model_s1, gen$train, gen$train_label)

##########

gen <- generate_train_test(N=1000, p1=c(0, 2, 0.5))
model <- abcrlda::abcrlda(gen$train, gen$train_label, gamma = 1, cost_10 = 5, cost_01 = 50)
err(model, gen$test, gen$test_label)
# confusionMatrix(reference = as.factor(gen$test_label),
#                 data = as.factor((predict(model, gen$test))$raw))
predict(model, gen$test)

plott(model, gen$train, gen$train_label)
plott(model, gen$test, gen$test_label)


##################################################################

data(iris)

traindata = iris[ which(iris[,ncol(iris)]=='virginica' |
                              iris[,ncol(iris)]=="versicolor"), 1:4]
trainlabel = factor(iris[ which(iris[,ncol(iris)]=='virginica' |
                                iris[,ncol(iris)]=="versicolor"), 5])

rr <- abcrlda(traindata, trainlabel, gamma = 0.5, cost_10 = 0.75)
predict(rr, traindata)
err(rr, traindata, as.numeric(trainlabel)-1)

# ------------ multivariate test -----------------------------------
library(Matrix)
library(matrixcalc)
library(corpcor)

generate_multivariate_train_test <- function(N=100, N_test=10000,
                                m0=c(0, 0, 0.5), m1=c(0, 1, 0.5), s0=NULL, s1=NULL){

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
set.seed(9)
s0 <- as.matrix(abs(forceSymmetric(Matrix(rnorm(100),10))))
s0 <- make.positive.definite(s0, tol=1e-3)
s1 <- as.matrix(abs(forceSymmetric(Matrix(rnorm(100),10))))
s1 <- make.positive.definite(s1, tol=1e-3)
m0 <- c(0,0,0,0,0,0,0,0,0,0)
m1 <- c(1,1,1,1,1,1,1,1,1,1)

gen <- generate_multivariate_train_test(N=10000,N_test=1000000,m0=m0,m1=m1,s0=s0,s1=s1)

crange <- seq(0.05, 0.95, by=0.05)
grange <- (c(5,10) %o% 10^(-1:3))
gs0 <- grid_search(gen$train, gen$train_label,
                   range_gamma = grange,
                   range_C_10 = crange,
                   method = "estimator")

gs1 <- grid_search(gen$train, gen$train_label,
                   range_gamma = grange,
                   range_C_10 = crange,
                   method = "cross")

model_s0 = abcrlda(gen$train, gen$train_label, gamma = gs0$gamma[1], cost_10 = gs0$C_10[1])
model_s1 = abcrlda(gen$train, gen$train_label, gamma = gs1$gamma[1], cost_10 = gs1$C_10[1])

e0 = err(model_s0, gen$test, gen$test_label)
e1 = err(model_s1, gen$test, gen$test_label)

gs0
e0
gs1
e1
confusionMatrix(reference = as.factor(gen$test_label),
                data = as.factor((predict(model_s0, gen$test))$raw))
confusionMatrix(reference = as.factor(gen$test_label),
                data = as.factor((predict(model_s1, gen$test))$raw))

