library("abcrlda")
bayes_error = pnorm(2)
N = 100
N_test = 10000
p0 = c(0, 0, 0.5)  # mean x, mean y, standart deviation
p1 = c(0, 2, 0.5)
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

err <- function(model, test, grouping){
  test0 <- test[grouping == 0,]
  test1 <- test[grouping == 1,]
  res0 = predict_abcrlda(model, test0)
  res1 = predict_abcrlda(model, test1)
  nerr0 = sum(res0)
  nerr1 = sum(!res1)
  err0 = nerr0/length(res0)
  err1 = nerr1/length(res1)
  return(structure(list(err0 = err0, err1 = err1,
                        errTotal = (nerr0 + nerr1)/length(grouping),
                        errCost = err0*model$cost_10 + err1*(1-model$cost_10))))
}

plott <- function(model, x, grouping, start = -2, finish = 2){
  p1 <- -(start * model$a[1] / model$a[2]) - model$m / model$a[2]
  p2 <- -(finish * model$a[1] / model$a[2]) - model$m / model$a[2]
  separator <- matrix(c(start, finish, p1, p2), ncol = 2)
  index <- which(grouping == 0)
  grouping[index] <- rgb(red = 1, green = 0, blue = 0, alpha = 0.8)
  grouping[-index] <- rgb(red = 0, green = 0, blue = 1, alpha = 0.8)
  plot(x, col = grouping,
       xlim=c(-2, 2), ylim = c(-1.5,3.5), pch = 16)
  lines(separator, lwd=3)
  legend('topright', c("1", "0"), col=c('blue', 'red'), pch=c(16,16))
}

model = abcrlda(train0, train1, gamma = 0.5, cost_10 = 0.5)
err(model, test, test_label)

plott(model, train, train_label)
plott(model, test, test_label)

# ----------
gs = grid_search(train0, train1,
            range_gamma = seq(0.01, 10.1, by=0.01),
            range_C_10 = seq(0.2, 0.8, by=0.02))
gs
model_s = abcrlda(train0, train1, gamma = gs$gamma[1], cost_10 = gs$C_10[1])

err(model_s, test0, test1)
plott(model_s, train0, train1)
plott(model_s, test0, test1)


##########

