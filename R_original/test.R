
data(iris)
x <- iris[,1:4]
rr <- rrlda(x, grouping=(iris[,5]), lambda=0.2, hp=0.75) ## perform rrlda
pred <- predict(rr, x) ## predict
table(as.numeric(pred$class), as.numeric(iris[,5])) ## show errors

length(unique(as.numeric(iris[,5])))
#----------
library("abcrlda")


# ---------- Iris dataset test -------------
data(iris)

traindata0 = iris[ which(iris[,ncol(iris)]=='virginica'), 1:4]
traindata1 = iris[ which(iris[,ncol(iris)]=='versicolor'), 1:4]

rr <- abcrlda(traindata0, traindata1, gamma = 0.5, cost_10 = 0.75)

res0 = predict_abcrlda(rr, as.matrix(traindata0))
res1 = predict_abcrlda(rr, as.matrix(traindata1))

err0 = sum(res0)/length(res0)
err1 = 1 - sum(res1)/length(res1)
list(err0, err1)
# table(as.numeric(libpred$class), as.numeric(iris[,5])) ## show errors



# ---------- Gausian distribution test -------------

# library("mvtnorm")
# sigma <- matrix(c(1,0,0,1), ncol=2)
# x <- rmvnorm(n=500, mean=c(0,0), sigma=sigma)
# colMeans(x)
# var(x)
# plot(x)
# sigma

bayes_error = pnorm(2)
std_my = 0.5
N = 100
x0 <- matrix(rnorm(N, mean = 0, sd = std_my), ncol=2)
x1 <- matrix(c(rnorm(N/2, mean = 0, sd = std_my), rnorm(N/2, mean = 2, sd = std_my)), ncol=2)


rr <- abcrlda(x0, x1, gamma = 0.5, cost_10 = 0.8)
res0 = predict_abcrlda(rr, x0)
res1 = predict_abcrlda(rr, x1)

err0 = sum(res0)/length(res0)
err1 = 1 - sum(res1)/length(res1)

p1 = -2 * rr$a[1] / rr$a[2] - rr$m / rr$a[2]
p2 = 2 * rr$a[1] / rr$a[2] - rr$m / rr$a[2]
separator = matrix(c(-2, 2, p1, p2), ncol = 2)

plot(x0, col = rgb(red = 1, green = 0, blue = 0, alpha = 0.8),
     xlim=c(-2, 2), ylim = c(-1.5,3.5), pch = 16)
points(x1, col = rgb(red = 0, green = 0, blue = 1, alpha = 0.8), pch = 16)
lines(separator, lwd=3)
legend('topright', c("1", "0"), col=c('blue', 'red'), pch=c(16,16))

# ---------- estimate ------

rr <- abcrlda(x0, x1, gamma = 0.5, cost_10 = 0.8)
# res0 = predict_abcrlda(rr, x0)
# res1 = predict_abcrlda(rr, x1)

# error_estimate_29(rr,0)
# error_estimate_29(rr,1)
bayes_error - risk_estimate_20(rr)

# -------------------- GRID Search -----------------

std_my = 0.5
N = 1000
x0 <- matrix(rnorm(N, mean = 0, sd = std_my), ncol=2)
x1 <- matrix(c(rnorm(N/2, mean = 0, sd = std_my), rnorm(N/2, mean = 2, sd = std_my)), ncol=2)

grid_search(x0, x1, range_gamma = seq(0.1, 10.1, by=0.1), range_C_10 = seq(0.5, 0.5, by=0.1))

#######
data(iris)

traindata0 = iris[ which(iris[,ncol(iris)]=='virginica'), 1:4]
traindata1 = iris[ which(iris[,ncol(iris)]=='versicolor'), 1:4]

grid_search(traindata0, traindata1, range_gamma = seq(0.1, 10.1, by=0.1), range_C_10 = seq(0.3, 0.8, by=0.1))

