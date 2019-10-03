data(iris)
train_data <- iris[which(iris[, ncol(iris)] == "virginica" |
                        iris[, ncol(iris)] == "versicolor"), 1:4]
train_label <- factor(iris[which(iris[, ncol(iris)] == "virginica" |
                                iris[,ncol(iris)] == "versicolor"), 5])
cost_range <- seq(0.1, 0.9, by=0.2)
gamma_range <- c(0.1,1,10,100,1000)

gs <- grid_search(train_data, train_label,
                  range_gamma = gamma_range,
                  range_C_10 = cost_range,
                  method = "estimator")
model <- abcrlda(train_data, train_label, gamma = gs$gamma[1], cost_10 = gs$C_10[1])