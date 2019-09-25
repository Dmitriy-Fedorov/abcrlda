data(iris)
traindata <- iris[which(iris[, ncol(iris)] == "virginica" |
                        iris[, ncol(iris)] == "versicolor"), 1:4]
trainlabel <- factor(iris[which(iris[, ncol(iris)] == "virginica" |
                                iris[,ncol(iris)] == "versicolor"), 5])
model <- abcrlda(traindata, trainlabel, gamma = 0.5, cost_10 = 0.75)
predict(model, traindata)