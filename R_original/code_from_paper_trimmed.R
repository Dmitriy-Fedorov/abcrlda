###############################################
# Author: Amin Zollanvari
# This is the main function. Using this file,
# we can train both RLDA and ABc-RLDA classifiers
# and compute their errors vs. gamma
###############################################


library("abcrlda")

trainData = read.csv(file="s5/SigPar_Train.csv", header=FALSE, sep=",")
testData = read.csv(file="s5/SigPar_Test.csv", header=FALSE, sep=",")


#N is the majority class or class 0

gammav=c(1,10,100)
cost = c(0.2,0.9)        	       #cost values
N = 500       							   #number of MC repetitions
n_train = c(120, 300)					   #samples
features = dim(trainData)[2]
samples = dim(trainData)[1]
un0 = sum(trainData[,features]=='N')
un1 = sum(trainData[,features]=='P')
r0 = un0/samples


##################################### TRAIN DATA  ##########

data0 = trainData[ which(trainData[,ncol(trainData)]=='N'), ]
data1 = trainData[ which(trainData[,ncol(trainData)]=='P'), ]

data0 = data0[,1:(features-1)]
data1 = data1[,1:(features-1)]

################## TEST DATA #####################################

testdata0 = testData[ which(testData[,ncol(trainData)]=='N'), ]
testdata1 = testData[ which(testData[,ncol(trainData)]=='P'), ]

testdata0 = testdata0[,1:(features-1)]
testdata1 = testdata1[,1:(features-1)]
test = as.matrix(rbind(testdata0,testdata1))

nt0 = nrow(testdata0)
nt1 = nrow(testdata1)

true0 = rep(1,nt0)
true1 = rep(2,nt1)
true = c(true0,true1)

####################################

for (nn in 1:length(n_train)) {

  n = n_train[nn]
  n0 = floor(r0*n)
  n1 = n - n0


  for (cc in 1:length(cost)) {

    c = cost[cc]

    matrix.RLD = numeric()
    matrix.ABCRLD = numeric()
    matrix.REFRLD = numeric()
    matrix.MAIRLD = numeric()

    for (ll in 1:length(gammav)) {

      gamma = gammav[ll]
      kappa = gamma

      doubasytruevABCRLD = numeric()

      for (j in 1:N) {
        sample0 = sort(sample( 1:un0 , n0 ))
        sample1 = sort(sample( 1:un1 , n1 ))

        train0 = data0[sample0,]
        train1 = data1[sample1,]


        X1 = t(as.matrix(train0))
        X2 = t(as.matrix(train1))
        S1 = cov(t(X1))
        S2 = cov(t(X2))


        S = ((n0-1)*S1+(n1-1)*S2)/(n0+n1-2)
        Hinv = (diag(dim(X1)[1])+gamma*S)
        H = solve(Hinv)


        ####### compute the results for abc_rlda #####################

        a.vec.abc.rlda = abc_rlda(X1,X2,H,S,c,gamma,kappa)
        a=a.vec.abc.rlda[2:(length(a.vec.abc.rlda))]
        m=a.vec.abc.rlda[1]
        error1 = 0
        error2 = 0

        for (i in 1:nrow(test)){
          predict = test_lda(test[i,],a,m,0)
          if(i <= nt0){
            error1 = error1 + abs(true[i]-predict)
          }else{
            error2 = error2 + abs(true[i]-predict)
          }
        }
        doub_asy_true_ABCRLD = c * error1/(nt0) + (1-c) * error2/(nt1)

        ############################################################

        doubasytruevABCRLD = c(doubasytruevABCRLD,doubasytrueABCRLD)

      } #for (j )

      matrix.ABCRLD=cbind(matrix.ABCRLD,doubasytruevABCRLD)


    } # for (ll)

    save(matrix.ABCRLD, file = paste(paste(paste("s5/real_matrix_ABCRLD_delta2_2_c_",c, sep=""),"n_",n, sep=""),".RData", sep=""))
  } #for (cc)
} #for (nn)





