---
layout: post
title: "Solving bridge regression using local quadratic approximation (LQA)"
date: 2015-02-24 11:27:50 -0500
comments: true
categories: [Programming, R, Statistics] 
---

Bridge regression is a broad class of penalized regression, and can be used in high-dimensional regression problems.

It includes the ridge (q=2) and lasso (q  =1) as special cases.

More technical details can be found [here](http://www.sciencedirect.com/science/article/pii/S0378375811001960). Below R code demonstrates:

1. sovling bridge regression using local quadratic approximation (LQA) and Newton–Raphson algorithm.
2. simulation of tuning parameters using 50/50/200 observations (training/validation/testing).


``` R
# The function to solve bridge regression
bridge<-function(x, y, lambda, q=1, eta=0.001){
  library(glmnet)
  # use ridge coefficients as a starting value
  beta.start<-coef(glmnet(x, y, alpha=0, lambda=lambda, intercept=F))
  beta.mat<-matrix(NA,ncol=ncol(beta.start),nrow=nrow(beta.start)-1)
  
  # take each lambda value in the grid
  for(i in 1:length(grid)){
    # initial beta without intercept
    beta_prev<-as.vector(beta.start[-1,i])
    # initial converge
    converge<-10^10
    
    # iteration until converge or two many iterations
    iteration<-0
    while(converge>eta && iteration<=100){
      iteration<-iteration+1
      # judge whether some beta is small enough
      del<-which(abs(beta_prev)<eta)
      
      # if all coefficient are small enough then stop the iteration
      if(length(del)==length(beta_prev)){
        beta_prev[del]<-0
        converge<-0
        
        # else if we need to remove some but not all of the coefficients  
      }else if(length(del)>0){
        # set these beta to 0
        beta_prev[del]<-0
        
        # update design matrix x
        x.new<-x;x.new<-x.new[,-del]
        
        # calculate the diagonal matrix involving penalty terms
        # and the next beta
        if(length(beta_prev)-length(del)==1){
          diag<-grid[i]*q*(abs(beta_prev[-del])^(q-2))/2
        }else{
          diag<-diag(grid[i]*q*(abs(beta_prev[-del])^(q-2))/2)
        }
        
        # next beta
        beta_curr<-beta_prev
        beta_curr[-del]<-solve(t(x.new)%*%x.new+diag)%*%t(x.new)%*%y
        
        # new converge value
        converge<-sum((beta_curr-beta_prev)^2)
        # next iteration
        beta_prev<-beta_curr
        
        # if we don't need to remove the coefficients
      }else{
        x.new<-x
        diag<-diag(grid[i]*q*(abs(beta_prev)^(q-2))/2)
        # next beta
        beta_curr<-solve(t(x.new)%*%x.new+diag)%*%t(x.new)%*%y
        
        # new converge value
        converge<-sum((beta_curr-beta_prev)^2)
        # next iteration
        beta_prev<-beta_curr
        
      }
    }
    
    beta.mat[,i]<-beta_prev
  }
  colnames(beta.mat)<-grid
  return(beta.mat)
}

# The function to find the optimal coefficients set
# that minimize mse of validation set
bridge.opt<-function(coef.matrix, newx, newy){
  # calculate the mse
  mse.validate<-NULL
  for(i in 1:ncol(coef.matrix)){
    mse.validate<-c(mse.validate, mean((newx%*%coef.matrix[,i]-newy)^2))
  }
  # find the optimal coefficients set
  coef.opt<-coef.matrix[,which.min(mse.validate)]
  return(coef.opt)
}


##################
# start simulation
# For this simulation, create three data sets 
# consisting of 50/50/200 observations (training/validation/testing).
# Use validation data to select the tuning parameter (lambda).


# q settings
q_seq<-c(0.1,0.5,1,2)
# lambda grid
grid=10^seq(10,-2,length=100)
# training and testing set
ntrain<-50
ntest<-200
# validation set size
nvalidate<-50

# repeat number
repeatnum<-100
# initialize MSE
MSE<-matrix(NA, nrow=repeatnum,ncol=length(q_seq))
colnames(MSE)<-q_seq

library(MASS)
# the beta
beta<-matrix(c(3,1.5,0,0,2,0,0,0),ncol=1)
# covariance matrix of X
cov<-matrix(NA,8,8)
for (i in 1:8){
  for (j in 1:8){
    cov[i,j]<-0.5^(abs(i-j))
  }
}

for (i in 1:repeatnum){
  # generate X
  x.train<-mvrnorm(n=ntrain,rep(0,8),cov)
  x.test <-mvrnorm(n=ntest,rep(0,8),cov)
  x.validate <-mvrnorm(n=nvalidate,rep(0,8),cov)
  # generate error term
  err.train<-matrix(rnorm(n=ntrain,mean=0,sd=3),ncol=1)
  err.test<-matrix(rnorm(n=ntest,mean=0,sd=3),ncol=1)
  err.validate<-matrix(rnorm(n=nvalidate,mean=0,sd=3),ncol=1)
  # calculate Y
  y.train<-x.train%*%beta+err.train
  y.test<-x.test%*%beta+err.test
  y.validate<-x.validate%*%beta+err.validate
  # centralize Y and standardize X
  y.train<-scale(y.train,scale=F);y.test<-scale(y.test,scale=F)
  y.validate<-scale(y.validate,scale=F)
  x.train<-scale(x.train);x.test<-scale(x.test)
  x.validate<-scale(x.validate)
  
  # run for each q
  for(j in 1:length(q_seq)){
    
    # coefficients for trainning set
    bridge.train<-bridge(x.train, y.train, grid, q=q_seq[j], eta=0.001)
    # find the optimal coefficients set that minimize mse of validation set
    coef.opt<-bridge.opt(bridge.train, x.validate, y.validate)
    # MSE on the test set using the optimal model
    MSE[i,j]<-mean((x.test%*%coef.opt-y.test)^2)
  }
}

# mean of MSEs under different models
apply(MSE,2,mean)
# sd of MSEs under different models
apply(MSE,2,sd)
# boxplot of MSEs under different models
boxplot(MSE, ylab="MSE", xlab="q")
```

![bridge MSE plots](http://bioops.info/images/uploads/2015/bridge.png)
