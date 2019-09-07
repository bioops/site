---
layout: post
title: "Estimate Gamma distribution parmaters using MME and MLE"
date: 2015-01-01 07:13:59 -0500
comments: true
categories: [Programming, R, Statistics]
---

This post shows how to estimate gamma distribution parameters using (a) moment of estimation (MME) and (b) maximum likelihood estimate (MLE).

The probability density function of Gamma distribution is

$$
\frac{1}{\Gamma (\alpha) \beta ^{\alpha}} x^{\alpha - 1} e^{- \frac{x}{\beta}}
$$

The MME:

$$
\hat{\alpha}=\frac{n\bar{X} ^2}{\sum_{i=1}^{n} (X_i-\bar{X})^2}
$$

$$
\hat{\beta}=\frac{\sum_{i=1}^{n} (X_i-\bar{X})^2}{n \bar{X}}
$$

We can calculate the MLE of $ \alpha $ using the Newton-Raphson method.

For $ k =1,2,...,$

$$
\hat{\alpha} ^{(k)}=\hat{\alpha} ^{(k-1)} - \frac{\ell'(\hat{\alpha} ^{(k-1)})}{\ell'' (\hat{\alpha} ^{(k-1)})}
$$

where

$$
\ell' (\alpha) = n \log \left(\frac{\alpha}{\bar{X}}\right)-n \frac{\Gamma '(\alpha)}{\Gamma(\alpha)}+\sum_{i=1}^{n} \log X_i
$$

$$
\ell'' (\alpha) = \frac{n}{\alpha} - n \left(\frac{\Gamma '(\alpha)}{\Gamma (\alpha)}\right)'
$$

Use the MME for the initial value of $ \alpha^{(0)} $, and stop the approximation when $ \vert \hat{\alpha}^{(k)}-\hat{\alpha}^{(k-1)} \vert < 0.0000001 $. The MLE of $ \beta $ can be found by $ \hat{\beta} = \bar{X} / \hat{\alpha} $.

Below is the R code.

``` R gamma.R

# (a) MME
gamma_MME<-function(x){
  n<-length(x)
  mean_x<-mean(x)
  alpha<-n*(mean_x^2)/sum((x-mean_x)^2)
  beta<-sum((x-mean_x)^2)/n/mean_x
  estimate_MME<-data.frame(alpha,beta)
  return(estimate_MME)
}


# (b) MLE
gamma_MLE<-function(x){
  n<-length(x)
  mean_x<-mean(x)
  
  # initiate the convergence and alpha value
  converg<-1000
  alpha_prev<-gamma_MME(x)$alpha
  
  # initiate two vectors to store alpha and beta in each step
  alpha_est<-alpha_prev
  beta_est<-mean_x/alpha_prev
  
  # Newton-Raphson
  while(converg>0.0000001){
    #first derivative of alpha_k-1
    der1<-n*log(alpha_prev/mean_x)-n*digamma(alpha_prev)+sum(log(x))
    #second derivative of alpha_k-1
    der2<-n/alpha_prev-n*trigamma(alpha_prev)
    #calculate next alpha
    alpha_next<-alpha_prev-der1/der2
    # get the convergence value
    converg<-abs(alpha_next-alpha_prev)
    # store estimators in each step
    alpha_est<-c(alpha_est, alpha_next)
    beta_est<-c(beta_est, mean_x/alpha_next)
    # go to next alpha
    alpha_prev<-alpha_next
  }

  alpha<-alpha_next
  beta<-mean_x/alpha_next
  estimate_MLE<-data.frame(alpha,beta)

  return(estimate_MLE)
}

# apply
x<-rgamma(100,2,scale=5)
gammma_MME(x)
gamma_MLE(x)

```
 
