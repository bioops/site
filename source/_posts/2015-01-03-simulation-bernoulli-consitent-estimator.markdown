---
layout: post
title: "The consistent estimator of bernouli distribution"
date: 2015-01-03 20:53:14 -0500
comments: true
categories: [Programming, R, Statistics]
---

This is a simple post showing the basic knowledge of statistics, the consistency.

For Bernoulli distribution, $ Y \sim B(n,p) $, $ \hat{p}=Y/n $ is a consistent estimator of $ p $, because:

$$
\lim_{n \to \infty} \left(p-\epsilon<\frac{Y}{n}<p+\epsilon\right)=1,
$$

for any positive number $ \epsilon $.

Here is the simulation to show the estimator is consitent.

``` R

# set parameters
n<-1000;p<-0.5;
# n Bernoulli trails
obs<-rbinom(n,1,p)
# estimate p on different number of trials.
phat<-cumsum(obs)/cumsum(rep(1,n))
# the convergence plot
plot(phat, type="l", xlab="Trails")
abline(h=p)

```

![convergence plot](http://bioops.info/images/uploads/2015/Rplot01.png)

``` R

# then, 100 repetitions

# set parameters
n<-1000;p<-0.5;B<-100;
# n*B Bernoulli trails
obs<-rbinom(n*B,1,p)
# convert n*B observations to a n*B matrix
obs_mat<-matrix(obs, nrow=n, ncol=B)
# a function to estimate p on different number of trials
est_p<-function(x,n) cumsum(x)/cumsum(rep(1,n))
# estimate p on different number of trials for each repetition
phat_mat<-apply(obs_mat,2, est_p, n=n)
# the convergence plot with 100 repetitions
matplot(phat_mat,type="l",lty=1,xlab="Trials",ylab="phat")

```

![convergence plots](http://bioops.info/images/uploads/2015/Rplot02.png)

