---
layout: post
title: "R function of Monte Carlo simulation to get the p-value from the joint cumulative distribution of an n-dimensional order statistic"
date: 2013-08-08 21:21
comments: true
categories: [Programming, R, Statistics]
---
I want to compute the P-value from the joint cumulative distribution of an
n-dimensional [order statistic](http://en.wikipedia.org/wiki/Order_statistic).

$$
P(r_1 r_2 ,..., r_n)=n! \int\limits_0^{r_1} \int\limits_{s_1}^{r_2} ... \int\limits_{s_{n-1}}^{r_n} ds_1ds_2...ds_n
$$


One efficient way is using the following recursive formula.

$$
P(r_1 r_2 ,..., r_n)=\sum_{i=1}^{n} (r_{n-i+1}-r_{n-i}) P(r_i r_2 ,..., r_{n-i},r_{n-i+2} ,..., r_n)
$$

**However, the facts are (or would be):**

1. I am too stupid to write a recursive function.
2. I didn't find the efficient formula at first.
3. In other cases, the efficient formula have not been derived yet, or too complicated to derive.

In Statistics [Monte Carlo simulation](http://en.wikipedia.org/wiki/Monte_Carlo_method) is a "quick" way to compute some complicated formulas. By saying "quick", I mean I can see the results without knowing or deriving "ugly" Math formulas. It's actually a very "slow" method in computing aspect.

Anyway, the R function is here.

``` R P_order_stat.R
# sub function of monte carlo simulation to get the p-value
P_order_stat <- function(ranks) {
  NumRep <- 10000 # number of replicates
  newvec <- sort(ranks) # sort the rank ratios
  pvalue <- 0 # inital pvalue
  for (i in 1:NumRep){
 
      # generate random uniform distributed data,
      # and then sort the simulated rank ratios
      newx <- sort(runif(length(ranks), min=0, max=1))  
 
      # if all the simulated data is lower than the input,
      # then sucess+1
      judge <- sum(newvec >= newx)
      if (judge == length(ranks)) pvalue<-pvalue+1 
  }
  pvalue <- pvalue / NumRep  # get the p-value
  return(pvalue)
}
```
