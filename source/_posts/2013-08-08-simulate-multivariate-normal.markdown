---
layout: post
title: "Simulate multivariate normal distribution using R"
date: 2013-08-08 20:38
comments: true
categories: [Statistics, R, Programming] 
---
I have been doing some research about co-expression network. "co-expression" means that genes have similar expression profiles across different conditions or tissues. In the network, genes are nodes, and "co-expression" relationship between two genes can be reprensented as edges. The co-expressed genes may involve in similar pathways or biological process.

In a small part of my research, I am testing some algorithms to detect co-expression relationship. One way to test algorithm is simulation. In an ideal (simple) case, the expression values of two co-expressed genes can be considered as bivariate normal distributed. To generate expression values of such gene pair or a group of genes given a correlation coefficient, is just to simulate multivariate normal distribution. [MASS](http://cran.r-project.org/web/packages/MASS/index.html) library in R has an function, [mvrnorm](http://stat.ethz.ch/R-manual/R-patched/library/MASS/html/mvrnorm.html), to do that, but it requires a covariance matrix.

The function below is to firstly generate the covariance matrix in order to use the mvnorm function. Because we only know the correlation coefficient, i.e. co-expression relationship (degree), the mean and variance of each gene's expression profile are random generated in the function. Then the matrix can be calulated as follows.

$$
\mu=\left( 
\begin{matrix}
  \mu_x \\
  \mu_y
 \end{matrix}
\right), \Sigma=\left( 
\begin{matrix}
  \sigma_x^2 & \rho \sigma_x \sigma_y \\
  \sigma_x \sigma_y & \sigma_y^2
 \end{matrix}
\right)
$$

``` R multi_norm.R

# function to simulate multivariate normal distribution
# given gene number, sample size and correlation coefficient
multi_norm <- function(gene_num,sample_num,R) { 
  # initial covariance matrix
  V <- matrix(data=NA, nrow=gene_num, ncol=gene_num)
  
  # mean for each gene
  meansmodule <- runif(gene_num, min=-3, max=3)
  # variance for each gene
  varsmodule <- runif(gene_num, min=0, max=5) 
  
  for (i in 1:gene_num) {
  # a two-level nested loop to generate covariance matrix
    for (j in 1:gene_num) {
      if (i == j) {
        # covariances on the diagonal
        V[i,j] <- varsmodule[i]
      } else {
        # covariances
        V[i,j] <- R * sqrt(varsmodule[i]) * sqrt(varsmodule[j]) 
      }
    }
  }
 
  # simulate multivariate normal distribution
  # given means and covariance matrix
  X <- t(mvrnorm(n = sample_num, meansmodule, V))
 
  return(X)
}

```
