---
layout: post
title: "Permutation test for principal component analysis"
date: 2015-01-02 00:35:51 -0500
comments: true
categories: [Programming, R, Statistics] 
---

The procedure of permutation test for PCA is as follows:

For each replicate,

1. Individually permute each column of the data matrix.

2. Conduct the PCA and find the proportion of variance explained by each of the components 1 to s. Store this information.

3. Repeat 1 and 2 R times.

At the end of this we will have a matrix with R rows and s columns that contains the proportion of variance explained by each component for each replicate.

Finally, compare the observed values from the original data to the set of values from the permutations in order to determine the approximate p-value.

The R code:

``` R pca_perm.R

# the fuction to assess the significance of the principal components.
sign.pc<-function(x,R=1000,s=10, cor=T,...){
  # run PCA
  pc.out<-princomp(x,cor=cor,...)
  # the proportion of variance of each PC
  pve=(pc.out$sdev^2/sum(pc.out$sdev^2))[1:s]
  
  # a matrix with R rows and s columns that contains
  # the proportion of variance explained by each pc
  # for each randomization replicate.
  pve.perm<-matrix(NA,ncol=s,nrow=R)
  for(i in 1:R){
    # permutation each column
    x.perm<-apply(x,2,sample)
    # run PCA
    pc.perm.out<-princomp(x.perm,cor=cor,...)
    # the proportion of variance of each PC.perm
    pve.perm[i,]=(pc.perm.out$sdev^2/sum(pc.perm.out$sdev^2))[1:s] 
  }
  # calcalute the p-values
  pval<-apply(t(pve.perm)>pve,1,sum)/R
  return(list(pve=pve,pval=pval))
}


# apply the function
library(RCurl)
data <- getURL("https://raw.githubusercontent.com/bioops/mis_scripts/master/statistics/data/pca.txt")
OCRdata <- read.table(text = data, header=T,sep="\t")
OCRdat<-OCRdata[,-1] #leave out location id column
sign.pc(OCRdat,cor=T)

```

The result:

<pre>
$pve
    Comp.1     Comp.2     Comp.3     Comp.4     Comp.5     Comp.6     Comp.7     Comp.8 
0.23129378 0.14864525 0.11552865 0.06741744 0.06274641 0.05858431 0.05033795 0.04484122 
    Comp.9    Comp.10 
0.03873311 0.03431297 

$pval
 [1] 0.000 0.000 0.000 1.000 1.000 0.996 1.000 1.000 1.000 1.000
</pre>
