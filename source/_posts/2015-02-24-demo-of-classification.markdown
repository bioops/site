---
layout: post
title: "Demo of classification"
date: 2015-02-24 12:35:23 -0500
comments: true
categories: [programming, r, statistics] 
---

R code demo of

1. Linear discriminant analysis (LDA)
2. Quadratic discriminant analysis (QDA)
3. k-nearest neighbor (KNN)

``` R

# read data
library(RCurl)
data <- getURL("https://raw.githubusercontent.com/bioops/mis_scripts/master/statistics/data/admission.txt")
admission<-read.table(text=data,header=T)
admission$CLASS<-as.factor(admission$CLASS)

# scale the data
admission$GMAT<-scale(admission$GMAT)
admission$GPA<-scale(admission$GPA)

# (1) LDA
library(MASS)
lda.fit<-lda(CLASS~GMAT+GPA,data=admission)
lda.pred<-predict(lda.fit)
# plot
# different symbols represents true classifications
# different colors are predicted classifications
plot(admission[,1:2],col=lda.pred$class, pch=as.numeric(admission$CLASS),
     xlab="GPA",ylab="GMAT",main="LDA")

```

![LDA](http://bioops.info/images/uploads/2015/class1.png)

``` R

# (2) QDA
qda.fit<-qda(CLASS~GMAT+GPA,data=admission)
qda.pred<-predict(qda.fit)
# plot
# different symbols represents true classifications
# different colors are predicted classifications
plot(admission[,1:2],col=qda.pred$class, pch=as.numeric(admission$CLASS),
     xlab="GPA",ylab="GMAT",main="QDA")

```

![QDA](http://bioops.info/images/uploads/2015/class2.png)

``` R

library(class)

# tune parameter using cross-validation (CV)
k<-seq(1,15,1) # different k
cv.err<-NULL # cv error
for (ki in k){
  knn.pred.cv<-knn.cv(admission[,1:2],admission[,3],k=ki)
  cv.err<-c(cv.err, mean(knn.pred.cv!=admission[,3]))
}
plot(k,cv.err, main="CV error vs k")

# using the optimal parameter
knn.pred<-knn(admission[,1:2],admission[,1:2], admission[,3],k=k[which.min(cv.err)])

# plot
# different symbols represents true classifications
# different colors are predicted classifications
plot(admission[,1:2],col=knn.pred, pch=as.numeric(admission$CLASS),
     xlab="GPA",ylab="GMAT",main="KNN")

```

![KNN](http://bioops.info/images/uploads/2015/class3.png)
