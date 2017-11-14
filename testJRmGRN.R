source("JRmGRN.R")
X1=read.csv("net1.csv",header=TRUE)
X1=as.matrix(X1)
X2=read.csv("net2.csv",header=TRUE)
X2=as.matrix(X2)

S1=cov(X1)
S2=cov(X2)
S=list(S1,S2)
n1=nrow(X1)
n2=nrow(X2)
ns=list(n1,n2)
p=ncol(X1)
b=n1

lambda1=0.143845
lambda2=0.115076
lambda3=0.143845
lambda4=2.157675
ans=JRmGRN(S,ns,lambda1=lambda1*b,lambda2=lambda2*b,lambda3=lambda3*b,lambda4=lambda4*b,verbose=FALSE)
