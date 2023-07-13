rm(list=ls())
install.packages("MASS")
# package MASS is needed for generating sample from multivariate normal distribution
library(MASS)
fun1<-function(d1,d2,s1)
{
  D=matrix(c(d1,0,0,d2),nrow=2,byrow=TRUE)
  S=matrix(c(d1^2,-s1,-s1,d2^2),nrow=2,byrow=TRUE)
  P=solve(D)%*%S%*%solve(D)
  mu=c(0,0)
  g<-mvrnorm(1,mu,P)
  T1=g[1];T2=g[2]
  A=min(T1,T2,na.rm=FALSE)
  return(A)
}
fun2<-function(d1,d2,s1)
{
  x<-replicate(1000,fun1(d1,d2,s1))
  y<-sort(x,decreasing=FALSE)
  c<-y[950]
  return(c)
}
fun3<-function(n11,n12,n21,n22,n31,n32,mu11,mu12,mu21,mu22,mu31,mu32,v11,v12,v21,v22,v31,v32)
{
  g11<-rnorm(n11,mu11,sqrt(v11));g12<-rnorm(n12,mu12,sqrt(v12))
  g21<-rnorm(n21,mu21,sqrt(v21));g22<-rnorm(n22,mu22,sqrt(v22))
  g31<-rnorm(n31,mu31,sqrt(v31));g32<-rnorm(n32,mu32,sqrt(v32))
  s11=var(g11);s12=var(g12);s21=var(g21);s22=var(g22);s31=var(g31);s32=var(g32)
  S1=(s11/n11+s12/n12)/4;S2=(s21/n21+s22/n22)/4;S3=(s31/n31+s32/n32)/4
  X1=(mean(g11)+mean(g12))/2;X2=(mean(g21)+mean(g22))/2;X3=(mean(g31)+mean(g32))/2
  v1<-sqrt(S1+S2);v2<-sqrt(S2+S3)
  T1=(X2-X1)/v1;T2=(X3-X2)/v2
  N=n11+n12+n21+n22+n31+n32
  nu11=n11/N;nu12=n12/N;nu21=n21/N;nu22=n22/N
  nu31=n31/N;nu32=n32/N
  d1=(sqrt(s21/nu21+s22/nu22+s11/nu11+s12/nu12))/2
  d2=(sqrt(s31/nu31+s32/nu32+s21/nu21+s22/nu22))/2
  s1=(s21/nu21+s22/nu22)/4
  out=fun2(d1,d2,s1)
  A=min(T1,T2,na.rm = FALSE)
  a=0
  if(A>out)
    a=a+1
  return(a)
}
fun4<-function(n11,n12,n21,n22,n31,n32,mu11,mu12,mu21,mu22,mu31,mu32,v11,v12,v21,v22,v31,v32)
{
  # find number of times maxt-T value > critical value among 1000 values
  out<-replicate(1000,fun3(n11,n12,n21,n22,n31,n32,mu11,mu12,mu21,mu22,mu31,mu32,v11,v12,v21,v22,v31,v32))
  p<-sum(out)/1000
  return(p)
}

##Size values of Table 3.3


X=replicate(5,fun4(5,6,5,6,5,6,0,0,0,0,0,0,4,75,67,53,96,52))
X;mean(X)
X=replicate(5,fun4(8,8,8,8,8,8,0,0,0,0,0,0,4,75,67,53,96,52))
X;mean(X)
X=replicate(5,fun4(15,5,6,10,18,20,0,0,0,0,0,0,4,75,67,53,96,52))
X;mean(X)
X=replicate(5,fun4(10,5,6,60,50,100,0,0,0,0,0,0,4,75,67,53,96,52))
X;mean(X)
X=replicate(5,fun4(5,10,15,20,25,30,0,0,0,0,0,0,4,75,67,53,96,52))
X;mean(X)
X=replicate(5,fun4(30,25,20,15,10,5,0,0,0,0,0,0,4,75,67,53,96,52))
X;mean(X)


X=replicate(5,fun4(5,6,5,6,5,6,0,0,0,0,0,0,52,96,53,67,75,4))
X;mean(X)
X=replicate(5,fun4(8,8,8,8,8,8,0,0,0,0,0,0,52,96,53,67,75,4))
X;mean(X)
X=replicate(5,fun4(15,5,6,10,18,20,0,0,0,0,0,0,52,96,53,67,75,4))
X;mean(X)
X=replicate(5,fun4(10,5,6,60,50,100,0,0,0,0,0,0,52,96,53,67,75,4))
X;mean(X)
X=replicate(5,fun4(5,10,15,20,25,30,0,0,0,0,0,0,52,96,53,67,75,4))
X;mean(X)
X=replicate(5,fun4(30,25,20,15,10,5,0,0,0,0,0,0,52,96,53,67,75,4))
X;mean(X)

X=replicate(5,fun4(5,6,5,6,5,6,0,0,0,0,0,0,10^2,7^2,3^2,2^2,2^2,1))
X;mean(X)
X=replicate(5,fun4(8,8,8,8,8,8,0,0,0,0,0,0,10^2,7^2,3^2,2^2,2^2,1))
X;mean(X)
X=replicate(5,fun4(15,5,6,10,18,20,0,0,0,0,0,0,10^2,7^2,3^2,2^2,2^2,1))
X;mean(X)
X=replicate(5,fun4(10,5,6,60,50,100,0,0,0,0,0,0,10^2,7^2,3^2,2^2,2^2,1))
X;mean(X)
X=replicate(5,fun4(5,10,15,20,25,30,0,0,0,0,0,0,10^2,7^2,3^2,2^2,2^2,1))
X;mean(X)
X=replicate(5,fun4(30,25,20,15,10,5,0,0,0,0,0,0,10^2,7^2,3^2,2^2,2^2,1))
X;mean(X)

X=replicate(5,fun4(5,6,5,6,5,6,0,0,0,0,0,0,1,2^2,2^2,3^2,7^2,10^2))
X;mean(X)
X=replicate(5,fun4(8,8,8,8,8,8,0,0,0,0,0,0,1,2^2,2^2,3^2,7^2,10^2))
X;mean(X)
X=replicate(5,fun4(15,5,6,10,18,20,0,0,0,0,0,0,1,2^2,2^2,3^2,7^2,10^2))
X;mean(X)
X=replicate(5,fun4(10,5,6,60,50,100,0,0,0,0,0,0,1,2^2,2^2,3^2,7^2,10^2))
X;mean(X)
X=replicate(5,fun4(5,10,15,20,25,30,0,0,0,0,0,0,1,2^2,2^2,3^2,7^2,10^2))
X;mean(X)
X=replicate(5,fun4(30,25,20,15,10,5,0,0,0,0,0,0,1,2^2,2^2,3^2,7^2,10^2))
X;mean(X)


