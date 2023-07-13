rm(list=ls())
install.packages("Mass")
library(Mass)
fun1<-function(n11,n12,n21,n22,n31,n32,v11,v12,v21,v22,v31,v32)
{
  g11<-rnorm(n11,0,sqrt(v11));g12<-rnorm(n12,0,sqrt(v12))
  g21<-rnorm(n21,0,sqrt(v21));g22<-rnorm(n22,0,sqrt(v22))
  g31<-rnorm(n31,0,sqrt(v31));g32<-rnorm(n32,0,sqrt(v32))
  s11=var(g11);s12=var(g12);s21=var(g21);s22=var(g22);s31=var(g31);s32=var(g32)
  S1=(s11/n11+s12/n12)/4;S2=(s21/n21+s22/n22)/4;S3=(s31/n31+s32/n32)/4
  X1=(mean(g11)+mean(g12))/2;X2=(mean(g21)+mean(g22))/2;X3=(mean(g31)+mean(g32))/2
  v1<-sqrt(S1+S2);v2<-sqrt(S2+S3)
  T1=(X2-X1)/v1;T2=(X3-X2)/v2
  # min-T value
  A=min(T1,T2,na.rm = FALSE)
  return(A)
}
fun2<-function(n11,n12,n21,n22,n31,n32,v11,v12,v21,v22,v31,v32)
{
  x<-replicate(1000,fun1(n11,n12,n21,n22,n31,n32,v11,v12,v21,v22,v31,v32))
  y<-sort(x,decreasing=FALSE)
  c<-y[950]
  return(c)
}
#mu11=3;mu12=0;mu13=0;mu21=3;mu22=0;mu23=0;mu31=3;mu32=0;mu33=0;n11=20;n12=20;n13=20;n21=20;n22=20;n23=20;n31=20;n32=20;n33=20;v11=1;v12=1;v13=1;v21=1;v22=1;v23=1;v31=1;v32=1;v33=1
fun3<-function(n11,n12,n21,n22,n31,n32,mu11,mu12,mu21,mu22,mu31,mu32,v11,v12,v21,v22,v31,v32)
{
  ##For generating samples from Normal distributions
  g11<-rnorm(n11,mu11,sqrt(v11));g12<-rnorm(n12,mu12,sqrt(v12))
  g21<-rnorm(n21,mu21,sqrt(v21));g22<-rnorm(n22,mu22,sqrt(v22))
  g31<-rnorm(n31,mu31,sqrt(v31));g32<-rnorm(n32,mu32,sqrt(v32))
  ##Generating samples from t distribution
  #Y11=rt(n11,3,ncp=0)/sqrt(3);Y12=rt(n12,3,ncp=0)/sqrt(3);Y13=rt(n13,3,ncp=0)/sqrt(3)
  #Y21=rt(n21,3,ncp=0)/sqrt(3);Y22=rt(n22,3,ncp=0)/sqrt(3);Y23=rt(n23,3,ncp=0)/sqrt(3)
  #Y31=rt(n31,3,ncp=0)/sqrt(3);Y32=rt(n32,3,ncp=0)/sqrt(3);Y33=rt(n33,3,ncp=0)/sqrt(3)
  #g11=sqrt(v11)*Y11+mu11;g12=sqrt(v12)*Y12+mu12;g13=sqrt(v13)*Y13+mu13
  #g21=sqrt(v21)*Y21+mu21;g22=sqrt(v22)*Y22+mu22;g23=sqrt(v23)*Y23+mu23
  #g31=sqrt(v31)*Y31+mu31;g32=sqrt(v32)*Y32+mu32;g33=sqrt(v33)*Y33+mu33
  ##Generating samples from exponential distribution
  #Y11=rexp(n11,1);Y12=rexp(n12,1);Y13=rexp(n13,1);Y21=rexp(n21,1);Y22=rexp(n22,1);Y23=rexp(n23,1)
  #Y31=rexp(n31,1);Y32=rexp(n32,1);Y33=rexp(n33,1)
  #g11=mu11+sqrt(v11)*Y11;g12=mu12+sqrt(v12)*Y12;g13=mu13+sqrt(v13)*Y13;g21=mu21+sqrt(v21)*Y21;g22=mu22+sqrt(v22)*Y22
  #g23=mu23+sqrt(v23)*Y23;g31=mu31+sqrt(v31)*Y31;g32=mu32+sqrt(v32)*Y32;g33=mu33+sqrt(v33)*Y33
  ##For generating samples from Gamma distribution
  #Y11=(rgamma(n11,2,1)-2)/sqrt(2);Y12=(rgamma(n12,2,1)-2)/sqrt(2);Y13=(rgamma(n13,2,1)-2)/sqrt(2);Y21=(rgamma(n21,2,1)-2)/sqrt(2)
  #Y22=(rgamma(n22,2,1)-2)/sqrt(2);Y23=(rgamma(n23,2,1)-2)/sqrt(2);Y31=(rgamma(n31,2,1)-2)/sqrt(2);Y32=(rgamma(n32,2,1)-2)/sqrt(2)
  #Y33=(rgamma(n33,2,1)-2)/sqrt(2)
  #g11=mu11+sqrt(v11)*Y11;g12=mu12+sqrt(v12)*Y12;g13=mu13+sqrt(v13)*Y13;g21=mu21+sqrt(v21)*Y21;g22=mu22+sqrt(v22)*Y22
  #g23=mu23+sqrt(v23)*Y23;g31=mu31+sqrt(v31)*Y31;g32=mu32+sqrt(v32)*Y32;g33=mu33+sqrt(v33)*Y33
  ##Generating samples from Laplace distribution
  #g11=rdoublex(n11,mu11,sqrt(v11)/sqrt(2));g12=rdoublex(n12,mu12,sqrt(v12)/sqrt(2))
  #g13=rdoublex(n13,mu13,sqrt(v13)/sqrt(2))
  #g21=rdoublex(n21,mu21,sqrt(v21)/sqrt(2));g22=rdoublex(n22,mu22,sqrt(v22)/sqrt(2))
  #g23=rdoublex(n23,mu23,sqrt(v23)/sqrt(2))
  #g31=rdoublex(n31,mu31,sqrt(v31)/sqrt(2));g32=rdoublex(n32,mu32,sqrt(v32)/sqrt(2))
  #g33=rdoublex(n33,mu33,sqrt(v33)/sqrt(2))
  ##Generating samples from Weibull distributions
  #Y11=(rweibull(n11,2,1)-sqrt(pi)/2)/sqrt(1-pi/4);Y12=(rweibull(n12,2,1)-sqrt(pi)/2)/sqrt(1-pi/4)
  #Y13=(rweibull(n13,2,1)-sqrt(pi)/2)/sqrt(1-pi/4)
  #Y21=(rweibull(n21,2,1)-sqrt(pi)/2)/sqrt(1-pi/4);Y22=(rweibull(n22,2,1)-sqrt(pi)/2)/sqrt(1-pi/4)
  #Y23=(rweibull(n23,2,1)-sqrt(pi)/2)/sqrt(1-pi/4)
  #Y31=(rweibull(n31,2,1)-sqrt(pi)/2)/sqrt(1-pi/4);Y32=(rweibull(n32,2,1)-sqrt(pi)/2)/sqrt(1-pi/4)
  #Y33=(rweibull(n33,2,1)-sqrt(pi)/2)/sqrt(1-pi/4)
  #g11=mu11+sqrt(v11)*Y11;g12=mu12+sqrt(v12)*Y12;g13=mu13+sqrt(v13)*Y13
  #g21=mu21+sqrt(v21)*Y21;g22=mu22+sqrt(v22)*Y22;g23=mu23+sqrt(v23)*Y23
  #g31=mu31+sqrt(v31)*Y31;g32=mu32+sqrt(v32)*Y32;g33=mu33+sqrt(v33)*Y33
  s11=var(g11);s12=var(g12);s21=var(g21);s22=var(g22);s31=var(g31);s32=var(g32)
  S1=(s11/n11+s12/n12)/4;S2=(s21/n21+s22/n22)/4;S3=(s31/n31+s32/n32)/4
  X1=(mean(g11)+mean(g12))/2;X2=(mean(g21)+mean(g22))/2;X3=(mean(g31)+mean(g32))/2
  v1<-sqrt(S1+S2);v2<-sqrt(S2+S3)
  T1=(X2-X1)/v1;T2=(X3-X2)/v2
  out=fun2(n11,n12,n21,n22,n31,n32,s11,s12,s21,s22,s31,s32)
  # min-T value
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

