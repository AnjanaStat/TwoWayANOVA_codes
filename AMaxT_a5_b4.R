rm(list=ls())
install.packages("MASS")
# package MASS is needed for generating sample from multivariate normal distribution
library(MASS)
# by fun1 and fun2 we find critical value of asymptotic Max-T test
fun1<-function(d1,d2,d3,d4,s1,s2,s3)
{
  #find the estimate of the diagonal matrix D
  D=matrix(c(d1,0,0,0,0,d2,0,0,0,0,d3,0,0,0,0,d4),nrow=4,byrow=TRUE)
  # find the estimate of the dispersion matrix
  S=matrix(c(d1^2,-s1,0,0,-s1,d2^2,-s2,0,0,-s2,d3^2,-s3,0,0,-s3,d4^2),nrow=4,byrow=TRUE)
  P=solve(D)%*%S%*%solve(D)
  mu=c(0,0,0,0)
  #generate sample from multivariate normal with mean vector mu and covariance matrix P
  g<-mvrnorm(1,mu,P)
  T1=g[1];T2=g[2];T3=g[3];T4=g[4]
  #find the max-T value
  A=max(T1,T2,T3,T4,na.rm = FALSE)
  return(A)
}
fun2<-function(d1,d2,d3,d4,s1,s2,s3)
{
  # find 500 max-T statistic values
  x<-replicate(1000,fun1(d1,d2,d3,d4,s1,s2,s3))
  # arrange them in increasing order
  y<-sort(x,decreasing=FALSE)
  # find the critical value at level 0.05
  c<-y[950]
  return(c)
}
fun3<-function(n11,n12,n13,n14,n21,n22,n23,n24,n31,n32,n33,n34,n41,n42,n43,n44,n51,n52,n53,n54,
               mu11,mu12,mu13,mu14,mu21,mu22,mu23,mu24,mu31,mu32,mu33,mu34,mu41,mu42,mu43,mu44,
               mu51,mu52,mu53,mu54,
               v11,v12,v13,v14,v21,v22,v23,v24,v31,v32,v33,v34,v41,v42,v43,v44,v51,v52,v53,v54)
{
  ###For generating samples from normal distributions
  
  g11<-rnorm(n11,mu11,sqrt(v11));g12<-rnorm(n12,mu12,sqrt(v12))
  g13<-rnorm(n13,mu13,sqrt(v13));g14<-rnorm(n14,mu14,sqrt(v14))
  g21<-rnorm(n21,mu21,sqrt(v21));g22<-rnorm(n22,mu22,sqrt(v22))
  g23<-rnorm(n23,mu23,sqrt(v23));g24<-rnorm(n24,mu24,sqrt(v24))
  g31<-rnorm(n31,mu31,sqrt(v31));g32<-rnorm(n32,mu32,sqrt(v32))
  g33<-rnorm(n33,mu33,sqrt(v33));g34<-rnorm(n34,mu34,sqrt(v34))
  g41<-rnorm(n41,mu41,sqrt(v41));g42<-rnorm(n42,mu42,sqrt(v42))
  g43<-rnorm(n43,mu43,sqrt(v43));g44<-rnorm(n44,mu44,sqrt(v44))
  g51<-rnorm(n51,mu51,sqrt(v51));g52<-rnorm(n52,mu52,sqrt(v52))
  g53<-rnorm(n53,mu53,sqrt(v53));g54<-rnorm(n54,mu54,sqrt(v54))
  
  ##For generating samples from $t$ distributions
  
  #Y11=rt(n11,3,ncp=0)/sqrt(3);Y12=rt(n12,3,ncp=0)/sqrt(3);Y13=rt(n13,3,ncp=0)/sqrt(3)
  #Y14=rt(n14,3,ncp=0)/sqrt(3)
  #Y21=rt(n21,3,ncp=0)/sqrt(3);Y22=rt(n22,3,ncp=0)/sqrt(3);Y23=rt(n23,3,ncp=0)/sqrt(3)
  #Y24=rt(n24,3,ncp=0)/sqrt(3)
  #Y31=rt(n31,3,ncp=0)/sqrt(3);Y32=rt(n32,3,ncp=0)/sqrt(3);Y33=rt(n33,3,ncp=0)/sqrt(3)
  #Y34=rt(n34,3,ncp=0)/sqrt(3)
  #Y41=rt(n41,3,ncp=0)/sqrt(3);Y42=rt(n42,3,ncp=0)/sqrt(3);Y43=rt(n43,3,ncp=0)/sqrt(3)
  #Y44=rt(n44,3,ncp=0)/sqrt(3)
  #Y51=rt(n51,3,ncp=0)/sqrt(3);Y52=rt(n52,3,ncp=0)/sqrt(3);Y53=rt(n53,3,ncp=0)/sqrt(3)
  #Y54=rt(n54,3,ncp=0)/sqrt(3)
  #g11=sqrt(v11)*Y11+mu11;g12=sqrt(v12)*Y12+mu12;g13=sqrt(v13)*Y13+mu13;g14=sqrt(v14)*Y14+mu14
  #g21=sqrt(v21)*Y21+mu21;g22=sqrt(v22)*Y22+mu22;g23=sqrt(v23)*Y23+mu23;g24=sqrt(v24)*Y24+mu24
  #g31=sqrt(v31)*Y31+mu31;g32=sqrt(v32)*Y32+mu32;g33=sqrt(v33)*Y33+mu33;g34=sqrt(v34)*Y34+mu34
  #g41=sqrt(v41)*Y41+mu41;g42=sqrt(v42)*Y42+mu42;g43=sqrt(v43)*Y43+mu43;g44=sqrt(v44)*Y44+mu44
  #g51=sqrt(v51)*Y51+mu51;g52=sqrt(v52)*Y52+mu52;g53=sqrt(v53)*Y53+mu53;g54=sqrt(v54)*Y54+mu54
  
  ##For generating samples from exponential distributions
  
  #Y11=rexp(n11,1);Y12=rexp(n12,1);Y13=rexp(n13,1);Y14=rexp(n14,1)
  #Y21=rexp(n21,1);Y22=rexp(n22,1);Y23=rexp(n23,1);Y24=rexp(n24,1)
  #Y31=rexp(n31,1);Y32=rexp(n32,1);Y33=rexp(n33,1);Y34=rexp(n34,1)
  #Y41=rexp(n41,1);Y42=rexp(n42,1);Y43=rexp(n43,1);Y44=rexp(n44,1)
  #Y51=rexp(n51,1);Y52=rexp(n52,1);Y53=rexp(n53,1);Y54=rexp(n54,1)
  #g11=mu11+sqrt(v11)*(Y11-1);g12=mu12+sqrt(v12)*(Y12-1);g13=mu13+sqrt(v13)*(Y13-1);g14=mu14+sqrt(v14)*(Y14-1)
  #g21=mu21+sqrt(v21)*(Y21-1);g22=mu22+sqrt(v22)*(Y22-1);g23=mu23+sqrt(v23)*(Y23-1);g24=mu24+sqrt(v24)*(Y24-1)
  #g31=mu31+sqrt(v31)*(Y31-1);g32=mu32+sqrt(v32)*(Y32-1);g33=mu33+sqrt(v33)*(Y33-1);g34=mu34+sqrt(v34)*(Y34-1)
  #g41=mu41+sqrt(v41)*(Y41-1);g42=mu42+sqrt(v42)*(Y42-1);g43=mu43+sqrt(v43)*(Y43-1);g44=mu44+sqrt(v44)*(Y44-1)
  #g51=mu51+sqrt(v51)*(Y51-1);g52=mu52+sqrt(v52)*(Y52-1);g53=mu53+sqrt(v53)*(Y53-1);g54=mu54+sqrt(v54)*(Y54-1)
  
  ##For generating samples from Gamma distributions
  
  #Y11=(rgamma(n11,2,1)-2)/sqrt(2);Y12=(rgamma(n12,2,1)-2)/sqrt(2);Y13=(rgamma(n13,2,1)-2)/sqrt(2);Y14=(rgamma(n14,2,1)-2)/sqrt(2)
  #Y21=(rgamma(n21,2,1)-2)/sqrt(2);Y22=(rgamma(n22,2,1)-2)/sqrt(2);Y23=(rgamma(n23,2,1)-2)/sqrt(2);Y24=(rgamma(n24,2,1)-2)/sqrt(2)
  #Y31=(rgamma(n31,2,1)-2)/sqrt(2);Y32=(rgamma(n32,2,1)-2)/sqrt(2);Y33=(rgamma(n33,2,1)-2)/sqrt(2);Y34=(rgamma(n34,2,1)-2)/sqrt(2)
  #Y41=(rgamma(n41,2,1)-2)/sqrt(2);Y42=(rgamma(n42,2,1)-2)/sqrt(2);Y43=(rgamma(n43,2,1)-2)/sqrt(2);Y44=(rgamma(n44,2,1)-2)/sqrt(2)
  #Y51=(rgamma(n51,2,1)-2)/sqrt(2);Y52=(rgamma(n52,2,1)-2)/sqrt(2);Y53=(rgamma(n53,2,1)-2)/sqrt(2);Y54=(rgamma(n54,2,1)-2)/sqrt(2)
  #g11=mu11+sqrt(v11)*Y11;g12=mu12+sqrt(v12)*Y12;g13=mu13+sqrt(v13)*Y13;g14=mu14+sqrt(v14)*Y14
  #g21=mu21+sqrt(v21)*Y21;g22=mu22+sqrt(v22)*Y22;g23=mu23+sqrt(v23)*Y23;g24=mu24+sqrt(v24)*Y24
  #g31=mu31+sqrt(v31)*Y31;g32=mu32+sqrt(v32)*Y32;g33=mu33+sqrt(v33)*Y33;g34=mu34+sqrt(v34)*Y34
  #g41=mu41+sqrt(v41)*Y41;g42=mu42+sqrt(v42)*Y42;g43=mu43+sqrt(v43)*Y43;g44=mu44+sqrt(v44)*Y44
  #g51=mu51+sqrt(v51)*Y51;g52=mu52+sqrt(v52)*Y52;g53=mu53+sqrt(v53)*Y53;g54=mu54+sqrt(v54)*Y54
  
  ##For generating samples from lognormal distributions
  
  #Y11=(rlnorm(n11,0,1)-exp(.5))/sqrt(exp(2)-exp(1));Y12=(rlnorm(n12,0,1)-exp(.5))/sqrt(exp(2)-exp(1))
  #Y13=(rlnorm(n13,0,1)-exp(.5))/sqrt(exp(2)-exp(1));Y14=(rlnorm(n14,0,1)-exp(.5))/sqrt(exp(2)-exp(1))
  #Y21=(rlnorm(n21,0,1)-exp(.5))/sqrt(exp(2)-exp(1));Y22=(rlnorm(n22,0,1)-exp(.5))/sqrt(exp(2)-exp(1))
  #Y23=(rlnorm(n23,0,1)-exp(.5))/sqrt(exp(2)-exp(1));Y24=(rlnorm(n24,0,1)-exp(.5))/sqrt(exp(2)-exp(1))
  #Y31=(rlnorm(n31,0,1)-exp(.5))/sqrt(exp(2)-exp(1));Y32=(rlnorm(n32,0,1)-exp(.5))/sqrt(exp(2)-exp(1))
  #Y33=(rlnorm(n33,0,1)-exp(.5))/sqrt(exp(2)-exp(1));Y34=(rlnorm(n14,0,1)-exp(.5))/sqrt(exp(2)-exp(1))
  #Y41=(rlnorm(n41,0,1)-exp(.5))/sqrt(exp(2)-exp(1));Y42=(rlnorm(n42,0,1)-exp(.5))/sqrt(exp(2)-exp(1))
  #Y43=(rlnorm(n43,0,1)-exp(.5))/sqrt(exp(2)-exp(1));Y44=(rlnorm(n44,0,1)-exp(.5))/sqrt(exp(2)-exp(1))
  #Y51=(rlnorm(n51,0,1)-exp(.5))/sqrt(exp(2)-exp(1));Y52=(rlnorm(n12,0,1)-exp(.5))/sqrt(exp(2)-exp(1))
  #Y53=(rlnorm(n53,0,1)-exp(.5))/sqrt(exp(2)-exp(1));Y54=(rlnorm(n14,0,1)-exp(.5))/sqrt(exp(2)-exp(1))
  #g11=mu11+sqrt(v11)*Y11;g12=mu12+sqrt(v12)*Y12;g13=mu13+sqrt(v13)*Y13;g14=mu14+sqrt(v14)*Y14
  #g21=mu21+sqrt(v21)*Y21;g22=mu22+sqrt(v22)*Y22;g23=mu23+sqrt(v23)*Y23;g24=mu24+sqrt(v24)*Y24
  #g31=mu31+sqrt(v31)*Y31;g32=mu32+sqrt(v32)*Y32;g33=mu33+sqrt(v33)*Y33;g34=mu34+sqrt(v34)*Y34
  #g41=mu41+sqrt(v41)*Y41;g42=mu42+sqrt(v42)*Y42;g43=mu43+sqrt(v43)*Y43;g44=mu44+sqrt(v44)*Y44
  #g51=mu51+sqrt(v51)*Y51;g52=mu52+sqrt(v52)*Y52;g53=mu53+sqrt(v53)*Y53;g54=mu54+sqrt(v54)*Y54
  
  ##For generating samples from Weibull distributions
  
  #Y11=(rweibull(n11,2,1)-sqrt(pi)/2)/sqrt(1-pi/4);Y12=(rweibull(n12,2,1)-sqrt(pi)/2)/sqrt(1-pi/4)
  #Y13=(rweibull(n13,2,1)-sqrt(pi)/2)/sqrt(1-pi/4);Y14=(rweibull(n14,2,1)-sqrt(pi)/2)/sqrt(1-pi/4)
  #Y21=(rweibull(n21,2,1)-sqrt(pi)/2)/sqrt(1-pi/4);Y22=(rweibull(n22,2,1)-sqrt(pi)/2)/sqrt(1-pi/4)
  #Y23=(rweibull(n23,2,1)-sqrt(pi)/2)/sqrt(1-pi/4);Y24=(rweibull(n24,2,1)-sqrt(pi)/2)/sqrt(1-pi/4)
  #Y31=(rweibull(n31,2,1)-sqrt(pi)/2)/sqrt(1-pi/4);Y32=(rweibull(n32,2,1)-sqrt(pi)/2)/sqrt(1-pi/4)
  #Y33=(rweibull(n33,2,1)-sqrt(pi)/2)/sqrt(1-pi/4);Y34=(rweibull(n34,2,1)-sqrt(pi)/2)/sqrt(1-pi/4)
  #Y41=(rweibull(n41,2,1)-sqrt(pi)/2)/sqrt(1-pi/4);Y42=(rweibull(n42,2,1)-sqrt(pi)/2)/sqrt(1-pi/4)
  #Y43=(rweibull(n43,2,1)-sqrt(pi)/2)/sqrt(1-pi/4);Y44=(rweibull(n44,2,1)-sqrt(pi)/2)/sqrt(1-pi/4)
  #Y51=(rweibull(n51,2,1)-sqrt(pi)/2)/sqrt(1-pi/4);Y52=(rweibull(n52,2,1)-sqrt(pi)/2)/sqrt(1-pi/4)
  #Y53=(rweibull(n53,2,1)-sqrt(pi)/2)/sqrt(1-pi/4);Y54=(rweibull(n34,2,1)-sqrt(pi)/2)/sqrt(1-pi/4)
  #g11=mu11+sqrt(v11)*Y11;g12=mu12+sqrt(v12)*Y12;g13=mu13+sqrt(v13)*Y13;g14=mu14+sqrt(v14)*Y14
  #g21=mu21+sqrt(v21)*Y21;g22=mu22+sqrt(v22)*Y22;g23=mu23+sqrt(v23)*Y23;g24=mu24+sqrt(v24)*Y24
  #g31=mu31+sqrt(v31)*Y31;g32=mu32+sqrt(v32)*Y32;g33=mu33+sqrt(v33)*Y33;g34=mu34+sqrt(v34)*Y34
  #g41=mu41+sqrt(v41)*Y41;g42=mu42+sqrt(v42)*Y42;g43=mu43+sqrt(v43)*Y43;g44=mu44+sqrt(v44)*Y44
  #g51=mu51+sqrt(v51)*Y51;g52=mu52+sqrt(v52)*Y52;g53=mu53+sqrt(v53)*Y53;g54=mu54+sqrt(v54)*Y54
  s11=var(g11);s12=var(g12);s13=var(g13);s14=var(g14)
  s21=var(g21);s22=var(g22);s23=var(g23);s24=var(g24)
  s31=var(g31);s32=var(g32);s33=var(g33);s34=var(g34)
  s41=var(g41);s42=var(g42);s43=var(g43);s44=var(g44)
  s51=var(g51);s52=var(g52);s53=var(g53);s54=var(g54)
  S1=(s11/n11+s12/n12+s13/n13+s14/n14)/16;S2=(s21/n21+s22/n22+s23/n23+s24/n24)/16
  S3=(s31/n31+s32/n32+s33/n33+s34/n34)/16;S4=(s41/n41+s42/n42+s43/n43+s44/n44)/16
  S5=(s51/n51+s52/n52+s53/n53+s54/n54)/16
  X1=(mean(g11)+mean(g12)+mean(g13)+mean(g14))/4;X2=(mean(g21)+mean(g22)+mean(g23)+mean(g24))/4
  X3=(mean(g31)+mean(g32)+mean(g33)+mean(g24))/4;X4=(mean(g41)+mean(g42)+mean(g43)+mean(g44))/4
  X5=(mean(g51)+mean(g52)+mean(g53)+mean(g54))/4
  v1<-sqrt(S1+S2);v2<-sqrt(S2+S3);v3<-sqrt(S3+S4);v4<-sqrt(S4+S5)
  T1=(X2-X1)/v1;T2=(X3-X2)/v2;T3=(X4-X3)/v3;T4=(X5-X4)/v4
  # max-T value
  A=max(T1,T2,T3,T4,na.rm = FALSE)
  N=n11+n12+n13+n14+n21+n22+n23+n24+n31+n32+n33+n34+n41+n42+n43+n44+n51+n52+n53+n54
  nu11=n11/N;nu12=n12/N;nu13=n13/N;nu14=n14/N
  nu21=n21/N;nu22=n22/N;nu23=n23/N;nu24=n24/N
  nu31=n31/N;nu32=n32/N;nu33=n33/N;nu34=n34/N
  nu41=n41/N;nu42=n42/N;nu43=n43/N;nu44=n44/N
  nu51=n51/N;nu52=n52/N;nu53=n53/N;nu54=n54/N
  d1=(sqrt(s21/nu21+s22/nu22+s23/nu23+s24/nu24+s11/nu11+s12/nu12+s13/nu13+s14/nu14))/4
  d2=(sqrt(s31/nu31+s32/nu32+s33/nu33+s34/nu34+s21/nu21+s22/nu22+s23/nu23+s24/nu24))/4
  d3=(sqrt(s41/nu41+s42/nu42+s43/nu43+s44/nu44+s31/nu31+s32/nu32+s33/nu33+s34/nu34))/4
  d4=(sqrt(s51/nu51+s52/nu52+s53/nu53+s54/nu54+s41/nu41+s42/nu42+s43/nu43+s44/nu44))/4
  s1=(s21/nu21+s22/nu22+s23/nu23+s24/nu24)/16;s2=(s31/nu31+s32/nu32+s33/nu33+s34/nu34)/16
  s3=(s41/nu41+s42/nu42+s43/nu43+s44/nu44)/16
  out=fun2(d1,d2,d3,d4,s1,s2,s3)
  # count if the statistic value is greater than the crtical value
  a=0
  if(A>out)
    a=a+1
  return(a)
}
fun4<-function(n11,n12,n13,n14,n21,n22,n23,n24,n31,n32,n33,n34,n41,n42,n43,n44,n51,n52,n53,n54,
               mu11,mu12,mu13,mu14,mu21,mu22,mu23,mu24,mu31,mu32,mu33,mu34,mu41,mu42,mu43,mu44,
               mu51,mu52,mu53,mu54,
               v11,v12,v13,v14,v21,v22,v23,v24,v31,v32,v33,v34,v41,v42,v43,v44,v51,v52,v53,v54)
{
  # find number of times maxt-T value > critical value among 1000 values
  out<-replicate(1000,fun3(n11,n12,n13,n14,n21,n22,n23,n24,n31,n32,n33,n34,n41,n42,n43,n44,n51,n52,n53,n54,
                           mu11,mu12,mu13,mu14,mu21,mu22,mu23,mu24,mu31,mu32,mu33,mu34,mu41,mu42,mu43,mu44,
                           mu51,mu52,mu53,mu54,
                           v11,v12,v13,v14,v21,v22,v23,v24,v31,v32,v33,v34,v41,v42,v43,v44,v51,v52,v53,v54))
  p<-sum(out)/1000
  return(p)
}
X=replicate(4,fun4(10,10,10,10,8,8,8,8,12,12,12,12,14,14,14,14,15,15,15,15,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1,0.2,0.3,0.4,
                   0.1,0.2,0.3,0.4,0.1,0.2,0.3,0.4,0.1,0.2,0.3,0.4,0.1,0.2,0.3,0.4))
mean(X)
X=replicate(4,fun4(10,10,10,10,8,8,8,8,12,12,12,12,14,14,14,14,15,15,15,15,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1,0.1,0.1,0.1,
                   0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1))
mean(X)
X=replicate(4,fun4(10,10,10,10,8,8,8,8,12,12,12,12,14,14,14,14,15,15,15,15,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1,0.2,0.3,0.4,
                   0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2))
mean(X)
X=replicate(4,fun4(10,10,10,10,8,8,8,8,12,12,12,12,14,14,14,14,15,15,15,15,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.2,0.2,0.2,0.2,
                   0.6,0.6,0.6,0.6,0.8,0.8,0.8,0.8,2,2,2,2,1.5,1.5,1.5,1.5))
mean(X)

X=replicate(4,fun4(10,15,18,20,10,15,18,20,10,15,18,20,10,15,18,20,10,15,18,20,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1,0.2,0.3,0.4,
                   0.1,0.2,0.3,0.4,0.1,0.2,0.3,0.4,0.1,0.2,0.3,0.4,0.1,0.2,0.3,0.4))
mean(X)
X=replicate(4,fun4(10,15,18,20,10,15,18,20,10,15,18,20,10,15,18,20,10,15,18,20,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1,0.1,0.1,0.1,
                   0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1))
mean(X)
X=replicate(4,fun4(10,15,18,20,10,15,18,20,10,15,18,20,10,15,18,20,10,15,18,20,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1,0.2,0.3,0.4,
                   0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2))
mean(X)
X=replicate(4,fun4(10,15,18,20,10,15,18,20,10,15,18,20,10,15,18,20,10,15,18,20,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.2,0.2,0.2,0.2,
                   0.6,0.6,0.6,0.6,0.8,0.8,0.8,0.8,2,2,2,2,1.5,1.5,1.5,1.5))
mean(X)

X=replicate(4,fun4(10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1,0.2,0.3,0.4,
                   0.1,0.2,0.3,0.4,0.1,0.2,0.3,0.4,0.1,0.2,0.3,0.4,0.1,0.2,0.3,0.4))
mean(X)
X=replicate(4,fun4(10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1,0.1,0.1,0.1,
                   0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1))
mean(X)
X=replicate(4,fun4(10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1,0.2,0.3,0.4,
                   0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2))
mean(X)
X=replicate(4,fun4(10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.2,0.2,0.2,0.2,
                   0.6,0.6,0.6,0.6,0.8,0.8,0.8,0.8,2,2,2,2,1.5,1.5,1.5,1.5))
mean(X)

X=replicate(4,fun4(10,8,9,15,16,17,10,16,8,10,12,14,15,17,8,9,10,12,13,10,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1,0.2,0.3,0.4,
                   0.1,0.2,0.3,0.4,0.1,0.2,0.3,0.4,0.1,0.2,0.3,0.4,0.1,0.2,0.3,0.4))
mean(X)
X=replicate(4,fun4(10,8,9,15,16,17,10,16,8,10,12,14,15,17,8,9,10,12,13,10,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1,0.1,0.1,0.1,
                   0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1))
mean(X)
X=replicate(4,fun4(10,8,9,15,16,17,10,16,8,10,12,14,15,17,8,9,10,12,13,10,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1,0.2,0.3,0.4,
                   0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2))
mean(X)
X=replicate(4,fun4(10,8,9,15,16,17,10,16,8,10,12,14,15,17,8,9,10,12,13,10,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.2,0.2,0.2,0.2,
                   0.6,0.6,0.6,0.6,0.8,0.8,0.8,0.8,2,2,2,2,1.5,1.5,1.5,1.5))
mean(X)

