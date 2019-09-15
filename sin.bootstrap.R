#' Estimate the confidence interval and predict the arrvial and departure number
#'
#' This r code is to give the mean's confidence interval of the total
#' arrived number when the parameters is known.
#'
#' @param u the arrival rate funcion and the service time PDF's parameters
#' @param T the totol time to work
#' @param flag the r file name which is to solve the problem about the specific sevice time PDF
#' @param t1 the beginning of the prediction
#' @param t2 the end of the prediction
#'
#' @return the confidence interval of the total(or the prediction number from t1 to t2)
#'         arrival number and departure number
#' @export
#'
#' @examples
#' m<-sin.bootstrap(u=c(10,5,24,2),T=24,flag='sin.exp',t1=25,t2=26)

sin.bootstrap<-function(c,T,flag,t1,t2){
  M<-matrix(0,1000,T)
  Nn<-matrix(0,1000,T)
  mMa<-rep(0,T)
  mMd<-rep(0,T)
  predM<-rep(0,1000)
  predNn<-rep(0,1000)
  flag1<-substitute(flag)
  if(flag == 'sin.lognorm'){

    sinlognorm<-function(k,y,miu,sigm){
      exp(-(log(k-y)-miu)^2/(2*sigm^2))/((k-y)*sigm*(2*pi)^0.5)
    }
    flag1<-sinlognorm
    flag2<-sin.lognorm

  }
  else if(flag == 'sin.exp'){

    sinexp<-function(k,y,miu){
      miu*exp(-miu*(k-y))
    }
    flag1<-sinexp
    flag2<-sin.exp

  }
  Ma<-function(c1,c2,c3,t){
    c1*t+c2*c3/(2*pi)*(1-cos(2*pi*t/c3))
  }
  lamda<-c[1]
  A<-c[2]
  T0<-c[3]
  miu<-c[4]

  for (j in 1:1000) {

    N=rpois(1,Ma(lamda,A,T0,T))
    ta=rep(0,N)
    for (i in 1:N){
      u=runif(1,0,1)
      ta[i]=uniroot(function(t) lamda*t+A*T0/(2*pi)*(1-cos(2*pi*t/T0))-u*Ma(lamda,A,T0,T),c(0,T))$root
    }

    ser=rexp(N,2)
    td=ta+ser

    na=rep(0,T)
    nd=rep(0,T)
    for (i in 1:T){
      na[i]= length(ta[ta<=i])
    }

    for (i in 1:T){
      nd[i]= length(td[td<=i])
    }

    x1<-flag2(na,nd,c)$'par'

    for (k in 1:T) {
      M[j,k]<-Ma(x1[1],x1[2],x1[3],k)
      mMa[k]<-Ma(c[1],c[2],c[3],k)
      if(flag == 'sin.lognorm'){
        Md<-function(a,b,c,y,miu,sigm){
          integrate(function(x) Ma(a,b,c,x)*flag1(y,x,miu,sigm),0,y)$value
        }

        Nn[j,k]<-Md(x1[1],x1[2],x1[3],k,x1[4],x1[5])
        mMd[k]<-Md(c[1],c[2],c[3],k,c[4],c[5])
      }
      else if(flag == 'sin.exp'){
        Md<-function(a,b,c,y,miu){
          integrate(function(x) Ma(a,b,c,x)*flag1(y,x,miu),0,y)$value
        }
        Nn[j,k]<-Md(x1[1],x1[2],x1[3],k,x1[4])
        mMd[k]<-Md(c[1],c[2],c[3],k,c[4])
      }
    }
      predM[j]<-Ma(x1[1],x1[2],x1[3],t2)-Ma(x1[1],x1[2],x1[3],t1)
      if(flag == 'sin.lognorm'){
        Md<-function(a,b,c,y,miu,sigm){
          integrate(function(x) Ma(a,b,c,x)*flag1(y,x,miu,sigm),0,y)$value
        }
        predNn[j]<-Md(x1[1],x1[2],x1[3],t2,x1[4],x1[5])-Md(x1[1],x1[2],x1[3],t1,x1[4],x1[5])
      }
      else if(flag == 'sin.exp'){
        Md<-function(a,b,c,y,miu){
          integrate(function(x) Ma(a,b,c,x)*flag1(y,x,miu),0,y)$value
        }
        predNn[j]<-Md(x1[1],x1[2],x1[3],t2,x1[4])-Md(x1[1],x1[2],x1[3],t1,x1[4])
      }

    }
    z1<-matrix(0,T,2)
    z2<-matrix(0,T,2)
    for (k in 1:T) {
      z1[k,]<-c(quantile(M[,k],0.025),quantile(M[,k],0.975))
      z2[k,]<-c(quantile(Nn[,k],0.025),quantile(Nn[,k],0.975))
    }
    z3<-c(quantile(predM,0.025),quantile(predM,0.975))
    z4<-c(quantile(predNn,0.025),quantile(predNn,0.975))
    list_data<-list(z1,z2,z3,z4)
    names(list_data)<-c('NA.CI','ND.CI','NA.predict','ND.predict')
    print(list_data)
    tt<-1:T
    par(mfrow=c(1,2))
    plot(tt,z1[,2],col='red',xlab='',ylab='',type='l',xlim=c(0,T),ylim=c(0,max(z1[,2])))
    par(new=TRUE)
    plot(tt,z1[,1],col='blue',xlab='',ylab='',type='l',xlim=c(0,T),ylim=c(0,max(z1[,2])))
    par(new=TRUE)
    plot(tt,mMa,col='green',axes=FALSE,xlab='time(h)',ylab='arrival',type='l',xlim=c(0,T),ylim=c(0,max(z1[,2])))
    plot(tt,z2[,2],col='red',xlab='',ylab='',xlim=c(0,T),ylim=c(0,max(z2[,2])),type='l')
    par(new=TRUE)
    plot(tt,z2[,1],col='blue',xlab='',ylab='',type='l',xlim=c(0,T),ylim=c(0,max(z2[,2])))
    par(new=TRUE)
    plot(tt,mMd,col='green',axes=FALSE,xlab='time(h)',ylab='departure',type='l',xlim=c(0,T),ylim=c(0,max(z2[,2])))
  }

#' Estimate the parameters
#'
#' Use MLE method to estimate the parameters when the arrival rate function is
#' cyclic arrival and the service time's PDF is exponential distribution.And also uses funcion 'nlminb'
#' which is to calculate the minimum point of the ML function.
#'
#' @param na the number of arrival
#' @param nd the number of departure
#' @param c the starting value of the parameters which the 'nlminb' need
#'
#' @return the estimation of the parameters
#' @export
#'
#' @examples
#' m<-sin.exp(na=c(2,2,2,3),nd=c(2,2,2,2),c=c(100,10,1))
sin.exp<-function(na,nd,c) {
  lsine <- function(x) {
    n=length(nd)

    da=na[2:n]-na[1:n-1]
    da=c(na[1],da)
    dd=nd[2:n]-nd[1:n-1]
    dd=c(nd[1],dd)

    w=pmin(da,dd)

    e=pmax(0,nd[2:n]-na[1:n-1])
    e=c(nd[1],e)

    t=1:n
    t1=t
    t1[na!=nd]=0
    t0=rep(0,n)
    for(i in 2:n){
      t0[i]=max(t1[1:i-1])
    }
    lamda=x[1]
    A=x[2]
    T0=x[3]
    miu=x[4]

    msin=function(lamda,A,T0,y){
      lamda*y+A*T0/(2*pi)*(1-cos(2*pi*y/T0))
    }

    p1=c(1,rep(-1,n-1))
    p=rep(0,n)
    q=rep(0,n)

    int3=integrate(function(y) miu*exp(-miu*(t[1]-y))*msin(lamda,A,T0,y),0,t[1])$value
    p2=rep((pexp(0,miu)*msin(lamda,A,T0,t[1])+int3)/msin(lamda,A,T0,t[1]),n)
    for(i in 2:n) {
      int4=integrate(function(y) miu*exp(-miu*(t[i]-y))*msin(lamda,A,T0,y),t[i-1],t[i])$value
      p2[i]=(pexp(0,miu)*msin(lamda,A,T0,t[i])-pexp(t[i]-t[i-1],miu)*msin(lamda,A,T0,t[i-1])+int4)/(msin(lamda,A,T0,t[i])-msin(lamda,A,T0,t[i-1]))
      if (t[i-1]!=t0[i]){
        int1=integrate(function(y) miu*exp(-miu*(t[i-1]-y))*msin(lamda,A,T0,y),t0[i],t[i-1])$value
        p[i]=(pexp(0,miu)*msin(lamda,A,T0,t[i-1])-pexp(t[i-1]-t0[i],miu)*msin(lamda,A,T0,t0[i])+int1)/(msin(lamda,A,T0,t[i-1])-msin(lamda,A,T0,t0[i]))

        int2=integrate(function(y) miu*exp(-miu*(t[i]-y))*msin(lamda,A,T0,y),t0[i],t[i-1])$value
        q[i]=(pexp(t[i]-t[i-1],miu)*msin(lamda,A,T0,t[i-1])-pexp(t[i]-t0[i],miu)*msin(lamda,A,T0,t0[i])+int2)/(msin(lamda,A,T0,t[i-1])-msin(lamda,A,T0,t0[i]))-p[i]
        p1[i]=q[i]/(1-p[i])
      }
    }

    sum1=0;
    lsum=0;
    for(i in 2:n) {
      sum1=sum1+da[i]*log(msin(lamda,A,T0,t[i])-msin(lamda,A,T0,t[i-1]))-(msin(lamda,A,T0,t[i])-msin(lamda,A,T0,t[i-1]))
      sum=0;
      if (t[i-1]==t0[i]){
        sum=dbinom(dd[i],da[i],p2[i])
      }
      else{
        for (j in e[i]:w[i]) {
          sum=sum+dbinom(dd[i]-j,na[i-1]-nd[i-1],p1[i])*dbinom(j,da[i],p2[i])
        }
      }
      lsum=lsum+log(sum)
    }
    -sum1-lsum-(da[1]*log(msin(lamda,A,T0,t[1]))-msin(lamda,A,T0,t[1]))-log(dbinom(dd[1],da[1],p2[1]))

  }
  z<-nlminb(c,lsine)
  return(z)
}

#' Estimate the parameters
#'
#' Use MLE method to estimate the parameters when the arrival rate function is
#' cyclic arrival and the service time's PDF is lognormal distribution.And also uses funcion 'nlminb'
#' which is to calculate the minimum point of the ML function.
#'
#' @param na the number of arrival
#' @param nd the number of departure
#' @param c the starting value of the parameters which the 'nlminb' need
#'
#' @return the estimation of the parameters
#' @export
#'
#' @examples
#' m<-sin.lognorm(na=c(2,2,2,3,1),nd=c(2,2,2,2,1),c=c(100,10,1,1,1))
sin.lognorm<-function(na,nd,c) {
  lsine <- function(x) {
    n=length(nd)

    da=na[2:n]-na[1:n-1]
    da=c(na[1],da)
    dd=nd[2:n]-nd[1:n-1]
    dd=c(nd[1],dd)

    w=pmin(da,dd)

    e=pmax(0,nd[2:n]-na[1:n-1])
    e=c(nd[1],e)

    t=1:n
    t1=t
    t1[na!=nd]=0
    t0=rep(0,n)
    for(i in 2:n){
      t0[i]=max(t1[1:i-1])
    }
    lamda=x[1]
    A=x[2]
    T0=x[3]
    miu=x[4]
    sigm=x[5]

    msin=function(y){
      lamda*y+A*T0/(2*pi)*(1-cos(2*pi*y/T0))
    }

    p1=c(1,rep(-1,n-1))
    p=rep(0,n)
    q=rep(0,n)

    int3=integrate(function(y) exp(-(log(t[1]-y)-miu)^2/(2*sigm^2))/((t[1]-y)*sigm*(2*pi)^0.5)*msin(y),0,t[1])$value
    temp<-try(int3)
    if('try-error' %in% class(temp)){
      print('The initial parameters need to be modified so that it can be calculated.')
      break
    }
    p2=rep((plnorm(0,miu,sigm)*msin(t[1])+int3)/msin(t[1]),n)
    for(i in 2:n) {
      int4=integrate(function(y) exp(-(log(t[i]-y)-miu)^2/(2*sigm^2))/((t[i]-y)*sigm*(2*pi)^0.5)*msin(y),t[i-1],t[i])$value
      p2[i]=(plnorm(0,miu,sigm)*msin(t[i])-plnorm(t[i]-t[i-1],miu,sigm)*msin(t[i-1])+int4)/(msin(t[i])-msin(t[i-1]))
      if (t[i-1]!=t0[i]){
        int1=integrate(function(y) exp(-(log(t[i-1]-y)-miu)^2/(2*sigm^2))/((t[i-1]-y)*sigm*(2*pi)^0.5)*msin(y),t0[i],t[i-1])$value
        p[i]=(plnorm(0,miu,sigm)*msin(t[i-1])-plnorm(t[i-1]-t0[i],miu,sigm)*msin(t0[i])+int1)/(msin(t[i-1])-msin(t0[i]))

        int2=integrate(function(y) exp(-(log(t[i]-y)-miu)^2/(2*sigm^2))/((t[i]-y)*sigm*(2*pi)^0.5)*msin(y),t0[i],t[i-1])$value
        q[i]=(plnorm(t[i]-t[i-1],miu,sigm)*msin(t[i-1])-plnorm(t[i]-t0[i],miu,sigm)*msin(t0[i])+int2)/(msin(T0,t[i-1])-msin(t0[i]))-p[i]
        p1[i]=q[i]/(1-p[i])
      }
    }

    sum1=0;
    lsum=0;
    for(i in 2:n) {
      sum1=sum1+da[i]*log(msin(t[i])-msin(t[i-1]))-(msin(t[i])-msin(t[i-1]))
      sum=0;
      if (t[i-1]==t0[i]){
        sum=dbinom(dd[i],da[i],p2[i])
      }
      else{
        for (j in e[i]:w[i]) {
          sum=sum+dbinom(dd[i]-j,na[i-1]-nd[i-1],p1[i])*dbinom(j,da[i],p2[i])
        }
      }
      lsum=lsum+log(sum)
    }
    -sum1-lsum-(da[1]*log(msin(t[1]))-msin(t[1]))-log(dbinom(dd[1],da[1],p2[1]))

  }
  z<-nlminb(c,lsine)
  return(z)
}
