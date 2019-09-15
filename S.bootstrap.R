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
#' m<-S.bootstrap(u=c(4720,0.1,190,1,0.1),T=86,flag='S.exp',t1=87,t2=88)

S.bootstrap<-function(u,T,flag,t1,t2){
  M<-matrix(0,1000,T)
  N<-matrix(0,1000,T)
  mMa<-rep(0,T)
  mMd<-rep(0,T)
  predM<-rep(0,1000)
  predN<-rep(0,1000)
  flag1<-substitute(flag)
  if(flag == 'S.lognorm'){

    slognorm<-function(k,y,miu,sigm){
      exp(-(log(k-y)-miu)^2/(2*sigm^2))/((k-y)*sigm*(2*pi)^0.5)
    }
    flag1<-slognorm
    flag2<-S.lognorm
  }
  else if(flag == 'S.exp'){

    sexp<-function(k,y,miu){
      miu*exp(-miu*(k-y))
    }
    flag1<-sexp
    flag2<-S.exp
  }
  Ma<-function(a,b,c,y){
    a*(1-exp(-b*y))/(1+c*exp(-b*y))
  }

  for (j in 1:1000) {


    m<-S.nand(u,T)
    z<-length(m)*0.5

    na<-m[1:z]

    nd<-m[(z+1):(2*z)]

    x1<-flag2(na,nd,u)$'par'


    for (k in 1:T) {
      M[j,k]<-Ma(x1[1],x1[2],x1[3],k)
      mMa[k]<-Ma(u[1],u[2],u[3],k)
      if(flag == 'S.lognorm'){
        Md<-function(a,b,c,y,miu,sigm){
          integrate(function(x) Ma(a,b,c,x)*flag1(y,x,miu,sigm),0,y)$value
        }
        N[j,k]<-Md(x1[1],x1[2],x1[3],k,x1[4],x1[5])
        mMd[k]<-Md(u[1],u[2],u[3],k,u[4],u[5])
      }
      else if(flag == 'S.exp'){
          Md<-function(a,b,c,y,miu){
            integrate(function(x) Ma(a,b,c,x)*flag1(y,x,miu),0,y)$value
          }
        N[j,k]<-Md(x1[1],x1[2],x1[3],k,x1[4])
        mMd[k]<-Md(u[1],u[2],u[3],k,u[4])
      }
    }
    predM[j]<-Ma(x1[1],x1[2],x1[3],t2)-Ma(x1[1],x1[2],x1[3],t1)
    if(flag == 'S.lognorm'){
      Md<-function(a,b,c,y,miu,sigm){
        integrate(function(x) Ma(a,b,c,x)*flag1(y,x,miu,sigm),0,y)$value
      }
      predN[j]<-Md(x1[1],x1[2],x1[3],t2,x1[4],x1[5])-Md(x1[1],x1[2],x1[3],t1,x1[4],x1[5])
    }
    else if(flag == 'S.exp'){
      Md<-function(a,b,c,y,miu){
        integrate(function(x) Ma(a,b,c,x)*flag1(y,x,miu),0,y)$value
      }
      predN[j]<-Md(x1[1],x1[2],x1[3],t2,x1[4])-Md(x1[1],x1[2],x1[3],t1,x1[4])
    }


  }
  z1<-matrix(0,T,2)
  z2<-matrix(0,T,2)
  for (k in 1:T) {
    z1[k,]<-c(quantile(M[,k],0.025),quantile(M[,k],0.975))
    z2[k,]<-c(quantile(N[,k],0.025),quantile(N[,k],0.975))
  }
  z3<-c(quantile(predM,0.025),quantile(predM,0.975))
  z4<-c(quantile(predN,0.025),quantile(predN,0.975))

  list_data<-list(z1,z2,z3,z4)
  names(list_data)<-c('NA.CI','ND.CI','NA.predict','ND.predict')
  print(list_data)
  tt<-1:T
  par(mfrow=c(1,2))
  plot(tt,z1[,2],col='red',type='l',xlab='',ylab='',xlim=c(0,T),ylim=c(0,max(z1[,2])))
  par(new=TRUE)
  plot(tt,z1[,1],col='blue',type='l',xlab='',ylab='',xlim=c(0,T),ylim=c(0,max(z1[,2])))
  par(new=TRUE)
  plot(tt,mMa,col='green',axes=FALSE,xlab='time(h)',ylab='arrival',type='l',xlim=c(0,T),ylim=c(0,max(z1[,2])))
  plot(tt,z2[,2],col='red',xlab='',ylab='',xlim=c(0,T),ylim=c(0,max(z2[,2])),type='l')
  par(new=TRUE)
  plot(tt,z2[,1],col='blue',type='l',xlab='',ylab='',xlim=c(0,T),ylim=c(0,max(z2[,2])))
  par(new=TRUE)
  plot(tt,mMd,col='green',axes=FALSE,xlab='time(h)',ylab='departure',type='l',xlim=c(0,T),ylim=c(0,max(z2[,2])))
}

#' Estimate the parameters
#'
#' Use MLE method to estimate the parameters when the arrival rate function is
#' S-shaped and the service time's PDF is exponential distribution.And also uses funcion 'nlminb'
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
#' m<-S.exp(na=c(2,2,2,3),nd=c(2,2,2,2),c=c(100,10,1))

S.exp<-function(na,nd,c) {
  lSe2 <- function(x) {

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

    a=x[1]
    b=x[2]
    c=x[3]
    miu=x[4]

    mS=function(a,b,c,y){
      a*(1-exp(-b*y))/(1+c*exp(-b*y))
    }

    p1=c(1,rep(-1,n-1))
    p=rep(0,n)
    q=rep(0,n)

    int3=integrate(function(y) miu*exp(-miu*(t[1]-y))*mS(a,b,c,y),0,t[1])$value
    p2=rep((pexp(0,miu)*mS(a,b,c,t[1])+int3)/mS(a,b,c,t[1]),n)
    for(i in 2:n) {
      int4=integrate(function(y) miu*exp(-miu*(t[i]-y))*mS(a,b,c,y),t[i-1],t[i])$value
      p2[i]=(pexp(0,miu)*mS(a,b,c,t[i])-pexp(t[i]-t[i-1],miu)*mS(a,b,c,t[i-1])+int4)/(mS(a,b,c,t[i])-mS(a,b,c,t[i-1]))
      if (t[i-1]!=t0[i]){
        int1=integrate(function(y) miu*exp(-miu*(t[i-1]-y))*mS(a,b,c,y),t0[i],t[i-1])$value
        p[i]=(pexp(0,miu)*mS(a,b,c,t[i-1])-pexp(t[i-1]-t0[i],miu)*mS(a,b,c,t0[i])+int1)/(mS(a,b,c,t[i-1])-mS(a,b,c,t0[i]))

        int2=integrate(function(y) miu*exp(-miu*(t[i]-y))*mS(a,b,c,y),t0[i],t[i-1])$value
        q[i]=(pexp(t[i]-t[i-1],miu)*mS(a,b,c,t[i-1])-pexp(t[i]-t0[i],miu)*mS(a,b,c,t0[i])+int2)/(mS(a,b,c,t[i-1])-mS(a,b,c,t0[i]))-p[i]
        p1[i]=q[i]/(1-p[i])
      }
    }

    sum1=0;
    lsum=0;
    for(i in 2:n) {
      sum1=sum1+da[i]*log(mS(a,b,c,t[i])-mS(a,b,c,t[i-1]))-(mS(a,b,c,t[i])-mS(a,b,c,t[i-1]))
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
    -sum1-lsum-(da[1]*log(mS(a,b,c,t[1]))-mS(a,b,c,t[1]))-log(dbinom(dd[1],da[1],p2[1]))

  }

  z<-nlminb(c,lSe2)
  return(z)
}

#' Estimate the parameters
#'
#' Use MLE method to estimate the parameters when the arrival rate function is
#' S-shaped and the service time's PDF is lognormal distribution.And also uses funcion 'nlminb'
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
#' m<-S.lognorm(na=c(2,2,2,3),nd=c(2,2,2,2),c=c(100,10,1))

S.lognorm<-function(na,nd,c) {
  lSn <- function(x) {
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
    a=x[1]
    b=x[2]
    c=x[3]
    miu=x[4]
    sigm=x[5]

    mS=function(a,b,c,y){
      a*(1-exp(-b*y))/(1+c*exp(-b*y))
    }

    p1=c(1,rep(-1,n-1))
    p=rep(0,n)
    q=rep(0,n)

    int3=integrate(function(y) exp(-(log(t[1]-y)-miu)^2/(2*sigm^2))/((t[1]-y)*sigm*(2*pi)^0.5)*mS(a,b,c,y),0,t[1])$value
    temp<-try(int3)
    if('try-error' %in% class(temp)){
      print('The initial parameters need to be modified so that it can be calculated.')
      break
    }
    p2=rep((plnorm(0,miu,sigm)*mS(a,b,c,t[1])+int3)/mS(a,b,c,t[1]),n)
    for(i in 2:n) {
      int4=integrate(function(y) exp(-(log(t[i]-y)-miu)^2/(2*sigm^2))/((t[i]-y)*sigm*(2*pi)^0.5)*mS(a,b,c,y),t[i-1],t[i])$value
      p2[i]=(plnorm(0,miu,sigm)*mS(a,b,c,t[i])-plnorm(t[i]-t[i-1],miu,sigm)*mS(a,b,c,t[i-1])+int4)/(mS(a,b,c,t[i])-mS(a,b,c,t[i-1]))
      if (t[i-1]!=t0[i]){
        int1=integrate(function(y) exp(-(log(t[i-1]-y)-miu)^2/(2*sigm^2))/((t[i-1]-y)*sigm*(2*pi)^0.5)*mS(a,b,c,y),t0[i],t[i-1])$value
        p[i]=(plnorm(0,miu,sigm)*mS(a,b,c,t[i-1])-plnorm(t[i-1]-t0[i],miu,sigm)*mS(a,b,c,t0[i])+int1)/(mS(a,b,c,t[i-1])-mS(a,b,c,t0[i]))

        int2=integrate(function(y) exp(-(log(t[i]-y)-miu)^2/(2*sigm^2))/((t[i]-y)*sigm*(2*pi)^0.5)*mS(a,b,c,y),t0[i],t[i-1])$value
        q[i]=(plnorm(t[i]-t[i-1],miu,sigm)*mS(a,b,c,t[i-1])-plnorm(t[i]-t0[i],miu,sigm)*mS(a,b,c,t0[i])+int2)/(mS(a,b,c,t[i-1])-mS(a,b,c,t0[i]))-p[i]
        p1[i]=q[i]/(1-p[i])
      }
    }

    sum1=0;
    lsum=0;
    for(i in 2:n) {
      sum1=sum1+da[i]*log(mS(a,b,c,t[i])-mS(a,b,c,t[i-1]))-(mS(a,b,c,t[i])-mS(a,b,c,t[i-1]))
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
    -sum1-lsum-(da[1]*log(mS(a,b,c,t[1]))-mS(a,b,c,t[1]))-log(dbinom(dd[1],da[1],p2[1]))

  }
  z<-nlminb(c,lSn)
  return(z)
}

#' Produce 'na' and 'nd'
#'
#' This code is to produce a sequence of the arrival number and
#' the departure number with specific parameters when the arrival rate
#' function is S-shaped which is needed in other funcion of this package.
#'
#' @param u the value of the specific parameters that already known
#' @param T the totle time to work
#'
#' @return a sequence of the arrival number and the departure number
#' @export
#'
#' @examples
#' m<-S.nand(u=c(100,10,1),T=50)
S.nand<-function(u,T){
  a<-u[1]
  b<-u[2]
  c<-u[3]
  M<-rep(0,1000)
  mS=function(y){
    a*(1-exp(-b*y))/(1+c*exp(-b*y))
  }
  # for (j in 1:1000) {
  N=rpois(1,mS(T))
  ta=rep(0,N)
  for (i in 1:N){
    u=runif(1,0,1)
    ta[i]=uniroot(function(t) mS(t)-u*mS(T),c(0,T))$root
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
  m<-c(na,nd)

}
