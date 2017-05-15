setwd("~/Documents/UWashington/courses/stat 572")
data<-read.csv("bosdata.csv")
str(data)
data<-data[,-5]
dat<-split(data,data$PTNUM)
the.data<-list()
cov.data<-matrix(nrow=364,ncol=2)
for(i in 1:length(dat))
{
  if(dat[[i]]$state[length(dat[[i]]$state)]==99)
  {
    dat[[i]]<-dat[[i]][-length(dat[[i]]$state),];
  }
}
for(i in 1:length(dat))
{
  the.data[[i]]<-matrix(ncol=3,nrow=nrow(dat[[i]]));
  colnames(the.data[[i]])=c("obsdata","obstimes","exact.times.ranks")
  the.data[[i]]<-data.frame(the.data[[i]])
  
  the.data[[i]]$obsdata<-dat[[i]]$state;
  the.data[[i]]$obstimes<-dat[[i]]$time;
  if(dat[[i]]$state[length(dat[[i]]$state)] %in% c(3,99) )
  {the.data[[i]]$exact.times.ranks<-dat[[i]]$state;}
  cov.data[i,]<-as.numeric(dat[[i]][1,c(1,4)])
}
colnames(cov.data)<-c("PTNUM","X")
rownames(cov.data)<-1:364
library(test4)
library(sqldf)
logit<-function(p){
  log(p/(1-p))
}
tprob<-function(t,Lambda)
{
  return(expm(Lambda*t))
}
obsstate=c(1,1,2,2,3);
minuslogL<-function(pars)
{
  e12_0=pars[1];
  e21=pars[2];
  nu12=pars[3];
  logtau1=pars[4];
  logtau2=pars[5];
  mu12=pars[6];
  mu13=pars[7];
  mu15=pars[8];
  mu23=exp(logtau1)*mu13;
  mu25=exp(logtau1)*mu15;
  mu31=pars[9];
  mu34=pars[10];
  mu35=pars[11];
  mu41=exp(logtau2)*mu31;
  mu45=exp(logtau2)*mu35;
  Lambda=rbind(c(-(mu12+mu13+mu15),mu12,mu13,0,mu15),
               c(0,-(mu23+mu25),mu23,0,mu25),
               c(mu31,0,-(mu31+mu34+mu35),mu34,mu35),
               c(mu41,0,0,-(mu41+mu45),mu45),
               c(0,0,0,0,0));
  pi0=pars[12];
  pislope=pars[13];
  Ltemp=list();
  logL=numeric(364);
for(pat in 1:364)
{
  e12=expit(logit(e12_0)+nu12*as.numeric(cov.data[pat,2]==4))
  emissionstates<-rbind(c(1-e12,e12,0),c(e21,1-e21,0),c(0,0,1));
  pi2=expit(logit(pi0)+pislope*as.numeric(cov.data[pat,2]==4))
  M<-list();
  Pi=c(1-pi2,0,pi2,0,0);
  if(nrow(dat[[pat]])>1)
{
    for(i in 2:nrow(dat[[pat]]))
{
  M[[i-1]]<-matrix(nrow=5,ncol=5);
  P=tprob(the.data[[pat]]$obstimes[i]-the.data[[pat]]$obstimes[i-1],Lambda);
  for(r in 1:5)
  {
    for(s in 1:5)
    {
      M[[i-1]][r,s]<-emissionstates[obsstate[s],the.data[[pat]]$obsdata[i]]* P[r,s]
    }
  }
}

Ltemp[[pat]]<-t(Pi)%*%M[[1]];
if(nrow(dat[[pat]])>2)
{for(i in 3:nrow(dat[[pat]]))
{
  Ltemp[[pat]]<-Ltemp[[pat]]%*%M[[i-1]]
}
}
L=Ltemp[[pat]]%*%rep(1,5)
logL[pat]=log(L);
}else
{s=0;
  for(i in 1:5)
  {s=s+Pi[i]*emissionstates[obsstate[i],the.data[[pat]]$obsdata[1]];}
logL[pat]=log(s);
}
}
totminuslogL=-sum(logL)
return(totminuslogL)
}
mleout2<-optim(fn = minuslogL, par=c(0.025,0.006,1.314,-1.093,-1.158,0.207,0.339,0.016,0.113,1.369,0.5,0.079,1.019),method = "BFGS",control=list(maxit=2500,reltol=1e-8),hessian=T)
###################################################
library(stats4)
library(clusterPower)
library(expm)
tprob<-function(t,Lambda)
{
  return(expm(Lambda*t))
}
obsstate=c(1,1,2,2,3);
minuslogL<-function(e12_0,e21,nu12,logtau1,logtau2,mu12,mu13,mu15,mu31,mu34,mu35,pi0,pislope)
{
  mu23=exp(logtau1)*mu13;
  mu25=exp(logtau1)*mu15;
  mu41=exp(logtau2)*mu31;
  mu45=exp(logtau2)*mu35;
  Lambda=rbind(c(-(mu12+mu13+mu15),mu12,mu13,0,mu15),
               c(0,-(mu23+mu25),mu23,0,mu25),
               c(mu31,0,-(mu31+mu34+mu35),mu34,mu35),
               c(mu41,0,0,-(mu41+mu45),mu45),
               c(0,0,0,0,0));
  Ltemp=list();
  logL=numeric(364);
  for(pat in 1:364)
  {
    e12=expit(logit(e12_0)+nu12*as.numeric(cov.data[pat,2]==4))
    emissionstates<-rbind(c(1-e12,e12,0),c(e21,1-e21,0),c(0,0,1));
    pi2=expit(logit(pi0)+pislope*as.numeric(cov.data[pat,2]==4))
    M<-list();
    Pi=c(1-pi2,0,pi2,0,0);
    if(nrow(dat[[pat]])>1)
    {
      for(i in 2:nrow(dat[[pat]]))
      {
        M[[i-1]]<-matrix(nrow=5,ncol=5);
        P=tprob(the.data[[pat]]$obstimes[i]-the.data[[pat]]$obstimes[i-1],Lambda);
        for(r in 1:5)
        {
          for(s in 1:5)
          {
            M[[i-1]][r,s]<-emissionstates[obsstate[s],the.data[[pat]]$obsdata[i]]* P[r,s]
          }
        }
      }
      
      Ltemp[[pat]]<-t(Pi)%*%M[[1]];
      if(nrow(dat[[pat]])>2)
      {for(i in 3:nrow(dat[[pat]]))
      {
        Ltemp[[pat]]<-Ltemp[[pat]]%*%M[[i-1]]
      }
      }
      L=Ltemp[[pat]]%*%rep(1,5)
      logL[pat]=log(L);
    }else
    {s=0;
    for(i in 1:5)
    {s=s+Pi[i]*emissionstates[obsstate[i],the.data[[pat]]$obsdata[1]];}
    logL[pat]=log(s);
    }
  }
  totminuslogL=-sum(logL)
  return(totminuslogL)
}
mleout2<-mle(minuslogl = minuslogL, start=list(e12_0=0.025,e21=0.006,nu12=1.314,logtau1=-1.093,logtau2=-1.158,mu12=0.207,mu13=0.339,mu15=0.016,mu31=0.113,mu34=1.369,mu35=0.5,pi0=0.079,pislope=1.019),method = "BFGS",control=list(maxit=2500,reltol=1e-8))
minuslogLhmm<-function(e12_0,e21,nu12,mu12,mu13,mu21,mu23,pi0,pislope)
{
  Lambda=rbind(c(-(mu12+mu13),mu12,mu13),
               c(mu21,-(mu21+mu23),mu23),
               c(0,0,0));
  Ltemp=list();
  logL=numeric(364);
  for(pat in 1:364)
  {
    e12=expit(logit(e12_0)+nu12*as.numeric(cov.data[pat,2]==4))
    emissionstates<-rbind(c(1-e12,e12,0),c(e21,1-e21,0),c(0,0,1));
    pi2=expit(logit(pi0)+pislope*as.numeric(cov.data[pat,2]==4))
    M<-list();
    Pi=c(1-pi2,pi2,0);
    if(nrow(dat[[pat]])>1)
    {
      for(i in 2:nrow(dat[[pat]]))
      {
        M[[i-1]]<-matrix(nrow=3,ncol=3);
        P=tprob(the.data[[pat]]$obstimes[i]-the.data[[pat]]$obstimes[i-1],Lambda);
        for(r in 1:3)
        {
          for(s in 1:3)
          {
            M[[i-1]][r,s]<-emissionstates[s,the.data[[pat]]$obsdata[i]]* P[r,s]
          }
        }
      }
      
      Ltemp[[pat]]<-t(Pi)%*%M[[1]];
      if(nrow(dat[[pat]])>2)
      {for(i in 3:nrow(dat[[pat]]))
      {
        Ltemp[[pat]]<-Ltemp[[pat]]%*%M[[i-1]]
      }
      }
      L=Ltemp[[pat]]%*%rep(1,3)
      logL[pat]=log(L);
    }else
    {s=0;
    for(i in 1:3)
    {s=s+Pi[i]*emissionstates[i,the.data[[pat]]$obsdata[1]];}
    logL[pat]=log(s);
    }
  }
  totminuslogL=-sum(logL)
  return(totminuslogL)
}
mlehmm<-mle(minuslogl = minuslogLhmm, start=list(e12_0=0.028,e21=0.007,nu12=1.306,mu12=0.197,mu13=0.029,mu21=0.035,mu23=0.195,pi0=0.094,pislope=0.937),method = "BFGS",control=list(maxit=2500,reltol=1e-8))


Lambdahat=rbind(c(-(mleout2@coef['mu12']+mleout2@coef['mu13']+mleout2@coef['mu15']),mleout2@coef['mu12'],mleout2@coef['mu13'],0,mleout2@coef['mu15']),c(0,-(mleout2@coef['mu23']+mleout2@coef['mu25']),mleout2@coef['mu23'],0,mleout2@coef['mu25']),c(mleout2@coef['mu31'],0,-(mleout2@coef['mu31']+mleout2@coef['mu34']+mleout2@coef['mu35']),mleout2@coef['mu34'],mleout2@coef['mu35']),c(mleout2@coef['mu41'],0,0,-(mleout2@coef['mu41']+mleout2@coef['mu45']),mleout2@coef['mu45']),c(0,0,0,0,0));

hazarddeathGhealthy = numeric(100);
hazarddeathGbos = numeric(100);
t=seq(from=0.1,to=10,by=0.1)
for(i in 1:100)
{
  hazarddeathGhealthy[i]=expm(Lambdahat * t[i])[1,3];
  hazarddeathGbos[i]=expm(Lambdahat * t[i])[2,3];
}
