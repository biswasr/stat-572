library(stats4)
library(clusterPower)
library(expm)

setwd("~/Documents/UWashington/courses/stat 572")
data<-read.csv("bosdata.csv")
str(data)
data<-data[,-5]
dat<-split(data,data$PTNUM)

print(length(dat))#nof patients

s=0;for(pat in 1:364){
  if(nrow(dat[[pat]])==1)
  {s=s+1}
}
print(s)#nof patients with 1 observation

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
  the.data[[i]]<-matrix(ncol=2,nrow=nrow(dat[[i]]));
  colnames(the.data[[i]])=c("obsdata","obstimes")
  the.data[[i]]<-data.frame(the.data[[i]])
  
  the.data[[i]]$obsdata<-dat[[i]]$state;
  the.data[[i]]$obstimes<-dat[[i]]$time;
  cov.data[i,]<-as.numeric(dat[[i]][1,c(1,4)])
}
colnames(cov.data)<-c("PTNUM","X")
rownames(cov.data)<-1:364

logit<-function(p){
  log(p/(1-p))
}
tprob<-function(t,Lambda)
{
  return(expm(Lambda*t))
}

obsstate=c(1,1,2,2,3);

minuslogL<-function(e12_0,e12_DL,e21,tau1,tau2,mu12,mu13,mu15,mu31,mu34,mu35,pi_0,pi_DL)
{#e12_HL=exp(e12_HL);e12_DL=exp(e12_DL);e21=exp(e21);tau1=exp(tau1);tau2=exp(tau2);mu12=exp(mu12);mu13=exp(mu13);mu15=exp(mu15);mu31=exp(mu31);mu34=exp(mu34);mu35=exp(mu35);pi_HL=exp(pi_HL);pi_DL=exp(pi_DL)
  mu23=(tau1)*mu13;
  mu25=(tau1)*mu15;
  mu41=(tau2)*mu31;
  mu45=(tau2)*mu35;
  Lambda=rbind(c(-(mu12+mu13+mu15),mu12,mu13,0,mu15),
               c(0,-(mu23+mu25),mu23,0,mu25),
               c(mu31,0,-(mu31+mu34+mu35),mu34,mu35),
               c(mu41,0,0,-(mu41+mu45),mu45),
               c(0,0,0,0,0));
  Ltemp=list();
  logL=numeric(364);
  for(pat in 1:364)
  {
    e12=expit(logit(e12_0)+e12_DL*as.numeric(cov.data[pat,2]==4))
    emissionstates<-rbind(c(1-e12,e12,0),c(e21,1-e21,0),c(0,0,1));
    pi2=expit(logit(pi_0)+pi_DL*as.numeric(cov.data[pat,2]==4))
    Pi=c(1-pi2,0,pi2,0,0);
    M<-vector("list",nrow(dat[[pat]])-1);
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
      {for(i in 2:(nrow(dat[[pat]])-1))
      {
        Ltemp[[pat]]<-Ltemp[[pat]]%*%M[[i]]
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

e12_0=0.018
e12_DL=0.058
e21=0.011
tau1=0.3623437
tau2=0.2544685
mu12=0.3870182
mu13=0.3897474
mu15=0.01201972
mu31=0.06266244
mu34=3.117166
mu35=0.7275771
pi_0=0.061
pi_DL=-0.018

print(minuslogL(e12_0,e12_DL,e21,tau1,tau2,mu12,mu13,mu15,mu31,mu34,mu35,pi_0,pi_DL))
