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

minuslogL<-function(e12_HL,e12_DL,e21,tau1,tau2,mu12,mu13,mu15,mu31,mu34,mu35,pi_HL,pi_DL)
{e12_0=exp(e12_HL);e12_DL=(exp(e12_DL)-e12_0);e21=exp(e21);tau1=exp(tau1);tau2=exp(tau2);mu12=exp(mu12);mu13=exp(mu13);mu15=exp(mu15);mu31=exp(mu31);mu34=exp(mu34);mu35=exp(mu35);pi_0=exp(pi_HL);pi_DL=exp(pi_DL)-pi_0;
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
    }#else
    #{s=0;
    #for(i in 1:5)
    #{s=s+Pi[i]*emissionstates[obsstate[i],the.data[[pat]]$obsdata[1]];}
    #logL[pat]=log(s);
    #}
  }
  totminuslogL=-sum(logL)
  return(totminuslogL)
}

e12_HL=log(0.018)
e12_DL=log(0.076)
e21=log(0.011)
tau1=log(0.3623437)
tau2=log(0.2544685)
mu12=log(0.3870182)
mu13=log(0.3897474)
mu15=log(0.01201972)
mu31=log(0.06266244)
mu34=log(3.117166)
mu35=log(0.7275771)
pi_HL=log(0.061)
pi_DL=log(0.043)

print(minuslogL(e12_HL,e12_DL,e21,tau1,tau2,mu12,mu13,mu15,mu31,mu34,mu35,pi_HL,pi_DL))
mleout2<-mle(minuslogl = minuslogL, start=list(e12_HL=e12_HL,e12_DL=e12_DL,e21=e21,tau1=tau1,tau2=tau2,mu12=mu12,mu13=mu13,mu15=mu15,mu31=mu31,mu34=mu34,mu35=mu35,pi_HL=pi_HL,pi_DL=pi_DL),control=list(maxit=2500,reltol=1e-8),method="BFGS")

H<-matrix(data=0,nrow=13,ncol=13);
for(i in 1:13)
{
    H[i,i]<-exp(mleout2@coef[i])
}
covexp<-H%*%(mleout2@vcov)%*%t(H)
seexp=sqrt(diag(covexp))
CI<-cbind(exp(mleout2@coef)-1.96*seexp,exp(mleout2@coef)+1.96*seexp)
CI<-cbind((mleout2@coef)-1.96*sqrt(diag(mleout2@vcov)),(mleout2@coef)+1.96*sqrt(diag(mleout2@vcov)))
round(exp(CI),2)

minuslogLhmm<-function(e12_HL,e12_DL,e21,mu12,mu13,mu21,mu23,pi_HL,pi_DL)
{
  e12_0=exp(e12_HL);e21=exp(e21);nu12=exp(e12_DL)-e12_0;mu12=exp(mu12);mu13=exp(mu13);mu21=exp(mu21);mu23=exp(mu23);pi0=exp(pi_HL);pislope=exp(pi_DL)-pi_0;
  Lambda=rbind(c(-(mu12+mu13),mu12,mu13),
               c(mu21,-(mu21+mu23),mu23),
               c(0,0,0));
  Ltemp=list();
  logL=numeric(364);
  for(pat in 1:364)
  {
    e12=expit(e12_0+nu12*as.numeric(cov.data[pat,2]==4))
    emissionstates<-rbind(c(1-e12,e12,0),c(e21,1-e21,0),c(0,0,1));
    pi2=expit(pi0+pislope*as.numeric(cov.data[pat,2]==4))
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
mlehmm<-mle(minuslogl = minuslogLhmm, start=list(e12_HL=log(0.028),e12_DL=log(1.306+0.028),e21=log(0.007),mu12=log(0.197),mu13=log(0.029),mu21=log(0.035),mu23=log(0.195),pi_HL=log(0.094),pi_DL=log(0.937+0.094)),control=list(maxit=2500,reltol=1e-8))

CI<-cbind((mlehmm@coef)-1.96*sqrt(diag(mlehmm@vcov)),(mlehmm@coef)+1.96*sqrt(diag(mlehmm@vcov)))
round(exp(CI),2)

H2<-matrix(data=0,nrow=9,ncol=9);
for(i in 1:9)
{
  H2[i,i]<-exp(mlehmm@coef[i])
}
covexp2<-H2%*%(mlehmm@vcov)%*%t(H2)
seexp2=sqrt(diag(covexp2))
CI2<-cbind(exp(mlehmm@coef)-1.96*seexp2,exp(mlehmm@coef)+1.96*seexp2)

###
mu12hat=exp(mleout2@coef['mu12'])
mu13hat=exp(mleout2@coef['mu13'])
mu15hat=exp(mleout2@coef['mu15'])
mu31hat=exp(mleout2@coef['mu31']);
mu34hat=exp(mleout2@coef['mu34']);
mu35hat=exp(mleout2@coef['mu35']);
tau1hat=exp(mleout2@coef['tau1'])
tau2hat=exp(mleout2@coef['tau2'])
mu23hat=tau1hat*mu13hat;
mu25hat=tau1hat*mu15hat;
mu41hat=tau2hat*mu31hat;
mu45hat=tau2hat*mu35hat;
pi_HLhat=expit(exp(mleout2@coef['pi_HL']))
pi_DLhat=expit(exp(mleout2@coef['pi_HL'])+exp(mleout2@coef['pi_DL']))

#mu12hat=exp(mu12)
#mu13hat=exp(mu13)
#mu15hat=exp(mu15)
#mu31hat=exp(mu31)
#mu34hat=exp(mu34)
#mu35hat=exp(mu35)
#mu23hat=exp(mu23)
#mu25hat=exp(mu25)
#mu41hat=exp(mu41)
#mu45hat=exp(mu45)

Lambdahat=rbind(c(-(mu12hat+mu13hat+mu15hat),mu12hat,mu13hat,0,mu15hat),c(0,-(mu23hat+mu25hat),mu23hat,0,mu25hat),c(mu31hat,0,-(mu31hat+mu34hat+mu35hat),mu34hat,mu35hat),c(mu41hat,0,0,-(mu41hat+mu45hat),mu45hat),c(0,0,0,0,0));
Lambdahat2=rbind(c(-(mu12hat+mu13hat+mu15hat),mu12hat,mu13hat,0,mu15hat),c(0,-(mu23hat+mu25hat),mu23hat,0,mu25hat),c(0,0,0,0,0),c(mu41hat,0,0,-(mu41hat+mu45hat),mu45hat),c(0,0,0,0,0));

B<-list();
B[[1]]=rbind(c(-1,1,0,0,0),c(0,0,0,0,0),c(0,0,0,0,0),c(0,0,0,0,0),c(0,0,0,0,0))#12
B[[2]]=rbind(c(-1,0,1,0,0),c(0,0,0,0,0),c(0,0,0,0,0),c(0,0,0,0,0),c(0,0,0,0,0))#13
B[[3]]=rbind(c(-1,0,0,0,1),c(0,0,0,0,0),c(0,0,0,0,0),c(0,0,0,0,0),c(0,0,0,0,0))#15
B[[4]]=rbind(c(0,0,0,0,0),c(0,0,0,0,0),c(1,0,-1,0,0),c(tau1hat,0,0,-tau1hat,0),c(0,0,0,0,0))#31
B[[5]]=rbind(c(0,0,0,0,0),c(0,0,0,0,0),c(0,0,-1,1,0),c(0,0,0,0,0),c(0,0,0,0,0))#34
B[[6]]=rbind(c(0,0,0,0,0),c(0,0,0,0,0),c(0,0,-1,0,1),c(0,0,0,0,0),c(0,0,0,0,0))#35
B[[7]]=rbind(c(0,0,0,0,0),c(0,-(mu13hat+mu15hat),mu13hat,0,mu15hat),c(0,0,0,0,0),c(0,0,0,0,0),c(0,0,0,0,0))#tau1 present in mu23=tau1*mu13, mu25=tau1*mu15
B[[8]]=rbind(c(0,0,0,0,0),c(0,0,0,0,0),c(0,0,0,0,0),c(mu31hat,0,0,-(mu31hat+mu35hat),mu35hat),c(0,0,0,0,0))#tau2 present in mu41=tau2*mu31, mu45=tau2*mu35

C<-list();
C[[1]]=rbind(c(-1,1,0,0,0),c(0,0,0,0,0),c(0,0,0,0,0),c(0,0,0,0,0),c(0,0,0,0,0))#12
C[[2]]=rbind(c(-1,0,1,0,0),c(0,0,0,0,0),c(0,0,0,0,0),c(0,0,0,0,0),c(0,0,0,0,0))#13
C[[3]]=rbind(c(-1,0,0,0,1),c(0,0,0,0,0),c(0,0,0,0,0),c(0,0,0,0,0),c(0,0,0,0,0))#15
C[[4]]=rbind(c(0,0,0,0,0),c(0,0,0,0,0),c(0,0,0,0,0),c(tau1hat,0,0,-tau1hat,0),c(0,0,0,0,0))#31
C[[5]]=rbind(c(0,0,0,0,0),c(0,0,0,0,0),c(0,0,0,0,0),c(0,0,0,0,0),c(0,0,0,0,0))#34
C[[6]]=rbind(c(0,0,0,0,0),c(0,0,0,0,0),c(0,0,0,0,0),c(0,0,0,0,0),c(0,0,0,0,0))#35
C[[7]]=rbind(c(0,0,0,0,0),c(0,-(mu13hat+mu15hat),mu13hat,0,mu15hat),c(0,0,0,0,0),c(0,0,0,0,0),c(0,0,0,0,0))#tau1 present in mu23=tau1*mu13, mu25=tau1*mu15
C[[8]]=rbind(c(0,0,0,0,0),c(0,0,0,0,0),c(0,0,0,0,0),c(mu31hat,0,0,-(mu31hat+mu35hat),mu35hat),c(0,0,0,0,0))#tau2 present in mu41=tau2*mu31, mu45=tau2*mu35

#Lambdahat=rbind(c(0,0,0,0,0),c(mu15hat,-(mu12hat+mu13hat+mu15hat),mu12hat,mu13hat,0),c(mu25hat,0,-(mu23hat+mu25hat),mu23hat,0),c(mu35hat,mu31hat,0,-(mu31hat+mu34hat+mu35hat),mu34hat),c(mu45hat,mu41hat,0,0,-(mu41hat+mu45hat)));#1=3,2=1_1,3=1_2,4=2_1,5=2_2
#Lambdahat2=rbind(c(0,0,0,0,0),c(0,0,0,0,0),c(mu15hat,mu13hat,-(mu12hat+mu13hat+mu15hat),mu12hat,0),c(mu25hat,mu23hat,0,-(mu23hat+mu25hat),0),c(mu45hat,0,mu41hat,0,-(mu41hat+mu45hat)));#1=3,2=2_1,3=1_1,4=1_2,5=2_2

t=seq(from=0.01,to=10,by=0.01)

hazarddeathGhealthy = numeric(length(t));
hazarddeathGbos = numeric(length(t));
cprobbosGhealthy = numeric(length(t));
cprobdeathGhealthy= numeric(length(t));
dendeathGhealthy= numeric(length(t));
cprobdeathGbos= numeric(length(t));
denddeathGbos= numeric(length(t));
cprobbosGhealthy= numeric(length(t));
denbosGhealthy= numeric(length(t));
hazardbosGhealthy=numeric(length(t));


for(i in 1:length(t))
{
  #den1<-expm(Lambdahat[2:5,2:5] * t[i] ) %*%Lambdahat[,1];
  #dendeathGhealthy[i]= den1[2];
  #cprobdeathGhealthy[i]= expm(Lambdahat * t[i])[2,1];
  #e12_HL,e12_DL,e21,tau1,tau2,mu12,mu13,mu15,mu31,mu34,mu35,pi_HL,pi_DL
  
  cprobmat=expm(Lambdahat * t[i])
  denmat=expm(Lambdahat * t[i])%*%Lambdahat;
  cprobdeathGhealthy[i]= cprobmat[1,5];
  for(k in 1:8)
  {dcprobmat15_dpar[k]=integrate(function(tau){(expm(Lambdahat * tau)%*%B[[k]]%*%expm(Lambdahat * (t[i]-tau)))[1,5]},lower=0,upper=t[i])}
  var_cprobdeathGhealthy[i]=t(c(0,0,0,dcprobmat15_dpar[7],dcprobmat15_dpar[8],dcprobmat15_dpar[1:6],0,0))%*%mleout2@vcov%*%c(0,0,0,dcprobmat15_dpar[7],dcprobmat15_dpar[8],dcprobmat15_dpar[1:6],0,0);
  leftCI_cprobdeathGhealthy[i]=cprobdeathGhealthy[i]-1.96*sqrt(var_cprobdeathGhealthy[i]);
  righttCI_cprobdeathGhealthy[i]=cprobdeathGhealthy[i]+1.96*sqrt(var_cprobdeathGhealthy[i]);
  
  dendeathGhealthy[i]=denmat[1,5]
  hazarddeathGhealthy[i]=dendeathGhealthy[i]/(1-cprobdeathGhealthy[i])
  #denddeathGbos[i]=den1[4];
  #cprobdeathGbos[i]=expm(Lambdahat)[4,1];
  
  cprobdeathGbos[i]=cprobmat[3,5];
  denddeathGbos[i]=denmat[3,5];
  hazarddeathGbos[i]=denddeathGbos[i]/(1-cprobdeathGbos[i])
  
  #den2<-expm(Lambdahat2[3:5,3:5]) %*%Lambdahat[,1:2];
  #denoncebosGhealthy=den2[3,2];
  cprobmat2=expm(Lambdahat2 * t[i])
  denmat2=expm(Lambdahat2 * t[i])%*%Lambdahat2;
  cprobbosGhealthy[i]= cprobmat2[1,3];
  for(k in 1:8)
  {dcprobmat13_dpar[k]=integrate(function(tau){(expm(Lambdahat2 * tau)%*%C[[k]]%*%expm(Lambdahat2 * (t[i]-tau)))[1,5]},lower=0,upper=t[i])}
  var_cprobbosGhealthy[i]=t(c(0,0,0,dcprobmat13_dpar[7],dcprobmat13_dpar[8],dcprobmat13_dpar[1:6],0,0))%*%mleout2@vcov%*%c(0,0,0,dcprobmat13_dpar[7],dcprobmat13_dpar[8],dcprobmat13_dpar[1:6],0,0);
  leftCI_cprobbosGhealthy[i]=cprobbosGhealthy[i]-1.96*sqrt(var_cprobbosGhealthy[i]);
  righttCI_cprobbosGhealthy[i]=cprobbosGhealthy[i]+1.96*sqrt(var_cprobbosGhealthy[i]);
  
  denbosGhealthy[i]=denmat2[1,3]
  hazardbosGhealthy[i]=denbosGhealthy[i]/(1-cprobbosGhealthy[i]);
}

#hazarddeathGbos<-denddeathGbos/cprobdeathGbos

plot(y=cprobdeathGhealthy,x=t,type="l",col="red",xlab="t-t0 years",ylab="Cum. probability of death",ylim=c(0,1))
lines(y=cprobdeathGbos,x=t,type="l",col="blue")

plot(y=cprobbosGhealthy,x=t,type="l",col="red",xlab="t-t0 years",ylab="Cum. probability of BOS",ylim=c(0,1))

plot(y=hazardbosGhealthy,x=t,type="l",col="red",xlab="t-t0 years",ylab="BOS rate/year",ylim=c(0,1))

plot(y=hazarddeathGbos,x=t,type="l",col="blue",xlab="t-t0 years",ylab="Mortality rate/year",ylim=c(0,0.8))
lines(y=hazarddeathGhealthy,x=t,type="l",col="red",xlab="t-t0 years",ylab="Mortality rate/year")

1-cprobbosGhealthy[500]
hazardbosGhealthy[1]
hazardbosGhealthy[500]
cprobdeathGbos[200]
1-cprobdeathGbos[200]
cprobdeathGhealthy[200]
cprobdeathGbos[200]
cprobdeathGbos[100]
1-cprobdeathGbos[100]
1-cprobdeathGbos[300]
1-cprobdeathGbos[500]
hazarddeathGbos[1]
hazarddeathGbos[5]
hazarddeathGbos[200]

#######
mu12hat=exp(mlehmm@coef['mu12']);
mu13hat=exp(mlehmm@coef['mu13']);
mu21hat=exp(mlehmm@coef['mu21']);
mu23hat=exp(mlehmm@coef['mu23']);
#pi_HLhat=expit(exp(mlehmm@coef['pi_HL']))
#pi_DLhat=expit(exp(mlehmm@coef['pi_HL'])+exp(mlehmm@coef['pi_DL']))
#mu12hat=exp(mu12)
#mu13hat=exp(mu13)
#mu15hat=exp(mu15)
#mu31hat=exp(mu31)
#mu34hat=exp(mu34)
#mu35hat=exp(mu35)
#mu23hat=exp(mu23)
#mu25hat=exp(mu25)
#mu41hat=exp(mu41)
#mu45hat=exp(mu45)

Lambdahat=rbind(c(-(mu12hat+mu13hat),mu12hat,mu13hat),
                       c(mu21hat,-(mu21hat+mu23hat),mu23),
                       c(0,0,0));
Lambdahat2=rbind(c(-(mu12hat+mu13hat),mu12hat,mu13hat),
                c(0,0,0),
                c(0,0,0));

#Lambdahat=rbind(c(0,0,0,0,0),c(mu15hat,-(mu12hat+mu13hat+mu15hat),mu12hat,mu13hat,0),c(mu25hat,0,-(mu23hat+mu25hat),mu23hat,0),c(mu35hat,mu31hat,0,-(mu31hat+mu34hat+mu35hat),mu34hat),c(mu45hat,mu41hat,0,0,-(mu41hat+mu45hat)));#1=3,2=1_1,3=1_2,4=2_1,5=2_2
#Lambdahat2=rbind(c(0,0,0,0,0),c(0,0,0,0,0),c(mu15hat,mu13hat,-(mu12hat+mu13hat+mu15hat),mu12hat,0),c(mu25hat,mu23hat,0,-(mu23hat+mu25hat),0),c(mu45hat,0,mu41hat,0,-(mu41hat+mu45hat)));#1=3,2=2_1,3=1_1,4=1_2,5=2_2

t=seq(from=0.01,to=10,by=0.01)

hazarddeathGhealthyhmm = numeric(length(t));
hazarddeathGboshmm = numeric(length(t));
cprobbosGhealthyhmm = numeric(length(t));
cprobdeathGhealthyhmm= numeric(length(t));
dendeathGhealthyhmm= numeric(length(t));
cprobdeathGboshmm= numeric(length(t));
denddeathGboshmm= numeric(length(t));
cprobbosGhealthyhmm= numeric(length(t));
denbosGhealthyhmm= numeric(length(t));
hazardbosGhealthyhmm=numeric(length(t));


for(i in 1:length(t))
{
  #den1<-expm(Lambdahat[2:5,2:5] * t[i] ) %*%Lambdahat[,1];
  #dendeathGhealthy[i]= den1[2];
  #cprobdeathGhealthy[i]= expm(Lambdahat * t[i])[2,1];
  
  
  cprobmathmm=expm(Lambdahat * t[i])
  denmathmm=expm(Lambdahat * t[i])%*%Lambdahat;
  cprobdeathGhealthyhmm[i]= cprobmathmm[1,3];
  dendeathGhealthyhmm[i]=denmathmm[1,3]
  hazarddeathGhealthyhmm[i]=dendeathGhealthyhmm[i]/(1-cprobdeathGhealthyhmm[i])
  #denddeathGbos[i]=den1[4];
  #cprobdeathGbos[i]=expm(Lambdahat)[4,1];
  
  cprobdeathGboshmm[i]=cprobmathmm[2,3];
  denddeathGboshmm[i]=denmathmm[2,3];
  hazarddeathGboshmm[i]=denddeathGboshmm[i]/(1-cprobdeathGboshmm[i])
  
  #den2<-expm(Lambdahat2[3:5,3:5]) %*%Lambdahat[,1:2];
  #denoncebosGhealthy=den2[3,2];
  cprobmat2hmm=expm(Lambdahat2 * t[i])
  denmat2hmm=expm(Lambdahat2 * t[i])%*%Lambdahat2;
  cprobbosGhealthyhmm[i]= cprobmat2hmm[1,2];
  denbosGhealthyhmm[i]=denmat2hmm[1,2]
  hazardbosGhealthyhmm[i]=denbosGhealthyhmm[i]/(1-cprobbosGhealthyhmm[i]);
  
}
lines(y=cprobdeathGhealthyhmm,x=t,lty=3,type="l",col="red",xlab="t-t0 years",ylab="Cum. probability of death",ylim=c(0,1))
lines(y=cprobdeathGboshmm,x=t,lty=3,type="l",col="blue")

lines(y=cprobbosGhealthyhmm,x=t,lty=3,type="l",col="red",xlab="t-t0 years",ylab="Cum. probability of BOS",ylim=c(0,1))

lines(y=hazardbosGhealthyhmm,x=t,lty=3,type="l",col="red",xlab="t-t0 years",ylab="BOS rate/year",ylim=c(0,1))

lines(y=hazarddeathGboshmm,x=t,lty=3,type="l",col="blue",xlab="t-t0 years",ylab="Mortality rate/year",ylim=c(0,0.8))
lines(y=hazarddeathGhealthyhmm,x=t,lty=3,type="l",col="red",xlab="t-t0 years",ylab="Mortality rate/year")
####
plot(y=cprobdeathGhealthy,x=t,type="l",col="red",xlab="t-t0 years",ylab="Cum. probability of death",ylim=c(0,1))
lines(y=cprobdeathGbos,x=t,type="l",col="blue")
lines(y=cprobdeathGhealthyhmm,x=t,lty=3,type="l",col="red",xlab="t-t0 years",ylab="Cum. probability of death",ylim=c(0,1))
lines(y=cprobdeathGboshmm,x=t,lty=3,type="l",col="blue")

plot(y=cprobbosGhealthy,x=t,type="l",col="red",xlab="t-t0 years",ylab="Cum. probability of BOS",ylim=c(0,1))
lines(y=cprobbosGhealthyhmm,x=t,lty=3,type="l",col="red",xlab="t-t0 years",ylab="Cum. probability of BOS",ylim=c(0,1))

plot(y=hazardbosGhealthy,x=t,type="l",col="red",xlab="t-t0 years",ylab="BOS rate/year",ylim=c(0,1))
lines(y=hazardbosGhealthyhmm,x=t,lty=3,type="l",col="red",xlab="t-t0 years",ylab="BOS rate/year",ylim=c(0,1))

plot(y=hazarddeathGbos,x=t,type="l",col="blue",xlab="t-t0 years",ylab="Mortality rate/year",ylim=c(0,0.8))
lines(y=hazarddeathGhealthy,x=t,type="l",col="red",xlab="t-t0 years",ylab="Mortality rate/year")
lines(y=hazarddeathGboshmm,x=t,lty=3,type="l",col="blue",xlab="t-t0 years",ylab="Mortality rate/year",ylim=c(0,0.8))
lines(y=hazarddeathGhealthyhmm,x=t,lty=3,type="l",col="red",xlab="t-t0 years",ylab="Mortality rate/year")

plot(y=1-cprobdeathGbos,x=t,type="l",col="blue",xlab="t-t0 years",ylab="Survival probability",ylim=c(0,1))
lines(y=1-cprobdeathGhealthy,x=t,type="l",col="red",xlab="t-t0 years",ylab="Survival probability")
lines(y=1-cprobdeathGboshmm,x=t,lty=3,type="l",col="blue",xlab="t-t0 years",ylab="Survival probability",ylim=c(0,1))
lines(y=1-cprobdeathGhealthyhmm,x=t,lty=3,type="l",col="red",xlab="t-t0 years",ylab="Survival probability")