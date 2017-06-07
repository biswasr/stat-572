mleout_lrtrev<-mle(minuslogl = minuslogL_lrt, start=list(e12_HL=e12_HL,e12_DL=e12_DL,e21=e21,tau1=tau1,tau2=tau2,mu12=mu12,mu13=mu13,mu15=mu15,mu34=mu34,mu35=mu35,pi_HL=pi_HL,pi_DL=pi_DL),control=list(maxit=2500,reltol=1e-8),method="Nelder-Mead")
minuslogL_lrt<-function(e12_HL,e12_DL,e21,tau1,tau2,mu12,mu13,mu15,mu34,mu35,pi_HL,pi_DL)
{e12_0=exp(e12_HL);e12_DL=(exp(e12_DL)-e12_0);e21=exp(e21);tau1=exp(tau1);tau2=exp(tau2);mu12=exp(mu12);mu13=exp(mu13);mu15=exp(mu15);mu31=0;mu34=exp(mu34);mu35=exp(mu35);pi_0=exp(pi_HL);pi_DL=exp(pi_DL)-pi_0;
mu23=(tau1)*mu13;
mu25=(tau1)*mu15;
mu41=0;
mu45=(tau2)*mu35;
Lambda=rbind(c(-(mu12+mu13+mu15),mu12,mu13,0,mu15),
             c(0,-(mu23+mu25),mu23,0,mu25),
             c(mu31,0,-(mu31+mu34+mu35),mu34,mu35),
             c(mu41,0,0,-(mu41+mu45),mu45),
             c(0,0,0,0,0));
Ltemp=list();
logL=numeric(N);
for(pat in 1:N)
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
