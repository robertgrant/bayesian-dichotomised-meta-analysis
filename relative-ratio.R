# R and Stan code for the simulation study and example analysis in the paper 
# "A taxonomy of thresholds used to dichotomise outcomes, and their inclusion in Bayesian meta-analysis"



################################################################################################
##################################  Relative-ratio threshold  ##################################
################################################################################################


library(meta)
library(MASS)
library(rstan)
library(parallel)


##################### Parameters for all 3 functions ######################


t<-10 # Number of trials
true.mean0<-155 # baseline mean
true.time<-(-5) # time effect
true.sigma<-15 # constant sigma at baseline and endpoint (intra-trial)
true.tau<-4 # random effect SD (inter-trial)
true.corr<-0.7 # correlation of baseline and endpoint



##################### Relative-ratio threshold ######################


# data to get the models running
armsize<-10
threshold<-0.9 # this is the ratio
effectsize<-(-10)
staniter<-2000

true.delta<-(effectsize) # treatment effect - we want to estimate this

# initialise tk (trials known) and td (trials dichotomised) for the first test-run model
td<-3
tk<-t-td

# generate data and study stats
u0<-rnorm(t,mean=true.mean0,sd=true.tau) # random intercepts

# ready to store means, sds & risks
# see Stan model below for definintions
mc0<-mt0<-mc1<-mt1<-sc0<-st0<-sc1<-st1<-rc1<-rt1<-rep(NA,t)
for (i in 1:t) {
  temp.c<-mvrnorm(n=armsize,mu=c(u0[i],u0[i]+true.time),
                  Sigma=matrix(c(true.sigma^2,true.corr*(true.sigma^2),
                                 true.corr*(true.sigma^2),true.sigma^2),
                               nrow=2))
  temp.t<-mvrnorm(n=armsize,mu=c(u0[i],u0[i]+true.time+true.delta),
                  Sigma=matrix(c(true.sigma^2,true.corr*(true.sigma^2),
                                 true.corr*(true.sigma^2),true.sigma^2),
                               nrow=2))
  mc0[i]<-mean(temp.c[,1])
  sc0[i]<-sd(temp.c[,1])
  mt0[i]<-mean(temp.t[,1])
  st0[i]<-sd(temp.t[,1])
  mc1[i]<-mean(temp.c[,2])
  sc1[i]<-sd(temp.c[,2])
  rc1[i]<- sum((temp.c[,2]/temp.c[,1])<=0.9)
  mt1[i]<-mean(temp.t[,2])
  st1[i]<-sd(temp.t[,2])
  rt1[i]<- sum((temp.t[,2]/temp.t[,1])<=0.9)
}
rm(temp.c)
rm(temp.t)
true.dt<-mean(mt1-mc1) # because all trials have same sample size

# Now we will run the Stan models once to get them compiled and stored in an R list


################# Stan model ##################
ma2model<-'
data {
int<lower=0> npart; // the number of participants per arm
int<lower=0> ntrial; // total number of trials (nm+nnom)
int<lower=0> nm; // number of studies with means
int<lower=0> nnom; // number without means
real mc0[ntrial]; // baseline means in control
real mt0[ntrial]; // baseline means in treatment
real<lower=0> sc0[ntrial]; // SDs in control baseline
real<lower=0> st0[ntrial]; // SDs in treatment baseline
real mt1[nm]; // treatment means where known
real<lower=0> st1[ntrial]; // treatment SD known for all
int<lower=0> rt1[nnom]; // number of responders in tx group
real mc1[nm]; // control means where known
real<lower=0> sc1[ntrial]; // SDs in control known for all
int<lower=0> rc1[nnom]; // number of responders in control
}
transformed data {
real<lower=0> sec0[ntrial];
real<lower=0> set0[ntrial];
real<lower=0> set1[ntrial];
real<lower=0> sec1[ntrial];
for (i in 1:ntrial)
sec0[i]<-sc0[i]/sqrt(npart); // standard errors
for (i in 1:ntrial)
set0[i]<-st0[i]/sqrt(npart); // standard errors
for (i in 1:ntrial)
set1[i]<-st1[i]/sqrt(npart); // standard errors
for (i in 1:ntrial)
sec1[i]<-sc1[i]/sqrt(npart); // standard errors
}
parameters {
real mu0; // global mean for both groups at baseline
real<lower=0> tau; // heterogeneity precision
real u[ntrial]; // random effect
real dc; // global mean of time difference
real dt; // global mean of treatment difference
}
transformed parameters {
real<lower=0> sdtau; // heterogeneity SD
real<lower=0, upper=1> riskt1[nnom]; // risk of normotension tx
real<lower=0, upper=1> riskc1[nnom]; // risk of normo control
real diff; // dt-dc
sdtau<-1/sqrt(tau);
for (i in 1:nnom)
riskt1[i]<-Phi_approx((-(((-0.668965)*(mu0+u[nm+i]))+(0.743294*(mu0+u[nm+i]+dt))))/sqrt((0.447561*st0[nm+i]*st0[nm+i])+(2*(-0.668965)*0.743294*0.7*st0[nm+i]*st1[nm+i])+(0.5524862*st1[nm+i]*st1[nm+i])));
for (i in 1:nnom)
riskc1[i]<-Phi_approx((-(((-0.668965)*(mu0+u[nm+i]))+(0.743294*(mu0+u[nm+i]+dc))))/sqrt((0.447561*sc0[nm+i]*sc0[nm+i])+(2*(-0.668965)*0.743294*0.7*sc0[nm+i]*sc1[nm+i])+(0.5524862*sc1[nm+i]*sc1[nm+i])));
diff<-(dt-dc);
}
model {
for (i in 1:ntrial)
u[i] ~ normal(0,sdtau);
for (i in 1:ntrial)
mc0[i] ~ normal(mu0+u[i],sec0[i]);
for (i in 1:ntrial)
mt0[i] ~ normal(mu0+u[i],set0[i]);
for (i in 1:nm)
mc1[i] ~ normal(mu0+u[i]+dc,sec1[i]);
for (i in 1:nm)
mt1[i] ~ normal(mu0+u[i]+dt,set1[i]);
for (i in 1:nnom)
rc1[i] ~ binomial(npart,riskc1[i]);
for (i in 1:nnom)
rt1[i] ~ binomial(npart,riskt1[i]);
mu0 ~ normal(150,20);
dc ~ normal(-5,10);
dt ~ normal(-13,10);
tau ~ gamma(0.01,0.01);
}
'

# Initial values
ma2inits<-list(
  init1=list(mu0=154,dc=-5,dt=-11,tau=0.08,
             u=seq(from=-2/sqrt(0.08),
                   to=2/sqrt(0.08),
                   length.out=t)),
  init2=list(mu0=158,dc=-3,dt=-14,tau=0.1,
             u=seq(from=-2/sqrt(0.1),
                   to=2/sqrt(0.1),
                   length.out=t)))

# Data to feed into Stan
ma2data<-list(ntrial=t,npart=armsize,nm=tk,nnom=td,threshold=threshold,
              mc0=mc0,sc0=sc0,
              mt0=mt0,st0=st0,
              mc1=mc1[1:tk],sc1=sc1,rc1=rc1[(tk+1):t],
              mt1=mt1[1:tk],st1=st1,rt1=rt1[(tk+1):t])

# Run model once and get it compiled
fit2<-stan(model_code=ma2model,
           data=ma2data,
           iter=staniter,
           chains=2,
           init=ma2inits,
           seed=434)
print(fit2)

# Stan model for Bayesian meta-analysis with baseline data only 
ma1model<-'
data {
int<lower=0> armsize; // the number of participants per arm
int<lower=0> tk; // number of studies with means
real mc0[tk]; // baseline means in control
real mt0[tk]; // baseline means in treatment
real<lower=0> sc0[tk]; // SDs in control baseline
real<lower=0> st0[tk]; // SDs in treatment baseline
real mt1[tk]; // treatment means where known
real<lower=0> st1[tk]; // treatment SD known for all
real mc1[tk]; // control means where known
real<lower=0> sc1[tk]; // SDs in control known for all
}
transformed data {
real<lower=0> sec0[tk];
real<lower=0> set0[tk];
real<lower=0> set1[tk];
real<lower=0> sec1[tk];
for (i in 1:tk)
sec0[i] = sc0[i]/sqrt(armsize); // standard errors
for (i in 1:tk)
set0[i] = st0[i]/sqrt(armsize); // standard errors
for (i in 1:tk)
set1[i] = st1[i]/sqrt(armsize); // standard errors
for (i in 1:tk)
sec1[i] = sc1[i]/sqrt(armsize); // standard errors
}
parameters {
real mu0; // global mean for both groups at baseline
real<lower=0> tau; // heterogeneity precision
real u[tk]; // random effect
real dc; // global mean of time difference
real dt; // global mean of treatment difference
}
transformed parameters {
real<lower=0> sdtau; // heterogeneity SD
real diff; // dt-dc
sdtau = 1/sqrt(tau);
diff = (dt-dc);
}
model {
for (i in 1:tk)
u[i] ~ normal(0,sdtau);
for (i in 1:tk)
mc0[i] ~ normal(mu0+u[i],sec0[i]);
for (i in 1:tk)
mt0[i] ~ normal(mu0+u[i],set0[i]);
for (i in 1:tk)
mc1[i] ~ normal(mu0+u[i]+dc,sec1[i]);
for (i in 1:tk)
mt1[i] ~ normal(mu0+u[i]+dt,set1[i]);
mu0 ~ normal(150,20);
dc ~ normal(-5,10);
dt ~ normal(-13,10);
tau ~ gamma(0.01,0.01);
}
'
# Initial values NOTE THAT TAU ALSO FEEDS INTO U[i]
ma1inits<-list(
  init1=list(mu0=154,dc=-5,dt=-11,tau=0.08,
             u=seq(from=-2/sqrt(0.08),
                   to=2/sqrt(0.08),
                   length.out=tk)),
  init2=list(mu0=158,dc=-3,dt=-14,tau=0.1,
             u=seq(from=-2/sqrt(0.1),
                   to=2/sqrt(0.1),
                   length.out=tk)))
# Data to feed into Stan
ma1data<-list(armsize=armsize,tk=tk,
              mc0=mc0[1:tk],sc0=sc0[1:tk],
              mt0=mt0[1:tk],st0=st0[1:tk],
              mc1=mc1[1:tk],sc1=sc1[1:tk],
              mt1=mt1[1:tk],st1=st1[1:tk])
# Run model
fit1<-stan(model_code=ma1model,
           data=ma1data,
           iter=staniter,
           chains=2,
           init=ma1inits,
           seed=434)
print(fit1)

############ Classical meta-analysis for comparison ##############
complete.ma1<-metacont(n.e=rep(armsize,t),mean.e=mt1,sd.e=st1,
                       n.c=rep(armsize,t),mean.c=mc1,sd.c=sc1,
                       comb.fixed=FALSE,comb.random=TRUE)
summary(complete.ma1)
# and with only the trials reporting means:
ma1<-metacont(n.e=rep(armsize,tk),mean.e=mt1[1:tk],sd.e=st1[1:tk],
              n.c=rep(armsize,tk),mean.c=mc1[1:tk],sd.c=sc1[1:tk],
              comb.fixed=FALSE,comb.random=TRUE)
summary(ma1)








# compile Stan models for each of the Bayes scenarios
# rd2fits is a list of the R objects pointing to the compiled models
# rd1fits likewise, for the baseline-data-only models
rd1fits<-rd2fits<-list(NULL)

# make some data
u0<-rnorm(t,mean=true.mean0,sd=true.tau) # random intercepts
mc0<-mt0<-mc1<-mt1<-sc0<-st0<-sc1<-st1<-rc1<-rt1<-rep(NA,t)
for (i in 1:t) {
  temp.c<-mvrnorm(n=armsize,mu=c(u0[i],u0[i]+true.time),
                  Sigma=matrix(c(true.sigma^2,true.corr*(true.sigma^2),
                                 true.corr*(true.sigma^2),true.sigma^2),
                               nrow=2))
  temp.t<-mvrnorm(n=armsize,mu=c(u0[i],u0[i]+true.time+true.delta),
                  Sigma=matrix(c(true.sigma^2,true.corr*(true.sigma^2),
                                 true.corr*(true.sigma^2),true.sigma^2),
                               nrow=2))
  mc0[i]<-mean(temp.c[,1])
  sc0[i]<-sd(temp.c[,1])
  mt0[i]<-mean(temp.t[,1])
  st0[i]<-sd(temp.t[,1])
  mc1[i]<-mean(temp.c[,2])
  sc1[i]<-sd(temp.c[,2])
  rc1[i]<- sum((temp.c[,2]/temp.c[,1])<=0.9)
  mt1[i]<-mean(temp.t[,2])
  st1[i]<-sd(temp.t[,2])
  rt1[i]<- sum((temp.t[,2]/temp.t[,1])<=0.9)
}
rm(temp.c)
rm(temp.t)

td<-2
tk<-t-td
ma2data<-list(ntrial=t,npart=armsize,nm=tk,nnom=td,threshold=threshold,
              mc0=mc0,sc0=sc0,
              mt0=mt0,st0=st0,
              mc1=mc1[1:tk],sc1=sc1,rc1=rc1[(tk+1):t],
              mt1=mt1[1:tk],st1=st1,rt1=rt1[(tk+1):t])
rd2fits[[1]]<-stan(model_code=ma2model,data=ma2data,iter=staniter,
                   init=ma2inits,chains=2)
rd2fits[[2]]<-fit2 # This already exists

ma1data<-list(armsize=armsize,tk=tk,
              mc0=mc0[1:tk],sc0=sc0[1:tk],
              mt0=mt0[1:tk],st0=st0[1:tk],
              mc1=mc1[1:tk],sc1=sc1[1:tk],
              mt1=mt1[1:tk],st1=st1[1:tk])
ma1inits<-list(
  init1=list(mu0=154,dc=-5,dt=-11,tau=0.08,
             u=seq(from=-2/sqrt(0.08),
                   to=2/sqrt(0.08),
                   length.out=tk)),
  init2=list(mu0=158,dc=-3,dt=-14,tau=0.1,
             u=seq(from=-2/sqrt(0.1),
                   to=2/sqrt(0.1),
                   length.out=tk)))
rd1fits[[1]]<-stan(model_code=ma1model,data=ma1data,iter=staniter,
                   init=ma1inits,chains=2)
rd1fits[[2]]<-fit1 # This already exists

for (td in 4:8) {
  tk<-t-td
  ma2data<-list(ntrial=t,npart=armsize,nm=tk,nnom=td,threshold=threshold,
                mc0=mc0,sc0=sc0,
                mt0=mt0,st0=st0,
                mc1=mc1[1:tk],sc1=sc1,rc1=rc1[(tk+1):t],
                mt1=mt1[1:tk],st1=st1,rt1=rt1[(tk+1):t])
  rd2fits[[td-1]]<-stan(model_code=ma2model,data=ma2data,
                        iter=staniter,init=ma2inits,chains=2)
  
  ma1data<-list(armsize=armsize,tk=tk,
                mc0=mc0[1:tk],sc0=sc0[1:tk],
                mt0=mt0[1:tk],st0=st0[1:tk],
                mc1=mc1[1:tk],sc1=sc1[1:tk],
                mt1=mt1[1:tk],st1=st1[1:tk])
  ma1inits<-list(
    init1=list(mu0=154,dc=-5,dt=-11,tau=0.08,
               u=seq(from=-2/sqrt(0.08),
                     to=2/sqrt(0.08),
                     length.out=tk)),
    init2=list(mu0=158,dc=-3,dt=-14,tau=0.1,
               u=seq(from=-2/sqrt(0.1),
                     to=2/sqrt(0.1),
                     length.out=tk)))
  rd1fits[[td-1]]<-stan(model_code=ma1model,data=ma1data,
                        iter=staniter,init=ma1inits,chains=2)
}


# remove the preliminary data
rm(armsize)
rm(threshold)
rm(effectsize)
rm(staniter)


# ... and we're ready to start

simrr<-function(simID=1,nsims=1,armsize=30,threshold=0.9,effectsize=(-10),
                sim1fits=rd1fits,sim2fits=rd2fits,
                seed=740903,staniter=2000,stanchains=2) {
  #monitorfile<-file(monitorfilename)
  # store truedts
  truedts<-rep(NA,nsims)
  # and compdts (change in studies reporting complete means)
  compdts<-matrix(NA,nrow=nsims,ncol=7)
  
  # store estimate and CI for: D-L, complete D-L, Bayes baseline-only, Bayes latent, Cox-Snell
  # in these 8 scenarios (the first will not contain bayes.*)
  sims<-array(NA,dim=c(nsims,15,8))
  
  for (k in 1:nsims) {
    u0<-rnorm(t,mean=true.mean0,sd=true.tau) # random intercepts
    mc0<-mt0<-mc1<-mt1<-sc0<-st0<-
      sc1<-st1<-rc1<-rt1<-
      pooled_sd<-smd<-sesmd<-logor<-selogor<-cs<-secs<-rep(NA,2*t)
    for (i in 1:t) { 
      temp.c<-mvrnorm(n=armsize,mu=c(u0[i],u0[i]+true.time),
                      Sigma=matrix(c(true.sigma^2,true.corr*(true.sigma^2),
                                     true.corr*(true.sigma^2),true.sigma^2),
                                   nrow=2))
      temp.t<-mvrnorm(n=armsize,mu=c(u0[i],u0[i]+true.time+true.delta),
                      Sigma=matrix(c(true.sigma^2,true.corr*(true.sigma^2),
                                     true.corr*(true.sigma^2),true.sigma^2),
                                   nrow=2))
      mc0[i]<-mean(temp.c[,1])
      sc0[i]<-sd(temp.c[,1])
      mt0[i]<-mean(temp.t[,1])
      st0[i]<-sd(temp.t[,1])
      mc1[i]<-mean(temp.c[,2])
      sc1[i]<-sd(temp.c[,2])
      rc1[i]<- sum((temp.c[,2]/temp.c[,1])<=0.9)
      mt1[i]<-mean(temp.t[,2])
      st1[i]<-sd(temp.t[,2])
      rt1[i]<- sum((temp.t[,2]/temp.t[,1])<=0.9)
      
      pooled_sd[i]<-sqrt((((armsize-1)*(sc1[i]^2))+
                            ((armsize-1)*(st1[i]^2)))/
                           ((2*armsize)-2))
      smd[i]<-((mt1[i]-mc1[i])/pooled_sd[i])*
        (1-(3/((4*2*armsize)-9)))
      sesmd[i]<-sqrt(((armsize+armsize)/(armsize*armsize))+
                       ((smd[i]^2)/((2*2*armsize)/3.94)))
      
      # logor[i]<-log(((armsize-rt1[i])*rc1[i])/(rt1[i]*(armsize-rc1[i])))
      # # OR flipped to match sign of SMD
      # selogor[i]<-sqrt((1/rt1[i])+(1/rc1[i])+(1/(armsize-rt1[i]))+(1/(armsize-rc1[i])))
      # Cox-Snell approximate conversions
      # disabled: no longer used:
      # cs[i]<-0.6061*logor[i]
      # secs[i]<-0.6061*selogor[i]
    }
    # drop any infinite logors and reduce to t trials
    # keepit<-is.finite(logor)&is.finite(selogor)
    keepit <- 1:t
    # logor<-logor[keepit][1:t]
    # selogor<-selogor[keepit][1:t]
    smd<-smd[keepit][1:t]
    sesmd<-sesmd[keepit][1:t]
    mc0<-mc0[keepit][1:t]
    mt0<-mt0[keepit][1:t]
    mc1<-mc1[keepit][1:t]
    mt1<-mt1[keepit][1:t]
    sc0<-sc0[keepit][1:t]
    st0<-st0[keepit][1:t]
    sc1<-sc1[keepit][1:t]
    st1<-st1[keepit][1:t]
    rc1<-rc1[keepit][1:t]
    rt1<-rt1[keepit][1:t]
    rm(temp.c)
    rm(temp.t)
    true.dt<-mean(mt1-mc1)
    truedts[k]<-true.dt ####### write to truedts #######
    
    # First, scenario one where all trials report means
    complete.ma1<-metacont(n.e=rep(armsize,t),mean.e=mt1,sd.e=st1,
                           n.c=rep(armsize,t),mean.c=mc1,sd.c=sc1,
                           comb.fixed=FALSE,comb.random=TRUE)
    sims[k,1,1]<-complete.ma1$TE.random ####### write D-L to sims #######
    sims[k,2,1]<-complete.ma1$lower.random ####### write D-L to sims #######
    sims[k,3,1]<-complete.ma1$upper.random ####### write D-L to sims #######
    
    # Now, the Bayes scenarios
    for (td in 2:8) { 
      cat(paste("Now running block",simID,"and sim number",k,"with td=",td,"\n"))
      tk<-t-td
      compdts[k,(td-1)]<-mean(mt1[1:tk]-mc1[1:tk]) ####### write to compdts #######
      # those missing means are at the end
      ma2data<-list(ntrial=t,npart=armsize,nm=tk,nnom=td,threshold=threshold,
                    mc0=mc0,sc0=sc0,
                    mt0=mt0,st0=st0,
                    mc1=mc1[1:tk],sc1=sc1,rc1=rc1[(tk+1):t],
                    mt1=mt1[1:tk],st1=st1,rt1=rt1[(tk+1):t])
      ma1data<-list(armsize=armsize,tk=tk,
                    mc0=mc0[1:tk],sc0=sc0[1:tk],
                    mt0=mt0[1:tk],st0=st0[1:tk],
                    mc1=mc1[1:tk],sc1=sc1[1:tk],
                    mt1=mt1[1:tk],st1=st1[1:tk])
      ma1inits<-list(init1=list(mu0=154,dc=-5,dt=-11,tau=0.08,
                                u=seq(from=-2/sqrt(0.08),
                                      to=2/sqrt(0.08),
                                      length.out=tk)),
                     init2=list(mu0=158,dc=-3,dt=-14,tau=0.1,
                                u=seq(from=-2/sqrt(0.1),
                                      to=2/sqrt(0.1),
                                      length.out=tk)))
      new2fit<-invisible(stan(fit=sim2fits[[td-1]],data=ma2data,
                              iter=staniter,chains=stanchains,refresh=0))
      new1fit<-invisible(stan(fit=sim1fits[[td-1]],data=ma1data,
                              iter=staniter,chains=stanchains,refresh=0))
      newma<-metacont(n.e=rep(armsize,tk),mean.e=mt1[1:tk],sd.e=st1[1:tk],
                      n.c=rep(armsize,tk),mean.c=mc1[1:tk],sd.c=sc1[1:tk],
                      comb.fixed=FALSE,comb.random=TRUE)
      compma<-metacont(n.e=rep(armsize,t),mean.e=mt1,sd.e=st1,
                       n.c=rep(armsize,t),mean.c=mc1,sd.c=sc1,
                       comb.fixed=FALSE,comb.random=TRUE)
      
      # SMD approximate MA using Cox-Snell
      # disabled: no longer used
      # csma<-metacont(n.e=rep(armsize,t),
      #                mean.e=c(mt1[1:tk],cs[(tk+1):t]),
      #                sd.e=c(st1[1:tk],(sqrt(armsize)*secs[(tk+1):t])),
      #                n.c=rep(armsize,t),
      #                mean.c=c(mc1[1:tk],rep(0,(t-tk))),
      #                sd.c=c(sc1[1:tk],(sqrt(armsize)*secs[(tk+1):t])),
      #                sm="SMD",
      #                comb.fixed=FALSE,comb.random=TRUE) # note assumption equal SDs across arms
      
      sims[k,1,td]<-newma$TE.random ####### write to sims #######
      sims[k,2,td]<-newma$lower.random ####### write to sims #######
      sims[k,3,td]<-newma$upper.random ####### write to sims #######
      sims[k,4,td]<-compma$TE.random ####### write to sims #######
      sims[k,5,td]<-compma$lower.random ####### write to sims #######
      sims[k,6,td]<-compma$upper.random ####### write to sims #######
      diffrow<-dim(summary(new1fit)$summary)[1]-1
      sims[k,7,td]<-summary(new1fit)$summary[diffrow,1] ####### write to sims #######
      sims[k,8,td]<-summary(new1fit)$summary[diffrow,4] ####### write to sims #######
      sims[k,9,td]<-summary(new1fit)$summary[diffrow,8] ####### write to sims #######
      diffrow<-dim(summary(new2fit)$summary)[1]-1
      sims[k,10,td]<-summary(new2fit)$summary[diffrow,1] ####### write to sims #######
      sims[k,11,td]<-summary(new2fit)$summary[diffrow,4] ####### write to sims #######
      sims[k,12,td]<-summary(new2fit)$summary[diffrow,8] ####### write to sims #######
      # disabled: no longer used
      # sims[k,13,td]<-csma$TE.random ####### write to sims #######
      # sims[k,14,td]<-csma$lower.random ####### write to sims #######
      # sims[k,15,td]<-csma$upper.random ####### write to sims #######
      
    }
  }
  return(sims)
}


# test run
# rr_sims <- simrr(nsims=10)
# rr_sims <- lapply(1:3,simrr,nsims=10)
# system.time(mclapply(X=1:3, FUN=simrr, nsims=3))

# parallel run
start_time <- Sys.time()
rr_sims <- mclapply(X=1:6, FUN=simrr, nsims=2000, armsize=100, mc.cores=6)
end_time <- Sys.time()
save.image('relative-ratio-armsize100.RData')
