# R and Stan code for the simulation study and example analysis in the paper 
# "A taxonomy of thresholds used to dichotomise outcomes, and their inclusion in Bayesian meta-analysis"



################################################################################################
#############################  Relative-difference threshold  ##################################
################################################################################################


library(meta)
library(MASS)
library(rstan)
library(parallel)
#options(mc.cores = parallel::detectCores())
#library(foreach)
#library(doParallel)


##################### Parameters for all 3 functions ######################


t<-10 # Number of trials
true.mean0<-155 # baseline mean
true.time<-(-5) # time effect
true.sigma<-15 # constant sigma at baseline and endpoint (intra-trial)
true.tau<-4 # random effect SD (inter-trial)
true.corr<-0.7 # correlation of baseline and endpoint



##################### Relative-difference threshold ######################


# data to get the models running
armsize<-10
threshold<-(-10)
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
  rc1[i]<-sum((temp.c[,2]-temp.c[,1])<=(threshold))
  mt1[i]<-mean(temp.t[,2])
  st1[i]<-sd(temp.t[,2])
  rt1[i]<-sum((temp.t[,2]-temp.t[,1])<=(threshold))
}
rm(temp.c)
rm(temp.t)
true.dt<-mean(mt1-mc1) # because all trials have same sample size

# Now we will run the Stan models once to get them compiled and stored in an R list


################# Stan model ##################
ma2model<-'
data {
int<lower=0> armsize; // the number of participants per arm
real threshold; // relative-difference threshold
int<lower=0> t; // total number of trials (tk+td)
int<lower=0> tk; // number of studies with means
int<lower=0> td; // number without means
real mc0[t]; // baseline means in control
real mt0[t]; // baseline means in treatment
real<lower=0> sc0[t]; // SDs in control baseline
real<lower=0> st0[t]; // SDs in treatment baseline
real mt1[tk]; // treatment means where known
real<lower=0> st1[t]; // treatment SD known for all
int<lower=0> rt1[td]; // number normotensive in tx group
real mc1[tk]; // control means where known
real<lower=0> sc1[t]; // SDs in control known for all
int<lower=0> rc1[td]; // number normotensive in control
}
transformed data {
real<lower=0> sec0[t];
real<lower=0> set0[t];
real<lower=0> set1[t];
real<lower=0> sec1[t];
for (i in 1:t)
sec0[i] = sc0[i]/sqrt(armsize); // standard errors
for (i in 1:t)
set0[i] = st0[i]/sqrt(armsize); // standard errors
for (i in 1:t)
set1[i] = st1[i]/sqrt(armsize); // standard errors
for (i in 1:t)
sec1[i] = sc1[i]/sqrt(armsize); // standard errors
}
parameters {
real mu0; // global mean for both groups at baseline
real<lower=0> tau; // heterogeneity precision
real u[t]; // random effect
real dc; // global mean of time difference
real dt; // global mean of treatment plus time difference
}
transformed parameters {
real<lower=0> sdtau; // heterogeneity SD
real<lower=0, upper=1> riskt1[td]; // risk of normotension in treatment arm
real<lower=0, upper=1> riskc1[td]; // risk of normotension in control arm
real diff; // dt-dc
sdtau = 1/sqrt(tau);
for (i in 1:td)
riskt1[i] = Phi_approx((threshold-dt)/(sqrt((st0[tk+i]*st0[tk+i])+(st1[tk+i]*st1[tk+i])-(2*0.7*st0[tk+i]*st1[tk+i]))));
for (i in 1:td)
riskc1[i] = Phi_approx((threshold-dc)/(sqrt((sc0[tk+i]*sc0[tk+i])+(sc1[tk+i]*sc1[tk+i])-(2*0.7*sc0[tk+i]*sc1[tk+i]))));
diff<-(dt-dc);
}
model {
for (i in 1:t)
u[i] ~ normal(0,sdtau);
for (i in 1:t)
mc0[i] ~ normal(mu0+u[i],sec0[i]);
for (i in 1:t)
mt0[i] ~ normal(mu0+u[i],set0[i]);
for (i in 1:tk)
mc1[i] ~ normal(mu0+u[i]+dc,sec1[i]);
for (i in 1:tk)
mt1[i] ~ normal(mu0+u[i]+dt,set1[i]);
for (i in 1:td)
rc1[i] ~ binomial(armsize,riskc1[i]);
for (i in 1:td)
rt1[i] ~ binomial(armsize,riskt1[i]);
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
ma2data<-list(t=t,armsize=armsize,tk=tk,td=td,threshold=threshold,
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
  rc1[i]<-sum((temp.c[,2]-temp.c[,1])<=(-10))
  mt1[i]<-mean(temp.t[,2])
  st1[i]<-sd(temp.t[,2])
  rt1[i]<-sum((temp.t[,2]-temp.t[,1])<=(-10))
}
rm(temp.c)
rm(temp.t)

td<-2
tk<-t-td
ma2data<-list(t=t,armsize=armsize,tk=tk,td=td,
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
  ma2data<-list(t=t,armsize=armsize,tk=tk,td=td,
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

simrd<-function(simID=1,nsims=1,armsize=30,threshold=(-10),effectsize=(-10),
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
    for (i in 1:(2*t)) { # we make double so we can drop any logor = -Inf
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
      rc1[i]<-sum((temp.c[,2]-temp.c[,1])<=(-10))
      mt1[i]<-mean(temp.t[,2])
      st1[i]<-sd(temp.t[,2])
      rt1[i]<-sum((temp.t[,2]-temp.t[,1])<=(-10))
      
      pooled_sd[i]<-sqrt((((armsize-1)*(sc1[i]^2))+
                            ((armsize-1)*(st1[i]^2)))/
                           ((2*armsize)-2))
      smd[i]<-((mt1[i]-mc1[i])/pooled_sd[i])*
        (1-(3/((4*2*armsize)-9)))
      sesmd[i]<-sqrt(((armsize+armsize)/(armsize*armsize))+
                       ((smd[i]^2)/((2*2*armsize)/3.94)))
      
      logor[i]<-log(((armsize-rt1[i])*rc1[i])/(rt1[i]*(armsize-rc1[i])))
      # OR flipped to match sign of SMD
      selogor[i]<-sqrt((1/rt1[i])+(1/rc1[i])+(1/(armsize-rt1[i]))+(1/(armsize-rc1[i])))
      # Cox-Snell approximate conversions
      cs[i]<-0.6061*logor[i]
      secs[i]<-0.6061*selogor[i]
    }
    # drop any infinite logors and reduce to t trials
    keepit<-is.finite(logor)&is.finite(selogor)
    logor<-logor[keepit][1:t]
    selogor<-selogor[keepit][1:t]
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
      ma2data<-list(t=t,armsize=armsize,tk=tk,td=td,threshold=threshold,
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
      csma<-metacont(n.e=rep(armsize,t),
                     mean.e=c(mt1[1:tk],cs[(tk+1):t]),
                     sd.e=c(st1[1:tk],(sqrt(armsize)*secs[(tk+1):t])),
                     n.c=rep(armsize,t),
                     mean.c=c(mc1[1:tk],rep(0,(t-tk))),
                     sd.c=c(sc1[1:tk],(sqrt(armsize)*secs[(tk+1):t])),
                     sm="SMD",
                     comb.fixed=FALSE,comb.random=TRUE) # note assumption equal SDs across arms
      
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
      sims[k,13,td]<-csma$TE.random ####### write to sims #######
      sims[k,14,td]<-csma$lower.random ####### write to sims #######
      sims[k,15,td]<-csma$upper.random ####### write to sims #######
      
    }
  }
  #save(sims,file=paste0("relative-difference-sims-",parallelID,".RData"))  
  return(sims)
  # note this replaces the outfile argument to simrd()
  
  # write.csv(sims,file=outfile)
  # the columns are:
  # 1)    point estimate, complete D-L, tk=10
  # 2)    lower CI, complete D-L, tk=10
  # 3)    upper CI, complete D-L, tk=10
  # 4-12) blank
  # 13-15)    point estimate, CI, D-L, tk=8
  # 16-18)    point estimate, CI, D-L, tk=10
  # 19-21)    point estimate, CI, Bayes baseline-only, tk=8
  # 12-24)    point estimate, CI, Bayes baseline & latent vars, tk=8 & td=2
  # 25-27)		point estimate, CI, Cox-Snell approximation
  # then so on for tk=7, 6,..., 2
  
  
  
  # # this stuff gets returned:
  # # bias
  # mcbias<-t(apply(sims[,c(1,4,7,10,13),],mean,MARGIN=c(2,3)))-effectsize
  # # coverage
  # mccoverage<-matrix(c(
  # 	apply((sims[,2,]<effectsize & sims[,3,]>effectsize),sum,MARGIN=2)/dim(sims)[1],
  # 	apply((sims[,5,]<effectsize & sims[,6,]>effectsize),sum,MARGIN=2)/dim(sims)[1],
  # 	apply((sims[,8,]<effectsize & sims[,9,]>effectsize),sum,MARGIN=2)/dim(sims)[1],
  # 	apply((sims[,11,]<effectsize & sims[,12,]>effectsize),sum,MARGIN=2)/dim(sims)[1],
  # 	apply((sims[,14,]<effectsize & sims[,15,]>effectsize),sum,MARGIN=2)/dim(sims)[1]),
  # 	nrow=8,ncol=4)
  # # confidence interval width
  # mcciwidth<-matrix(c(
  # 	apply(sims[,3,]-sims[,2,],mean,MARGIN=2),
  # 	apply(sims[,6,]-sims[,5,],mean,MARGIN=2),
  # 	apply(sims[,9,]-sims[,8,],mean,MARGIN=2),
  # 	apply(sims[,12,]-sims[,11,],mean,MARGIN=2),
  # 	apply(sims[,15,]-sims[,14,],mean,MARGIN=2)),
  # 	nrow=8,ncol=4)
  
#  close(monitorfile)
  
  # return(list(mcbias=mcbias,
  # 						mccoverage=mccoverage,
  # 						mcciwidth=mcciwidth))
}


# test run
#system.time(simrd(nsims=3))
#system.time(lapply(1:3,simrd,nsims=10))
#system.time(mclapply(X=1:3, FUN=simrd, nsims=10))

# parallel run
start_time <- Sys.time()
rd_sims <- mclapply(X=1:6, FUN=simrd, nsims=2000, armsize=100, mc.cores=6)
end_time <- Sys.time()
save.image('relative-difference-armsize100.RData')
