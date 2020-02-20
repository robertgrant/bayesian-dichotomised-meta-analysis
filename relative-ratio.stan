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
  for (i in 1:ntrial) sec0[i]<-sc0[i]/sqrt(npart); // standard errors
  for (i in 1:ntrial) set0[i]<-st0[i]/sqrt(npart); // standard errors
  for (i in 1:ntrial) set1[i]<-st1[i]/sqrt(npart); // standard errors
  for (i in 1:ntrial) sec1[i]<-sc1[i]/sqrt(npart); // standard errors
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
  for (i in 1:nnom) riskt1[i]<-Phi_approx((-(((-0.668965)*(mu0+u[nm+i]))+(0.743294*(mu0+u[nm+i]+dt))))/sqrt((0.447561*st0[nm+i]*st0[nm+i])+(2*(-0.668965)*0.743294*0.7*st0[nm+i]*st1[nm+i])+(0.5524862*st1[nm+i]*st1[nm+i])));
  for (i in 1:nnom) riskc1[i]<-Phi_approx((-(((-0.668965)*(mu0+u[nm+i]))+(0.743294*(mu0+u[nm+i]+dc))))/sqrt((0.447561*sc0[nm+i]*sc0[nm+i])+(2*(-0.668965)*0.743294*0.7*sc0[nm+i]*sc1[nm+i])+(0.5524862*sc1[nm+i]*sc1[nm+i])));
  diff<-(dt-dc);
}
model {
  for (i in 1:ntrial) u[i] ~ normal(0,sdtau);
  for (i in 1:ntrial) mc0[i] ~ normal(mu0+u[i],sec0[i]);
  for (i in 1:ntrial) mt0[i] ~ normal(mu0+u[i],set0[i]);
  for (i in 1:nm) mc1[i] ~ normal(mu0+u[i]+dc,sec1[i]);
  for (i in 1:nm) mt1[i] ~ normal(mu0+u[i]+dt,set1[i]);
  for (i in 1:nnom) rc1[i] ~ binomial(npart,riskc1[i]);
  for (i in 1:nnom) rt1[i] ~ binomial(npart,riskt1[i]);
  mu0 ~ normal(150,20);
  dc ~ normal(-5,10);
  dt ~ normal(-13,10);
  tau ~ gamma(0.01,0.01);
}
