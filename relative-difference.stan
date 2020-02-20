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
  for (i in 1:t) sec0[i] = sc0[i]/sqrt(armsize); // standard errors
  for (i in 1:t) set0[i] = st0[i]/sqrt(armsize); // standard errors
  for (i in 1:t) set1[i] = st1[i]/sqrt(armsize); // standard errors
  for (i in 1:t) sec1[i] = sc1[i]/sqrt(armsize); // standard errors
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
  for (i in 1:td) riskt1[i] = Phi_approx((threshold-dt)/(sqrt((st0[tk+i]*st0[tk+i])+(st1[tk+i]*st1[tk+i])-(2*0.7*st0[tk+i]*st1[tk+i]))));
  for (i in 1:td) riskc1[i] = Phi_approx((threshold-dc)/(sqrt((sc0[tk+i]*sc0[tk+i])+(sc1[tk+i]*sc1[tk+i])-(2*0.7*sc0[tk+i]*sc1[tk+i]))));
  diff<-(dt-dc);
}
model {
  for (i in 1:t) u[i] ~ normal(0,sdtau);
  for (i in 1:t) mc0[i] ~ normal(mu0+u[i],sec0[i]);
  for (i in 1:t) mt0[i] ~ normal(mu0+u[i],set0[i]);
  for (i in 1:tk) mc1[i] ~ normal(mu0+u[i]+dc,sec1[i]);
  for (i in 1:tk) mt1[i] ~ normal(mu0+u[i]+dt,set1[i]);
  for (i in 1:td) rc1[i] ~ binomial(armsize,riskc1[i]);
  for (i in 1:td) rt1[i] ~ binomial(armsize,riskt1[i]);
  mu0 ~ normal(150,20);
  dc ~ normal(-5,10);
  dt ~ normal(-13,10);
  tau ~ gamma(0.01,0.01);
}
