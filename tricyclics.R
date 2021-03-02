# R and Stan code accompanying the paper "A taxonomy of thresholds 
# used to dichotomise outcomes, and their inclusion in Bayesian meta-analysis"
#
# Case study application to the Cochrane review of tricyclic
# anti-depressants in children
#

library(meta)
library(rstan)
options(mc.cores=4)
rstan_options(auto_write=TRUE)

# All stats
label<-c("Bernstein1990",
         "Bernstein2000",
         "Birmaher1998BDI",
         "Birmaher1998HDRS",
         "Geller8992CDRS",
         "Geller8992KSADSP",
         "Geller1990",
         "Hughes1990",
         "Kashani1984",
         "Keller2001",
         "Klein1998",
         "Kramer1981",
         "Kutcher1994",
         "Kye1996",
         "Petti1982",
         "PuigAntich1987")
study_id<-c(1,2,3,3,4,4,5,6,7,8,9,10,11,12,13,14)
# N used in mean in each group
n_mean_t<-c(9,31,13,13,26,26,12,13,5,88,18,10,17,12,3,16)
n_mean_c<-c(7,32,14,14,24,24,19,14,4,85,18,10,25,10,3,22)
# baseline mean in placebo group:
m0_c<-c(36.5,52.5,24.1,21.9,
        49.6,3.89,51.4,NA,
        27.75,18.97,21.33,25.2,
        23.77,15.3,20.67,3.0)
# baseline SD in placebo group:
s0_c<-c(3.4,10.8,9.8,9.11,
        4.6,0.5,3.7,NA,
        4.35,4.1,5.2,0.8,
        5.31,6.0,22.3,0.66)
# baseline mean in tricyclics group:
m0_t<-c(41.0,46.8,29.3,22.9,
        49.9,3.98,51.3,NA,
        25.6,18.11,21.44,16.2,
        22.63,12.3,41.33,3.1)
# baseline SD in tricyclics group:
s0_t<-c(12.2,9.5,12.6,5.11,
        4.2,0.42,4.4,NA,
        2.41,4.17,3.7,1.2,
        5.17,5.6,19.3,0.43)
# endpoint mean in placebo group:
m_c<-c(30.1,45.7,10.1,8.6,
       32.0,2.2,37.8,NA,
       23.75,9.88,14.61,22.7,
       13.42,7.8,15.33,1.9)
# endpoint SD in placebo group:
s_c<-c(NA,16.5,11.1,11.5,
       9.8,0.78,9.1,NA,
       5.56,7.74,2.1,0.7,
       8.43,5.5,15.04,0.86)
# endpoint mean in tricyclics group:
m_t<-c(29.5,34.6,10.1,7.7,
       32.9,2.41,34.7,NA,
       18.4,9.2,10.83,10.3,
       12.68,4.7,7.0,1.9)
# endpoint SD in tricyclics group:
s_t<-c(NA,8.9,11.8,8.0,
       11.4,0.81,7.8,NA,
       3.29,7.85,2.1,0.4,
       8.68,6.3,8.19,0.68)
# N used in responders (risk) in placebo group:
n_risk_c<-c(NA,NA,14,14,
            24,NA,19,14,
            NA,87,18,NA,
            25,10,NA,NA)
# responders in placebo group:
r_c<-c(NA,NA,5,5,
       4,NA,4,7,
       NA,40,9,NA,
       9,9,NA,NA)
# N used in responders (risk) in tricyclics group:
n_risk_t<-c(NA,NA,13,13,
            26,NA,11,13,
            NA,94,18,NA,
            17,12,NA,NA)
# responders in tricyclics groups:
r_t<-c(NA,NA,5,5,
       8,NA,1,6,
       NA,47,13,NA,
       8,11,NA,NA)
# row numbers for certain conditions of dichotomisation:
dichot_only<-8
dichot_and_mean<-c(3,4,5,7,10,11,13,14)
# row numbers for missing stats
missing_s_endpoint<-1
# outcomes (1=CDRS, 2=HDRS, 3=BID, 4=BDI, 5=K-SADS-P, 6=DACL)
outcomes<-c(1,1,4,2,
            1,5,1,1,
            3,2,2,6,
            2,2,3,5)
# assemble data frame
datamat<-data.frame(label,
                    study_id,
                    model.matrix(~0+as.factor(outcomes)),
                    n_mean_c,n_mean_t,
                    m0_c,s0_c,m0_t,s0_t,
                    m_c,s_c,m_t,s_t,
                    n_risk_c,r_c,
                    n_risk_t,r_t)
colnames(datamat)[3:8]<-c("CDRS","HDRS","BID","BDI","KSADSP","DACL")

# imputation of Bernstein 1990's endpoint SDs.
# plot(c(log(as.numeric(datamat[-c(1,8),'s0_c'])),
#        log(as.numeric(datamat[-c(1,8),'s0_t']))),
#      c(log(as.numeric(datamat[-c(1,8),'s_c'])),
#        log(as.numeric(datamat[-c(1,8),'s_t']))),
#      xlab='baseline log-SDs',ylab='endpoint log-SDs')
# abline(a=0,b=1)
# abline(v=log(as.numeric(datamat[1,'s0_c'])))
# abline(v=log(as.numeric(datamat[1,'s0_t'])))
# lmsd<-lm(c(log(as.numeric(datamat[-c(1,8),'s_c'])),
#            log(as.numeric(datamat[-c(1,8),'s_t']))) ~
#            c(log(as.numeric(datamat[-c(1,8),'s0_c'])),
#              log(as.numeric(datamat[-c(1,8),'s0_t']))))
# abline(a=lmsd$coefficients[1],b=lmsd$coefficients[2],lty=3)
# so we single-impute with the baseline SDs
datamat[missing_s_endpoint,'s_c']<-
  datamat[missing_s_endpoint,'s0_c']
datamat[missing_s_endpoint,'s_t']<-
  datamat[missing_s_endpoint,'s0_t']

# exploratory dataviz
# viz<-datamat[order(outcomes),]
# plot(3*(1:NROW(viz)),viz$m0_c,ylim=c(0,80))
# points(3*(1:NROW(viz)),viz$m_c,pch=19)
# points(3*(1:NROW(viz))+0.6,viz$m0_t,col='#b62525')
# points(3*(1:NROW(viz))+0.6,viz$m_t,pch=19,col='#b62525')

stancode<-
'
data {
  int n_m; // number of rows contributing means to likelihood
  int n_p; // number of rows contributing proportions to likelihood
  int n; // n_p+n_m, total rows
  int m; // number of studies
  int study_m[n_m]; // study ID (for random effect)
  int study_p[n_p]; // study ID (for random effect)
  real phi; // rotation angle for relative ratio dichotomisation (radians)
  int n_mean_c[n_m]; // n in control group
  int n_mean_t[n_m]; // n in treatment group
  real m0_c[n_m]; // baseline mean in control
  real m_c[n_m]; // endpoint mean in control
  real m0_t[n_m]; // baseline mean in treatment
  real m_t[n_m]; // endpoint mean in treatment
  real s0_c[n_m]; // baseline SD in control
  real s_c[n_m]; // endpoint SD in control
  real s0_t[n_m]; // baseline SD in treatment
  real s_t[n_m]; // endpoint SD in treatment
  int n_risk_c[n_p]; // n in control group
  int n_risk_t[n_p]; // n in treatment group
  int r_c[n_p]; // responders in control group
  int r_t[n_p]; // responders in treatment group
  real cdrs_m[n_m]; // indicator variable for CDRS
  real hdrs_m[n_m]; // indicator variable for HDRS
  real bdi_m[n_m]; // indicator variable for BDI
  real bid_m[n_m]; // indicator variable for BID
  real ksadsp_m[n_m]; // indicator variable for K-SADS-P
  real dacl_m[n_m]; // indicator variable for DACL
  real cdrs_p[n_p]; // indicator variable for CDRS
  real hdrs_p[n_p]; // indicator variable for HDRS
  real bdi_p[n_p]; // indicator variable for BDI
  real bid_p[n_p]; // indicator variable for BID
  real ksadsp_p[n_p]; // indicator variable for K-SADS-P
  real dacl_p[n_p]; // indicator variable for DACL
}

transformed data {
  real se0_c[n_m];
  real se0_t[n_m];
  real se_c[n_m];
  real se_t[n_m];

  // make standard errors
  for(i in 1:n_m){
    se0_c[i] = s0_c[i] / sqrt(n_mean_c[i]);
    se0_t[i] = s0_t[i] / sqrt(n_mean_t[i]);
    se_c[i] = s_c[i] / sqrt(n_mean_c[i]);
    se_t[i] = s_t[i] / sqrt(n_mean_t[i]);
  }
}

parameters {
  real<lower=0,upper=1> rho; // correlation
  real mu0; // mean at baseline, CRDS scale
  real mu_c; // mean at endpoint in control group
  real mu_t; // mean at endpoint in treatment group
  real<lower=0> tau; // heterogeneity SD
  real gamma_hdrs; // scaling from CRDS
  real gamma_bid; // scaling from CRDS
  real gamma_bdi; // scaling from CRDS
  real gamma_ksadsp; // scaling from CRDS
  real gamma_dacl; // scaling from CRDS
  real u[m]; // random effect (CRDS scale)
  real<lower=0> s_p[n_p]; // unknown SD in dichotomised studies
}

transformed parameters {
  real delta; // time effect
  real theta; // treatment effect
  // means for each row (trials reporting means):
  real mu_m_0[n_m]; 
  real mu_m_c[n_m]; 
  real mu_m_t[n_m];
  // means for each row (trials reporting proportions):
  real mu_p_0[n_p];
  real mu_p_c[n_p];
  real mu_p_t[n_p];
  // probability of being a responder:
  real p_c[n_p];
  real p_t[n_p];

  delta = mu_c-mu0;
  theta = mu_t-mu_c;
  for(i in 1:n_m) {
    mu_m_0[i] = (cdrs_m[i] +
                 (hdrs_m[i] * gamma_hdrs) +
                 (bid_m[i] * gamma_bid) +
                 (bdi_m[i] * gamma_bdi) +
                 (ksadsp_m[i] * gamma_ksadsp) +
                 (dacl_m[i] * gamma_dacl)) *
                (mu0 + u[study_m[i]]);
    mu_m_c[i] = (cdrs_m[i] + 
                  (hdrs_m[i] * gamma_hdrs) +
               (bid_m[i] * gamma_bid) +
               (bdi_m[i] * gamma_bdi) +
               (ksadsp_m[i] * gamma_ksadsp) +
               (dacl_m[i] * gamma_dacl)) *
              (mu_c + u[study_m[i]]);
    mu_m_t[i] = (cdrs_m[i] + 
                  (hdrs_m[i] * gamma_hdrs) +
               (bid_m[i] * gamma_bid) +
               (bdi_m[i] * gamma_bdi) +
               (ksadsp_m[i] * gamma_ksadsp) +
               (dacl_m[i] * gamma_dacl)) *
              (mu_t + u[study_m[i]]);
  }

  for(i in 1:n_p) {
    // predicted means for dichotomised trials
    mu_p_0[i] = (cdrs_p[i] + 
                  (hdrs_p[i] * gamma_hdrs) +
               (bid_p[i] * gamma_bid) +
               (bdi_p[i] * gamma_bdi) +
               (ksadsp_p[i] * gamma_ksadsp) +
               (dacl_p[i] * gamma_dacl)) *
              (mu0 + u[study_p[i]]);
    mu_p_c[i] = (cdrs_p[i] + 
                  (hdrs_p[i] * gamma_hdrs) +
               (bid_p[i] * gamma_bid) +
               (bdi_p[i] * gamma_bdi) +
               (ksadsp_p[i] * gamma_ksadsp) +
               (dacl_p[i] * gamma_dacl)) *
              (mu_c + u[study_p[i]]);
    mu_p_t[i] = (cdrs_p[i] + 
                  (hdrs_p[i] * gamma_hdrs) +
               (bid_p[i] * gamma_bid) +
               (bdi_p[i] * gamma_bdi) +
               (ksadsp_p[i] * gamma_ksadsp) +
               (dacl_p[i] * gamma_dacl)) *
              (mu_t + u[study_p[i]]);

    // probability of being a "responder" in each row (relative-ratio):
    p_c[i] = Phi_approx(
              (-cos(phi)*mu_p_c[i] - sin(phi)*mu_p_0[i]) / 
              sqrt((cos(phi)*s_p[i])^2 + 
                    (sin(phi)*s_p[i])^2 - 
              2.0*sin(phi)*cos(phi)*rho*s_p[i]*s_p[i]));
    p_t[i] = Phi_approx(
              (-cos(phi)*mu_p_t[i] - sin(phi)*mu_p_0[i]) / 
              sqrt((cos(phi)*s_p[i])^2 + 
                    (sin(phi)*s_p[i])^2 - 
              2.0*sin(phi)*cos(phi)*rho*s_p[i]*s_p[i]));
  }
}

model {
  // model priors:
  rho ~ beta(4,4);
  mu0 ~ normal(45,10);
  mu_c ~ normal(40,10);
  mu_t ~ normal(32,10);
  tau ~ cauchy(0,5);
  s_p ~ chi_square(8);
  gamma_hdrs ~ uniform(0,5);
  gamma_bid ~ uniform(0,5);
  gamma_bdi ~ uniform(0,5);
  gamma_ksadsp ~ uniform(0,5);
  gamma_dacl ~ uniform(0,5);
  
  // random effect:
  u ~ normal(0,tau);

  // likelihoods for means
  for(i in 1:n_m) {
    m0_c[i] ~ normal(mu_m_0[i],se0_c[i]);
    m0_t[i] ~ normal(mu_m_0[i],se0_t[i]);
    m_c[i] ~ normal(mu_m_c[i],se_c[i]);
    m_t[i] ~ normal(mu_m_t[i],se_t[i]);
  }
  // likelihoods for proportions
  for(i in 1:n_p){
    r_c[i] ~ binomial(n_risk_c[i],p_c[i]);
    r_t[i] ~ binomial(n_risk_t[i],p_t[i]);
  }
  /*
  For simplicity:
  * we do not subtract the bottom of each scale
  * we treat SDs as known perfectly
  * we single-impute the missing endpoint SDs in Bernstein 1990 by LOCF
  * we impute one SD for all times and arms in each dichotomised study
  */
}
'

stanmodel <- stan_model(model_code=stancode)
# for checking the mean part of the model
#datamat<-dplyr::filter(datamat,label!="Hughes1990")
mm<-dplyr::filter(datamat,label!="Hughes1990")
pm<-dplyr::filter(datamat,label=="Hughes1990")
stanfit1 <- sampling(stanmodel,
                     data=list(n_m=NROW(mm),
                               n_p=NROW(pm),
                               n=NROW(mm)+NROW(pm),
                               m=14,
                               study_m=mm$study_id,
                               study_p=as.array(pm$study_id),
                               n_mean_c=mm$n_mean_c,
                               n_mean_t=mm$n_mean_t,
                               m0_c=mm$m0_c,
                               m_c=mm$m_c,
                               m0_t=mm$m0_t,
                               m_t=mm$m_t,
                               s0_c=mm$s0_c,
                               s_c=mm$s_c,
                               s0_t=mm$s0_t,
                               s_t=mm$s_t,
                               cdrs_m=mm$CDRS,
                               hdrs_m=mm$HDRS,
                               bdi_m=mm$BDI,
                               bid_m=mm$BID,
                               ksadsp_m=mm$KSADSP,
                               dacl_m=mm$DACL,
                               cdrs_p=as.array(pm$CDRS),
                               hdrs_p=as.array(pm$HDRS),
                               bdi_p=as.array(pm$BDI),
                               bid_p=as.array(pm$BID),
                               ksadsp_p=as.array(pm$KSADSP),
                               dacl_p=as.array(pm$DACL),
                               phi=atan(-0.5),
                               n_risk_c=as.array(pm$n_risk_c),
                               n_risk_t=as.array(pm$n_risk_t),
                               r_c=as.array(pm$r_c),
                               r_t=as.array(pm$r_t)),
                      chains=3,
                      cores=3,
                      iter=5000,
                      warmup=1000,
                      seed=313116)

stansumm1<-summary(stanfit1)$summary
# checking outputs:
traceplot(stanfit1,pars=c('mu0',
                          'delta',
                          'theta'))
traceplot(stanfit1,pars='u')
traceplot(stanfit1,pars='s_p')
traceplot(stanfit1,pars='p_c')
traceplot(stanfit1,pars='p_t')
traceplot(stanfit1,pars=c('gamma_hdrs',
                          'gamma_bdi',
                          'gamma_bid',
                          'gamma_ksadsp',
                          'gamma_dacl'))
temp<-extract(stanfit1,permuted=FALSE)
pairs(temp[,1,c(1:10,25,26)],cex=0.2,col='#00000020')

# output stats: 
stansumm1[c(1,2,5,6,7,8,9,10,25,26,27),c(1,4,8)]





################################################
# analysis 2, including Geller I & II as dichotomised
mm<-dplyr::filter(datamat,
                  label!="Hughes1990" &
                    label!="Geller8992CDRS" &
                    label!="Geller1990")
pm<-dplyr::filter(datamat,label=="Hughes1990" |
                    label=="Geller8992CDRS" |
                    label=="Geller1990")
rinits<-function(n_p) {
  list(rho=runif(1,0.1,0.9),
       mu0=runif(1,40,50),
       mu_c=runif(1,35,45),
       mu_t=runif(1,30,40),
       tau=runif(1,5,10),
       gamma_hdrs=runif(1,0.2,0.8),
       gamma_bid=runif(1,0.2,0.8),
       gamma_bdi=runif(1,0.2,0.8),
       gamma_ksadsp=runif(1,0.2,0.8),
       gamma_dacl=runif(1,0.2,0.8),
       s_p=runif(n_p,5,15))
}
set.seed(718)
stanfit2 <- sampling(stanmodel,
                     data=list(n_m=NROW(mm),
                               n_p=NROW(pm),
                               n=NROW(mm)+NROW(pm),
                               m=14,
                               study_m=mm$study_id,
                               study_p=as.array(pm$study_id),
                               n_mean_c=mm$n_mean_c,
                               n_mean_t=mm$n_mean_t,
                               m0_c=mm$m0_c,
                               m_c=mm$m_c,
                               m0_t=mm$m0_t,
                               m_t=mm$m_t,
                               s0_c=mm$s0_c,
                               s_c=mm$s_c,
                               s0_t=mm$s0_t,
                               s_t=mm$s_t,
                               cdrs_m=mm$CDRS,
                               hdrs_m=mm$HDRS,
                               bdi_m=mm$BDI,
                               bid_m=mm$BID,
                               ksadsp_m=mm$KSADSP,
                               dacl_m=mm$DACL,
                               cdrs_p=as.array(pm$CDRS),
                               hdrs_p=as.array(pm$HDRS),
                               bdi_p=as.array(pm$BDI),
                               bid_p=as.array(pm$BID),
                               ksadsp_p=as.array(pm$KSADSP),
                               dacl_p=as.array(pm$DACL),
                               phi=atan(-0.5),
                               n_risk_c=as.array(pm$n_risk_c),
                               n_risk_t=as.array(pm$n_risk_t),
                               r_c=as.array(pm$r_c),
                               r_t=as.array(pm$r_t)),
                     chains=3,
                     cores=3,
                     iter=5000,
                     warmup=1000,
                     seed=313118,
                     init=list(rinits(NROW(pm)),
                               rinits(NROW(pm)),
                               rinits(NROW(pm))))
# real u[m]; // random effect (CRDS scale)
# real<lower=0> s_p[n_p]; // unknown SD in dichotomised studies

stansumm2<-summary(stanfit2)$summary
# checking outputs:
traceplot(stanfit2,pars=c('mu0',
                          'delta',
                          'theta'))
traceplot(stanfit2,pars='u')
traceplot(stanfit2,pars='s_p')
traceplot(stanfit2,pars='p_c')
traceplot(stanfit2,pars='p_t')
traceplot(stanfit2,pars=c('gamma_hdrs',
                          'gamma_bdi',
                          'gamma_bid',
                          'gamma_ksadsp',
                          'gamma_dacl'))
temp<-extract(stanfit2,permuted=FALSE)
pairs(temp[,1,c(1:10,25,26)],cex=0.2,col='#00000020')

# output stats: 
stansumm2[c(1,2,5,6,7,8,9,10,25,26),c(1,4,8)]

# Geller posterior quantiles
gellerpars<-c("mu_p_0[1]",
              "mu_p_c[1]",
              "mu_p_t[1]",
              "mu_p_0[2]",
              "mu_p_c[2]",
              "mu_p_t[2]")
temp<-extract(stanfit2,
              permuted=TRUE,
              pars=gellerpars)

mean(temp[["mu_p_0[1]"]]<49.6)
mean(temp[["mu_p_0[1]"]]<49.9)
mean(temp[["mu_p_c[1]"]]<32.0)
mean(temp[["mu_p_t[1]"]]<32.9)

mean(temp[["mu_p_0[2]"]]<51.4)
mean(temp[["mu_p_0[2]"]]<51.3)
mean(temp[["mu_p_c[2]"]]<37.8)
mean(temp[["mu_p_t[2]"]]<34.7)


