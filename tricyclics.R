# Application of Bayesian meta-analysis
# to the Cochrane review of tricyclics for children
#
#
# this version includes all dichotomised trials so that we can compare
# their parameters against their true means
#
# No adjustment for skewed / constrained outcomes is applied.
# Assumptions in the dichotomoised trial Hughes (1990):
#	correlation=0.3, SD baseline=6.6, endpoint=10.6


library(meta)
library(rstan)

md<-2 # dichotomised trials

# starting with only the CDRS trials (4):
nc_cdrs<-c(7,32,24)
nt_cdrs<-c(9,31,26)
mc0_cdrs<-c(36.5,52.5,49.6)
mt0_cdrs<-c(41,46.8,49.9)
sc0_cdrs<-c(3.4,10.8,4.6)
st0_cdrs<-c(12.2,9.5,4.2)
mc1_cdrs<-c(30.1,45.7,32)
mt1_cdrs<-c(29.5,34.6,32.9)
sc1_cdrs<-c(NA,16.5,9.8)
st1_cdrs<-c(NA,8.9,11.4)
rc1<-c(4,7)
rt1<-c(1,6)
nrc1<-c(19,14)
nrt1<-c(11,13)

# single imputation! 
# later, try data synthesis
sc1_cdrs[1]<-mean(sc1_cdrs[2:3])
st1_cdrs[1]<-mean(st1_cdrs[2:3])

# HDRS trials (5)
nc_hdrs<-c(14,85,18,25,10)
nt_hdrs<-c(13,88,18,17,12)
mc0_hdrs<-c(21.9,18.97,21.33,23.77,15.3)
mt0_hdrs<-c(22.9,18.11,21.44,22.63,12.3)
sc0_hdrs<-c(9.11,4.104,5.2,5.31,6.0)
st0_hdrs<-c(5.11,4.169,3.7,5.17,5.6)
mc1_hdrs<-c(8.6,9.88,14.61,13.42,7.8)
mt1_hdrs<-c(7.7,9.20,10.83,12.68,4.7)
sc1_hdrs<-c(11.5,7.742,2.1,8.43,5.5)
st1_hdrs<-c( 8.0,7.853,2.1,8.68,6.3)

# BID trials (2)
nc_bid<-c(4,3)
nt_bid<-c(5,3)
mc0_bid<-c(27.75,20.67)
mt0_bid<-c(25.6,41.33)
sc0_bid<-c(4.35,22.3)
st0_bid<-c(2.41,19.3)
mc1_bid<-c(23.75,15.33)
mt1_bid<-c(18.4,7)
sc1_bid<-c(5.56,15.04)
st1_bid<-c(3.29,8.19)

# K-SADS-P trial (1)
nc_ksadsp<-22
nt_ksadsp<-16
mc0_ksadsp<-3
mt0_ksadsp<-3.1
sc0_ksadsp<-0.66
st0_ksadsp<-0.43
mc1_ksadsp<-1.9
mt1_ksadsp<-1.9
sc1_ksadsp<-0.86
st1_ksadsp<-0.68

madata<-list(nc_cdrs=nc_cdrs,nt_cdrs=nt_cdrs,
		mc0_cdrs=mc0_cdrs,mt0_cdrs=mt0_cdrs,
		mc1_cdrs=mc1_cdrs,mt1_cdrs=mt1_cdrs,
		sc0_cdrs=sc0_cdrs,st0_cdrs=st0_cdrs,
		sc1_cdrs=sc1_cdrs,st1_cdrs=st1_cdrs,
		rc1=rc1,rt1=rt1,nrc1=nrc1,nrt1=nrt1,
		nc_hdrs=nc_hdrs,nt_hdrs=nt_hdrs,
		mc0_hdrs=mc0_hdrs,mt0_hdrs=mt0_hdrs,
		mc1_hdrs=mc1_hdrs,mt1_hdrs=mt1_hdrs,
		sc0_hdrs=sc0_hdrs,st0_hdrs=st0_hdrs,
		sc1_hdrs=sc1_hdrs,st1_hdrs=st1_hdrs,
		nc_bid=nc_bid,nt_bid=nt_bid,
		mc0_bid=mc0_bid,mt0_bid=mt0_bid,
		mc1_bid=mc1_bid,mt1_bid=mt1_bid,
		sc0_bid=sc0_bid,st0_bid=st0_bid,
		sc1_bid=sc1_bid,st1_bid=st1_bid,
		nc_ksadsp=nc_ksadsp,nt_ksadsp=nt_ksadsp,
		mc0_ksadsp=mc0_ksadsp,mt0_ksadsp=mt0_ksadsp,
		mc1_ksadsp=mc1_ksadsp,mt1_ksadsp=mt1_ksadsp,
		sc0_ksadsp=sc0_ksadsp,st0_ksadsp=st0_ksadsp,
		sc1_ksadsp=sc1_ksadsp,st1_ksadsp=st1_ksadsp,
		md=md)

mainits<-list(
		init1=list(mu0=50,dc=-8,dt=-15,sdu=12,
			     	u_cdrs=c(2,-5,5),
				ux=c(-1,0.4),
				u_hdrs=c(-5,-3,-1,4,7),
				u_bid=c(-4,4),
				u_ksadsp=2,
				hdrs=0.7,bid=0.6,ksadsp=0.07),
		init2=list(mu0=40,dc=-11,dt=-14,sdu=8,
			     	u_cdrs=c(-2,5,-5),
				ux=c(1,1),
				u_hdrs=c(-1,4,-3,2,5),
				u_bid=c(4,-4),
				u_ksadsp=-2,
				hdrs=0.6,bid=0.7,ksadsp=0.05))

# Stan model
mamodel<-'
	data {
		real mc0_cdrs[3];
		real mt0_cdrs[3];
		real mc1_cdrs[3];
		real mt1_cdrs[3];
		real<lower=0> sc0_cdrs[3];
		real<lower=0> st0_cdrs[3];
		real<lower=0> sc1_cdrs[3];
		real<lower=0> st1_cdrs[3];
		int<lower=0> nc_cdrs[3];
		int<lower=0> nt_cdrs[3];
		real mc0_hdrs[5];
		real mt0_hdrs[5];
		real mc1_hdrs[5];
		real mt1_hdrs[5];
		real<lower=0> sc0_hdrs[5];
		real<lower=0> st0_hdrs[5];
		real<lower=0> sc1_hdrs[5];
		real<lower=0> st1_hdrs[5];
		int<lower=0> nc_hdrs[5];
		int<lower=0> nt_hdrs[5];
		real mc0_bid[2];
		real mt0_bid[2];
		real mc1_bid[2];
		real mt1_bid[2];
		real<lower=0> sc0_bid[2];
		real<lower=0> st0_bid[2];
		real<lower=0> sc1_bid[2];
		real<lower=0> st1_bid[2];
		int<lower=0> nc_bid[2];
		int<lower=0> nt_bid[2];
		real mc0_ksadsp;
		real mt0_ksadsp;
		real mc1_ksadsp;
		real mt1_ksadsp;
		real<lower=0> sc0_ksadsp;
		real<lower=0> st0_ksadsp;
		real<lower=0> sc1_ksadsp;
		real<lower=0> st1_ksadsp;
		int<lower=0> nc_ksadsp;
		int<lower=0> nt_ksadsp;
		int<lower=0> rc1[2];
		int<lower=0> rt1[2];
		int<lower=0> nrc1[2];
		int<lower=0> nrt1[2];
	}
	transformed data {
		real<lower=0> sec0_cdrs[3];
		real<lower=0> set0_cdrs[3];
		real<lower=0> set1_cdrs[3];
		real<lower=0> sec1_cdrs[3];
		real<lower=0> sec0_hdrs[5];
		real<lower=0> set0_hdrs[5];
		real<lower=0> set1_hdrs[5];
		real<lower=0> sec1_hdrs[5];
		real<lower=0> sec0_bid[2];
		real<lower=0> set0_bid[2];
		real<lower=0> set1_bid[2];
		real<lower=0> sec1_bid[2];
		real<lower=0> sec0_ksadsp;
		real<lower=0> set0_ksadsp;
		real<lower=0> set1_ksadsp;
		real<lower=0> sec1_ksadsp;
		for (i in 1:3) sec0_cdrs[i] = sc0_cdrs[i]/sqrt(nc_cdrs[i]); 
		for (i in 1:3) set0_cdrs[i] = st0_cdrs[i]/sqrt(nt_cdrs[i]); 
		for (i in 1:3) set1_cdrs[i] = st1_cdrs[i]/sqrt(nt_cdrs[i]); 
		for (i in 1:3) sec1_cdrs[i] = sc1_cdrs[i]/sqrt(nc_cdrs[i]); 
		for (i in 1:5) sec0_hdrs[i] = sc0_hdrs[i]/sqrt(nc_hdrs[i]); 
		for (i in 1:5) set0_hdrs[i] = st0_hdrs[i]/sqrt(nt_hdrs[i]); 
		for (i in 1:5) set1_hdrs[i] = st1_hdrs[i]/sqrt(nt_hdrs[i]); 
		for (i in 1:5) sec1_hdrs[i] = sc1_hdrs[i]/sqrt(nc_hdrs[i]); 
		for (i in 1:2) sec0_bid[i] = sc0_bid[i]/sqrt(nc_bid[i]); 
		for (i in 1:2) set0_bid[i] = st0_bid[i]/sqrt(nt_bid[i]); 
		for (i in 1:2) set1_bid[i] = st1_bid[i]/sqrt(nt_bid[i]); 
		for (i in 1:2) sec1_bid[i] = sc1_bid[i]/sqrt(nc_bid[i]); 
		sec0_ksadsp = sc0_ksadsp/sqrt(nc_ksadsp); 
		set0_ksadsp = st0_ksadsp/sqrt(nt_ksadsp); 
		set1_ksadsp = st1_ksadsp/sqrt(nt_ksadsp); 
		sec1_ksadsp = sc1_ksadsp/sqrt(nc_ksadsp); 
	}
	parameters {
		real ux[2]; // random effect for the dichotomised trials
		real<lower=0> sdu;
		real mu0; // global mean for both groups at baseline
		real u_cdrs[3]; 
		real u_hdrs[5]; 
		real u_bid[2]; 
		real u_ksadsp; 
		real dc; // global mean of control change
		real dt; // global mean of treatment change
		real<lower=0> hdrs;
		real<lower=0> bid;
		real<lower=0> ksadsp;
	}
	transformed parameters {
//		real<lower=0> sdu; // heterogeneity SD
		real<lower=0, upper=1> riskt[2]; // risk of response tx
		real<lower=0, upper=1> riskc[2]; // risk of response control
		real diff; // dt-dc
		real mu_cdrs[3]; // study means for monitoring
		real mu_hdrs[5];
		real mu_bid[2];
		real mu_ksadsp;
		real mu_x[2];
		for (i in 1:2) riskt[i] = Phi_approx(-((-0.4472*(mu0+ux[i]))+(0.8944*(mu0+ux[i]+dt)))/9.0273);
		for (i in 1:2) riskc[i] = Phi_approx(-((-0.4472*(mu0+ux[i]))+(0.8944*(mu0+ux[i]+dc)))/9.0273);
		diff = (dt-dc);
		for (i in 1:3) mu_cdrs[i] = mu0+u_cdrs[i];
		for (i in 1:5) mu_hdrs[i] = hdrs*(mu0+u_hdrs[i]);
		for (i in 1:2) mu_bid[i] = bid*(mu0+u_bid[i]);
		mu_ksadsp = ksadsp*(mu0+u_ksadsp);
		for (i in 1:2) mu_x[i] = mu0+ux[i];
	}
	model {
		for (i in 1:3) u_cdrs[i] ~ normal(0,sdu);
		for (i in 1:5) u_hdrs[i] ~ normal(0,sdu);
		for (i in 1:2) u_bid[i] ~ normal(0,sdu);
		u_ksadsp ~ normal(0,sdu);
		for (i in 1:2) ux[i] ~ normal(0,sdu);
		for (i in 1:3) mc0_cdrs[i] ~ normal(mu0+u_cdrs[i],sec0_cdrs[i]);
		for (i in 1:3) mt0_cdrs[i] ~ normal(mu0+u_cdrs[i],set0_cdrs[i]);
		for (i in 1:3) mc1_cdrs[i] ~ normal(mu0+u_cdrs[i]+dc,sec1_cdrs[i]);
		for (i in 1:3) mt1_cdrs[i] ~ normal(mu0+u_cdrs[i]+dt,set1_cdrs[i]);
		for (i in 1:5) mc0_hdrs[i] ~ normal(hdrs*(mu0+u_hdrs[i]),sec0_hdrs[i]);
		for (i in 1:5) mt0_hdrs[i] ~ normal(hdrs*(mu0+u_hdrs[i]),set0_hdrs[i]);
		for (i in 1:5) mc1_hdrs[i] ~ normal(hdrs*(mu0+u_hdrs[i]+dc),sec1_hdrs[i]);
		for (i in 1:5) mt1_hdrs[i] ~ normal(hdrs*(mu0+u_hdrs[i]+dt),set1_hdrs[i]);
		for (i in 1:2) mc0_bid[i] ~ normal(bid*(mu0+u_bid[i]),sec0_bid[i]);
		for (i in 1:2) mt0_bid[i] ~ normal(bid*(mu0+u_bid[i]),set0_bid[i]);
		for (i in 1:2) mc1_bid[i] ~ normal(bid*(mu0+u_bid[i]+dc),sec1_bid[i]);
		for (i in 1:2) mt1_bid[i] ~ normal(bid*(mu0+u_bid[i]+dt),set1_bid[i]);
		mc0_ksadsp ~ normal(ksadsp*(mu0+u_ksadsp),sec0_ksadsp);
		mt0_ksadsp ~ normal(ksadsp*(mu0+u_ksadsp),set0_ksadsp);
		mc1_ksadsp ~ normal(ksadsp*(mu0+u_ksadsp+dc),sec1_ksadsp);
		mt1_ksadsp ~ normal(ksadsp*(mu0+u_ksadsp+dt),set1_ksadsp);
		for (i in 1:2) rc1[i] ~ binomial(nrc1[i],riskc[i]);
		for (i in 1:2) rt1[i] ~ binomial(nrt1[i],riskt[i]);
		mu0 ~ normal(45,9);
		dc ~ normal(-5,10);
		dt ~ normal(-13,10);
		sdu ~ cauchy(0,5);
		hdrs ~ uniform(0,5);
		bid ~ uniform(0,5);
		ksadsp ~ uniform(0,5);
	}
'
	
# work on basis that SD at baseline is 6.6 and at endpoint is 10.6
# and correlation is 0.3

# fit Stan model
fit1<-stan(model_code=mamodel,
		data=madata,
		iter=10000,
		chains=2,
		init=mainits,
		seed=434)
print(fit1)

# fit classical model
#complete.ma1<-metacont(n.e=nt,mean.e=mt1,sd.e=st1,
#			    n.c=nc,mean.c=mc1,sd.c=sc1,
#				comb.fixed=FALSE,comb.random=TRUE)
#summary(complete.ma1)

