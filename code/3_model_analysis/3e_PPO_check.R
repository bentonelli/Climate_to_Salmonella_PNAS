
library(MCMCvis)
library(dplyr)
library(scoringRules)
library(ggplot2)
library(tidyr)

full_data <- readRDS("data/3_model_in/model_full_data_2024.rds")
full_model <- readRDS("data/4_model_out/mdl_fit.rds")

set.seed(1234)
iter_length <- 12000

mu_alpha <- rnorm(iter_length,0,4)
sigma_alpha <- rgamma(iter_length,shape=2,rate=2)

nu <- rnorm(iter_length,0,1)
theta <- rnorm(iter_length,0,1)

upsilon <- rnorm(iter_length,1,1)

gamma1 <- rnorm(iter_length,0,1);
mu_omega1 <- rnorm(iter_length,0,3);
mu_eta1 <- rnorm(iter_length,0,3);

gamma2 <- rnorm(iter_length,0,1);
mu_omega2 <- rnorm(iter_length,0,3);
mu_eta2 <- rnorm(iter_length,0,3);

sigma_omega1 <- rgamma(iter_length,3,3);
sigma_eta1 <- rgamma(iter_length,3,3);

sigma_omega2 <- rgamma(iter_length,3,3);
sigma_eta2 <- rgamma(iter_length,3,3);

sigma_phi <- rnorm(iter_length,1,0.5); 
zeta1 <- rnorm(iter_length,1,2);
zeta2 <- rnorm(iter_length,1,2);

psi1 <- rnorm(iter_length,0,5);
psi2 <- rnorm(iter_length,0,5);

kappa1 <- rnorm(iter_length,0,3);
kappa2 <- rnorm(iter_length,0,3);

beta1 <- rnorm(iter_length,0,3);
beta2 <- rnorm(iter_length,0,3);

lambda1 <- rexp(iter_length,0.05); 
lambda2 <- rexp(iter_length,0.05);

#species level sigmas (8 for each region)
sigmas <- rgamma(iter_length,shape=3,rate=2)
sigma_rep <- c()
for (ii in 1:16){
  sigma_rep <- cbind(sigma_rep,sigmas)
}

PR <- cbind(mu_alpha,sigma_alpha,nu,theta,
            upsilon,gamma1,mu_omega1,mu_eta1,
            gamma2,mu_omega2,mu_eta2,
            sigma_omega1,sigma_eta1,sigma_omega2,sigma_eta2,
            sigma_phi,
            zeta1,zeta2,
            psi1,psi2,
            kappa1,kappa2,beta1,beta2,lambda1,lambda2,sigma_rep)

MCMCtrace(full_model,params=c("mu_alpha","sigma_alpha","nu","theta",
                              "upsilon","gamma1","mu_omega1","mu_eta1",
                              "gamma2","mu_omega2","mu_eta2",
                              "sigma_omega1","sigma_eta1","sigma_omega2","sigma_eta2",
                              'sigma_phi',"zeta1",'zeta2',
                              'psi1','psi2',
                              "kappa1","kappa2",
                              "beta1",'beta2',"lambda1","lambda2",
                              "sigma1","sigma2"),priors = PR,PPO_out = TRUE,post_zm = FALSE)
