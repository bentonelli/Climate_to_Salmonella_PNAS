// Model code to test hypotheses for connections between climate, cone production,
// irruptions, and salmonellosis outbreaks.

data {
  // Shared data
  int<lower=1> Nyrs;
  int<lower=1> Ncls;
  int<lower=1> Nsp;

  //Spatial data
  int<lower=0> N_edges;
  int<lower=1, upper=Ncls> node1 [N_edges]; // node1[i] adjacent to node2[i]
  int<lower=1, upper=Ncls> node2 [N_edges]; // and node1[i] < node2[i]

  real mean_wf_val;
  real mean_ef_val;

  //Delta-t data
  int<lower=0> N_dt_yrs_obs_else;

  int<lower=0> dt_obs_ind [N_dt_yrs_obs_else];
  int<lower=0> dt_miss_ind [(Nyrs-N_dt_yrs_obs_else)];

  // Masting data
  int<lower = 0> N_obs_mast[Nyrs];    // number of non-missing for each year
  int<lower = 0> N_miss_mast[Nyrs];   // number missing for each year

  vector[Ncls] yy_obs_mast[Nyrs];    //Observed masting data

  real mean_mast_val;
  real sd_mean_mast;

  vector[Ncls] dt_obs[Nyrs]; //Delta-t matrix

  vector[Ncls] fc; //Forest cover

  vector[Ncls] er_identity [2]; // Regional identity of cells

  int<lower=0> Ncells_wf; // Number of cells in each region
  int<lower=0> Ncells_nf;

  //Irruption data
  int<lower=1> N_irr_obs_year;

  int<lower=1> irr_obs_yrs[N_irr_obs_year];
  int<lower=1> irr_missing_yrs;

  vector[Nsp] irr_ind_e[Nyrs];
  vector[Nsp] irr_ind_w[Nyrs];

  int<lower = 0> ii_obs_mast[Nyrs, Ncls];    // indices of observed data
  int<lower = 0> ii_miss_mast[Nyrs, Ncls];   // indices of missing data

  //Disease data
  int <lower=0,upper=1> outbreak_yn_west[Nyrs];
  int <lower=0,upper=1> outbreak_yn_east[Nyrs];

  vector[Nyrs] west_disease_ts;
  vector[Nyrs] east_disease_ts;

  int<lower=0> N_dis_unobserved;

  int<lower=0> disease_ts_observed[Nyrs - N_dis_unobserved];
  int<lower=0> disease_ts_unobserved[N_dis_unobserved];

  int<lower=0> dis_obs_vect[Nyrs];
}
parameters {

  // Masting (cone production) parameters
  vector<lower = 0.001>[Ncls] yy_miss_else[Nyrs];
  vector<lower=0.001>[Ncls] pre_data_mast;

  vector[Nyrs] alpha;

  real nu;
  real theta;
  vector[Ncls] phi [Nyrs];
  real<lower = 0> sigma_phi;
  real<lower = 0> upsilon;

  real mu_alpha;
  real<lower=0> sigma_alpha;

  //Irruption parameters
  vector[Nsp] irr_ind_miss_e[Nyrs];
  vector[Nsp] irr_ind_miss_w[Nyrs];

  vector[Nsp] pre_data_irr_e;
  vector[Nsp] pre_data_irr_w;

  real gamma1;
  real mu_omega1;
  real mu_eta1;
  real gamma2;
  real mu_omega2;
  real mu_eta2;
  real<lower = 0> sigma_omega1;
  real<lower = 0> sigma_eta1;
  real<lower = 0> sigma_omega2;
  real<lower = 0> sigma_eta2;
  vector[Nsp] omega1;
  vector[Nsp] omega2;

  vector[Nsp] eta1;
  vector[Nsp] eta2;

  vector<lower = 0> [Nsp] sigma1;
  vector<lower = 0> [Nsp] sigma2;

  //Correlation matrix for MVN
  corr_matrix[Nsp] Rho1;
  corr_matrix[Nsp] Rho2;

  // Disease parameters
  vector[Nyrs] dis_miss_e;
  vector[Nyrs] dis_miss_w;

  real zeta1;
  real zeta2;

  real kappa1;
  real kappa2;

  real beta1;
  real beta2;

  real<lower = 0> lambda1;
  real<lower = 0> lambda2;

  real psi1;
  real psi2;

}
transformed parameters {

  real mast_total_mean;
  real mast_year_sd;
  
  
  // Delta T matrices
  vector[Ncls] dt_else[Nyrs];
  
  // New matrix for masting
  vector<lower = 0.001> [Ncls] yy_else[Nyrs];

  vector[Nyrs] year_sum_x_fc_ef_else;
  vector[Nyrs] year_sum_x_fc_wf_else;
  vector[Nyrs] ERM;
  vector[Nyrs] WRM;

  vector[Nyrs] mast_yr_means;
  
  vector[Nsp] irr_e[Nyrs];
  vector[Nsp] irr_w[Nyrs];

  vector[Nyrs] disease_e;
  vector[Nyrs] disease_w;

  //Fill in dt from observed data
  //This is all known but could be adapted to incorporate more missing years
  // (e.g. hindcasting or further forwards forecasting)
  dt_else[dt_obs_ind,] = dt_obs[dt_obs_ind,];

  for (i in 1:Nyrs){

    //Fill in known masting data
    yy_else[i, ii_obs_mast[i, 1:N_obs_mast[i]]] = yy_obs_mast[i, ii_obs_mast[i, 1:N_obs_mast[i]]];
    yy_else[i, ii_miss_mast[i, 1:N_miss_mast[i]]] = yy_miss_else[i, ii_miss_mast[i, 1:N_miss_mast[i]]];

    // Below lines get sums for two regions - western forests and northern forests,
    // Standardize by max values
    // Subtract 0.001 to un-bias regional masting estimates
    // (small value was added to model local masting as gamma-distributed)
    year_sum_x_fc_wf_else[i] = sum((yy_else[i]-0.001) .*fc .* er_identity[1])/mean_wf_val;
    year_sum_x_fc_ef_else[i] = sum((yy_else[i]-0.001) .*fc .* er_identity[2])/mean_ef_val;

    mast_yr_means[i] = mean(yy_else[i]);
  }
  
  //Average local masting value across all years, sites
  mast_total_mean = mean(mast_yr_means);
  
  //Record variation in yearly masting averages 
  mast_year_sd = sd(mast_yr_means);

  //Center
  ERM = (year_sum_x_fc_ef_else - mean(year_sum_x_fc_ef_else));
  WRM = (year_sum_x_fc_wf_else - mean(year_sum_x_fc_wf_else));

  //Fill in known irruption data
  irr_w[irr_missing_yrs,] = irr_ind_miss_w[irr_missing_yrs,];
  irr_w[irr_obs_yrs,] = irr_ind_w[irr_obs_yrs,];
  
  irr_e[irr_missing_yrs,] = irr_ind_miss_e[irr_missing_yrs,];
  irr_e[irr_obs_yrs,] = irr_ind_e[irr_obs_yrs,];
  
  //Fill in known disease data
  disease_e[disease_ts_observed] = east_disease_ts[disease_ts_observed];
  disease_e[disease_ts_unobserved] = dis_miss_e[disease_ts_unobserved];
  
  disease_w[disease_ts_observed] = west_disease_ts[disease_ts_observed];
  disease_w[disease_ts_unobserved] = dis_miss_w[disease_ts_unobserved];
}
model {
  vector[Ncls] mu_mast[Nyrs]; // Linear predictor for masting
  
  alpha ~ normal(mu_alpha,sigma_alpha);
  mu_alpha ~ normal(0,4);
  sigma_alpha ~ gamma(2,2);
  
  // Soft constraint on mean masting across all sitess, years = for identifiability
  mast_total_mean ~ normal(1,0.025);
  
  nu ~ normal(0,1);
  theta ~ normal(0,1);
  
  upsilon ~ normal(1,1); 
  
  gamma1 ~ normal(0,1);
  mu_omega1 ~ normal(0,3);
  mu_eta1 ~ normal(0,3);

  gamma2 ~ normal(0,1);
  mu_omega2 ~ normal(0,3);
  mu_eta2 ~ normal(0,3);
  
  sigma_omega1 ~ gamma(3,3);
  sigma_eta1 ~ gamma(3,3);
  
  sigma_omega2 ~ gamma(3,3);
  sigma_eta2 ~ gamma(3,3);

  omega1 ~ normal(mu_omega1,sigma_omega1);
  omega2 ~ normal(mu_omega2,sigma_omega2);

  eta1 ~ normal(mu_eta1,sigma_eta1);
  eta2 ~ normal(mu_eta2,sigma_eta2);

  sigma_phi ~ normal(1,0.5); 
  zeta1 ~ normal(1,2);
  zeta2 ~ normal(1,2);

  psi1 ~ normal(0,5);
  psi2 ~ normal(0,5);

  kappa1 ~ normal(0,3);
  kappa2 ~ normal(0,3);

  beta1 ~ normal(0,3);
  beta2 ~ normal(0,3);

  lambda1 ~ exponential(0.05); 
  lambda2 ~ exponential(0.05);

  // Broad prior from observed data, represents unknown masting in year before data start
  pre_data_mast ~ gamma(upsilon,upsilon/mast_total_mean); 

  // Set a prior on the distribution of total masting data based on other years
  // This allows for realistic regional-level masting patterns
  mean(pre_data_mast) ~ normal(mast_total_mean,mast_year_sd);

  // Prior from observed data, represents unknown irr in year before data start
  pre_data_irr_e ~ normal(0,1); 
  pre_data_irr_w ~ normal(0,1);
  
  // For MVN irruptions, covariance matrices
  target += lkj_corr_lpdf(Rho1 | 2);
  target += lkj_corr_lpdf(Rho2 | 2);

  //Species-level irruption sigmas
  sigma1 ~ gamma(3,2); 
  sigma2 ~ gamma(3,2);
  
  for (i in 1:Nyrs){

    // To avoid unknown parameters taking on crazy (e.g. Inf, or 10^-100) values, 
    // fill in the missing matrix with ones where yy_obs will be used
    yy_miss_else[i, ii_obs_mast[i, 1:N_obs_mast[i]]] ~ normal(0,1); 

    // Prior for phis to total around 0, for each year
    target += -0.5 * dot_self(phi[i, node1] - phi[i, node2]);
    sum(phi[i]) ~ normal(0, .001 * Ncls);

    // To account for missing data in the year immeadediately preceeding the study period,
    // the first year of the data draws from dummy variables representing delta-t, masting, and irruptions
    if(i == 1){
      mu_mast[i] = exp(alpha[i] + nu * pre_data_mast + theta * dt_else[i]+phi[i]*sigma_phi);
      //Masting, in the first year
      yy_else[i] ~ gamma(upsilon,upsilon/mu_mast[i]);

      //Irruptions
      target += multi_normal_lpdf( irr_w[i] | (gamma1 + omega1 * WRM[i] + eta1 .* pre_data_irr_w), quad_form_diag(Rho1 , sigma1));
      target += multi_normal_lpdf( irr_e[i] | (gamma2 + omega2 * ERM[i] + eta2 .* pre_data_irr_e), quad_form_diag(Rho2 , sigma2));
    } else {

      //Masting
      mu_mast[i] = exp(alpha[i] + nu * yy_else[i-1] + theta * dt_else[i]+phi[i]*sigma_phi);
      yy_else[i] ~ gamma(upsilon,upsilon/mu_mast[i]);

      //Irruptions
      target += multi_normal_lpdf( irr_w[i] | (gamma1 + omega1 * WRM[i] + eta1 .* irr_w[i-1]), quad_form_diag(Rho1 , sigma1));
      target += multi_normal_lpdf( irr_e[i] | (gamma2 + omega2 * ERM[i] + eta2 .* irr_e[i-1]), quad_form_diag(Rho2 , sigma2));
    }

    //Disease
    if (dis_obs_vect[i] == 1){

       outbreak_yn_west[i] ~ bernoulli_logit(psi1 + beta1 * irr_w[i,1]);
       outbreak_yn_east[i] ~ bernoulli_logit(psi2 + beta2 * irr_e[i,1]);

       if (outbreak_yn_west[i] == 1){
         disease_w[i] ~ gamma(lambda1,(lambda1/exp(zeta1 + kappa1 * irr_w[i,1])));
       }

       if (outbreak_yn_east[i] == 1){
         //Model outbreak liklihood
         disease_e[i] ~ gamma(lambda2,(lambda2/exp(zeta2 + kappa2 * irr_e[i,1])));
       }
     }
  }
}

// Predict future disease liklihood and size. Irruption predictions are modeled above.
 generated quantities {
   real dis_w_pred[Nyrs];
   real dis_e_pred[Nyrs];

   int<lower=0,upper=1> outbreak_w_pred[Nyrs];
   int<lower=0,upper=1> outbreak_e_pred[Nyrs];

   for (l in 1:Nyrs){

     // Outbreak, y/n
     outbreak_w_pred[l] = bernoulli_logit_rng(psi1 + beta1 * irr_w[l,1]);
     outbreak_e_pred[l] = bernoulli_logit_rng(psi2 + beta2 * irr_e[l,1]);

     // Size of outbreak
     if (outbreak_w_pred[l] == 1){
       dis_w_pred[l] = gamma_rng(lambda1,(lambda1/exp(zeta1 + kappa1 * irr_w[l,1])));
     } else {
       dis_w_pred[l] = 0;
     }

     if (outbreak_e_pred[l] == 1){
      dis_e_pred[l] = gamma_rng(lambda2,(lambda2/exp(zeta2 + kappa2 * irr_e[l,1])));
     } else {
       dis_e_pred[l] = 0;
     }
   }
 }
