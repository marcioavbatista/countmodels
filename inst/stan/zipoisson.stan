
functions{
#include chunks/links.stan
#include chunks/mylib.stan
}

data {
  int<lower=0> n;
  int<lower=0> p;
  int<lower=0> q;
  array[n] int<lower=0> y;
  matrix[n, p] X;
  matrix[n, q] Z;
  int link1;
  int link2;
  row_vector[p] x_mean;
  vector<lower=0>[p] x_sd;
  row_vector[q] z_mean;
  vector<lower=0>[q] z_sd;
  int<lower=0, upper=1> approach;
  real mu_beta;
  real<lower=0> sigma_beta;
  real mu_psi;
  real<lower=0> sigma_psi;
}

parameters {
  vector[q] psi_std;
  vector[p] beta_std;
}

transformed parameters{
  vector[q] psi;
  vector[p] beta;

  if(p==1){
    beta[1] = beta_std[1]/x_sd[1];
  }else{
    beta[2:p] = beta_std[2:p] ./ x_sd[2:p];
    beta[1] = beta_std[1]/x_sd[1] - x_mean[2:p]*beta[2:p];
  }

  if(q==1){
    psi[1] = psi_std[1]/z_sd[1];
  }else{
    psi[2:q] = psi_std[2:q] ./ z_sd[2:q];
    psi[1] = psi_std[1]/z_sd[1] - z_mean[2:q]*psi[2:q];
  }
}

model{
    // likelihood:
    vector[n] loglik = loglik_zipoisson(y, X, Z, beta_std, psi_std, link1, link2);
    target += sum(loglik);
    if(approach==1){
      // prior distributions:
      beta_std ~ normal(mu_beta, sigma_beta);
      psi_std ~ normal(mu_psi, sigma_psi);
    }
}


generated quantities{
  vector[n] log_lik = loglik_zipoisson(y, X, Z, beta, psi, link1, link2);
}

