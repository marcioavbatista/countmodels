
functions{
#include chunks/links.stan
#include chunks/mylib.stan
}

data {
  int<lower=0> n;
  int<lower=0> p;
  array[n] int<lower=0> y;
  matrix[n, p] X;
  int link;
  row_vector[p] x_mean;
  vector<lower=0>[p] x_sd;
  int<lower=0, upper=1> approach;
  real mu_beta;
  real<lower=0> sigma_beta;
  real<lower=0> alpha_psi;
  real<lower=0> beta_psi;
}

parameters {
  real<lower=0, upper=1> psi;
  vector[p] beta_std;
}

transformed parameters{
  vector[p] beta;

  if(p==1){
    beta[1] = beta_std[1]/x_sd[1];
  }else{
    beta[2:p] = beta_std[2:p] ./ x_sd[2:p];
    beta[1] = beta_std[1]/x_sd[1] - x_mean[2:p]*beta[2:p];
  }
}

model{
    // likelihood:
    vector[n] loglik = loglik_zapoisson(y, X, beta_std, psi, link);
    target += sum(loglik);
    if(approach==1){
      // prior distributions:
      beta_std ~ normal(mu_beta, sigma_beta);
      psi ~ beta(alpha_psi, beta_psi);
    }
}


generated quantities{
  vector[n] log_lik = loglik_zapoisson(y, X, beta_std, psi, link);
}

