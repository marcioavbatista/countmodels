
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
  int<lower=1, upper = 2> dist;
  int<lower=1, upper = 2> case_id;
  row_vector[p] x_mean;
  vector<lower=0>[p] x_sd;
  row_vector[q] z_mean;
  vector<lower=0>[q] z_sd;
  int<lower=0, upper=1> approach;
  vector[p] mu_beta;
  matrix[p,p]sigma_beta;
  vector[q] mu_psi;
  matrix[q,q] sigma_psi;
  real<lower = 0> a_theta;
  real<lower = 0> b_theta;
}

parameters {
  vector[q] psi_std;
  vector[p] beta_std;
  real <lower=0> theta; 
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


    if(approach==1){
        beta_std ~ multi_normal(mu_beta, sigma_beta);
        psi_std ~ multi_normal(mu_psi, sigma_psi);
        if(dist == 2){
            theta ~ gamma(a_theta,b_theta);
        }
    }




    if(dist == 1){
        if(case_id == 1){
            target += sum(loglik_zipoisson(y, X, Z, beta_std, psi_std, link1, link2));
        }else{
            target += sum(loglik_zapoisson(y, X, Z, beta_std, psi_std, link1, link2));
        }
        
        
        
    } else{
        if(case_id == 1){
            target += sum(loglik_zinegbin(y, X, Z, beta_std, psi_std,theta, link1, link2));
        }else{
            target += sum(loglik_zanegbin(y, X, Z, beta_std, psi_std,theta, link1, link2));
        }
        
    }

}


generated quantities{

    vector[n] loglik;

    if(dist == 1){
        if(case_id == 1){
            loglik = loglik_zipoisson(y, X, Z, beta_std, psi_std, link1, link2);
        }else{
            loglik = loglik_zapoisson(y, X, Z, beta_std, psi_std, link1, link2);
        }
        
        
        
    } else{
        if(case_id == 1){
            loglik = loglik_zinegbin(y, X, Z, beta_std, psi_std,theta, link1, link2);
        }else{
            loglik = loglik_zanegbin(y, X, Z, beta_std, psi_std,theta, link1, link2);
        }
        
    }
}

