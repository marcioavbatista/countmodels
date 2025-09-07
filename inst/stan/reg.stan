
functions{
#include chunks/links.stan
#include chunks/mylib.stan
}

data{
  int<lower=1> n;
  int<lower=1> p;
  array[n] int<lower=0> y;
  matrix[n, p] X;
  int<lower = 1, upper = 3> link;
  int<lower = 1, upper = 2> dist;
  row_vector[p] x_mean;
  vector<lower=0>[p] x_sd;
  int<lower=0, upper=1> approach;
  vector[p] mu_beta;
  matrix[p,p] sigma_beta;
  real<lower = 0> a_theta;
  real<lower = 0> b_theta;
}

parameters{
  vector[p] beta_std;
  real <lower=0> theta;
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

    if(approach==1){
        beta_std ~ multi_normal(mu_beta, sigma_beta);
        if(dist == 2){
            theta ~ gamma(a_theta,b_theta);
        }
    }




    if(dist == 1){
        target += sum(loglik_poisson(y, X, beta_std, link));
        
    } else{
        target += sum(loglik_negbin(y, X, beta_std, theta, link));
    }

}

