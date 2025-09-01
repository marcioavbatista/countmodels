vector loglik_poisson(array[] int y, matrix X, vector beta, int link1){
  int n = num_elements(y);
  vector[n] lprob;
  vector[n] lp1 = X*beta;
  vector[n] mu = linkinv_poisson(lp1, link1);
  for(i in 1:n){
    lprob[i] = poisson_lpmf(y[i] | mu[i]);
  }
  
  
return lprob;
}

vector loglik_zipoisson(array[] int y, matrix X, matrix Z, vector beta, vector psi, int link1,int link2){
  int n = num_elements(y);
  vector[n] lprob;
  vector[n] lp1 = X*beta;
  vector[n] lp2 = Z*psi;
  vector[n] mu = linkinv_poisson(lp1, link1);
  vector[n] omega = linkinv_bern(lp2, link2);

  for(i in 1:n){
    if(y[i] == 0){
      lprob[i] = log_sum_exp(bernoulli_lpmf(0|omega[i]) , (bernoulli_lpmf(1|omega[i])) + poisson_lpmf(y[i] | mu[i]));
    }else{
      lprob[i] = bernoulli_lpmf(1|omega[i]) + poisson_lpmf(y[i] | mu[i]);
    }
  }

  
  
return lprob;
}

vector loglik_zapoisson(array[] int y, matrix X, vector beta, real psi, int link){
  int n = num_elements(y);
  vector[n] lprob;
  vector[n] lp1 = X*beta;
  vector[n] mu = linkinv_poisson(lp1, link);

  for(i in 1:n){
    if(y[i] == 0){
      lprob[i] = bernoulli_lpmf(0|psi);
    }else{
      lprob[i] = bernoulli_lpmf(1|psi) + ( poisson_lpmf(y[i] | mu[i]) - log( 1 - exp(poisson_lpmf(0|mu[i])) ) );
    }
  }

  
  
return lprob;
}



