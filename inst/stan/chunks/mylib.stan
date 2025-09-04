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

vector loglik_negbin(array[] int y, matrix X, vector beta,real theta, int link1){
  int n = num_elements(y);
  vector[n] lprob;
  vector[n] lp1 = X*beta;
  vector[n] mu = linkinv_poisson(lp1, link1);
  vector[n] alpha = mu * theta;

  for(i in 1:n){
    lprob[i] = neg_binomial_lpmf(y[i] | alpha[i], theta);
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
      lprob[i] = log_sum_exp(log(omega[i]) , (bernoulli_lpmf(1|omega[i])) + poisson_lpmf(y[i] | mu[i]));
    }else{
      lprob[i] = log1m(omega[i]) + poisson_lpmf(y[i] | mu[i]);
    }
  }

  
  
return lprob;
}

vector loglik_zapoisson(array[] int y, matrix X, matrix Z, vector beta, vector psi, int link1,int link2){
  int n = num_elements(y);
  vector[n] lprob;
  vector[n] lp1 = X*beta;
  vector[n] lp2 = Z*psi;
  vector[n] mu = linkinv_poisson(lp1, link1);
  vector[n] omega = linkinv_bern(lp2, link2);

  for(i in 1:n){
    if(y[i] == 0){
      lprob[i] = log(omega[i]);
    }else{
      lprob[i] = log1m(omega[i]) + poisson_lpmf(y[i] | mu[i]) - log1m_exp(-omega[i]);
    }
  }

  
  
return lprob;
}


vector loglik_zinegbin(array[] int y, matrix X, matrix Z, vector beta, vector psi, real theta ,int link1,int link2){
  int n = num_elements(y); 
  vector[n] lprob;
  vector[n] lp1 = X*beta;
  vector[n] lp2 = Z*psi;
  vector[n] mu = linkinv_poisson(lp1, link1);
  vector[n] omega = linkinv_bern(lp2, link2);
  vector[n] alpha = mu * theta;

  for(i in 1:n){
    if(y[i] == 0){
      lprob[i] = log_sum_exp(log(omega[i]) , (bernoulli_lpmf(1|omega[i])) + neg_binomial_lpmf(y[i]| alpha[i], theta));
    }else{
      lprob[i] = log1m(omega[i]) + neg_binomial_lpmf(y[i] | alpha[i], theta);
    }
  }

  
  
return lprob;
}

vector loglik_zanegbin(array[] int y, matrix X, matrix Z, vector beta, vector psi, real theta , int link1,int link2){
  int n = num_elements(y); 
  vector[n] lprob;
  vector[n] lp1 = X*beta;
  vector[n] lp2 = Z*psi;
  vector[n] mu = linkinv_poisson(lp1, link1);
  vector[n] omega = linkinv_bern(lp2, link2);
  vector[n] alpha = mu * theta;

  for(i in 1:n){
    if(y[i] == 0){
      lprob[i] = log(omega[i]);
    }else{
      lprob[i] = log1m(omega[i]) + neg_binomial_lpmf(y[i] | alpha[i], theta) - log1m_exp(-omega[i]);
    }
  }

  
  
return lprob;
}



