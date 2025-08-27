functions{

vector inv_logit2(vector lp){
  int n = num_elements(lp);
  vector[n] x;
  for(i in 1:n){
     x[i] = 1/(1+exp(-lp[i]));
  }
  return(x);
}

/**
 * Code extracted/adapted from R package rstanarm
 *
 * @param lp Linear predictor vector
 * @param link An integer indicating the link function
 * @return A vector, i.e. inverse-link(lp)
 */
vector linkinv_poisson(vector lp, int link) {
  if (link == 1)      return exp(lp);     // log
  else if (link == 2) return lp;          // identity
  else if (link == 3) return(square(lp)); // sqrt
  else reject("Invalid link");
  return lp; // never reached
}


/**
 * * Code extracted/adapted from R package rstanarm
 *
 * @param lp Linear predictor vector
 * @param link An integer indicating the link function
 * @return A vector, i.e. inverse-link(lp)
 */
vector linkinv_bern(vector lp, int link) {
  //if (link == 1)      return(inv_logit(lp)); // logit
  if (link == 1)      return(inv_logit2(lp)); // logit
  else if (link == 2) return(Phi(lp)); // probit
  else if (link == 3) return(inv_cloglog(lp)); // cloglog
  else if (link == 4) return(atan(lp) / pi() + 0.5);  // Cauchy
  else reject("Invalid link");
  return lp; // never reached
}

}
