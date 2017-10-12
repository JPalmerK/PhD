// static occupancy model
data {
  int<lower=0> nsite;
  int<lower=0> ndata;
  int<lower=0,upper=1> Y[ndata];
  //real x[ndata];
  int<lower=0, upper=30> site[ndata];
  real<lower=0, upper=1> p;
  //real<lower=0, upper=15000> shape1[ndata];
  //real<lower=0, upper=850000> shape2[ndata];
}
parameters {
  //real lo_p;       // log-odds detectibility
  real lo_psi1;    // log-odds initial occupancy (Intercept)
  real alpha_occ[nsite]; // Random effect due to site
  //real bx;         // covariate influencing occupancy
  //real<lower=0,upper=1> p[ndata];
  
}

transformed parameters {

  real psi[ndata] ; 
  for(ii in 1:ndata){
    psi[ii] <- inv_logit(lo_psi1 + alpha_occ[site[ii]]);
  }
}



model {
  // local variables
  
  // priors
  lo_psi1 ~ normal(0,20);
  alpha_occ ~ normal(0,20);
  

  
  // likelihood
  for(ii in 1:ndata){
    
    if (Y[ii] > 0){
      // likelihood occupied and see Y detections
      increment_log_prob(log(psi[site[ii]]) + bernoulli_log(Y[ii], p));
    } else {
      // marginalize over likelihood:
        // (1) occupied and see zero detections
      // (2) unoccupied
      increment_log_prob(
        log_sum_exp(
          log(psi[site[ii]]) + bernoulli_log(Y[ii],p), //(1)
          log(1-psi[site[ii]])                           //(2)
        )
      );
    }
  }
  

}