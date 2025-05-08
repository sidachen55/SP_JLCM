// paquid data analysis
// G=2 clusters with non-centered parameterization 

functions{
// ------------------------------------------------------
//     LINEAR PREDICTOR FOR THE LONGITUDINAL SUBMODEL                
// ------------------------------------------------------ 
    vector linear_predictor(vector x, vector times, int[] ID, vector theta, matrix bi){
         int N = num_elements(times);
         vector[N] out;
         out = theta[1] + bi[ID,1] + theta[2]*times + bi[ID,2].*times + theta[3]*square(times) + bi[ID,3].*square(times) + theta[4]*x[ID];
         return out;
    } 
// ------------------------------------------------------ 
}


data{
  int N; 
  int n;
  vector[N] y;
  vector[N] times;
  int<lower=1,upper=n> ID[N];
  vector[n] Time;
  vector[n] status;
  vector[n] x1;
  vector[n] x2;
  int K;
  vector[K] xk;
  vector[K] wk;
  int<lower=1,upper=N> start[n];
  int<lower=1,upper=N> stop[n];
}

transformed data {
  // Precompute adjusted times for quadrature
  matrix[n, K] adjTime;
  for (j in 1:K) {
    adjTime[, j] = Time / 2 * (xk[j] + 1);
  }
}

parameters{
  vector[4] theta1;
  vector[4] theta2;
  vector[4] beta1;
  vector[4] beta2;
  real alpha1;
  real alpha2;
  real<lower=0> lambda1;
  real<lower=0> lambda2;
  real<lower=0> var_e1;
  real<lower=0> var_e2;
  
  // Non-centered parameterization: standard normal params
  matrix[n,3] z1i;
  matrix[n,3] z2i;
  
  // Standard deviations instead of variances
  vector<lower=0>[3] var_b1;
  vector<lower=0>[3] var_b2;
  
  vector[3] psi1;
}


transformed parameters {
  // Transform z to b using the standard deviations
  matrix[n,3] b1i;
  matrix[n,3] b2i;
  
  // Non-centered parameterization transformation
  b1i[,1] = sqrt(var_b1[1]) * z1i[,1];
  b1i[,2] = sqrt(var_b1[2]) * z1i[,2];
  b1i[,3] = sqrt(var_b1[3]) * z1i[,3];
  b2i[,1] = sqrt(var_b2[1]) * z2i[,1];
  b2i[,2] = sqrt(var_b2[2]) * z2i[,2];
  b2i[,3] = sqrt(var_b2[3]) * z2i[,3];
}

model{
{
   vector[N] mu1;
   vector[N] mu2;
   vector[n] Haz1;
   vector[n] Haz2;
   matrix[n,K] cumHazK1;
   matrix[n,K] cumHazK2;
   vector[n] cumHaz1;
   vector[n] cumHaz2;
   vector[n] loglik_JM1;
   vector[n] loglik_JM2;
   
// ------------------------------------------------------
//        LOG-LIKELIHOOD FOR LONGITUDINAL SUBMODEL                
// ------------------------------------------------------
   
   // Linear predictors
   mu1 = linear_predictor(x1, times, ID, theta1, b1i);
   mu2 = linear_predictor(x1, times, ID, theta2, b2i);

// ------------------------------------------------------
//        LOG-LIKELIHOOD FOR SURVIVAL SUBMODEL                
// ------------------------------------------------------
  Haz1 = lambda1 * pow(Time, lambda1 - 1) .* 
                 exp(beta1[1] + beta1[2] * x1 + beta1[3] * x2 + beta1[4] * (x1 .* x2) + alpha1 * (theta1[1] + b1i[,1] + (theta1[2] + b1i[,2]) .* Time + (theta1[3] + b1i[,3]) .* square(Time) + theta1[4] * x1));
  Haz2 = lambda2 * pow(Time, lambda2 - 1) .* 
                 exp(beta2[1] + beta2[2] * x1 + beta2[3] * x2 + beta2[4] * (x1 .* x2) + alpha2 * (theta2[1] + b2i[,1] + (theta2[2] + b2i[,2]) .* Time + (theta2[3] + b2i[,3]) .* square(Time) + theta2[4] * x1));

   // Optimized Gauss-Legendre quadrature
   for (j in 1:K) {
      // Compute linear predictors for quadrature points
      vector[n] adjLP1 = theta1[1] + b1i[,1] + (theta1[2] + b1i[,2]) .* adjTime[, j] + (theta1[3] + b1i[,3]) .* square(adjTime[, j]) + theta1[4] * x1;
      vector[n] adjLP2 = theta2[1] + b2i[,1] + (theta2[2] + b2i[,2]) .* adjTime[, j] + (theta2[3] + b2i[,3]) .* square(adjTime[, j]) + theta2[4] * x1;
      
      // Compute hazards at quadrature points
      cumHazK1[, j] = lambda1 * pow(adjTime[, j], lambda1 - 1) .* 
                     exp(beta1[1] + beta1[2] * x1 + beta1[3] * x2 + beta1[4] * (x1 .* x2) + alpha1 * adjLP1);
      cumHazK2[, j] = lambda2 * pow(adjTime[, j], lambda2 - 1) .* 
                     exp(beta2[1] + beta2[2] * x1 + beta2[3] * x2 + beta2[4] * (x1 .* x2) + alpha2 * adjLP2);
   }

   // Compute cumulative hazards
   for (i in 1:n) {
      cumHaz1[i] = Time[i] / 2 * dot_product(wk, cumHazK1[i,]);
      cumHaz2[i] = Time[i] / 2 * dot_product(wk, cumHazK2[i,]);
   }
 
// ------------------------------------------------------
//                       log-mixture of JM                       
// ------------------------------------------------------

   // Vectorized log-likelihood calculation
   for (i in 1:n) {
      // Compute longitudinal log-likelihood for subject i directly
      real longit1_i = normal_lpdf(y[start[i]:stop[i]] | mu1[start[i]:stop[i]], sqrt(var_e1));
      real longit2_i = normal_lpdf(y[start[i]:stop[i]] | mu2[start[i]:stop[i]], sqrt(var_e2));
      
      loglik_JM1[i] = longit1_i + status[i]*log(Haz1[i]) - cumHaz1[i];
      loglik_JM2[i] = longit2_i + status[i]*log(Haz2[i]) - cumHaz2[i];
   }

   // Mixture model likelihood
   for (i in 1:n) {
      vector[2] log_terms;
      vector[2] logits;
      vector[2] pi;
      logits[1] = psi1[1] + x1[i] * psi1[2] + x2[i] * psi1[3];
      logits[2] = 0; 
      pi = softmax(logits);
      log_terms[1] = log(pi[1]) + loglik_JM1[i];
      log_terms[2] = log(pi[2]) + loglik_JM2[i];
      target += log_sum_exp(log_terms); 
   }
}

// ------------------------------------------------------
//                       LOG-PRIORS                       
// ------------------------------------------------------
  // Longitudinal fixed effects
  target += normal_lpdf(theta1 | 0, 2);
  target += normal_lpdf(theta2 | 0, 2);

  // Survival fixed effects
  target += normal_lpdf(beta1 | 0, 3); 
  target += normal_lpdf(beta2 | 0, 3);

  // Association parameters
  target += normal_lpdf(alpha1 | 0, 3);
  target += normal_lpdf(alpha2 | 0, 3);

  // Shape parameter (Weibull hazard)
  target += gamma_lpdf(lambda1 | 2, 0.5);
  target += gamma_lpdf(lambda2 | 2, 0.5);
  
  // Residual error variance
  target += normal_lpdf(var_e1 | 0, 0.2);
  target += normal_lpdf(var_e2 | 0, 0.2);
  
  // Random-effects variances
  target += gamma_lpdf(var_b1 | 1.5, 1.5);
  target += gamma_lpdf(var_b2 | 1.5, 1.5);
  target += -1.0 / pow(var_b1[1] + var_b1[2] + var_b1[3], 0.5);
  target += -1.0 / pow(var_b2[1] + var_b2[2] + var_b2[3], 0.5);
  
  // Standard normal priors for z (non-centered parameterization)
  target += std_normal_lpdf(to_vector(z1i));
  target += std_normal_lpdf(to_vector(z2i));
  
  // mixing proportion
  target += normal_lpdf(psi1 | 0, 1);
}


generated quantities {
   int<lower=1, upper=2> comp[n];
   vector[n] log_lik;
   vector[N] mu1;
   vector[N] mu2;
   vector[n] Haz1;
   vector[n] Haz2;
   matrix[n,K] cumHazK1;
   matrix[n,K] cumHazK2;
   vector[n] cumHaz1;
   vector[n] cumHaz2;
   vector[n] loglik_JM1;
   vector[n] loglik_JM2;
   vector[n] log_lik_cond;
   real log_posterior_original; // NEW: to store the CP log posterior
   
   mu1 = linear_predictor(x1, times, ID, theta1, b1i);
   mu2 = linear_predictor(x1, times, ID, theta2, b2i);

   Haz1 = lambda1 * pow(Time, lambda1 - 1) .* 
                 exp(beta1[1] + beta1[2] * x1 + beta1[3] * x2 + beta1[4] * (x1 .* x2) + alpha1 * (theta1[1] + b1i[,1] + (theta1[2] + b1i[,2]) .* Time + (theta1[3] + b1i[,3]) .* square(Time) + theta1[4] * x1));
   Haz2 = lambda2 * pow(Time, lambda2 - 1) .* 
                 exp(beta2[1] + beta2[2] * x1 + beta2[3] * x2 + beta2[4] * (x1 .* x2) + alpha2 * (theta2[1] + b2i[,1] + (theta2[2] + b2i[,2]) .* Time + (theta2[3] + b2i[,3]) .* square(Time) + theta2[4] * x1));

   // Optimized Gauss-Legendre quadrature (repeat from model block)
   for (j in 1:K) {
      // Compute linear predictors for quadrature points
      vector[n] adjLP1 = theta1[1] + b1i[,1] + (theta1[2] + b1i[,2]) .* adjTime[, j] + (theta1[3] + b1i[,3]) .* square(adjTime[, j]) + theta1[4] * x1;
      vector[n] adjLP2 = theta2[1] + b2i[,1] + (theta2[2] + b2i[,2]) .* adjTime[, j] + (theta2[3] + b2i[,3]) .* square(adjTime[, j]) + theta2[4] * x1;
      
      // Compute hazards at quadrature points
      cumHazK1[, j] = lambda1 * pow(adjTime[, j], lambda1 - 1) .* 
                     exp(beta1[1] + beta1[2] * x1 + beta1[3] * x2 + beta1[4] * (x1 .* x2) + alpha1 * adjLP1);
      cumHazK2[, j] = lambda2 * pow(adjTime[, j], lambda2 - 1) .* 
                     exp(beta2[1] + beta2[2] * x1 + beta2[3] * x2 + beta2[4] * (x1 .* x2) + alpha2 * adjLP2);
   }

   // Compute cumulative hazards
   for (i in 1:n) {
      cumHaz1[i] = Time[i] / 2 * dot_product(wk, cumHazK1[i,]);
      cumHaz2[i] = Time[i] / 2 * dot_product(wk, cumHazK2[i,]);
   }

   // Vectorized log-likelihood calculation (repeat for generated quantities)
   for (i in 1:n) {
      // Compute longitudinal log-likelihood for subject i directly
      real longit1_i = normal_lpdf(y[start[i]:stop[i]] | mu1[start[i]:stop[i]], sqrt(var_e1));
      real longit2_i = normal_lpdf(y[start[i]:stop[i]] | mu2[start[i]:stop[i]], sqrt(var_e2));
      
      loglik_JM1[i] = longit1_i + status[i]*log(Haz1[i]) - cumHaz1[i];
      loglik_JM2[i] = longit2_i + status[i]*log(Haz2[i]) - cumHaz2[i];
   }

   for (i in 1:n) {
      vector[2] log_terms;
      vector[2] logits;
      vector[2] pi;
      logits[1] = psi1[1] + x1[i] * psi1[2] + x2[i] * psi1[3];
      logits[2] = 0; 
      pi = softmax(logits);
      log_terms[1] = log(pi[1]) + loglik_JM1[i];
      log_terms[2] = log(pi[2]) + loglik_JM2[i];
      log_lik[i] = log_sum_exp(log_terms); 
      vector[2] safe_terms = log_terms;
      for (j in 1:2) if (is_inf(safe_terms[j])) safe_terms[j] = max(log_terms) - 1000;
      comp[i] = categorical_logit_rng(safe_terms);
      if (comp[i] == 1) {
          log_lik_cond[i] = loglik_JM1[i];
      } else if (comp[i] == 2) {
          log_lik_cond[i] = loglik_JM2[i];
      }
   }
   
   // Compute log posterior in original parameter space (centered parameterization)
   log_posterior_original = 0;
   
   // Add prior terms
   // Longitudinal fixed effects
   log_posterior_original += normal_lpdf(theta1 | 0, 2);
   log_posterior_original += normal_lpdf(theta2 | 0, 2);

   // Survival fixed effects
   log_posterior_original += normal_lpdf(beta1 | 0, 3); 
   log_posterior_original += normal_lpdf(beta2 | 0, 3);

   // Association parameters
   log_posterior_original += normal_lpdf(alpha1 | 0, 3);
   log_posterior_original += normal_lpdf(alpha2 | 0, 3);

   // Shape parameter (Weibull hazard)
   log_posterior_original += gamma_lpdf(lambda1 | 2, 0.5);
   log_posterior_original += gamma_lpdf(lambda2 | 2, 0.5);
  
   // Residual error variance
   log_posterior_original += normal_lpdf(var_e1 | 0, 0.2);
   log_posterior_original += normal_lpdf(var_e2 | 0, 0.2);
  
   // Random-effects variances
   log_posterior_original += gamma_lpdf(var_b1 | 1.5, 1.5);
   log_posterior_original += gamma_lpdf(var_b2 | 1.5, 1.5);
   log_posterior_original += -1.0 / pow(var_b1[1] + var_b1[2] + var_b1[3], 0.5);
   log_posterior_original += -1.0 / pow(var_b2[1] + var_b2[2] + var_b2[3], 0.5);
   
   // Mixing proportion
   log_posterior_original += normal_lpdf(psi1 | 0, 1);
   
   // Add original CP random effects priors (instead of NCP priors)
   for (i in 1:n) {
     log_posterior_original += normal_lpdf(b1i[i,1] | 0, sqrt(var_b1[1]));
     log_posterior_original += normal_lpdf(b1i[i,2] | 0, sqrt(var_b1[2]));
     log_posterior_original += normal_lpdf(b1i[i,3] | 0, sqrt(var_b1[3]));
     log_posterior_original += normal_lpdf(b2i[i,1] | 0, sqrt(var_b2[1]));
     log_posterior_original += normal_lpdf(b2i[i,2] | 0, sqrt(var_b2[2]));
     log_posterior_original += normal_lpdf(b2i[i,3] | 0, sqrt(var_b2[3]));
   }
   
   // Add likelihood (already calculated)
   log_posterior_original += sum(log_lik);
}
