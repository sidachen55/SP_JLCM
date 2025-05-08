// paquid data analysis
// G=4 clusters with non-centered parameterization 

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
  vector[4] theta3;
  vector[4] theta4;
  vector[4] beta1;
  vector[4] beta2;
  vector[4] beta3;
  vector[4] beta4;
  real alpha1;
  real alpha2;
  real alpha3;
  real alpha4;
  real<lower=0> lambda1;
  real<lower=0> lambda2;
  real<lower=0> lambda3;
  real<lower=0> lambda4;
  real<lower=0> var_e1;
  real<lower=0> var_e2;
  real<lower=0> var_e3;
  real<lower=0> var_e4;
  
  // Non-centered parameterization: standard normal params
  matrix[n,3] z1i;
  matrix[n,3] z2i;
  matrix[n,3] z3i;
  matrix[n,3] z4i;
  
  // Standard deviations instead of variances
  vector<lower=0>[3] var_b1;
  vector<lower=0>[3] var_b2;
  vector<lower=0>[3] var_b3;
  vector<lower=0>[3] var_b4;
  
  vector[3] psi1;
  vector[3] psi2;
  vector[3] psi3;
}


transformed parameters {
  // Transform z to b using the standard deviations
  matrix[n,3] b1i;
  matrix[n,3] b2i;
  matrix[n,3] b3i;
  matrix[n,3] b4i;
  
  // Non-centered parameterization transformation
  b1i[,1] = sqrt(var_b1[1]) * z1i[,1];
  b1i[,2] = sqrt(var_b1[2]) * z1i[,2];
  b1i[,3] = sqrt(var_b1[3]) * z1i[,3];
  b2i[,1] = sqrt(var_b2[1]) * z2i[,1];
  b2i[,2] = sqrt(var_b2[2]) * z2i[,2];
  b2i[,3] = sqrt(var_b2[3]) * z2i[,3];
  b3i[,1] = sqrt(var_b3[1]) * z3i[,1];
  b3i[,2] = sqrt(var_b3[2]) * z3i[,2];
  b3i[,3] = sqrt(var_b3[3]) * z3i[,3];
  b4i[,1] = sqrt(var_b4[1]) * z4i[,1];
  b4i[,2] = sqrt(var_b4[2]) * z4i[,2];
  b4i[,3] = sqrt(var_b4[3]) * z4i[,3];
}

model{
{
   vector[N] mu1;
   vector[N] mu2;
   vector[N] mu3;
   vector[N] mu4;
   vector[n] Haz1;
   vector[n] Haz2;
   vector[n] Haz3;
   vector[n] Haz4;
   matrix[n,K] cumHazK1;
   matrix[n,K] cumHazK2;
   matrix[n,K] cumHazK3;
   matrix[n,K] cumHazK4;
   vector[n] cumHaz1;
   vector[n] cumHaz2;
   vector[n] cumHaz3;
   vector[n] cumHaz4;
   vector[n] loglik_JM1;
   vector[n] loglik_JM2;
   vector[n] loglik_JM3;
   vector[n] loglik_JM4;
   
// ------------------------------------------------------
//        LOG-LIKELIHOOD FOR LONGITUDINAL SUBMODEL                
// ------------------------------------------------------
   
   // Linear predictors
   mu1 = linear_predictor(x1, times, ID, theta1, b1i);
   mu2 = linear_predictor(x1, times, ID, theta2, b2i);
   mu3 = linear_predictor(x1, times, ID, theta3, b3i);
   mu4 = linear_predictor(x1, times, ID, theta4, b4i);

// ------------------------------------------------------
//        LOG-LIKELIHOOD FOR SURVIVAL SUBMODEL                
// ------------------------------------------------------
  Haz1 = lambda1 * pow(Time, lambda1 - 1) .* 
                 exp(beta1[1] + beta1[2] * x1 + beta1[3] * x2 + beta1[4] * (x1 .* x2) + alpha1 * (theta1[2] + b1i[,2] + (2*(theta1[3] + b1i[,3])) .* Time));
  Haz2 = lambda2 * pow(Time, lambda2 - 1) .* 
                 exp(beta2[1] + beta2[2] * x1 + beta2[3] * x2 + beta2[4] * (x1 .* x2) + alpha2 * (theta2[2] + b2i[,2] + (2*(theta2[3] + b2i[,3])) .* Time));
  Haz3 = lambda3 * pow(Time, lambda3 - 1) .* 
                 exp(beta3[1] + beta3[2] * x1 + beta3[3] * x2 + beta3[4] * (x1 .* x2) + alpha3 * (theta3[2] + b3i[,2] + (2*(theta3[3] + b3i[,3])) .* Time));
  Haz4 = lambda4 * pow(Time, lambda4 - 1) .* 
                 exp(beta4[1] + beta4[2] * x1 + beta4[3] * x2 + beta4[4] * (x1 .* x2) + alpha4 * (theta4[2] + b4i[,2] + (2*(theta4[3] + b4i[,3])) .* Time));

   // Optimized Gauss-Legendre quadrature
   for (j in 1:K) {
      // Compute linear predictors for quadrature points
      vector[n] adjLP1 = theta1[2] + b1i[,2] + (2*(theta1[3] + b1i[,3])) .* adjTime[, j];
      vector[n] adjLP2 = theta2[2] + b2i[,2] + (2*(theta2[3] + b2i[,3])) .* adjTime[, j]; 
      vector[n] adjLP3 = theta3[2] + b3i[,2] + (2*(theta3[3] + b3i[,3])) .* adjTime[, j];
      vector[n] adjLP4 = theta4[2] + b4i[,2] + (2*(theta4[3] + b4i[,3])) .* adjTime[, j];
      // Compute hazards at quadrature points
      cumHazK1[, j] = lambda1 * pow(adjTime[, j], lambda1 - 1) .* 
                     exp(beta1[1] + beta1[2] * x1 + beta1[3] * x2 + beta1[4] * (x1 .* x2) + alpha1 * adjLP1);
      cumHazK2[, j] = lambda2 * pow(adjTime[, j], lambda2 - 1) .* 
                     exp(beta2[1] + beta2[2] * x1 + beta2[3] * x2 + beta2[4] * (x1 .* x2) + alpha2 * adjLP2);
      cumHazK3[, j] = lambda3 * pow(adjTime[, j], lambda3 - 1) .* 
                     exp(beta3[1] + beta3[2] * x1 + beta3[3] * x2 + beta3[4] * (x1 .* x2) + alpha3 * adjLP3);
      cumHazK4[, j] = lambda4 * pow(adjTime[, j], lambda4 - 1) .* 
                     exp(beta4[1] + beta4[2] * x1 + beta4[3] * x2 + beta4[4] * (x1 .* x2) + alpha4 * adjLP4);
   }

   // Compute cumulative hazards
   for (i in 1:n) {
      cumHaz1[i] = Time[i] / 2 * dot_product(wk, cumHazK1[i,]);
      cumHaz2[i] = Time[i] / 2 * dot_product(wk, cumHazK2[i,]);
      cumHaz3[i] = Time[i] / 2 * dot_product(wk, cumHazK3[i,]);
      cumHaz4[i] = Time[i] / 2 * dot_product(wk, cumHazK4[i,]);
   }
 
// ------------------------------------------------------
//                       log-mixture of JM                       
// ------------------------------------------------------

   // Vectorized log-likelihood calculation
   for (i in 1:n) {
      // Compute longitudinal log-likelihood for subject i directly
      real longit1_i = normal_lpdf(y[start[i]:stop[i]] | mu1[start[i]:stop[i]], sqrt(var_e1));
      real longit2_i = normal_lpdf(y[start[i]:stop[i]] | mu2[start[i]:stop[i]], sqrt(var_e2));
      real longit3_i = normal_lpdf(y[start[i]:stop[i]] | mu3[start[i]:stop[i]], sqrt(var_e3));
      real longit4_i = normal_lpdf(y[start[i]:stop[i]] | mu4[start[i]:stop[i]], sqrt(var_e4));
      
      loglik_JM1[i] = longit1_i + status[i]*log(Haz1[i]) - cumHaz1[i];
      loglik_JM2[i] = longit2_i + status[i]*log(Haz2[i]) - cumHaz2[i];
      loglik_JM3[i] = longit3_i + status[i]*log(Haz3[i]) - cumHaz3[i];
      loglik_JM4[i] = longit4_i + status[i]*log(Haz4[i]) - cumHaz4[i];
   }

   // Mixture model likelihood
   for (i in 1:n) {
      vector[4] log_terms;
      vector[4] logits;
      vector[4] pi;
      logits[1] = psi1[1] + x1[i] * psi1[2] + x2[i] * psi1[3];
      logits[2] = psi2[1] + x1[i] * psi2[2] + x2[i] * psi2[3];
      logits[3] = psi3[1] + x1[i] * psi3[2] + x2[i] * psi3[3];
      logits[4] = 0; 
      pi = softmax(logits);
      log_terms[1] = log(pi[1]) + loglik_JM1[i];
      log_terms[2] = log(pi[2]) + loglik_JM2[i];
      log_terms[3] = log(pi[3]) + loglik_JM3[i];
      log_terms[4] = log(pi[4]) + loglik_JM4[i];
      target += log_sum_exp(log_terms); 
   }
}

// ------------------------------------------------------
//                       LOG-PRIORS                       
// ------------------------------------------------------
  // Longitudinal fixed effects
  target += normal_lpdf(theta1 | 0, 2);
  target += normal_lpdf(theta2 | 0, 2);
  target += normal_lpdf(theta3 | 0, 2);
  target += normal_lpdf(theta4 | 0, 2);

  // Survival fixed effects
  target += normal_lpdf(beta1 | 0, 3); 
  target += normal_lpdf(beta2 | 0, 3);
  target += normal_lpdf(beta3 | 0, 3);
  target += normal_lpdf(beta4 | 0, 3);

  // Association parameters
  target += normal_lpdf(alpha1 | 0, 3);
  target += normal_lpdf(alpha2 | 0, 3);
  target += normal_lpdf(alpha3 | 0, 3);
  target += normal_lpdf(alpha4 | 0, 3);

  // Shape parameter (Weibull hazard)
  target += gamma_lpdf(lambda1 | 2, 0.5);
  target += gamma_lpdf(lambda2 | 2, 0.5);
  target += gamma_lpdf(lambda3 | 2, 0.5);
  target += gamma_lpdf(lambda4 | 2, 0.5);
  
  // Residual error variance
  target += normal_lpdf(var_e1 | 0, 0.2);
  target += normal_lpdf(var_e2 | 0, 0.2);
  target += normal_lpdf(var_e3 | 0, 0.2);
  target += normal_lpdf(var_e4 | 0, 0.2);
  
  // Random-effects variances
  target += gamma_lpdf(var_b1 | 1.5, 1.5);
  target += gamma_lpdf(var_b2 | 1.5, 1.5);
  target += gamma_lpdf(var_b3 | 1.5, 1.5);
  target += gamma_lpdf(var_b4 | 1.5, 1.5);
  target += -1.0 / pow(var_b1[1] + var_b1[2] + var_b1[3], 0.5);
  target += -1.0 / pow(var_b2[1] + var_b2[2] + var_b2[3], 0.5);
  target += -1.0 / pow(var_b3[1] + var_b3[2] + var_b3[3], 0.5);
  target += -1.0 / pow(var_b4[1] + var_b4[2] + var_b4[3], 0.5);
  
  // Standard normal priors for z (non-centered parameterization)
  target += std_normal_lpdf(to_vector(z1i));
  target += std_normal_lpdf(to_vector(z2i));
  target += std_normal_lpdf(to_vector(z3i));
  target += std_normal_lpdf(to_vector(z4i));
  
  // mixing proportion
  target += normal_lpdf(psi1 | 0, 1);
  target += normal_lpdf(psi2 | 0, 1);
  target += normal_lpdf(psi3 | 0, 1);
}

generated quantities {
   int<lower=1, upper=4> comp[n];
   vector[n] log_lik;
   vector[N] mu1;
   vector[N] mu2;
   vector[N] mu3;
   vector[N] mu4;
   vector[n] Haz1;
   vector[n] Haz2;
   vector[n] Haz3;
   vector[n] Haz4;
   matrix[n,K] cumHazK1;
   matrix[n,K] cumHazK2;
   matrix[n,K] cumHazK3;
   matrix[n,K] cumHazK4;
   vector[n] cumHaz1;
   vector[n] cumHaz2;
   vector[n] cumHaz3;
   vector[n] cumHaz4;
   vector[n] loglik_JM1;
   vector[n] loglik_JM2;
   vector[n] loglik_JM3;
   vector[n] loglik_JM4;
   vector[n] log_lik_cond;
   real log_posterior_original; // To store the CP log posterior
   
   mu1 = linear_predictor(x1, times, ID, theta1, b1i);
   mu2 = linear_predictor(x1, times, ID, theta2, b2i);
   mu3 = linear_predictor(x1, times, ID, theta3, b3i);
   mu4 = linear_predictor(x1, times, ID, theta4, b4i);

   Haz1 = lambda1 * pow(Time, lambda1 - 1) .* 
                 exp(beta1[1] + beta1[2] * x1 + beta1[3] * x2 + beta1[4] * (x1 .* x2) + alpha1 * (theta1[2] + b1i[,2] + (2*(theta1[3] + b1i[,3])) .* Time));
   Haz2 = lambda2 * pow(Time, lambda2 - 1) .* 
                 exp(beta2[1] + beta2[2] * x1 + beta2[3] * x2 + beta2[4] * (x1 .* x2) + alpha2 * (theta2[2] + b2i[,2] + (2*(theta2[3] + b2i[,3])) .* Time));
   Haz3 = lambda3 * pow(Time, lambda3 - 1) .* 
                 exp(beta3[1] + beta3[2] * x1 + beta3[3] * x2 + beta3[4] * (x1 .* x2) + alpha3 * (theta3[2] + b3i[,2] + (2*(theta3[3] + b3i[,3])) .* Time));
   Haz4 = lambda4 * pow(Time, lambda4 - 1) .* 
                 exp(beta4[1] + beta4[2] * x1 + beta4[3] * x2 + beta4[4] * (x1 .* x2) + alpha4 * (theta4[2] + b4i[,2] + (2*(theta4[3] + b4i[,3])) .* Time));

   // Optimized Gauss-Legendre quadrature (repeat from model block)
   for (j in 1:K) {
      // Compute linear predictors for quadrature points
      vector[n] adjLP1 = theta1[2] + b1i[,2] + (2*(theta1[3] + b1i[,3])) .* adjTime[, j];
      vector[n] adjLP2 = theta2[2] + b2i[,2] + (2*(theta2[3] + b2i[,3])) .* adjTime[, j]; 
      vector[n] adjLP3 = theta3[2] + b3i[,2] + (2*(theta3[3] + b3i[,3])) .* adjTime[, j];
      vector[n] adjLP4 = theta4[2] + b4i[,2] + (2*(theta4[3] + b4i[,3])) .* adjTime[, j];
      // Compute hazards at quadrature points
      cumHazK1[, j] = lambda1 * pow(adjTime[, j], lambda1 - 1) .* 
                     exp(beta1[1] + beta1[2] * x1 + beta1[3] * x2 + beta1[4] * (x1 .* x2) + alpha1 * adjLP1);
      cumHazK2[, j] = lambda2 * pow(adjTime[, j], lambda2 - 1) .* 
                     exp(beta2[1] + beta2[2] * x1 + beta2[3] * x2 + beta2[4] * (x1 .* x2) + alpha2 * adjLP2);
      cumHazK3[, j] = lambda3 * pow(adjTime[, j], lambda3 - 1) .* 
                     exp(beta3[1] + beta3[2] * x1 + beta3[3] * x2 + beta3[4] * (x1 .* x2) + alpha3 * adjLP3);
      cumHazK4[, j] = lambda4 * pow(adjTime[, j], lambda4 - 1) .* 
                     exp(beta4[1] + beta4[2] * x1 + beta4[3] * x2 + beta4[4] * (x1 .* x2) + alpha4 * adjLP4);
   }

   // Compute cumulative hazards
   for (i in 1:n) {
      cumHaz1[i] = Time[i] / 2 * dot_product(wk, cumHazK1[i,]);
      cumHaz2[i] = Time[i] / 2 * dot_product(wk, cumHazK2[i,]);
      cumHaz3[i] = Time[i] / 2 * dot_product(wk, cumHazK3[i,]);
      cumHaz4[i] = Time[i] / 2 * dot_product(wk, cumHazK4[i,]);
   }

   // Vectorized log-likelihood calculation (repeat for generated quantities)
   for (i in 1:n) {
      // Compute longitudinal log-likelihood for subject i directly
      real longit1_i = normal_lpdf(y[start[i]:stop[i]] | mu1[start[i]:stop[i]], sqrt(var_e1));
      real longit2_i = normal_lpdf(y[start[i]:stop[i]] | mu2[start[i]:stop[i]], sqrt(var_e2));
      real longit3_i = normal_lpdf(y[start[i]:stop[i]] | mu3[start[i]:stop[i]], sqrt(var_e3));
      real longit4_i = normal_lpdf(y[start[i]:stop[i]] | mu4[start[i]:stop[i]], sqrt(var_e4));
      
      loglik_JM1[i] = longit1_i + status[i]*log(Haz1[i]) - cumHaz1[i];
      loglik_JM2[i] = longit2_i + status[i]*log(Haz2[i]) - cumHaz2[i];
      loglik_JM3[i] = longit3_i + status[i]*log(Haz3[i]) - cumHaz3[i];
      loglik_JM4[i] = longit4_i + status[i]*log(Haz4[i]) - cumHaz4[i];
   }

   for (i in 1:n) {
      vector[4] log_terms;
      vector[4] logits;
      vector[4] pi;
      logits[1] = psi1[1] + x1[i] * psi1[2] + x2[i] * psi1[3];
      logits[2] = psi2[1] + x1[i] * psi2[2] + x2[i] * psi2[3];
      logits[3] = psi3[1] + x1[i] * psi3[2] + x2[i] * psi3[3];
      logits[4] = 0; 
      pi = softmax(logits);
      log_terms[1] = log(pi[1]) + loglik_JM1[i];
      log_terms[2] = log(pi[2]) + loglik_JM2[i];
      log_terms[3] = log(pi[3]) + loglik_JM3[i];
      log_terms[4] = log(pi[4]) + loglik_JM4[i];
      log_lik[i] = log_sum_exp(log_terms); 
      vector[4] safe_terms = log_terms;
      for (j in 1:4) if (is_inf(safe_terms[j])) safe_terms[j] = max(log_terms) - 1000;
      comp[i] = categorical_logit_rng(safe_terms);
      if (comp[i] == 1) {
          log_lik_cond[i] = loglik_JM1[i];
      } else if (comp[i] == 2) {
          log_lik_cond[i] = loglik_JM2[i];
      } else if (comp[i] == 3) {
          log_lik_cond[i] = loglik_JM3[i];
      } else if (comp[i] == 4) {
          log_lik_cond[i] = loglik_JM4[i];
      }
   }
   
   // Compute log posterior in original parameter space (centered parameterization)
   log_posterior_original = 0;
   
   // Add prior terms
   // Longitudinal fixed effects
   log_posterior_original += normal_lpdf(theta1 | 0, 2);
   log_posterior_original += normal_lpdf(theta2 | 0, 2);
   log_posterior_original += normal_lpdf(theta3 | 0, 2);
   log_posterior_original += normal_lpdf(theta4 | 0, 2);

   // Survival fixed effects
   log_posterior_original += normal_lpdf(beta1 | 0, 3); 
   log_posterior_original += normal_lpdf(beta2 | 0, 3);
   log_posterior_original += normal_lpdf(beta3 | 0, 3);
   log_posterior_original += normal_lpdf(beta4 | 0, 3);

   // Association parameters
   log_posterior_original += normal_lpdf(alpha1 | 0, 3);
   log_posterior_original += normal_lpdf(alpha2 | 0, 3);
   log_posterior_original += normal_lpdf(alpha3 | 0, 3);
   log_posterior_original += normal_lpdf(alpha4 | 0, 3);

   // Shape parameter (Weibull hazard)
   log_posterior_original += gamma_lpdf(lambda1 | 2, 0.5);
   log_posterior_original += gamma_lpdf(lambda2 | 2, 0.5);
   log_posterior_original += gamma_lpdf(lambda3 | 2, 0.5);
   log_posterior_original += gamma_lpdf(lambda4 | 2, 0.5);
  
   // Residual error variance
   log_posterior_original += normal_lpdf(var_e1 | 0, 0.2);
   log_posterior_original += normal_lpdf(var_e2 | 0, 0.2);
   log_posterior_original += normal_lpdf(var_e3 | 0, 0.2);
   log_posterior_original += normal_lpdf(var_e4 | 0, 0.2);
  
   // Random-effects variances
   log_posterior_original += gamma_lpdf(var_b1 | 1.5, 1.5);
   log_posterior_original += gamma_lpdf(var_b2 | 1.5, 1.5);
   log_posterior_original += gamma_lpdf(var_b3 | 1.5, 1.5);
   log_posterior_original += gamma_lpdf(var_b4 | 1.5, 1.5);
   log_posterior_original += -1.0 / pow(var_b1[1] + var_b1[2] + var_b1[3], 0.5);
   log_posterior_original += -1.0 / pow(var_b2[1] + var_b2[2] + var_b2[3], 0.5);
   log_posterior_original += -1.0 / pow(var_b3[1] + var_b3[2] + var_b3[3], 0.5);
   log_posterior_original += -1.0 / pow(var_b4[1] + var_b4[2] + var_b4[3], 0.5);
   
   // Mixing proportion
   log_posterior_original += normal_lpdf(psi1 | 0, 1);
   log_posterior_original += normal_lpdf(psi2 | 0, 1);
   log_posterior_original += normal_lpdf(psi3 | 0, 1);
   
   // Add original CP random effects priors (instead of NCP priors)
   for (i in 1:n) {
     log_posterior_original += normal_lpdf(b1i[i,1] | 0, sqrt(var_b1[1]));
     log_posterior_original += normal_lpdf(b1i[i,2] | 0, sqrt(var_b1[2]));
     log_posterior_original += normal_lpdf(b1i[i,3] | 0, sqrt(var_b1[3]));
     log_posterior_original += normal_lpdf(b2i[i,1] | 0, sqrt(var_b2[1]));
     log_posterior_original += normal_lpdf(b2i[i,2] | 0, sqrt(var_b2[2]));
     log_posterior_original += normal_lpdf(b2i[i,3] | 0, sqrt(var_b2[3]));
     log_posterior_original += normal_lpdf(b3i[i,1] | 0, sqrt(var_b3[1]));
     log_posterior_original += normal_lpdf(b3i[i,2] | 0, sqrt(var_b3[2]));
     log_posterior_original += normal_lpdf(b3i[i,3] | 0, sqrt(var_b3[3]));
     log_posterior_original += normal_lpdf(b4i[i,1] | 0, sqrt(var_b4[1]));
     log_posterior_original += normal_lpdf(b4i[i,2] | 0, sqrt(var_b4[2]));
     log_posterior_original += normal_lpdf(b4i[i,3] | 0, sqrt(var_b4[3]));
   }
   
   // Add likelihood (already calculated)
   log_posterior_original += sum(log_lik);
}
