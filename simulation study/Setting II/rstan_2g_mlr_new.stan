// stan code for Setting II (G=2)
// Model 2
// use NCP for random effects

functions{
// ------------------------------------------------------
//     LINEAR PREDICTOR FOR THE LONGITUDINAL SUBMODEL                
// ------------------------------------------------------ 
    vector linear_predictor(vector x, vector times, int[] ID, vector theta, matrix bi){
         int N = num_elements(times);
         vector[N] out;

         out = theta[1] + bi[ID,1] + theta[2]*times + bi[ID,2].*times + theta[3]*x[ID];

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


parameters{
  vector[3] theta1;
  vector[3] theta2;
  vector[2] beta1;
  vector[2] beta2;
  real alpha1;
  real alpha2;
  real<lower=0> lambda1;
  real<lower=0> lambda2;
  real<lower=0> var_e1;
  real<lower=0> var_e2;
  matrix[n,2] z1i;
  matrix[n,2] z2i;
  vector<lower=0>[2] var_b1;
  vector<lower=0>[2] var_b2;
  vector[3] psi1;
}


transformed parameters {
 // Transform z to b using the standard deviations
  matrix[n,2] b1i;
  matrix[n,2] b2i;
  
  // Non-centered parameterization transformation
  b1i[,1] = sqrt(var_b1[1]) * z1i[,1];
  b1i[,2] = sqrt(var_b1[2]) * z1i[,2];
  b2i[,1] = sqrt(var_b2[1]) * z2i[,1];
  b2i[,2] = sqrt(var_b2[2]) * z2i[,2];
}

model{
{
   vector[N] mu1;
   vector[N] mu2;
   vector[N] longit1;
   vector[N] longit2;
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

   // Longitudinal Normal log-likelihood
   
   for(i in 1:N){
      longit1[i] = normal_lpdf(y[i] | mu1[i], sqrt(var_e1));
      longit2[i] = normal_lpdf(y[i] | mu2[i], sqrt(var_e2));
   }
  
// ------------------------------------------------------
//        LOG-LIKELIHOOD FOR SURVIVAL SUBMODEL                
// ------------------------------------------------------
Haz1 = lambda1 * pow(Time, lambda1 - 1) .* 
                 exp(beta1[1] + beta1[2] * x2 + alpha1 * (theta1[1] + b1i[,1] + (theta1[2] + b1i[,2]) .* Time + theta1[3] * x1));
Haz2 = lambda2 * pow(Time, lambda2 - 1) .* 
                 exp(beta2[1] + beta2[2] * x2 + alpha2 * (theta2[1] + b2i[,1] + (theta2[2] + b2i[,2]) .* Time + theta2[3] * x1));


for (j in 1:K) {
    vector[n] adjTime = Time / 2 * (xk[j] + 1);
    cumHazK1[, j] = lambda1 * pow(adjTime, lambda1 - 1) .* 
                    exp(beta1[1] + beta1[2] * x2 + alpha1 * (theta1[1] + b1i[,1] + (theta1[2] + b1i[,2]) .* adjTime + theta1[3] * x1));
    cumHazK2[, j] = lambda2 * pow(adjTime, lambda2 - 1) .* 
                    exp(beta2[1] + beta2[2] * x2 + alpha2 * (theta2[1] + b2i[,1] + (theta2[2] + b2i[,2]) .* adjTime + theta2[3] * x1));
}

for (i in 1:n) {
    cumHaz1[i] = Time[i] / 2 * dot_product(wk, cumHazK1[i,]);
    cumHaz2[i] = Time[i] / 2 * dot_product(wk, cumHazK2[i,]);
}
 
// ------------------------------------------------------
//                       log-mixture of JM                       
// ------------------------------------------------------


 for (i in 1:n) {
  loglik_JM1[i] = sum(longit1[start[i]:stop[i]]) + status[i]*log(Haz1[i]) - cumHaz1[i];
  loglik_JM2[i] = sum(longit2[start[i]:stop[i]]) + status[i]*log(Haz2[i]) - cumHaz2[i];
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
    target += log_sum_exp(log_terms); 
}
}

// ------------------------------------------------------
//                       LOG-PRIORS                       
// ------------------------------------------------------
  // Longitudinal fixed effects
  target += normal_lpdf(theta1 | 0, 5);
  target += normal_lpdf(theta2 | 0, 5);

  // Survival fixed effects
  target += normal_lpdf(beta1 | 0, 5); 
  target += normal_lpdf(beta2 | 0, 5);

  // Association parameters
  target += normal_lpdf(alpha1 | 0, 5);
  target += normal_lpdf(alpha2 | 0, 5);

  // Shape parameter (Weibull hazard)
  target += gamma_lpdf(lambda1 | 2, 0.5);
  target += gamma_lpdf(lambda2 | 2, 0.5);
  
  // Residual error variance
  target += normal_lpdf(var_e1 | 0, 0.5);
  target += normal_lpdf(var_e2 | 0, 0.5);
  
  // Random-effects variances
  target += gamma_lpdf(var_b1 | 1.5, 1.5);
  target += gamma_lpdf(var_b2 | 1.5, 1.5);
  target += -1.0 / pow(var_b1[1] + var_b1[2], 0.5);
  target += -1.0 / pow(var_b2[1] + var_b2[2], 0.5);
  
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
   vector[N] longit1;
   vector[N] longit2;
   vector[n] Haz1;
   vector[n] Haz2;
   matrix[n,K] cumHazK1;
   matrix[n,K] cumHazK2;
   vector[n] cumHaz1;
   vector[n] cumHaz2;
   vector[n] loglik_JM1;
   vector[n] loglik_JM2;
   real log_posterior_original;
   
   mu1 = linear_predictor(x1, times, ID, theta1, b1i);
   mu2 = linear_predictor(x1, times, ID, theta2, b2i);

   for(i in 1:N){
      longit1[i] = normal_lpdf(y[i] | mu1[i], sqrt(var_e1));
      longit2[i] = normal_lpdf(y[i] | mu2[i], sqrt(var_e2));
   }
Haz1 = lambda1 * pow(Time, lambda1 - 1) .* 
                 exp(beta1[1] + beta1[2] * x2 + alpha1 * (theta1[1] + b1i[,1] + (theta1[2] + b1i[,2]) .* Time + theta1[3] * x1));
Haz2 = lambda2 * pow(Time, lambda2 - 1) .* 
                 exp(beta2[1] + beta2[2] * x2 + alpha2 * (theta2[1] + b2i[,1] + (theta2[2] + b2i[,2]) .* Time + theta2[3] * x1));


for (j in 1:K) {
    vector[n] adjTime = Time / 2 * (xk[j] + 1);

    cumHazK1[, j] = lambda1 * pow(adjTime, lambda1 - 1) .* 
                    exp(beta1[1] + beta1[2] * x2 + alpha1 * (theta1[1] + b1i[,1] + (theta1[2] + b1i[,2]) .* adjTime + theta1[3] * x1));
    cumHazK2[, j] = lambda2 * pow(adjTime, lambda2 - 1) .* 
                    exp(beta2[1] + beta2[2] * x2 + alpha2 * (theta2[1] + b2i[,1] + (theta2[2] + b2i[,2]) .* adjTime + theta2[3] * x1));
}

for (i in 1:n) {
    cumHaz1[i] = Time[i] / 2 * dot_product(wk, cumHazK1[i,]);
    cumHaz2[i] = Time[i] / 2 * dot_product(wk, cumHazK2[i,]);
}
for (i in 1:n) {
  loglik_JM1[i] = sum(longit1[start[i]:stop[i]]) + status[i]*log(Haz1[i]) - cumHaz1[i];
  loglik_JM2[i] = sum(longit2[start[i]:stop[i]]) + status[i]*log(Haz2[i]) - cumHaz2[i];
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
    comp[i] = categorical_logit_rng(log_terms);
}
// Compute log posterior in original parameter space (centered parameterization)
   log_posterior_original = 0;
   
   // Add prior terms
   // Longitudinal fixed effects
   log_posterior_original += normal_lpdf(theta1 | 0, 5);
   log_posterior_original += normal_lpdf(theta2 | 0, 5);

   // Survival fixed effects
   log_posterior_original += normal_lpdf(beta1 | 0, 5); 
   log_posterior_original += normal_lpdf(beta2 | 0, 5);

   // Association parameters
   log_posterior_original += normal_lpdf(alpha1 | 0, 5);
   log_posterior_original += normal_lpdf(alpha2 | 0, 5);

   // Shape parameter (Weibull hazard)
   log_posterior_original += gamma_lpdf(lambda1 | 2, 0.5);
   log_posterior_original += gamma_lpdf(lambda2 | 2, 0.5);
  
   // Residual error variance
   log_posterior_original += normal_lpdf(var_e1 | 0, 0.5);
   log_posterior_original += normal_lpdf(var_e2 | 0, 0.5);
  
   // Random-effects variances
   log_posterior_original += gamma_lpdf(var_b1 | 1.5, 1.5);
   log_posterior_original += gamma_lpdf(var_b2 | 1.5, 1.5);
   log_posterior_original += -1.0 / pow(var_b1[1] + var_b1[2], 0.5);
   log_posterior_original += -1.0 / pow(var_b2[1] + var_b2[2], 0.5);
   
   // Mixing proportion
   log_posterior_original += normal_lpdf(psi1 | 0, 1);
   
   // Add original CP random effects priors (instead of NCP priors)
   for (i in 1:n) {
     log_posterior_original += normal_lpdf(b1i[i,1] | 0, sqrt(var_b1[1]));
     log_posterior_original += normal_lpdf(b1i[i,2] | 0, sqrt(var_b1[2]));
     log_posterior_original += normal_lpdf(b2i[i,1] | 0, sqrt(var_b2[1]));
     log_posterior_original += normal_lpdf(b2i[i,2] | 0, sqrt(var_b2[2]));
   }
   
   // Add likelihood (already calculated)
   log_posterior_original += sum(log_lik);
}
