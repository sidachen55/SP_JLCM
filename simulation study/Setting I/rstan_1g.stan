// final rstan version (homogeneous mixture weights)
// G=1 clusters

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
  vector[2] beta1;
  real alpha1;
  real<lower=0> lambda1;
  real<lower=0> var_e1;
  matrix[n,2] b1i;
  vector<lower=0>[2] var_b1;
}


transformed parameters {
  matrix[2,2] Sigma1;               // Covariance matrix for component 1 random effect
  Sigma1 = diag_matrix(var_b1);
}

model{
{
   vector[N] mu1;
   vector[N] longit1;
   vector[n] Haz1;
   matrix[n,K] cumHazK1;
   vector[n] cumHaz1;
   
// ------------------------------------------------------
//        LOG-LIKELIHOOD FOR LONGITUDINAL SUBMODEL                
// ------------------------------------------------------
   
   // Linear predictors
   mu1 = linear_predictor(x1, times, ID, theta1, b1i);

   // Longitudinal Normal log-likelihood
   
   for(i in 1:N){
      longit1[i] = normal_lpdf(y[i] | mu1[i], sqrt(var_e1));
   }
  
// ------------------------------------------------------
//        LOG-LIKELIHOOD FOR SURVIVAL SUBMODEL                
// ------------------------------------------------------
Haz1 = lambda1 * pow(Time, lambda1 - 1) .* 
                 exp(beta1[1] + beta1[2] * x2 + alpha1 * (theta1[1] + b1i[,1] + (theta1[2] + b1i[,2]) .* Time + theta1[3] * x1));

for (j in 1:K) {
    vector[n] adjTime = Time / 2 * (xk[j] + 1);

    cumHazK1[, j] = lambda1 * pow(adjTime, lambda1 - 1) .* 
                    exp(beta1[1] + beta1[2] * x2 + alpha1 * (theta1[1] + b1i[,1] + (theta1[2] + b1i[,2]) .* adjTime + theta1[3] * x1));
}

for (i in 1:n) {
    cumHaz1[i] = Time[i] / 2 * dot_product(wk, cumHazK1[i,]);
}
 
// ------------------------------------------------------
//                       log-mixture of JM                       
// ------------------------------------------------------


 for (i in 1:n) {
  target += sum(longit1[start[i]:stop[i]]) + status[i]*log(Haz1[i]) - cumHaz1[i];
}
}

// ------------------------------------------------------
//                       LOG-PRIORS                       
// ------------------------------------------------------
  // Longitudinal fixed effects
  target += normal_lpdf(theta1 | 0, 5);

  // Survival fixed effects
  target += normal_lpdf(beta1 | 0, 5); 

  // Association parameters
  target += normal_lpdf(alpha1 | 0, 5);

  // Shape parameter (Weibull hazard)
  target += gamma_lpdf(lambda1 | 2, 0.5);
  
  // Residual error variance
  target += normal_lpdf(var_e1 | 0, 0.5);
  
  // Random-effects variances
  target += inv_gamma_lpdf(var_b1 | 0.01, 0.01); 
  
  // log-density for random effects
  for(i in 1:n){ 
    target += multi_normal_lpdf(b1i[i,1:2] | rep_vector(0,2), Sigma1); 
  }
}


generated quantities {
   vector[n] log_lik;
   vector[N] mu1;
   vector[N] longit1;
   vector[n] Haz1;
   matrix[n,K] cumHazK1;
   vector[n] cumHaz1;
   
   mu1 = linear_predictor(x1, times, ID, theta1, b1i);

   for(i in 1:N){
      longit1[i] = normal_lpdf(y[i] | mu1[i], sqrt(var_e1));
   }
Haz1 = lambda1 * pow(Time, lambda1 - 1) .* 
                 exp(beta1[1] + beta1[2] * x2 + alpha1 * (theta1[1] + b1i[,1] + (theta1[2] + b1i[,2]) .* Time + theta1[3] * x1));

for (j in 1:K) {
    vector[n] adjTime = Time / 2 * (xk[j] + 1);

    cumHazK1[, j] = lambda1 * pow(adjTime, lambda1 - 1) .* 
                    exp(beta1[1] + beta1[2] * x2 + alpha1 * (theta1[1] + b1i[,1] + (theta1[2] + b1i[,2]) .* adjTime + theta1[3] * x1));
}

for (i in 1:n) {
    cumHaz1[i] = Time[i] / 2 * dot_product(wk, cumHazK1[i,]);
}

for (i in 1:n) {
  log_lik[i] = sum(longit1[start[i]:stop[i]]) + status[i]*log(Haz1[i]) - cumHaz1[i];
}
}
