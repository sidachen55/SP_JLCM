################################################################################
# for sim setting I, S4 (G=5 clusters)
################################################################################
library(rstan)
library(loo)
library(statmod)
library(lcmm)
library(e1071)
library(combinat)

task_id_string <- Sys.getenv("SLURM_ARRAY_TASK_ID")
task_id <- as.numeric(task_id_string)
set.seed(task_id)
input_file <- paste0("data_", task_id, ".RData")
load(input_file)

################################################################################
# Gauss-Legendre quadrature (15 points)
glq <- gauss.quad(15, kind = "legendre")
xk <- glq$nodes   # nodes
wk <- glq$weights # weights
K <- length(xk)   # K-points

################################################################################

################################################################################
# initialization for jags
# fit JLCM with G=5 with jointlcmm
################################################################################
# data preparation
data_combined <- merge(data$longitudinal, data$survival, by = "id")
data_combined$x1 <- data_combined$x1.x
data_combined <- data_combined[, !(names(data_combined) %in% c("x1.x", "x1.y"))]

# Fit the initial model with one latent class for initialization
init_model <- Jointlcmm(
  fixed = y ~ times + x1,
  random = ~ times,
  subject = "id",
  survival = Surv(time, status) ~ x2,
  hazard = "Weibull",
  ng = 1,
  data = data_combined
)

# Fit the joint latent class model with 5 latent classes
final_model <- gridsearch(rep = 30, maxiter = 15, minit = init_model,
                          Jointlcmm(
                            fixed = y ~ times + x1,
                            mixture = ~ times + x1,  # Class-specific effects
                            random = ~ times,
                            subject = "id",
                            survival = Surv(time, status) ~ mixture(x2),
                            hazard = "Weibull",
                            hazardtype = "Specific",  
                            ng = 5,
                            data = data_combined,
                            verbose = FALSE
                          ))


# extract results from the fitted model
coefficients <- coef(final_model)
fit_summary <- summary(final_model)


if(is.null(fit_summary[,"coef"])) {
  init_function <- function() {
    list(
      alpha1 = rnorm(1,0,1), 
      alpha2 = rnorm(1,0,1), 
      alpha3 = rnorm(1,0,1),
      alpha4 = rnorm(1,0,1),
      alpha5 = rnorm(1,0,1),
      psi1 = rnorm(1,0,1),
      psi2 = rnorm(1,0,1),
      psi3 = rnorm(1,0,1),
      psi4 = rnorm(1,0,1),
      var_b1 = c(runif(1, 0.1, 5), runif(1, 0.1, 2)),
      var_b2 = c(runif(1, 0.1, 5), runif(1, 0.1, 2)),
      var_b3 = c(runif(1, 0.1, 5), runif(1, 0.1, 2)),
      var_b4 = c(runif(1, 0.1, 5), runif(1, 0.1, 2)),
      var_b5 = c(runif(1, 0.1, 5), runif(1, 0.1, 2))
    )
  }
} else {
init_function <- function() {
  list(
    var_e1 = unname((coefficients["stderr"]^2) * runif(1,0.8,1.2)),
    var_e2 = unname((coefficients["stderr"]^2) * runif(1,0.8,1.2)),
    var_e3 = unname((coefficients["stderr"]^2) * runif(1,0.8,1.2)),
    var_e4 = unname((coefficients["stderr"]^2) * runif(1,0.8,1.2)),
    var_e5 = unname((coefficients["stderr"]^2) * runif(1,0.8,1.2)),
    beta1 = pmin(pmax(unname(c((coefficients["event1 +/-sqrt(Weibull2) class 1"]^2)*log(coefficients["event1 +/-sqrt(Weibull1) class 1"]^2) * runif(1,0.8,1.2), coefficients["x2 class1"] * runif(1,0.8,1.2))),-5), 5),
    beta2 = pmin(pmax(unname(c((coefficients["event1 +/-sqrt(Weibull2) class 2"]^2)*log(coefficients["event1 +/-sqrt(Weibull1) class 2"]^2) * runif(1,0.8,1.2), coefficients["x2 class2"] * runif(1,0.8,1.2))),-5), 5),
    beta3 = pmin(pmax(unname(c((coefficients["event1 +/-sqrt(Weibull2) class 3"]^2)*log(coefficients["event1 +/-sqrt(Weibull1) class 3"]^2) * runif(1,0.8,1.2), coefficients["x2 class3"] * runif(1,0.8,1.2))),-5), 5),
    beta4 = pmin(pmax(unname(c((coefficients["event1 +/-sqrt(Weibull2) class 4"]^2)*log(coefficients["event1 +/-sqrt(Weibull1) class 4"]^2) * runif(1,0.8,1.2), coefficients["x2 class4"] * runif(1,0.8,1.2))),-5), 5),
    beta5 = pmin(pmax(unname(c((coefficients["event1 +/-sqrt(Weibull2) class 5"]^2)*log(coefficients["event1 +/-sqrt(Weibull1) class 5"]^2) * runif(1,0.8,1.2), coefficients["x2 class5"] * runif(1,0.8,1.2))),-5), 5),
    alpha1 = rnorm(1,0,1), 
    alpha2 = rnorm(1,0,1), 
    alpha3 = rnorm(1,0,1),
    alpha4 = rnorm(1,0,1),
    alpha5 = rnorm(1,0,1),
    theta1 = (matrix(fit_summary[,"coef"], 5, 3)[1,])* runif(1,0.7,1.4),
    theta2 = (matrix(fit_summary[,"coef"], 5, 3)[2,])* runif(1,0.7,1.4),
    theta3 = (matrix(fit_summary[,"coef"], 5, 3)[3,])* runif(1,0.7,1.4),
    theta4 = (matrix(fit_summary[,"coef"], 5, 3)[4,])* runif(1,0.7,1.4),
    theta5 = (matrix(fit_summary[,"coef"], 5, 3)[5,])* runif(1,0.7,1.4),
    psi1 = rnorm(1,0,1),
    psi2 = rnorm(1,0,1),
    psi3 = rnorm(1,0,1),
    psi4 = rnorm(1,0,1),
    var_b1 = c(runif(1, 0.1, 5), runif(1, 0.1, 2)),
    var_b2 = c(runif(1, 0.1, 5), runif(1, 0.1, 2)),
    var_b3 = c(runif(1, 0.1, 5), runif(1, 0.1, 2)),
    var_b4 = c(runif(1, 0.1, 5), runif(1, 0.1, 2)),
    var_b5 = c(runif(1, 0.1, 5), runif(1, 0.1, 2))
  )
}
}


init_list <- lapply(1:6, function(x) init_function())

##########################################
fit <- stan(file = "rstan_5g_new.stan", 
            data = list(N=nrow(data$longitudinal), n=nrow(data$survival), 
                        y=data$longitudinal$y, times=data$longitudinal$times, 
                        ID=rep(1:nrow(data$survival),table(data$longitudinal$id)),
                        Time=data$survival$time, 
                        status=data$survival$status, x1=data$survival$x1, x2=data$survival$x2,
                        K=K, xk=xk, wk=wk, 
                        start=data$survival$start, stop=data$survival$stop),
            init = init_list,
            warmup = 3000,                 
            iter = 6000,
            thin = 3,
            chains = 6,
            seed = 2024,
            cores = 6,
            save_warmup=F) 
##########################################

sampling_time <- get_elapsed_time(fit)
chain_total_times <- rowSums(sampling_time)
total_sampling_time <- max(chain_total_times)/3600

########################
# chain selection
########################
lp_allchain <- extract(fit, pars = "log_posterior_original", permuted = FALSE)
log_weights <- numeric(6)

for (chain in 1:6) {
  lp_chain <- lp_allchain[, chain, ]
  HPD_ind <- sort.list(lp_chain,decreasing = TRUE)[1:(0.8*(1000))]
  lp_hpd_chain <- lp_chain[HPD_ind]
  lp_hpd_chain_min <- min(lp_hpd_chain)
  log_weights[chain] <- lp_hpd_chain_min - log(sum(exp(-lp_hpd_chain+lp_hpd_chain_min)))
}
best_chain <- which.max(log_weights)

################################################################
# Extract the `comp` parameter from the fit object
samples_comp <- extract(fit, pars = "comp", permuted = FALSE)
comp_extract <- unname(as.matrix(samples_comp[, best_chain, ]))
nsample <- 1000
n <- nrow(data$survival)
classprob<-matrix(NA,5,n)
for (i in 1:5) {
  for (j in 1:n){
    sp<-rep(NA,nsample)
    for (k in 1:(nsample)) {sp[k]<-1*(comp_extract[k,j]==i)}
    classprob[i,j]<-sum(sp)/(nsample)
  }
}
ld<-rep(NA,n)
for (i in 1:n){ld[i]<-which.max(classprob[,i])}
################################################################################
# loo, waic, and margllk estimates for model comparison 
################################################################################

loglik_allchain <- extract(fit, pars = "log_lik_cond", permuted = FALSE)
loglik_extract <- as.matrix(loglik_allchain[,best_chain,])

loo_result <- loo(loglik_extract,cores=6)
waic_result <- waic(loglik_extract,cores=6)



fit_smy <- list( waic_condllk=waic_result,
                 loo_condllk=loo_result,
                 time = total_sampling_time,
                 decoding=table(ld)
)


save(fit_smy, file=paste0("stan_ms_5g_",task_id,".rdata"))

