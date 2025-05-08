################################################################################
# for sim setting II (G=2 clusters)
# cov dependent membership submodel (Model 2)
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
#######################################
# initialization for jags
# fit JLCM with G=2 with jointlcmm
#######################################
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

# Fit the joint latent class model with 3 latent classes
final_model <- gridsearch(rep = 30, maxiter = 15, minit = init_model,
                          Jointlcmm(
                            fixed = y ~ times + x1,
                            mixture = ~ times + x1,  # Class-specific effects
                            random = ~ times,
                            subject = "id",
                            survival = Surv(time, status) ~ mixture(x2),
                            hazard = "Weibull",
                            hazardtype = "Specific",  
                            ng = 2,
                            data = data_combined,
                            verbose = FALSE
                          ))


# extract results from the fitted model
coefficients <- coef(final_model)
fit_summary <- summary(final_model)
posterior_probs <- final_model$pprob[, c("probYT1", "probYT2")]

init_function <- function() {
  list(
    var_e1 = unname((coefficients["stderr"]^2) * runif(1,0.8,1.2)),
    var_e2 = unname((coefficients["stderr"]^2) * runif(1,0.8,1.2)),
    beta1 = pmax(unname(c((coefficients["event1 +/-sqrt(Weibull2) class 1"]^2)*log(coefficients["event1 +/-sqrt(Weibull1) class 1"]^2) * runif(1,0.8,1.2), coefficients["x2 class1"] * runif(1,0.8,1.2))),-5),
    beta2 = pmax(unname(c((coefficients["event1 +/-sqrt(Weibull2) class 2"]^2)*log(coefficients["event1 +/-sqrt(Weibull1) class 2"]^2) * runif(1,0.8,1.2), coefficients["x2 class2"] * runif(1,0.8,1.2))),-5),
    alpha1 = rnorm(1,0,1), 
    alpha2 = rnorm(1,0,1), 
    theta1 = (matrix(fit_summary[,"coef"], 2, 3)[1,])* runif(1,0.8,1.4),
    theta2 = (matrix(fit_summary[,"coef"], 2, 3)[2,])* runif(1,0.8,1.4)
  )
}

init_list <- lapply(1:6, function(x) init_function())

##########################################
fit <- stan(file = "rstan_2g_mlr_new.stan", 
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


samples <- extract(fit, pars=c("theta1[1]","theta1[2]", "theta1[3]","theta2[1]","theta2[2]", "theta2[3]", "beta1[1]", "beta1[2]","beta2[1]", "beta2[2]", "alpha1", "alpha2", 
                               "lambda1", "lambda2",  "var_e1", "var_e2", "var_b1[1]", "var_b1[2]", "var_b2[1]", "var_b2[2]"), permuted=FALSE)
posterior_samples <- samples[,best_chain,]
calc_summary <- function(samples) {
  mean_val <- mean(samples)
  quantiles <- quantile(samples, probs=c(0.025, 0.975),names = F)
  return(c(mean = mean_val, `lower_2.5` = quantiles[1], `upper_97.5` = quantiles[2]))
}
#######################
#
######################

# Extract the `comp` parameter from the fit object
samples_comp <- extract(fit, pars = "comp", permuted = FALSE)
comp_extract <- unname(as.matrix(samples_comp[, best_chain, ]))


nsample <- 1000
n <- nrow(data$survival)
classprob<-matrix(NA,2,n)
for (i in 1:2) {
  for (j in 1:n){
    sp<-rep(NA,nsample)
    for (k in 1:(nsample)) {sp[k]<-1*(comp_extract[k,j]==i)}
    classprob[i,j]<-sum(sp)/(nsample)
  }
}
ld<-rep(NA,n)
for (i in 1:n){ld[i]<-which.max(classprob[,i])}
# s <- c(rep(1,n1),rep(2,n2)) # for S1 only

hd<-rep(NA,factorial(2))
for (j in 1:factorial(2)) {
  pm<-permn(c(1:2))[[j]]
  ld_pm<-rep(NA,n)
  for (k in 1:2) {ld_pm[which(ld==k)]<-pm[k]}
  hd[j]<-hamming.distance(ld_pm,s)
}
pm<-permn(c(1:2))[[which.min(hd)]]
spm<-order(pm)

post_class <- rep(NA,n)
for (k in 1:2) {post_class[which(ld==k)]<-pm[k]}

class_acc <- 1-hamming.distance(post_class,s)/n



fit_smy <- list( theta1_3 = calc_summary(posterior_samples[, paste0("theta", spm[1], "[3]")]),
                 theta2_3 = calc_summary(posterior_samples[, paste0("theta", spm[2], "[3]")]),
                 
                 theta1_1 = calc_summary(posterior_samples[, paste0("theta", spm[1], "[1]")]),
                 theta2_1 = calc_summary(posterior_samples[, paste0("theta", spm[2], "[1]")]),
                 
                 theta1_2 = calc_summary(posterior_samples[, paste0("theta", spm[1], "[2]")]),
                 theta2_2 = calc_summary(posterior_samples[, paste0("theta", spm[2], "[2]")]),
                 
                 beta1_2 = calc_summary(posterior_samples[, paste0("beta", spm[1], "[2]")]),
                 beta2_2 = calc_summary(posterior_samples[, paste0("beta", spm[2], "[2]")]),
                 
                 beta1_1 = calc_summary(posterior_samples[, paste0("beta", spm[1], "[1]")]),
                 beta2_1 = calc_summary(posterior_samples[, paste0("beta", spm[2], "[1]")]),
                 
                 alpha1 = calc_summary(posterior_samples[, paste0("alpha", spm[1], "")]),
                 alpha2 = calc_summary(posterior_samples[, paste0("alpha", spm[2], "")]),
                 
                 var_e1 = calc_summary(posterior_samples[, paste0("var_e", spm[1], "")]),
                 var_e2 = calc_summary(posterior_samples[, paste0("var_e", spm[2], "")]),
                 
                 lambda1 = calc_summary(posterior_samples[, paste0("lambda", spm[1], "")]),
                 lambda2 = calc_summary(posterior_samples[, paste0("lambda", spm[2], "")]),
                 
                 D1_11 = calc_summary(posterior_samples[, paste0("var_b", spm[1], "[1]")]),
                 D1_22 = calc_summary(posterior_samples[, paste0("var_b", spm[1], "[2]")]),
                 D2_11 = calc_summary(posterior_samples[, paste0("var_b", spm[2], "[1]")]),
                 D2_22 = calc_summary(posterior_samples[, paste0("var_b", spm[2], "[2]")]),
                 
                 time = total_sampling_time,
                 class_acc = class_acc
)


save(fit_smy, file=paste0("stan_2g_mlr_new_",task_id,".rdata"))


#############################################
# Create trace plots for specified parameters
#############################################
pdf(paste0("stan_2g_mlr_new_tr_",task_id,".pdf"))
traceplot(fit, pars = c("alpha1", "alpha2", "beta1", "beta2", "var_e1", "var_e2", "var_b1", "var_b2" ,"lp__"), inc_warmup = F) 
dev.off()