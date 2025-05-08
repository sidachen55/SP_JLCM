################################################################################
# for sim model1
# G=1
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
# initialization 
# fit JLCM with G=1 with jointlcmm
################################################################################
# data preparation
data_combined <- merge(data$longitudinal, data$survival, by = "id")
data_combined$x1 <- data_combined$x1.x
data_combined <- data_combined[, !(names(data_combined) %in% c("x1.x", "x1.y"))]

# Fit the initial model with one latent class for initialization
final_model <- Jointlcmm(
  fixed = y ~ times + x1,
  random = ~ times,
  subject = "id",
  survival = Surv(time, status) ~ x2,
  hazard = "Weibull",
  ng = 1,
  data = data_combined
)


# extract results from the fitted model
coefficients <- coef(final_model)
fit_summary <- summary(final_model)


init_function <- function() {
  list(
    var_e1 = unname((coefficients["stderr"]^2) * runif(1,0.8,1.2)),
    beta1 = pmax(unname(c((coefficients["event1 +/-sqrt(Weibull2)"]^2)*log(coefficients["event1 +/-sqrt(Weibull1)"]^2) * runif(1,0.8,1.2), coefficients["x2"] * runif(1,0.8,1.2))),-5),
    alpha1 = rnorm(1,0,1), 
    theta1 = as.numeric(fit_summary[,"coef"] * runif(1,0.8,1.2)),
    var_b1 = c(runif(1, 0.1, 5), runif(1, 0.1, 2))
  )
}


init_list <- lapply(1:6, function(x) init_function())



################################################################################
fit <- stan(file = "rstan_1g.stan", 
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


sampling_time <- get_elapsed_time(fit)
chain_total_times <- rowSums(sampling_time)
total_sampling_time <- max(chain_total_times)/3600

################################################################################
# chain selection
################################################################################
lp_allchain <- extract(fit, pars = "lp__", permuted = FALSE)
log_weights <- numeric(6)

for (chain in 1:6) {
  lp_chain <- lp_allchain[, chain, ]
  HPD_ind <- sort.list(lp_chain,decreasing = TRUE)[1:(0.8*(1000))]
  lp_hpd_chain <- lp_chain[HPD_ind]
  lp_hpd_chain_min <- min(lp_hpd_chain)
  log_weights[chain] <- lp_hpd_chain_min - log(sum(exp(-lp_hpd_chain+lp_hpd_chain_min)))
}
best_chain <- which.max(log_weights)

################################################################################
# loo, waic, and margllk estimates for model comparison 
################################################################################

loglik_allchain <- extract(fit, pars = "log_lik", permuted = FALSE)
loglik_extract <- as.matrix(loglik_allchain[,best_chain,])

loo_result <- loo(loglik_extract,cores=6)
waic_result <- waic(loglik_extract,cores=6)



fit_smy <- list( waic_condllk=waic_result,
                 loo_condllk=loo_result,
                 time = total_sampling_time
)


save(fit_smy, file=paste0("stan_ms_1g_",task_id,".rdata"))

#############################################
# Create trace plots for specified parameters
#############################################
pdf(paste0("stan_ms_1g_tr_",task_id,".pdf"))
traceplot(fit, pars = c("var_e1", "theta1", "var_b1", "alpha1", "lp__"), inc_warmup = F) 
dev.off()