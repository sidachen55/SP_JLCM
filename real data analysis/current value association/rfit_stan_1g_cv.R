###################################################
# application to the paquid data (stan version)
# G=1
###################################################
library(rstan)
library(loo)
library(statmod)
library(lcmm)
library(NormPsy)
library(bayesplot)
set.seed(2024)

######################################################
# pre-process
# time scale for both processes: (age-65)/10  
######################################################
paquid$normMMSE <- normMMSE(paquid$MMSE) # define normMMSE
paquid$age65 <- (paquid$age - 65) / 10 # define age65
paquidS <- paquid[paquid$agedem > paquid$age_init, ] # retain 499 subj
summary(as.numeric(table(paquidS$ID))) 

longdata <- data.frame(id=paquidS$ID, y=paquidS$normMMSE, times=paquidS$age65, x1=paquidS$CEP)
longdata <- na.omit(longdata) # NA in y is removed;  n=499 [no subj removed]
longdata$y <- as.vector(scale(longdata$y))

survdata <- data.frame(id=paquidS$ID, time=(paquidS$agedem-65)/10, status=paquidS$dem, x1=paquidS$CEP, x2=paquidS$male)
survdata <- survdata[!duplicated(survdata$id), ]
data <- list(longitudinal=longdata, survival=survdata)
summary(as.numeric(table(data$longitudinal$id)))


data$survival$start <- data$survival$stop <- NA
aux <- 0
for(i in 1:nrow(data$survival)){
  pos <- which(data$longitudinal$id == data$survival$id[i])
  data$survival$start[i] <- aux + 1
  data$survival$stop[i] <- aux + length(pos)
  aux <- data$survival$stop[i]
}

################################################################################
# Gauss-Legendre quadrature (15 points)
glq <- gauss.quad(15, kind = "legendre")
xk <- glq$nodes   # nodes
wk <- glq$weights # weights
K <- length(xk)   # K-points

################################################################################

################################################################################
# initialization for rstan
# fit JLCM with G=1 with jointlcmm
################################################################################
paquidS$agedem65 <- (paquidS$agedem-65)/10
paquidS$std_normMMSE <- as.vector(scale(paquidS$normMMSE))

final_model <- Jointlcmm(
  fixed =std_normMMSE ~ poly(age65, degree = 2, raw = TRUE) + CEP,
  random = ~ poly(age65, degree = 2, raw = TRUE),
  subject = "ID",
  survival = Surv(agedem65, dem) ~ CEP + male,
  hazard = "Weibull", 
  ng = 1,
  data = paquidS)


# extract results from the fitted JLCM
coefficients <- coef(final_model)
fit_summary <- summary(final_model)


init_function <- function() {
  list(
    var_e1 = unname((coefficients["stderr"]^2) * runif(1,0.8,1.2)),
    beta1 = pmax(unname(c((coefficients["event1 +/-sqrt(Weibull2)"]^2)*log(coefficients["event1 +/-sqrt(Weibull1)"]^2) * runif(1,0.8,1.2), coefficients["CEP"] * runif(1,0.8,1.2), coefficients["male"] * runif(1,0.8,1.2),rnorm(1,0,1))),-5),
    alpha1 = runif(1,-5,0), 
    theta1 = as.numeric(fit_summary[,"coef"] * runif(1,0.8,1.2))
  )
}


init_list <- lapply(1:6, function(x) init_function())



##########################################
fit <- stan(file = "stan_1g.stan", 
            data = list(N=nrow(data$longitudinal), n=nrow(data$survival), 
                        y=data$longitudinal$y, times=data$longitudinal$times, 
                        ID=rep(1:nrow(data$survival),table(data$longitudinal$id)),
                        Time=data$survival$time, 
                        status=data$survival$status, x1=data$survival$x1, x2=data$survival$x2,
                        K=K, xk=xk, wk=wk, 
                        start=data$survival$start, stop=data$survival$stop),
            init = init_list,
            warmup = 4000,                 
            iter = 12000,
            thin = 4,
            chains = 6,
            seed = 2024,
            cores = 6,
            save_warmup=F) 


################################################################################
# chain selection
################################################################################
lp_allchain <- extract(fit, pars = "lp__", permuted = FALSE)
log_weights <- numeric(6)

for (chain in 1:6) {
  lp_chain <- lp_allchain[, chain, ]
  HPD_ind <- sort.list(lp_chain,decreasing = TRUE)[1:(0.8*(length(lp_chain)))]
  lp_hpd_chain <- lp_chain[HPD_ind]
  lp_hpd_chain_min <- min(lp_hpd_chain)
  log_weights[chain] <- lp_hpd_chain_min - log(sum(exp(-lp_hpd_chain+lp_hpd_chain_min)))
}
best_chain <- which.max(log_weights)

################################################################################
# loo, waic, estimates for model comparison 
################################################################################

loglik_allchain <- extract(fit, pars = "log_lik", permuted = FALSE)
loglik_extract <- as.matrix(loglik_allchain[,best_chain,])

loo_result <- loo(loglik_extract,cores=6)
waic_result <- waic(loglik_extract,cores=6)

IC_smy <- list(
  waic_condllk=waic_result,
  loo_condllk=loo_result)


save(IC_smy, file=paste0("IC_smy_1g.rdata"))


