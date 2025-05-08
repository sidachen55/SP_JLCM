###################################################
# application to the paquid data (stan version)
# G=4
# without delay entry, use age65 as time scale
###################################################
library(rstan)
library(loo)
library(statmod)
library(lcmm)
library(NormPsy)
library(ggplot2)
library(survival)
library(survminer)
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
# initialization for jags
# fit JLCM with G=4 with jointlcmm
################################################################################
paquidS$agedem65 <- (paquidS$agedem-65)/10
paquidS$std_normMMSE <- as.vector(scale(paquidS$normMMSE))


init_model <- Jointlcmm(
  std_normMMSE ~ poly(age65, degree = 2, raw = TRUE) + CEP,
  random = ~ poly(age65, degree = 2, raw = TRUE),
  subject = "ID",
  survival = Surv(agedem65, dem) ~ CEP + male,
  hazard = "Weibull", 
  ng = 1,
  data = paquidS)

final_model <- gridsearch(rep = 30, maxiter = 15, minit = init_model,
                          Jointlcmm(
                            fixed = std_normMMSE ~ poly(age65, degree = 2, raw = TRUE) + CEP,
                            mixture = ~ poly(age65, degree = 2, raw = TRUE) + CEP,
                            random = ~ poly(age65, degree = 2, raw = TRUE),
                            subject = "ID",
                            survival = Surv(agedem65, dem) ~ mixture(CEP) + mixture(male),
                            hazard = "Weibull",
                            ng = 4,
                            data = paquidS,
                            verbose = FALSE
                          ))

# extract results from the fitted JLCM
coefficients <- coef(final_model)
fit_summary <- summary(final_model)
posterior_probs <- final_model$pprob[, c("probYT1", "probYT2", "probYT3", "probYT4")]


init_function <- function() {
  list(
    var_e1 = unname((coefficients["stderr"]^2) * runif(1,0.8,1.2)),
    var_e2 = unname((coefficients["stderr"]^2) * runif(1,0.8,1.2)),
    var_e3 = unname((coefficients["stderr"]^2) * runif(1,0.8,1.2)),
    var_e4 = unname((coefficients["stderr"]^2) * runif(1,0.8,1.2)),
    beta1 = pmax(unname(c((coefficients["event1 +/-sqrt(Weibull2) class 1"]^2)*log(coefficients["event1 +/-sqrt(Weibull1) class 1"]^2) * runif(1,0.8,1.2), coefficients["CEP class1"] * runif(1,0.8,1.2), coefficients["male class1"] * runif(1,0.8,1.2), rnorm(1,0,1))),-5),
    beta2 = pmax(unname(c((coefficients["event1 +/-sqrt(Weibull2) class 2"]^2)*log(coefficients["event1 +/-sqrt(Weibull1) class 2"]^2) * runif(1,0.8,1.2), coefficients["CEP class2"] * runif(1,0.8,1.2), coefficients["male class2"] * runif(1,0.8,1.2), rnorm(1,0,1))),-5),
    beta3 = pmax(unname(c((coefficients["event1 +/-sqrt(Weibull2) class 3"]^2)*log(coefficients["event1 +/-sqrt(Weibull1) class 3"]^2) * runif(1,0.8,1.2), coefficients["CEP class3"] * runif(1,0.8,1.2), coefficients["male class3"] * runif(1,0.8,1.2), rnorm(1,0,1))),-5),
    beta4 = pmax(unname(c((coefficients["event1 +/-sqrt(Weibull2) class 4"]^2)*log(coefficients["event1 +/-sqrt(Weibull1) class 4"]^2) * runif(1,0.8,1.2), coefficients["CEP class4"] * runif(1,0.8,1.2), coefficients["male class4"] * runif(1,0.8,1.2), rnorm(1,0,1))),-5),
    alpha1 = runif(1,-5,0), 
    alpha2 = runif(1,-5,0), 
    alpha3 = runif(1,-5,0),
    alpha4 = runif(1,-5,0),
    theta1 = (matrix(fit_summary[,"coef"], 4, 4)[1,])* runif(1,0.7,1.4),
    theta2 = (matrix(fit_summary[,"coef"], 4, 4)[2,])* runif(1,0.7,1.4),
    theta3 = (matrix(fit_summary[,"coef"], 4, 4)[3,])* runif(1,0.7,1.4),
    theta4 = (matrix(fit_summary[,"coef"], 4, 4)[4,])* runif(1,0.7,1.4),
    psi1 = rnorm(3,0,1),
    psi2 = rnorm(3,0,1),
    psi3 = rnorm(3,0,1)
  )
}


init_list <- lapply(1:6, function(x) init_function())

##########################################
fit <- stan(file = "stan_4g_cs.stan", 
            data = list(N=nrow(data$longitudinal), n=nrow(data$survival), 
                        y=data$longitudinal$y, times=data$longitudinal$times, 
                        ID=rep(1:nrow(data$survival),table(data$longitudinal$id)),
                        Time=data$survival$time, 
                        status=data$survival$status, x1=data$survival$x1, x2=data$survival$x2,
                        K=K, xk=xk, wk=wk, 
                        start=data$survival$start, stop=data$survival$stop),
            init = init_list,
            warmup = 6000,                 
            iter = 14000,
            thin = 4,
            chains = 6,
            seed = 2024,
            cores = 6,
            save_warmup=F) 
##########################################


################################################################################
# chain selection
################################################################################
lp_allchain <- extract(fit, pars = "log_posterior_original", permuted = FALSE)
log_weights <- numeric(6)

for (chain in 1:6) {
  lp_chain <- lp_allchain[, chain, ]
  HPD_ind <- sort.list(lp_chain,decreasing = TRUE)[1:(0.8*(length(lp_chain)))]
  lp_hpd_chain <- lp_chain[HPD_ind]
  lp_hpd_chain_min <- min(lp_hpd_chain)
  log_weights[chain] <- lp_hpd_chain_min - log(sum(exp(-lp_hpd_chain+lp_hpd_chain_min)))
}
best_chain <- which.max(log_weights)


posterior <- as.array(fit)
posterior_chain <- posterior[,best_chain, ]  
mcmc_trace(posterior_chain, pars = c("var_e1", "var_e2", "var_e3", "var_e4", "alpha1", "alpha2", "alpha3", "alpha4","beta1[1]","beta1[2]","beta2[1]","beta3[1]", "beta4[1]", "lambda1", "lambda2", "lambda3"))
mcmc_trace(posterior_chain, pars = c( "theta1[1]","theta1[2]","theta1[3]","theta1[4]", "theta2[1]","theta2[2]","theta2[3]","theta2[4]", "theta3[1]","theta3[2]","theta3[3]","theta3[4]", "psi1[1]","psi1[2]","psi1[3]", "psi2[1]","psi2[2]","psi2[3]"))
# as expected, mixing for parameters associated with the 'tiny' clusters is poor due to weak identifiability

################################################################################
# loo, waic, and margllk estimates for model comparison 
################################################################################

loglik_allchain <- extract(fit, pars = "log_lik_cond", permuted = FALSE)
loglik_extract <- as.matrix(loglik_allchain[,best_chain,])

loo_result <- loo(loglik_extract,cores=6)
waic_result <- waic(loglik_extract,cores=6)


IC_smy <- list(
  waic_condllk=waic_result,
  loo_condllk=loo_result)


save(IC_smy, file=paste0("IC_smy_4g_cs.rdata"))


################################################################################

# Extract the `comp` parameter from the fit object
samples_comp <- extract(fit, pars = "comp", permuted = FALSE)
comp_extract <- unname(as.matrix(samples_comp[, best_chain, ]))

nsample <- dim(comp_extract)[1]
n <- dim(comp_extract)[2]
classprob<-matrix(NA,4,n)
for (i in 1:4) {
  for (j in 1:n){
    sp<-rep(NA,nsample)
    for (k in 1:(nsample)) {sp[k]<-1*(comp_extract[k,j]==i)}
    classprob[i,j]<-sum(sp)/(nsample)
  }
}
ld<-rep(NA,n)
for (i in 1:n){ld[i]<-which.max(classprob[,i])} # MAP class allocation

################################################################################
# plots of longitudinal data by cluster
################################################################################
paquid$normMMSE <- normMMSE(paquid$MMSE) # define normMMSE
paquid$age65 <- (paquid$age - 65) / 10 # define age65
paquidS <- paquid[paquid$agedem > paquid$age_init, ] # retain 499 subj

longdata <- data.frame(id=paquidS$ID, y=paquidS$normMMSE, times=paquidS$age65, x1=paquidS$CEP)
longdata <- na.omit(longdata) # NA in y is removed;  n=499 [no subj removed]
longdata$y <- as.vector(scale(longdata$y))

survdata <- data.frame(id=paquidS$ID, time=(paquidS$agedem-65)/10, status=paquidS$dem, x1=paquidS$CEP, x2=paquidS$male)
survdata <- survdata[!duplicated(survdata$id), ]
data <- list(longitudinal=longdata, survival=survdata)

# 
unique_ids <- unique(longdata$id)
ld_df <- data.frame(id = unique_ids, ld = ld)
merged_longdata <- merge(longdata, ld_df, by = "id")

ggplot(merged_longdata, aes(x = times, y = y, group = id, color = as.factor(ld))) +
  geom_line() +  
  geom_point() +
  labs(x = "age65", y = "normalized MMSE", title = "Longitudinal trajectories by estimated subgroups",
       color = "") +  # Name the legend as "Cluster"
  scale_color_manual(values = c("blue", "green", "red", "orange"), 
                     labels = c("subgroup 1", "subgroup 2", "subgroup 3", "subgroup 4")) + 
  theme_minimal() +
  theme(legend.position = "top")  


survival_data <- data$survival
survival_data$cluster <- ld

km_fit_cluster <- survfit(Surv(time, status) ~ cluster, data = survival_data)

# Generate KM plot for each cluster
km_plot_clusters <- ggsurvplot(km_fit_cluster,
                               data = survival_data,
                               xlab = "age65",
                               ylab = "Survival probability",
                               title = "Kaplan-Meier curves by estimated subgroups",
                               palette = c("blue","green", "red", "orange" ),  
                               conf.int = FALSE,  
                               legend.title = "",
                               legend.labs = c("subgroup 1", "subgroup 2", "subgroup 3", "subgroup 4"))

# Print the plot
print(km_plot_clusters)