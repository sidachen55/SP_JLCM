################################################################################
# Model 2 scenario 3: 2 sub-populations (G=2)
# cov: age, gender
################################################################################
task_id_string <- Sys.getenv("SLURM_ARRAY_TASK_ID")
task_id <- as.numeric(task_id_string)

set.seed(task_id)


library(MASS)
###################################
# function to simulate from a JM [1 subject]
###################################
sim_JM <- function(x1, x2, theta, sigma.e, Sigma.b, beta, alpha, lambda, tmax, ind){

  # Random effects
  b <- mvrnorm(1, rep(0, nrow(Sigma.b)), Sigma.b)
  
  # Simulating survival data using inverse method [Weibull baseline]
  eta <- as.vector(cbind(1, x2) %*% beta)
  u <- runif(1)
  
  invS <- function(t,u,...){
    h <- function(s,...){
      XX <- cbind(1,s,x1)
      ZZ <- cbind(1,s)
      fval <- as.vector(XX%*%theta + ZZ%*%b)
      return(exp(log(lambda)+(lambda-1)*log(s)+eta+alpha*fval))
    }
    return(integrate(h,lower=0,upper=t)$value + log(u))
  }
  
  trueTime <- NULL
  Up <- 50
  tries <- 5
  Root <- try(uniroot(invS,interval=c(1e-05,Up),u=u)$root,TRUE)
  while(inherits(Root,"try-error") && tries>0){
    tries <- tries - 1 
    Up <- Up + 500
    Root <- try(uniroot(invS,interval=c(1e-05,Up),u=u)$root,TRUE)	
  }
  trueTime <- ifelse(!inherits(Root,"try-error"), Root, NA)	
  
  
  # True time, censoring, observed time and event indicator
  trueTime <- ifelse(is.na(trueTime),tmax,trueTime)
  Ctimes <- runif(1,0,tmax) # can be changed
  Time <- min(trueTime,Ctimes)
  status <- as.numeric(trueTime<=Ctimes)
  
  # Longitudinal time points
  nmaxobs <- ceiling(tmax) # 
  nobs <- ceiling(Time) + 1 # 
  
  times <- c(seq(0,Time,len=nobs)[-nobs],rep(NA,nmaxobs-nobs+1))
  
  
  # Simulating longitudinal data
  repX <- rep(x1,times=nmaxobs)
  repb0 <- rep(b[1],times=nmaxobs)
  repb1 <- rep(b[2],times=nmaxobs)
  e <- rnorm(nmaxobs,0,sigma.e)
  Y <- theta[1] + repb0 + (theta[2]+repb1)*times + theta[3]*repX + e
  id <- rep(ind,times=nmaxobs)
  dta.long <- data.frame(id = id, y = Y, times = times, x1=repX)
  
  # Longitudinal and survival data
  dta.long <- dta.long[which(times <= rep(Time,nmaxobs)),]
  dta.surv <- data.frame(id = ind, time = Time, status = status, x1=x1, x2 = x2)
  
  return(list(longitudinal=dta.long,survival=dta.surv))
}



# Longitudinal parameters
theta1 <- c(8.03, -0.16, -5.86)
theta2 <- c(-8.03, 0.46, 12.20)
theta3 <- c(0.03, -0.01, -1.96)
sigma1.e <- sigma2.e <- sigma3.e <- 0.69
Sigma1.b <- diag(c(0.87, 0.02))
Sigma2.b <- diag(c(0.02, 0.91))
Sigma3.b <- diag(c(0.28, 0.31))

# Survival parameters
beta1 <- c(-4.85, -0.02)
beta2 <- c(-4.85, 0.09)
beta3 <- c(2.85, -0.12)
alpha1 <- 0.38
alpha2 <- 0.08
alpha3 <- 0.58
lambda1 <- 1.8
lambda2 <- 1.4
lambda3 <- 1.8
tmax <- 19.5

# cluster membership model parameters [class 2 is ref cat]
phi_10 <- 2  # Intercept for class 1
phi_11 <- 4   # class 1 coefficient for x1
phi_12 <- -0.1  # class 2 coefficient for x2

################################################################################
# data generation
################################################################################
n <- 900 # total sample size

# baseline covariates for each subject
x1 <- rep(NA,n)
for (i in 1:n) {x1[i] <- sample(c(0,1), 1, replace = TRUE)} # in the longitudinal submodel
x2 <- rnorm(n, 45, 15.69578) # in the survival submodel

# generate membership prob matrix
cluster_prob <- matrix(NA, n, 2)
for (i in 1:n) {
  linear_predictor <- phi_10 + x1[i] * phi_11 + x2[i] * phi_12
  exp_linear_predictor <- exp(linear_predictor)
  cluster_prob[i,1] <- exp_linear_predictor / (1 + exp_linear_predictor)
  cluster_prob[i,2] <- 1 - cluster_prob[i,1]
}


# generate cluster indicator
s<-rep(NA,n)
for (k in 1:n){s[k]<-sample(1:2,1,prob=cluster_prob[k,])}


table(s)



#
for (k in 1:1){
  if (s[k]==1) {
    data1 <- sim_JM(x1[k], x2[k], theta=theta1, sigma.e=sigma1.e, Sigma.b=Sigma1.b,
                    beta=beta1, alpha=alpha1, lambda=lambda1, tmax=tmax, ind=k)
    longdata <- data1$longitudinal
    survdata <- data1$survival
  } else if (s[k]==2) {
    data1 <- sim_JM(x1[k], x2[k], theta=theta2, sigma.e=sigma2.e, Sigma.b=Sigma2.b,
                    beta=beta2, alpha=alpha2, lambda=lambda2, tmax=tmax, ind=k)
    longdata <- data1$longitudinal
    survdata <- data1$survival
  }
}

for (k in 2:n){
  if (s[k]==1) {
    data1 <- sim_JM(x1[k], x2[k], theta=theta1, sigma.e=sigma1.e, Sigma.b=Sigma1.b,
                    beta=beta1, alpha=alpha1, lambda=lambda1, tmax=tmax, ind=k)
    longdata <- rbind(longdata, data1$longitudinal)
    survdata <- rbind(survdata, data1$survival)
  } else if (s[k]==2) {
    data1 <- sim_JM(x1[k], x2[k], theta=theta2, sigma.e=sigma2.e, Sigma.b=Sigma2.b,
                    beta=beta2, alpha=alpha2, lambda=lambda2, tmax=tmax, ind=k)
    longdata <- rbind(longdata, data1$longitudinal)
    survdata <- rbind(survdata, data1$survival)
  }
}

data <- list(longitudinal=longdata, survival=survdata)




# Start and stop rows for each individual 
data$survival$start <- data$survival$stop <- NA
aux <- 0
for(i in 1:nrow(data$survival)){
  pos <- which(data$longitudinal$id == data$survival$id[i])
  data$survival$start[i] <- aux + 1
  data$survival$stop[i] <- aux + length(pos)
  aux <- data$survival$stop[i]
}

summary(as.numeric(table(data$longitudinal$id)))
table(data$survival$status)
table(data$survival$status)[1]/length(data$survival$id) # censoring rate
summary(data$survival$time)

save.image(paste0("data_",task_id,".RData"))