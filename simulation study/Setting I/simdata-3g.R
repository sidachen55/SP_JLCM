################################################################################
# scenario: 3 sub-populations (G=3)
################################################################################
task_id_string <- Sys.getenv("SLURM_ARRAY_TASK_ID")
task_id <- as.numeric(task_id_string)

set.seed(task_id)


library(MASS)
###################################
# function to simulate from a JM 
###################################
sim_JM <- function(n, theta, sigma.e, Sigma.b, beta, alpha, lambda, tmax, init=1){
  
  # baseline covariates for each subject
  x1 <- sample(c(0,1), n, replace = TRUE) # in the longitudinal submodel
  x2 <- rnorm(n, 45, 15.69578) # in the survival submodel
  
  # Random effects
  b <- mvrnorm(n, rep(0, nrow(Sigma.b)), Sigma.b)
  
  # Simulating survival data using inverse method [Weibull baseline]
  eta <- as.vector(cbind(1, x2) %*% beta)
  u <- runif(n)
  
  invS <- function(t,u,i,...){
    h <- function(s,...){
      XX <- cbind(1,s,x1[i])
      ZZ <- cbind(1,s)
      fval <- as.vector(XX%*%theta + rowSums(ZZ*b[rep(i,nrow(ZZ)),]))
      return(exp(log(lambda)+(lambda-1)*log(s)+eta[i]+alpha*fval))
    }
    return(integrate(h,lower=0,upper=t)$value + log(u))
  }
  
  trueTime <- NULL
  for(i in 1:n){
    Up <- 50
    tries <- 5
    Root <- try(uniroot(invS,interval=c(1e-05,Up),u=u[i],i=i)$root,TRUE)
    while(inherits(Root,"try-error") && tries>0){
      tries <- tries - 1 
      Up <- Up + 500
      Root <- try(uniroot(invS,interval=c(1e-05,Up),u=u[i],i=i)$root,TRUE)	
    }
    trueTime[i] <- ifelse(!inherits(Root,"try-error"), Root, NA)	
  }
  
  # True time, censoring, observed time and event indicator
  trueTime <- ifelse(is.na(trueTime),tmax,trueTime)
  Ctimes <- runif(n,0,tmax) # can be changed
  Time <- pmin(trueTime,Ctimes)
  status <- as.numeric(trueTime<=Ctimes)
  
  # Longitudinal time points
  nmaxobs <- ceiling(tmax) # 
  nobs <- ceiling(Time) + 1 # 
  times <- NULL
  for(i in 1:n){
    times <- c(times,c(seq(0,Time[i],len=nobs[i])[-nobs[i]],rep(NA,nmaxobs-nobs[i]+1)))
  }
  
  # Simulating longitudinal data
  repX <- rep(x1,times=rep(nmaxobs,n))
  repb0 <- rep(b[,1],times=rep(nmaxobs,n))
  repb1 <- rep(b[,2],times=rep(nmaxobs,n))
  e <- rnorm(n*nmaxobs,0,sigma.e)
  Y <- theta[1] + repb0 + (theta[2]+repb1)*times + theta[3]*repX + e
  id <- rep(init:(init+n-1),times=rep(nmaxobs,n))
  dta.long <- data.frame(id = id, y = Y, times = times,x1=repX)
  
  # Longitudinal and survival data
  dta.long <- dta.long[which(times <= rep(Time,rep(nmaxobs,n))),]
  dta.surv <- data.frame(id = unique(id), time = Time, status = status, x1=x1,x2 = x2)
  
  return(list(longitudinal=dta.long,survival=dta.surv, b=b))
}



#################################################################
# parameter setting (motivated from Andrinopoulou et al (2020))
#################################################################
n1 <- 100
n2 <- 300
n3 <- 500

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








# Simulating from JM1

data1 <- sim_JM(n=n1, theta=theta1, sigma.e=sigma1.e, Sigma.b=Sigma1.b,
                beta=beta1, alpha=alpha1, lambda=lambda1, tmax=tmax)

# Simulating from JM2

data2 <- sim_JM(n=n2, theta=theta2, sigma.e=sigma2.e, Sigma.b=Sigma2.b,
                beta=beta2, alpha=alpha2, lambda=lambda2, tmax=tmax, init=n1+1)

# Simulating from JM3
data3 <- sim_JM(n=n3, theta=theta3, sigma.e=sigma3.e, Sigma.b=Sigma3.b,
                beta=beta3, alpha=alpha3, lambda=lambda3, tmax=tmax, init=n1+n2+1)


#### 
data <- list(longitudinal=rbind(data1$longitudinal, data2$longitudinal, data3$longitudinal), survival=rbind(data1$survival, data2$survival, data3$survival))


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