library(lcmm)
library(NormPsy)
library(JMbayes)
library(rjags)
library(xtable)
library(splines)

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


surv_data <- data.frame(id=data$survival$id, time=data$survival$time, status=data$survival$status, x1=data$survival$x1, x2=data$survival$x2, x12=data$survival$x1*data$survival$x2)
surv_data_subset <- surv_data[, c("id", "time", "status", "x2" ,"x12")]
long_data <- merge(data$longitudinal, surv_data_subset, by = "id")

data <- data.frame(echotime=long_data$times, x1=long_data$x1, IDnr=long_data$id, y0=long_data$y, years=long_data$time, status=long_data$status, x2=long_data$x2, x12=long_data$x12)
data.id <- data.frame(IDnr=surv_data$id, years=surv_data$time, status=surv_data$status, x1=surv_data$x1, x2=surv_data$x2, x12=surv_data$x12)

##################################################################################
# Longitudinal submodel
fm1 <- lme(y0 ~ x1 + echotime + I(echotime^2), data = data,
           na.action = na.omit,
           random =~ echotime + I(echotime^2) | IDnr)
lmeObject <- fm1

timeVar <- "echotime"
lag <- 0
survMod <- "spline-PH"

id <- data$IDnr 
offset <- as.vector(c(1, 1 + cumsum(tapply(id, id, length))))

Time <- data.id$years

# Survival submodel
W <- model.matrix(~ -1 + x1 + x2 + x12, data.id)

# Design matrices
formYx <- formula(lmeObject)
TermsX <- lmeObject$terms
mfX <- model.frame(TermsX, data = data)
X <- model.matrix(formYx, mfX)

formYz <- formula(lmeObject$modelStruct$reStruct[[1]])
mfZ <- model.frame(terms(formYz), data = data)
TermsZ <- attr(mfZ, "terms")
Z <- model.matrix(formYz, mfZ)

#
data.id <- data[!duplicated(id), ]
data.id[[timeVar]] <- pmax(Time - 0, 0)


mfX.id <- model.frame(TermsX, data = data.id)  
mfZ.id <- model.frame(TermsZ, data = data.id)  
Xtime <- model.matrix(formYx, mfX.id)
Ztime <- model.matrix(formYz, mfZ.id)

###
####################################################
eventD <- data.id$status
nT <- length(Time)
zeros <- numeric(nT)

y.long <- model.response(mfX, "numeric")
y <- list(y = y.long, offset = offset, logT = log(Time),
          eventD = eventD, zeros = zeros, lag = lag)

# Time and data for approximation
gaussKronrod <- JMbayes:::gaussKronrod
wk <- gaussKronrod()$wk
sk <- gaussKronrod()$sk

ordsk <- order(sk)
sk <- sk[ordsk]
wk <- wk[ordsk]

K <- length(sk)

P <- Time/2
st <- outer(P, sk + 1)
id.GK <- rep(seq_along(Time), each = K)

data.id2 <- data.id[id.GK, ]
data.id2[[timeVar]] <- c(t(st))

# Design matrices for approximation
mfX <- model.frame(TermsX, data = data.id2)   
mfZ <- model.frame(TermsZ, data = data.id2)   
Xs <- model.matrix(formYx, mfX)
Zs <- model.matrix(formYz, mfZ)

################################################################################

con <- list(n.chains = 1, n.iter = 50000,
            n.burnin = 25000, n.thin = 10, n.adapt = 500, K = 1000,
            C = 5000,knots = NULL, ObsTimes.knots = TRUE, lng.in.kn = 5, ordSpline = 3, bugs.seed = 1, quiet = FALSE)

x <- list(X = X, Z = Z,  W = if (survMod == "weibull-PH") {
  if (is.null(W)) cbind(rep(1, nT), rep(0, nT)) else cbind(1,
                                                           W)
} else {
  if (is.null(W)) cbind(rep(0, nT), rep(0, nT)) else {
    if (ncol(W) == 1) cbind(W, rep(0, nT)) else W
  }
}
) 

# Baseline hazard
kn <- if (is.null(con$knots)) {
  pp <- seq(0, 1, length.out = con$lng.in.kn + 2)
  pp <- tail(head(pp, -1), -1)
  tt <- if (con$ObsTimes.knots) {
    Time
  } else {  Time[event == 1]    }
  quantile(tt, pp, names = FALSE)
} else {
  con$knots
}
kn <- kn[kn < max(Time)]
rr <- sort(c(rep(range(Time, st), con$ordSpline), kn))
con$knots <- rr

WBH <- splineDesign(rr, Time, ord = con$ordSpline)
if (any(colSums(WBH) == 0))
  stop("\nsome of the knots of the B-splines basis are set outside the range",
       "\n   of the observed event times for one of the strata; refit the model",
       "\n   setting the control argument 'equal.strata.knots' to FALSE.")

# design matrices for the baseline hazard for the 15-point Gauss-Kronrod quadrature rule approximation
WBHs <- splineDesign(rr, c(t(st)), ord = con$ordSpline)

x <- c(x, list(WBH = WBH, WBHs = WBHs))

ncX <- ncol(X)
ncZ <- ncol(Z)
ncW <- ncol(x$W)
ncWBH <- ncol(x$WBH)
ncWBHs <- ncol(x$WBHs)
C <- con$C
nb <- ncZ 


##
betas1 <- rep(0, ncX)
betas2 <- rep(0, ncX)
betas3 <- rep(0, ncX)
betas4 <- rep(0, ncX)
betas5 <- rep(0, ncX)
betas6 <- rep(0, ncX)

var.betas1 <- rep(con$K, ncX)
var.betas2 <- rep(con$K, ncX)
var.betas3 <- rep(con$K, ncX)
var.betas4 <- rep(con$K, ncX)
var.betas5 <- rep(con$K, ncX)
var.betas6 <- rep(con$K, ncX)

alphas1 <- 0
alphas2 <- 0
alphas3 <- 0
alphas4 <- 0
alphas5 <- 0
alphas6 <- 0

var.alphas1 <- 100
var.alphas2 <- 100
var.alphas3 <- 100
var.alphas4 <- 100
var.alphas5 <- 100
var.alphas6 <- 100

gammas1 <- rep(0,(ncW))
gammas2 <- rep(0,(ncW))
gammas3 <- rep(0,(ncW))
gammas4 <- rep(0,(ncW))
gammas5 <- rep(0,(ncW))
gammas6 <- rep(0,(ncW))

var.gammas1 <- rep(con$K, (ncW))
var.gammas2 <- rep(con$K, (ncW))
var.gammas3 <- rep(con$K, (ncW))
var.gammas4 <- rep(con$K, (ncW))
var.gammas5 <- rep(con$K, (ncW))
var.gammas6 <- rep(con$K, (ncW))

Bs.gammas1 <- rep(0, (ncWBH))
Bs.gammas2 <- rep(0, (ncWBH))
Bs.gammas3 <- rep(0, (ncWBH))
Bs.gammas4 <- rep(0, (ncWBH))
Bs.gammas5 <- rep(0, (ncWBH))
Bs.gammas6 <- rep(0, (ncWBH))

var.Bs.gammas1 <- rep(con$K, (ncWBH))
var.Bs.gammas2 <- rep(con$K, (ncWBH))
var.Bs.gammas3 <- rep(con$K, (ncWBH))
var.Bs.gammas4 <- rep(con$K, (ncWBH))
var.Bs.gammas5 <- rep(con$K, (ncWBH))
var.Bs.gammas6 <- rep(con$K, (ncWBH))

mu01 <- rep(0,(ncZ))
mu02 <- rep(0,(ncZ))
mu03 <- rep(0,(ncZ))
mu04 <- rep(0,(ncZ))
mu05 <- rep(0,(ncZ))
mu06 <- rep(0,(ncZ))

b <- cbind(data.matrix(ranef(lmeObject)))
nY <- nrow(b)


Data <- list(N = nY, K = K, offset = offset, X = X, Xtime = Xtime, ncX = ncol(X), ncZ = ncol(Z), 
             y = y$y, 
             Xs = Xs, 
             Z = Z, Ztime = Ztime,  
             Zs = Zs, 
             event = eventD, zeros = zeros,  
             C = C, P = P,
             wk = wk, nb = nb, 
             W = x$W,
             ncW = ncol(x$W), 
             WBH = WBH, 
             ncWBH = ncol(x$WBH),
             WBHs = WBHs, 
             ncWBHs = ncol(x$WBH),priorA.tau = 0.01,
             priorB.tau = 0.01, prior.cl = rep(6.9,6),
             mu01 = mu01, mu02 = mu02, mu03 = mu03, mu04 = mu04, mu05 = mu05, mu06 = mu06,
             priorMean.betas1 = betas1, priorMean.betas2 = betas2, priorMean.betas3 = betas3, priorMean.betas4 = betas4, priorMean.betas5 = betas5, priorMean.betas6 = betas6,
             priorTau.betas1 = diag(1/var.betas1), priorTau.betas2 = diag(1/var.betas2), priorTau.betas3 = diag(1/var.betas3), priorTau.betas4 = diag(1/var.betas4), priorTau.betas5 = diag(1/var.betas5), priorTau.betas6 = diag(1/var.betas6),
             priorMean.gammas1 = gammas1, priorMean.gammas2 = gammas2, priorMean.gammas3 = gammas3, priorMean.gammas4 = gammas4, priorMean.gammas5 = gammas5, priorMean.gammas6 = gammas6,
             priorTau.gammas1 = diag(1/var.gammas1), priorTau.gammas2 = diag(1/var.gammas2), priorTau.gammas3 = diag(1/var.gammas3), priorTau.gammas4 = diag(1/var.gammas4), priorTau.gammas5 = diag(1/var.gammas5), priorTau.gammas6 = diag(1/var.gammas6),
             priorMean.alphas1 = alphas1, priorMean.alphas2 = alphas2, priorMean.alphas3 = alphas3, priorMean.alphas4 = alphas4, priorMean.alphas5 = alphas5, priorMean.alphas6 = alphas6,
             priorTau.alphas1 = 1/var.alphas1, priorTau.alphas2 = 1/var.alphas2, priorTau.alphas3 = 1/var.alphas3, priorTau.alphas4 = 1/var.alphas4, priorTau.alphas5 = 1/var.alphas5, priorTau.alphas6 = 1/var.alphas6,
             priorMean.Bs.gammas1 = Bs.gammas1, priorMean.Bs.gammas2 = Bs.gammas2, priorMean.Bs.gammas3 = Bs.gammas3, priorMean.Bs.gammas4 = Bs.gammas4, priorMean.Bs.gammas5 = Bs.gammas5, priorMean.Bs.gammas6 = Bs.gammas6,
             priorTau.Bs.gammas1 = diag(1/var.Bs.gammas1), priorTau.Bs.gammas2 = diag(1/var.Bs.gammas2), priorTau.Bs.gammas3 = diag(1/var.Bs.gammas3), priorTau.Bs.gammas4 = diag(1/var.Bs.gammas4), priorTau.Bs.gammas5 = diag(1/var.Bs.gammas5), priorTau.Bs.gammas6 = diag(1/var.Bs.gammas6))

Data$inv.D <- matrix(0, nb, nb)
diag(Data$inv.D) <- NA

parms <- c("pr", "v", "tau", "alphas1", "alphas2", "inv.D", "gammas1", "betas3")


# Model in JAGS
model.fit <- jags.model(file = "JMsuper.txt", data = Data, n.chains = con$n.chains, 
                        n.adapt = con$n.adapt, quiet = con$quiet)

update(model.fit, con$n.burnin)
res <- coda.samples(model.fit, parms,  n.iter = con$n.iter - con$n.burnin, thin = con$n.thin)

####
# Extract posterior samples of v from your run
bss <- do.call(rbind, res)
v_index <- grep("^v\\[", colnames(bss))
W <- as.matrix(bss[, v_index])  # Rows = posterior samples, Columns = subjects

# Number of subjects and classes
N <- ncol(W)
Class <- 6

# Function from author: mode
Mod <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# Compute how many classes are "tiny" (used by ≤ 1% of subjects)
EPTC_0.01 <- apply(W, 1, function(x){
  sum(sum(x == 1) <= N/100,
      sum(x == 2) <= N/100,
      sum(x == 3) <= N/100,
      sum(x == 4) <= N/100,
      sum(x == 5) <= N/100,
      sum(x == 6) <= N/100)
})
cl_0.01 <- Class - Mod(EPTC_0.01)

EPTC_0.05 <- apply(W, 1, function(x){
  sum(sum(x == 1) <= 5*N/100,
      sum(x == 2) <= 5*N/100,
      sum(x == 3) <= 5*N/100,
      sum(x == 4) <= 5*N/100,
      sum(x == 5) <= 5*N/100,
      sum(x == 6) <= 5*N/100)
})
cl_0.05 <- Class - Mod(EPTC_0.05)

EPTC_0.10 <- apply(W, 1, function(x){
  sum(sum(x == 1) <= 10*N/100,
      sum(x == 2) <= 10*N/100,
      sum(x == 3) <= 10*N/100,
      sum(x == 4) <= 10*N/100,
      sum(x == 5) <= 10*N/100,
      sum(x == 6) <= 10*N/100)
})
cl_0.10 <- Class - Mod(EPTC_0.10)

EPTC_0.15 <- apply(W, 1, function(x){
  sum(sum(x == 1) <= 15*N/100,
      sum(x == 2) <= 15*N/100,
      sum(x == 3) <= 15*N/100,
      sum(x == 4) <= 15*N/100,
      sum(x == 5) <= 15*N/100,
      sum(x == 6) <= 15*N/100)
})
cl_0.15 <- Class - Mod(EPTC_0.15)

fit_smy <- list(cl_0.01=cl_0.01, cl_0.05=cl_0.05, cl_0.10=cl_0.10, cl_0.15=cl_0.15)













