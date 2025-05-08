## ----setup, include=FALSE-----------------------------------------------------
rmarkdown::find_pandoc(version = '2.9.2.1')
knitr::opts_chunk$set(echo = TRUE)

## ----IV1, echo = TRUE, eval = TRUE--------------------------------------------
library(PartialNetwork)
set.seed(123)
# Number of groups
M             <- 30
# size of each group
N             <- rep(50,M)
# individual effects
beta          <- c(2,1,1.5) 
# endogenous effects
alpha         <- 0.4
# std-dev errors
se            <- 1
# network distribution
distr         <- runif(sum(N*(N-1)))
distr         <- vec.to.mat(distr, N, normalise = FALSE)
# covariates
X             <- cbind(rnorm(sum(N),0,5),rpois(sum(N),7))
# true network
G0            <- sim.network(distr)
# normalise 
G0norm        <- norm.network(G0)
# simulate dependent variable use an external package
y             <- CDatanet::simsar(~ X, Glist = G0norm, theta = c(alpha, beta, se))
y             <- y$y
# generate instruments 
instr         <- sim.IV(distr, X, y, replication = 1, power = 1)
GY1c1         <- instr[[1]]$G1y       # proxy for Gy (draw 1)
GXc1          <- instr[[1]]$G1X[,,1]  # proxy for GX (draw 1)
GXc2          <- instr[[1]]$G2X[,,1]  # proxy for GX (draw 2)
# build dataset
# keep only instrument constructed using a different draw than the one used to proxy Gy
dataset           <- as.data.frame(cbind(y, X, GY1c1, GXc1, GXc2)) 
colnames(dataset) <- c("y","X1","X2","G1y", "G1X1", "G1X2", "G2X1", "G2X2") 

## ----IV2, echo = TRUE, eval = TRUE, message=FALSE-----------------------------
library(AER)
# Same draws
out.iv1           <- ivreg(y ~ X1 + X2 + G1y | X1 + X2 + G1X1 + G1X2, data = dataset)
summary(out.iv1)

## ----IV3, echo = TRUE, eval = TRUE, message=FALSE-----------------------------
# Different draws
out.iv2           <- ivreg(y ~ X1 + X2 + G1y | X1 + X2 + G2X1 + G2X2, data = dataset)
summary(out.iv2)

## ----IV3bias------------------------------------------------------------------
ddS     <- as.matrix(cbind(1, dataset[,c("X1", "X2", "G1y")]))         #\ddot{S} 
dZ      <- as.matrix(cbind(1, dataset[,c("X1", "X2", "G2X1", "G2X2")]))#\dot{Z} 
dZddS   <- crossprod(dZ, ddS)/sum(N)                                  
W       <- solve(crossprod(dZ)/sum(N))                                 
matM    <- solve(crossprod(dZddS, W%*%dZddS), crossprod(dZddS, W))    
maxbias <- apply(sapply(1:1000, function(...){
  dddGy <- peer.avg(sim.network(distr, normalise = TRUE) , y)
  abs(matM%*%crossprod(dZ, dddGy - dataset$G1y)/sum(N))
}), 1, max); names(maxbias) <- c("(Intercept)", "X1", "X2", "G1y")
{cat("Maximal absolute bias\n"); print(maxbias)}

## ----IV4, echo = TRUE, eval = TRUE--------------------------------------------
rm(list = ls())
library(PartialNetwork)
set.seed(123)
# Number of groups
M             <- 30
# size of each group
N             <- rep(50,M)
# individual effects
beta          <- c(2,1,1.5) 
# contextual effects
gamma         <- c(5, -3)
# endogenous effects
alpha         <- 0.4
# std-dev errors
se            <- 1
# network distribution
distr         <- runif(sum(N*(N-1)))
distr         <- vec.to.mat(distr, N, normalise = FALSE)
# covariates
X             <- cbind(rnorm(sum(N),0,5),rpois(sum(N),7))
# true network
G0            <- sim.network(distr)
# normalise 
G0norm        <- norm.network(G0)
# GX
GX            <- peer.avg(G0norm, X)
# simulate dependent variable use an external package
y             <- CDatanet::simsar(~ X + GX, Glist = G0norm, 
                                  theta = c(alpha, beta, gamma, se))
y             <- y$y
# generate instruments 
# we need power = 2 for models with contextual effetcs
instr         <- sim.IV(distr, X, y, replication = 1, power = 2)
GY1c1         <- instr[[1]]$G1y       # proxy for Gy (draw 1)
GXc1          <- instr[[1]]$G1X[,,1]  # proxy for GX (draw 1)
GXc2          <- instr[[1]]$G2X[,,1]  # proxy for GX (draw 2)
GXc2sq        <- instr[[1]]$G2X[,,2]  # proxy for G^2X (draw 2)
# build dataset
# keep only instrument constructed using a different draw than the one used to proxy Gy
dataset           <- as.data.frame(cbind(y, X, GX, GY1c1, GXc1, GXc2, GXc2sq)) 
colnames(dataset) <- c("y","X1","X2", "GX1", "GX2", "G1y", "G1X1", "G1X2", "G2X1", "G2X2",
                       "G2X1sq", "G2X2sq") 

## ----IV5, echo = TRUE, eval = TRUE, message=FALSE-----------------------------
# Different draws
out.iv2           <- ivreg(y ~ X1 + X2 + GX1 + GX2 + G1X1 + G1X2 + G1y | X1 + X2 + GX1 + 
                             GX2 + G1X1 + G1X2 + G2X1 + G2X2 + G2X1sq + G2X2sq, 
                           data = dataset)
summary(out.iv2)

## ----IV5bias------------------------------------------------------------------
ddS     <- as.matrix(cbind(1, dataset[,c("X1", "X2", "GX1", "GX2", "G1X1", "G1X2", 
                                         "G1y")]))                
dZ      <- as.matrix(cbind(1, dataset[,c("X1", "X2", "GX1", "GX2", "G1X1", 
                                         "G1X2", "G2X1", "G2X2", "G2X1sq", "G2X2sq")]))
dZddS   <- crossprod(dZ, ddS)/sum(N)                                  
W       <- solve(crossprod(dZ)/sum(N))                                 
matM    <- solve(crossprod(dZddS, W%*%dZddS), crossprod(dZddS, W))    
maxbias <- apply(sapply(1:1000, function(...){
  dddGy <- peer.avg(sim.network(distr, normalise = TRUE) , y)
  abs(matM%*%crossprod(dZ, dddGy - dataset$G1y)/sum(N))
}), 1, max); names(maxbias) <- c("(Intercept)", "X1", "X2", "GX1", "GX2", "G1X1", 
                                 "G1X2", "G1y")
{cat("Maximal absolute bias\n"); print(maxbias)}

## ----smm1, echo = TRUE, eval = TRUE-------------------------------------------
rm(list = ls())
library(PartialNetwork)
set.seed(123)
# Number of groups
M             <- 100
# size of each group
N             <- rep(30,M)
# individual effects
beta          <- c(2, 1, 1.5, 5, -3) 
# endogenous effects
alpha         <- 0.4
# std-dev errors
se            <- 1 
# network distribution
distr         <- runif(sum(N*(N-1)))
distr         <- vec.to.mat(distr, N, normalise = FALSE)
# covariates
X             <- cbind(rnorm(sum(N),0,5),rpois(sum(N),7))
# true network
G0            <- sim.network(distr)
# normalise 
G0norm        <- norm.network(G0)
# Matrix GX
GX            <- peer.avg(G0norm, X)
# simulate dependent variable use an external package
y             <- CDatanet::simsar(~ X + GX, Glist = G0norm, theta = c(alpha, beta, se))
Gy            <- y$Gy
y             <- y$y
# build dataset
dataset           <- as.data.frame(cbind(y, X, Gy, GX)) 
colnames(dataset) <- c("y","X1","X2", "Gy", "GX1", "GX2") 

## ----Smmload, echo = FALSE, eval = TRUE, message=FALSE------------------------
load(url('https://raw.githubusercontent.com/ahoundetoungan/PartialNetwork/master/datavignettes/smm.rda', "rb"))

## ----Smm2, echo = TRUE, eval = FALSE, message=FALSE---------------------------
# out.smm1       <- smmSAR(y ~ X1 + X2 | Gy | GX1 + GX2, dnetwork = distr, contextual = T,
#                          smm.ctr  = list(R = 1, print = F), data = dataset)
# summary(out.smm1)

## ----Smm2a, echo = FALSE, eval = TRUE, message=FALSE--------------------------
print(smm$smm1)

## ----Smm3, echo = TRUE, eval = FALSE, message=FALSE---------------------------
# out.smm2       <- smmSAR(y ~ X1 + X2 || GX1 + GX2, dnetwork = distr, contextual = T,
#                          smm.ctr  = list(R = 1, print = F), data = dataset)
# summary(out.smm2)

## ----Smm3a, echo = FALSE, eval = TRUE, message=FALSE--------------------------
print(smm$smm2)

## ----Smm4, echo = TRUE, eval = FALSE, message=FALSE---------------------------
# out.smm3       <- smmSAR(y ~ X1 + X2 | Gy, dnetwork = distr, contextual = T,
#                          smm.ctr  = list(R = 100, print = F), data = dataset)
# summary(out.smm3)

## ----Smm4a, echo = FALSE, eval = TRUE, message=FALSE--------------------------
print(smm$smm3)

## ----Smm5, echo = TRUE, eval = FALSE, message=FALSE---------------------------
# out.smm4       <- smmSAR(y ~ X1 + X2, dnetwork = distr, contextual = T,
#                          smm.ctr  = list(R = 100, print = F), data = dataset)
# summary(out.smm4)

## ----Smm5a, echo = FALSE, eval = TRUE, message=FALSE--------------------------
print(smm$smm4)

## ----smm6a, echo = TRUE, eval = FALSE-----------------------------------------
# rm(list = ls())

## ----smm6b, echo = TRUE, eval = TRUE------------------------------------------
library(PartialNetwork)
set.seed(123)
# Number of groups
M             <- 200
# size of each group
N             <- rep(30,M)
# individual effects
beta          <- c(1, 1, 1.5, 5, -3)
# endogenous effects
alpha         <- 0.4
# std-dev errors
se            <- 1
# network distribution
distr         <- runif(sum(N*(N-1)))
distr         <- vec.to.mat(distr, N, normalise = FALSE)
# covariates
X             <- cbind(rnorm(sum(N),0,5), rpois(sum(N),7))
# Groups' fixed effects
# In order to have groups' heterogeneity correlated to X (fixed effects),
# We consider the quantile of X2 at 25% in the group
eff           <- unlist(lapply(1:M, function(x) 
  rep(quantile(X[(c(0, cumsum(N))[x]+1):(cumsum(N)[x]),2], probs = 0.25), each = N[x])))
print(c("cor(eff, X1)" = cor(eff, X[,1]), "cor(eff, X2)" = cor(eff, X[,2]))) 
# We can see that eff is correlated to X2. We can confirm that the correlation is 
# strongly significant.
print(c("p.value.cor(eff, X1)" = cor.test(eff, X[,1])$p.value,
        "p.value.cor(eff, X2)" = cor.test(eff, X[,2])$p.value))
# true network
G0            <- sim.network(distr)
# normalise 
G0norm        <- norm.network(G0)
# Matrix GX
GX            <- peer.avg(G0norm, X)
# simulate dependent variable use an external package
y             <- CDatanet::simsar(~ -1 + eff + X + GX, Glist = G0norm, 
                                  theta = c(alpha, beta, se))
Gy            <- y$Gy
y             <- y$y
# build dataset
dataset           <- as.data.frame(cbind(y, X, Gy, GX)) 
colnames(dataset) <- c("y","X1","X2", "Gy", "GX1", "GX2") 

## ----Smm7, echo = TRUE, eval = FALSE, message=FALSE---------------------------
# out.smmeff1 <- smmSAR(y ~ X1 + X2 || GX1 + GX2, dnetwork = distr, contextual = T,
#                       fixed.effects = T, smm.ctr  = list(R = 1, print = F),
#                       data = dataset)
# summary(out.smmeff1)

## ----Smm7a, echo = FALSE, eval = TRUE, message=FALSE--------------------------
print(smm$smme1)

## ----Smm8, echo = TRUE, eval = FALSE, message=FALSE---------------------------
# out.smmeff2 <- smmSAR(y ~ X1 + X2 | Gy, dnetwork = distr, contextual = T,
#                       fixed.effects = T, smm.ctr  = list(R = 100, print = F),
#                       data = dataset)
# summary(out.smmeff2)

## ----Smm8a, echo = FALSE, eval = TRUE, message=FALSE--------------------------
print(smm$smme2)

## ----Smm9, echo = TRUE, eval = FALSE, message=FALSE---------------------------
# out.smmeff3 <- smmSAR(y ~ X1 + X2, dnetwork = distr, contextual = T, fixed.effects = T,
#                       smm.ctr  = list(R = 100, print = F), data = dataset)
# summary(out.smmeff3)

## ----Smm9a, echo = FALSE, eval = TRUE, message=FALSE--------------------------
print(smm$smme3)

## ----Smm10, echo = TRUE, eval = FALSE, message=FALSE--------------------------
# fMC <- function(...){
#   # Number of groups
#   M             <- 200
#   # size of each group
#   N             <- rep(30,M)
#   # individual effects
#   beta          <- c(1, 1, 1.5, 5, -3)
#   # endogenous effects
#   alpha         <- 0.4
#   # std-dev errors
#   se            <- 1
#   # network distribution
#   distr         <- runif(sum(N*(N-1)))
#   distr         <- vec.to.mat(distr, N, normalise = FALSE)
#   # covariates
#   X             <- cbind(rnorm(sum(N),0,5), rpois(sum(N),7))
#   # Groups' fixed effects
#   # We defined the groups' fixed effect as the quantile at 25% of X2 in the group
#   # This implies that the effects are correlated with X
#   eff           <- unlist(lapply(1:M, function(x)
#     rep(quantile(X[(c(0, cumsum(N))[x]+1):(cumsum(N)[x]),2], probs = 0.25), each = N[x])))
#   # true network
#   G0            <- sim.network(distr)
#   # normalise
#   G0norm        <- norm.network(G0)
#   # Matrix GX
#   GX            <- peer.avg(G0norm, X)
#   # simulate dependent variable use an external package
#   y             <- CDatanet::simsar(~ -1 + eff + X + GX, Glist = G0norm,
#                                     theta = c(alpha, beta, se))
#   Gy            <- y$Gy
#   y             <- y$y
#   # build dataset
#   dataset           <- as.data.frame(cbind(y, X, Gy, GX))
#   colnames(dataset) <- c("y","X1","X2", "Gy", "GX1", "GX2")
#   out.smmeff1 <- smmSAR(y ~ X1 + X2 || GX1 + GX2, dnetwork = distr, contextual = T,
#                         fixed.effects = T, smm.ctr  = list(R = 1, print = F),
#                         data = dataset)
#   out.smmeff2 <- smmSAR(y ~ X1 + X2 | Gy, dnetwork = distr, contextual = T,
#                         fixed.effects = T, smm.ctr  = list(R = 100, print = F),
#                         data = dataset)
#   out.smmeff3 <- smmSAR(y ~ X1 + X2, dnetwork = distr, contextual = T, fixed.effects = T,
#                       smm.ctr  = list(R = 100, print = F), data = dataset)
#   out         <- data.frame("GX.observed"   = out.smmeff1$estimates,
#                             "Gy.observed"   = out.smmeff2$estimates,
#                             "None.observed" = out.smmeff3$estimates)
#   out
# }

## ----Smm10a, echo = TRUE, eval = FALSE, message=FALSE-------------------------
# smm.Monte.C   <- lapply(1:250, fMC)

## ----Smm10b, echo = TRUE, eval = FALSE, message=FALSE-------------------------
# Reduce('+', smm.Monte.C)/250

## ----Smm10c, echo = FALSE, eval = TRUE, message=FALSE-------------------------
print(smm$smmmc)

## ----smmsave, echo = FALSE, eval = FALSE--------------------------------------
# smm <- list(smm1 = out.smm1, smm2 = out.smm2,
#             smm3 = out.smm3, smm4 = out.smm4,
#             smme1 = out.smmeff1, smme2 = out.smmeff2,
#             smme3 = out.smmeff3, smmmc = smm$smmmc)
# save(smm, file = "~/Dropbox/Papers - In progress/Partial Network/Package/AH/PartialNetwork/datavignettes/smm.rda")

## ----Smmlogitload, echo = FALSE, eval = TRUE, message=FALSE-------------------
load(url('https://raw.githubusercontent.com/ahoundetoungan/PartialNetwork/master/datavignettes/smmlogit.rda', "rb"))

## ----smmp1, echo = TRUE, eval = FALSE-----------------------------------------
# library(PartialNetwork)
# rm(list = ls())
# set.seed(123)
# # Number of groups
# M        <- 100
# # size of each group
# N        <- rep(30,M)
# # covariates
# X        <- cbind(rnorm(sum(N),0,5),rpois(sum(N),7))
# # network formation model parameter
# rho      <- c(-0.8, 0.2, -0.1)
# # individual effects
# beta     <- c(2, 1, 1.5, 5, -3)
# # endogenous effects
# alpha    <- 0.4
# # std-dev errors
# se       <- 1
# # network
# tmp      <- c(0, cumsum(N))
# X1l      <- lapply(1:M, function(x) X[c(tmp[x] + 1):tmp[x+1],1])
# X2l      <- lapply(1:M, function(x) X[c(tmp[x] + 1):tmp[x+1],2])
# dist.net <- function(x, y) abs(x - y)
# X1.mat   <- lapply(1:M, function(m) {
#   matrix(kronecker(X1l[[m]], X1l[[m]], FUN = dist.net), N[m])})
# X2.mat   <- lapply(1:M, function(m) {
#   matrix(kronecker(X2l[[m]], X2l[[m]], FUN = dist.net), N[m])})
# Xnet     <- as.matrix(cbind("Const" = 1,
#                             "dX1"   = mat.to.vec(X1.mat),
#                             "dX2"   = mat.to.vec(X2.mat)))
# ynet     <- Xnet %*% rho
# ynet     <- c(1*((ynet + rlogis(length(ynet))) > 0))
# G0       <- vec.to.mat(ynet, N, normalise = FALSE)
# # normalise
# G0norm   <- norm.network(G0)
# # Matrix GX
# GX       <- peer.avg(G0norm, X)
# # simulate dependent variable use an external package
# y        <- CDatanet::simsar(~ X + GX, Glist = G0norm, theta = c(alpha, beta, se))
# Gy       <- y$Gy
# y        <- y$y
# # build dataset
# dataset           <- as.data.frame(cbind(y, X, Gy, GX))
# colnames(dataset) <- c("y","X1","X2", "Gy", "GX1", "GX2")

## ----smmp2, echo = TRUE, eval = FALSE-----------------------------------------
# nNet      <- nrow(Xnet) # network formation model sample size
# Aobs      <- sample(1:nNet, round(0.3*nNet)) # We observed 30%
# # We can estimate rho using the gml function from the stats package
# logestim  <- glm(ynet[Aobs] ~ -1 + Xnet[Aobs,], family = binomial(link = "logit"))
# slogestim <- summary(logestim)
# rho.est   <- logestim$coefficients
# rho.var   <- slogestim$cov.unscaled # we also need the covariance of the estimator

## ----smmp3, echo = TRUE, eval = FALSE-----------------------------------------
# d.logit     <- lapply(1:M, function(x) {
#   out       <- 1/(1 + exp(-rho.est[1] - rho.est[2]*X1.mat[[x]] -
#                             rho.est[3]*X2.mat[[x]]))
#   diag(out) <- 0
#   out})

## ----Smmp4, echo = TRUE, eval = FALSE, message=FALSE--------------------------
# smm.logit   <- smmSAR(y ~ X1 + X2, dnetwork = d.logit, contextual = T,
#                       smm.ctr  = list(R = 100, print = F), data = dataset)
# summary(smm.logit)

## ----Smmp4b, echo = FALSE, eval = TRUE, message=FALSE-------------------------
print(smmlo$smm1)

## ----Smmp5, echo = TRUE, eval = FALSE, message=FALSE--------------------------
# fdist        <- function(rho.est, rho.var, M, X1.mat, X2.mat){
#   rho.est1   <- MASS::mvrnorm(mu = rho.est, Sigma = rho.var)
#   lapply(1:M, function(x) {
#   out        <- 1/(1 + exp(-rho.est1[1] - rho.est1[2]*X1.mat[[x]] -
#                              rho.est1[3]*X2.mat[[x]]))
#   diag(out)  <- 0
#   out})
# }

## ----Smmp6a, echo = TRUE, eval = FALSE, message=FALSE-------------------------
# fdist_args  <- list(rho.est = rho.est, rho.var = rho.var, M = M, X1.mat = X1.mat,
#                     X2.mat = X2.mat)
# summary(smm.logit, dnetwork = d.logit, data = dataset, .fun = fdist, .args = fdist_args,
#         sim = 500, ncores = 8) # ncores performs simulations in parallel

## ----Smmp6c, echo = FALSE, eval = TRUE, message=FALSE-------------------------
print(smmlo$smm2)

## ----smmlogitsave, echo = FALSE, eval = FALSE---------------------------------
# smm2  <- summary(smm.logit, dnetwork = d.logit, data = dataset, .fun = fdist, .args = fdist_args, sim = 500, ncores = 8)
# smmlo <- list(smm1 = smm.logit, smm2 = smm2)
# save(smmlo, file = "~/Dropbox/Papers - In progress/Partial Network/Package/AH/PartialNetwork/datavignettes/smmlogit.rda")

## ----BayesNone0, echo=FALSE---------------------------------------------------
load(url('https://raw.githubusercontent.com/ahoundetoungan/PartialNetwork/master/datavignettes/out.none.rda', "rb"))
#save(out.none1, out.none2.1, out.none2.2, out.none3.1, out.none3.2, file = "~/Dropbox/Papers - In progress/Partial Network/Package/AH/PartialNetwork/datavignettes/out.none.rda")

## ----BayesNone1as, eval=FALSE-------------------------------------------------
# rm(list = ls())
# library(PartialNetwork)
# set.seed(123)
# # EXAMPLE I: WITHOUT NETWORK FORMATION MODEL
# # Number of groups
# M             <- 50
# # size of each group
# N             <- rep(30,M)
# # individual effects
# beta          <- c(2,1,1.5)
# # contextual effects
# gamma         <- c(5,-3)
# # endogenous effects
# alpha         <- 0.4
# # std-dev errors
# se            <- 1
# # network distribution
# distr         <- runif(sum(N*(N-1)))
# distr         <- vec.to.mat(distr, N, normalise = FALSE)
# # covariates
# X             <- cbind(rnorm(sum(N),0,5),rpois(sum(N),7))
# # true network
# G0            <- sim.network(distr)
# # normalize
# G0norm        <- norm.network(G0)
# # GX
# GX            <- peer.avg(G0norm, X)
# # simulate dependent variable use an external package
# y             <- CDatanet::simsar(~ X + GX, Glist = G0norm,
#                                   theta = c(alpha, beta, gamma, se))
# y             <- y$y
# # dataset
# dataset       <- as.data.frame(cbind(y, X1 = X[,1], X2 = X[,2]))

## ----BayesNone1ae, eval=FALSE-------------------------------------------------
# # Example I-1: When the network is fully observed
# out.none1     <- mcmcSAR(formula = y ~ X1 + X2, contextual = TRUE, G0.obs = "all",
#                          G0 = G0, data = dataset, iteration = 2e4)
# summary(out.none1)

## ----BayesNone1b, echo=FALSE--------------------------------------------------
summary(out.none1)

## ----BayesNone21s, eval=FALSE-------------------------------------------------
# # Example I-2: When a part of the network is observed
# # 60% of the network data is observed
# G0.obs       <- lapply(N, function(x) matrix(rbinom(x^2, 1, 0.6), x))

## ----BayesNone21e, eval=FALSE-------------------------------------------------
# # replace the non-observed part of the network by 0 (missing links)
# G0.start     <- lapply(1:M, function(x) G0[[x]]*G0.obs[[x]])
# # Use network with missing data as the true network
# out.none2.1  <- mcmcSAR(formula = y ~ X1 + X2, contextual = TRUE, G0.obs = "all",
#                         G0 = G0.start,   data = dataset, iteration = 2e4)
# summary(out.none2.1) # the peer effets seem overestimated

## ----BayesNone21b, echo=FALSE-------------------------------------------------
summary(out.none2.1) # the peer effets seem overestimated

## ----BayesNone22a, eval=FALSE-------------------------------------------------
# out.none2.2  <- mcmcSAR(formula = y ~ X1 + X2, contextual = TRUE, G0.obs = G0.obs,
#                         G0 = G0.start, data = dataset,
#                         mlinks = list(dnetwork = distr), iteration = 2e4)
# summary(out.none2.2)

## ----BayesNone22b, echo=FALSE-------------------------------------------------
summary(out.none2.2)

## ----BayesNone31as, eval=FALSE------------------------------------------------
# # Example I-3: When only the network distribution is available
# # Simulate a fictitious network and use as true network
# G0.tmp       <- sim.network(distr)
# out.none3.1  <- mcmcSAR(formula = y ~ X1 + X2, contextual = TRUE, G0.obs = "all",
#                         G0 = G0.tmp, data = dataset, iteration = 2e4)
# summary(out.none3.1)  # the peer effets seem overestimated

## ----BayesNone31b, echo=FALSE-------------------------------------------------
summary(out.none3.1)  # the peer effets seem overestimated

## ----BayesNone32a, eval=FALSE-------------------------------------------------
# out.none3.2  <- mcmcSAR(formula = y ~ X1 + X2, contextual = TRUE, G0.obs = "none",
#                         data = dataset, mlinks = list(dnetwork = distr), iteration = 2e4)
# summary(out.none3.2)

## ----BayesNone32b, echo=FALSE-------------------------------------------------
summary(out.none3.2)  

## ----BayesLog0, echo=FALSE----------------------------------------------------
load(url('https://raw.githubusercontent.com/ahoundetoungan/PartialNetwork/master/datavignettes/out.logi.rda', "rb"))
# save(out.logi2.2, out.logi3.2, slogestim, sout.logi4.1, sout.logi4.2, sout.logi4.3, sout.selb1, sout.selb2, file = "~/Dropbox/Papers - In progress/Partial Network/Package/AH/PartialNetwork/datavignettes/out.logi.rda")

## ----BayesLog1as, eval=FALSE--------------------------------------------------
# # EXAMPLE II: NETWORK FORMATION MODEL: LOGIT
# rm(list = ls())
# library(PartialNetwork)
# set.seed(123)
# # Number of groups
# M             <- 50
# # size of each group
# N             <- rep(30,M)
# # individual effects
# beta          <- c(2,1,1.5)
# # contextual effects
# gamma         <- c(5,-3)
# # endogenous effects
# alpha         <- 0.4
# # std-dev errors
# se            <- 2
# # parameters of the network formation model
# rho           <- c(-2, -.5, .2)
# # covariates
# X             <- cbind(rnorm(sum(N),0,5),rpois(sum(N),7))
# # compute distance between individuals
# tmp           <- c(0, cumsum(N))
# X1l           <- lapply(1:M, function(x) X[c(tmp[x] + 1):tmp[x+1],1])
# X2l           <- lapply(1:M, function(x) X[c(tmp[x] + 1):tmp[x+1],2])
# dist.net      <- function(x, y) abs(x - y)
# X1.mat        <- lapply(1:M, function(m) {
#   matrix(kronecker(X1l[[m]], X1l[[m]], FUN = dist.net), N[m])})
# X2.mat        <- lapply(1:M, function(m) {
#   matrix(kronecker(X2l[[m]], X2l[[m]], FUN = dist.net), N[m])})
# # true network
# Xnet          <- as.matrix(cbind("Const" = 1,
#                                  "dX1"   = mat.to.vec(X1.mat),
#                                  "dX2"   = mat.to.vec(X2.mat)))
# ynet          <- Xnet %*% rho
# ynet          <- 1*((ynet + rlogis(length(ynet))) > 0)
# G0            <- vec.to.mat(ynet, N, normalise = FALSE)
# G0norm        <- norm.network(G0)
# # GX
# GX            <- peer.avg(G0norm, X)
# # simulate dependent variable use an external package
# y             <- CDatanet::simsar(~ X + GX, Glist = G0norm,
#                                      theta = c(alpha, beta, gamma, se))
# y             <- y$y
# # dataset
# dataset       <- as.data.frame(cbind(y, X1 = X[,1], X2 = X[,2]))

## ----BayesLog21as, eval=FALSE-------------------------------------------------
# # Example II-1: When a part of the network is observed
# # 60% of the network data is observed
# G0.obs       <- lapply(N, function(x) matrix(rbinom(x^2, 1, 0.6), x))
# # replace the non-observed part of the network by 0
# G0.start     <- lapply(1:M, function(x) G0[[x]]*G0.obs[[x]])
# # Infer the missing links in the network data
# mlinks       <- list(model = "logit", mlinks.formula = ~ dX1 + dX2,
#                      mlinks.data = as.data.frame(Xnet))

## ----BayesLog21ae, eval=FALSE-------------------------------------------------
# out.logi2.2  <- mcmcSAR(formula = y ~ X1 + X2, contextual = TRUE, G0.obs = G0.obs,
#                         G0 = G0.start, data = dataset, mlinks = mlinks,
#                         iteration = 2e4)
# summary(out.logi2.2)

## ----BayesLog22b, echo=FALSE--------------------------------------------------
summary(out.logi2.2) 

## ----BayesLog22c, fig.height = 3, fig.align = "center"------------------------
plot(out.logi2.2, plot.type = "sim", mar = c(3, 2.1, 1, 1))

## ----BayesLog22cc, fig.height = 3, fig.align = "center"-----------------------
plot(out.logi2.2, plot.type = "sim", which.parms = "rho", mar = c(3, 2.1, 1, 1))

## ----BayesLog31as, eval=FALSE-------------------------------------------------
# # Example II-2: When only the network distribution is available
# # Infer the network data
# # We only provide estimate of rho and its variance
# Gvec         <- mat.to.vec(G0, ceiled = TRUE)
# logestim     <- glm(Gvec ~ -1 + Xnet, family = binomial(link = "logit"))
# slogestim    <- summary(logestim)
# estimates    <- list("rho"     = logestim$coefficients,
#                      "var.rho" = slogestim$cov.unscaled,
#                      "N"       = N)
# mlinks       <- list(model = "logit", mlinks.formula = ~ dX1 + dX2,
#                      mlinks.data = as.data.frame(Xnet), estimates = estimates)
# out.logi3.2  <- mcmcSAR(formula = y ~ X1 + X2, contextual = TRUE, G0.obs = "none",
#                         data = dataset, mlinks = mlinks, iteration = 2e4)
# summary(out.logi3.2)

## ----BayesLog32b, echo=FALSE--------------------------------------------------
summary(out.logi3.2) 

## ----BayesLog32c, fig.height = 3, fig.align = "center"------------------------
plot(out.logi3.2, plot.type = "sim", mar = c(3, 2.1, 1, 1))

## ----BayesLog32cc, fig.height = 3, fig.align = "center"-----------------------
plot(out.logi3.2, plot.type = "sim", which.parms = "rho", mar = c(3, 2.1, 1, 1))

## ----BayesLog21ae3, eval=FALSE------------------------------------------------
# estimates    <- list("rho"     = logestim$coefficients,
#                      "var.rho" = slogestim$cov.unscaled)
# mlinks       <- list(model = "logit", mlinks.formula = ~ dX1 + dX2,
#                      mlinks.data = as.data.frame(Xnet), estimates = estimates)
# out.logi4.1  <- mcmcSAR(formula = y ~ X1 + X2, contextual = TRUE, G0.obs = G0.obs,
#                         G0 = G0.start, data = dataset, mlinks = mlinks,
#                         iteration = 2e4)
# summary(out.logi4.1)

## ----BayesLog21ae4, eval=TRUE, echo=FALSE-------------------------------------
print(sout.logi4.1)

## ----BayesLog21ae1, eval=TRUE-------------------------------------------------
print(slogestim)

## ----BayesLog21ae5, eval=FALSE------------------------------------------------
# prior        <- list("rho"     = c(0, 0, 0),
#                      "var.rho" = diag(3))
# mlinks       <- list(model = "logit", mlinks.formula = ~ dX1 + dX2,
#                      mlinks.data = as.data.frame(Xnet), prior = prior)
# out.logi4.2  <- mcmcSAR(formula = y ~ X1 + X2, contextual = TRUE, G0.obs = G0.obs,
#                         G0 = G0.start, data = dataset, mlinks = mlinks,
#                         iteration = 2e4)
# summary(out.logi4.2)

## ----BayesLog21ae6, echo=FALSE------------------------------------------------
print(sout.logi4.2)

## ----BayesLog21ae8, echo=FALSE------------------------------------------------
print(sout.logi4.3)

## ----Bayesard0, echo=FALSE----------------------------------------------------
load(url('https://raw.githubusercontent.com/ahoundetoungan/PartialNetwork/master/datavignettes/out.ard.rda', "rb"))
# save(data.plot1, data.plot2, genz, genv, zi, vk, file = "~/Dropbox/Papers - In progress/Partial Network/Package/AH/PartialNetwork/datavignettes/out.ard.rda")

## ----ard1, eval=FALSE---------------------------------------------------------
# rm(list = ls())
# library(PartialNetwork)
# set.seed(123)
# # LATENT SPACE MODEL
# N           <- 500
# genzeta     <- 1
# mu          <- -1.35
# sigma       <- 0.37
# K           <- 12    # number of traits
# P           <- 3     # Sphere dimension
# # ARD parameters
# # Generate z (spherical coordinates)
# genz        <- rvMF(N, rep(0,P))
# # Generate nu  from a Normal(mu, sigma^2) (The gregariousness)
# gennu       <- rnorm(N, mu, sigma)
# # compute degrees
# gend        <- N*exp(gennu)*exp(mu+0.5*sigma^2)*exp(logCpvMF(P,0) - logCpvMF(P,genzeta))
# # Link probabilities
# distr       <- sim.dnetwork(gennu, gend, genzeta, genz)
# # Adjacency matrix
# G           <- sim.network(distr)
# # Generate vk, the trait location
# genv        <- rvMF(K, rep(0, P))
# # set fixed some vk  distant
# genv[1,]    <- c(1, 0, 0)
# genv[2,]    <- c(0, 1, 0)
# genv[3,]    <- c(0, 0, 1)
# # eta, the intensity parameter
# geneta      <- abs(rnorm(K, 2, 1))
# # Build traits matrix
# densityatz  <- matrix(0, N, K)
# for(k in 1:K){
#   densityatz[,k] <- dvMF(genz, genv[k,]*geneta[k])
# }
# trait       <- matrix(0, N, K)
# NK          <- floor(runif(K, 0.8, 0.95)*colSums(densityatz)/apply(densityatz, 2, max))
# for (k in 1:K) {
#   trait[,k] <- rbinom(N, 1, NK[k]*densityatz[,k]/sum(densityatz[,k]))
# }
# # Build ADR
# ARD         <- G %*% trait
# # generate b
# genb        <- numeric(K)
# for(k in 1:K){
#   genb[k]   <- sum(G[,trait[,k]==1])/sum(G)
# }

## ----ard1as, eval=FALSE-------------------------------------------------------
# # Example1: ARD is observed for the whole population
# # initialization
# d0     <- exp(rnorm(N)); b0 <- exp(rnorm(K)); eta0 <- rep(1,K)
# zeta0  <- 2; z0 <- matrix(rvMF(N, rep(0,P)), N); v0 <- matrix(rvMF(K,rep(0, P)), K)
# # We should fix some vk and bk
# vfixcolumn      <- 1:5
# bfixcolumn      <- c(3, 7, 9)
# b0[bfixcolumn]  <- genb[bfixcolumn]
# v0[vfixcolumn,] <- genv[vfixcolumn,]
# start           <- list("z" = z0, "v" = v0, "d" = d0, "b" = b0, "eta" = eta0,
#                         "zeta" = zeta0)

## ----ard1e, eval=FALSE--------------------------------------------------------
# # MCMC
# estim.ard1      <- mcmcARD(Y = ARD, traitARD = trait, start = start, fixv = vfixcolumn,
#                            consb = bfixcolumn, iteration = 5000)

## ----ardaa1, eval=FALSE-------------------------------------------------------
# # plot coordinates of individual 123
# i     <- 123
# zi    <- estim.ard1$simulations$z[i,,]
# par(mfrow = c(3, 1), mar = c(2.1, 2.1, 1, 1))
# invisible(lapply(1:3, function(x) {
#   plot(zi[x,], type = "l", ylab = "", col = "blue", ylim = c(-1, 1))
#   abline(h = genz[i, x], col = "red")
# }))

## ----ardaa1a, echo=FALSE, fig.height = 4, fig.align = "center"----------------
i     <- 123
par(mfrow = c(3, 1), mar = c(2.1, 2.1, 1, 1))
invisible(lapply(1:3, function(x) {
  plot(zi[x,], type = "l", ylab = "", col = "blue", ylim = c(-1, 1))
  abline(h = genz[i, x], col = "red")
}))

## ----ardaa2, eval=FALSE-------------------------------------------------------
# # plot coordinates of the trait 8
# k     <- 8
# vk    <- estim.ard1$simulations$v[k,,]
# par(mfrow = c(3, 1), mar = c(2.1, 2.1, 1, 1))
# invisible(lapply(1:3, function(x) {
#   plot(vk[x,], type = "l", ylab = "", col = "blue", ylim = c(-1, 1))
#   abline(h = genv[k, x], col = "red")
# }))

## ----ardaa2a, echo=FALSE, fig.height = 4, fig.align = "center"----------------
k     <- 8
par(mfrow = c(3, 1), mar = c(2.1, 2.1, 1, 1))
invisible(lapply(1:3, function(x) {
  plot(vk[x,], type = "l", ylab = "", col = "blue", ylim = c(-1, 1))
  abline(h = genv[k, x], col = "red")
}))

## ----ardb, message=FALSE, eval=FALSE------------------------------------------
# # plot degree
# library(ggplot2)
# data.plot1 <- data.frame(True_degree  = gend,
#                          Estim_degree = colMeans(tail(estim.ard1$simulations$d, 2500)))
# ggplot(data = data.plot1, aes(x = True_degree, y = Estim_degree)) +
#    geom_abline(col = "red") + geom_point(col = "blue")

## ----ardbb, message=FALSE, echo=FALSE, fig.height = 3, fig.align = "center"----
library(ggplot2)
ggplot(data = data.plot1, aes(x = True_degree, y = Estim_degree)) + 
   geom_abline(col = "red") + geom_point(col = "blue")

## ----ard2as, eval=FALSE-------------------------------------------------------
# # Example2: ARD is observed for 70% population
# # sample with ARD
# n          <- round(0.7*N)
# # individual with ARD
# iselect    <- sort(sample(1:N, n, replace = FALSE))
# ARDs       <- ARD[iselect,]
# traits     <- trait[iselect,]
# # initialization
# d0         <- d0[iselect]; z0 <- z0[iselect,]
# start      <- list("z" = z0, "v" = v0, "d" = d0, "b" = b0, "eta" = eta0, "zeta" = zeta0)

## ----ard2ae, eval=FALSE-------------------------------------------------------
# # MCMC
# estim.ard2 <- mcmcARD(Y = ARDs, traitARD = traits, start = start, fixv = vfixcolumn,
#                       consb = bfixcolumn, iteration = 5000)

## ----ard2b, eval=FALSE--------------------------------------------------------
# # estimation for non ARD
# # we need a logical vector indicating if the i-th element has ARD
# hasARD      <- (1:N) %in% iselect
# # we use the matrix of traits to estimate distance between individuals
# estim.nard2 <- fit.dnetwork(estim.ard2, X = trait, obsARD = hasARD, m = 1)

## ----arde, eval=FALSE---------------------------------------------------------
# # estimated degree
# estd          <- estim.nard2$degree
# data.plot2    <- data.frame(True_degree  = gend,
#                             Estim_degree = estd,
#                             Has_ARD      = ifelse(hasARD, "YES", "NO"))
# ggplot(data = data.plot2, aes(x = True_degree, y = Estim_degree, colour = Has_ARD)) +
#   geom_abline(col = "red") + geom_point()

## ----ardee, message=FALSE, echo=FALSE, fig.height = 3, fig.align = "center"----
ggplot(data = data.plot2, aes(x = True_degree, y = Estim_degree, colour = Has_ARD)) + 
  geom_abline(col = "red") + geom_point() 

## ----Bayesard00, echo=FALSE---------------------------------------------------
load(url('https://raw.githubusercontent.com/ahoundetoungan/PartialNetwork/master/datavignettes/out.lsp.rda', "rb"))
#save(out.lspa1, out.lspa2, file = "~/Dropbox/Papers - In progress/Partial Network/Package/AH/PartialNetwork/datavignettes/out.lsp.rda")

## ----bayesard1, eval=FALSE----------------------------------------------------
# rm(list = ls())
# library(PartialNetwork)
# set.seed(123)
# M             <- 30
# N             <- rep(60, M)
# genzeta       <- 3
# mu            <- -1.35
# sigma         <- 0.37
# K             <- 12    # number of traits
# P             <- 3     # Sphere dimension
# 
# # IN THIS LOOP, WE GENERATE DATA FOLLOWING BREZA ET AL. (2020) AND
# # ESTIMATE THEIR LATENT SPACE MODEL FOR EACH SUB-NETWORK.
# estimates     <- list()
# list.trait    <- list()
# G0            <- list()
# for (m in 1:M) {
#   #######################################################################################
#   #######                             SIMULATION STAGE                           ########
#   #######################################################################################
#   # ARD parameters
#   # Generate z (spherical coordinates)
#   genz    <- rvMF(N[m], rep(0,P))
#   # Generate nu  from a Normal(mu, sigma^2) (The gregariousness)
#   gennu   <- rnorm(N[m],mu,sigma)
#   # compute degrees
#   gend    <- N[m]*exp(gennu)*exp(mu+0.5*sigma^2)*exp(logCpvMF(P,0) - logCpvMF(P,genzeta))
#   # Link probabilities
#   distr   <- sim.dnetwork(gennu, gend, genzeta, genz)
#   # Adjacency matrix
#   G        <- sim.network(distr)
#   G0[[m]]  <- G
#   # Generate vk, the trait location
#   genv     <- rvMF(K, rep(0, P))
#   # set fixed some vk  distant
#   genv[1,] <- c(1, 0, 0)
#   genv[2,] <- c(0, 1, 0)
#   genv[3,] <- c(0, 0, 1)
#   # eta, the intensity parameter
#   geneta   <-abs(rnorm(K, 2, 1))
#   # Build traits matrix
#   densityatz       <- matrix(0, N[m], K)
#   for(k in 1:K){
#     densityatz[,k] <- dvMF(genz, genv[k,]*geneta[k])
#   }
#   trait       <- matrix(0, N[m], K)
#   NK          <- floor(runif(K, .8, .95)*colSums(densityatz)/apply(densityatz, 2, max))
#   for (k in 1:K) {
#     trait[,k] <- rbinom(N[m], 1, NK[k]*densityatz[,k]/sum(densityatz[,k]))
#   }
#   list.trait[[m]]  <- trait
#   # Build ADR
#   ARD         <- G %*% trait
#   # generate b
#   genb        <- numeric(K)
#   for(k in 1:K){
#     genb[k]  <- sum(G[,trait[,k]==1])/sum(G) + 1e-8
#   }
# 
#   #######################################################################################
#   #######                             ESTIMATION STAGE                           ########
#   #######################################################################################
#   # initialization
#   d0     <- gend; b0 <- exp(rnorm(K)); eta0 <- rep(1,K); zeta0 <- genzeta
#   z0     <- matrix(rvMF(N[m], rep(0,P)), N[m]); v0 <- matrix(rvMF(K,rep(0, P)), K)
#   # We should fix some vk and bk
#   vfixcolumn      <- 1:5
#   bfixcolumn      <- c(1, 3, 5, 7, 9, 11)
#   b0[bfixcolumn]  <- genb[bfixcolumn]
#   v0[vfixcolumn,] <- genv[vfixcolumn,]
#   start           <- list("z" = z0, "v" = v0, "d" = d0, "b" = b0, "eta" = eta0,
#                           "zeta" = zeta0)
#   estimates[[m]]  <- mcmcARD(Y = ARD, traitARD = trait, start = start, fixv = vfixcolumn,
#                              consb = bfixcolumn, sim.d = FALSE, sim.zeta = FALSE,
#                              iteration = 5000, ctrl.mcmc = list(print = FALSE))
# }
# 
# # SIMULATE DATA FOR THE OUTCOME MODEL
# # individual effects
# beta          <- c(2,1,1.5)
# # contextual effects
# gamma         <- c(5,-3)
# # endogenous effects
# alpha         <- 0.4
# # std-dev errors
# se            <- 1
# # covariates
# X             <- cbind(rnorm(sum(N),0,5),rpois(sum(N),7))
# # Normalise G0
# G0norm        <- norm.network(G0)
# # GX
# GX            <- peer.avg(G0norm, X)
# # simulate dependent variable use an external package
# y             <- CDatanet::simsar(~ X + GX, Glist = G0norm,
#                                   theta = c(alpha, beta, gamma, se))
# y             <- y$y
# # dataset
# dataset       <- as.data.frame(cbind(y, X1 = X[,1], X2 = X[,2]))

## ----bayesard2s, eval=FALSE---------------------------------------------------
# mlinks       <- list(model = "latent space", estimates = estimates)
# out.lspa1    <- mcmcSAR(formula = y ~ X1 + X2, contextual = TRUE, G0.obs = "none",
#                         data = dataset, mlinks = mlinks, iteration = 2e4)
# summary(out.lspa1)

## ----Bayesard2b, echo=FALSE---------------------------------------------------
summary(out.lspa1) 

## ----Bayesard2c, fig.height = 3, fig.align = "center"-------------------------
plot(out.lspa1, plot.type = "sim", mar = c(3, 2.1, 1, 1))

## ----bayesard3, eval=FALSE----------------------------------------------------
# rm(list = ls())
# library(PartialNetwork)
# set.seed(123)
# M             <- 30
# N             <- rep(60, M)
# genzeta       <- 3
# mu            <- -1.35
# sigma         <- 0.37
# K             <- 12
# P             <- 3
# 
# # IN THIS LOOP, WE GENERATE DATA FOLLOWING BREZA ET AL. (2020) AND
# # ESTIMATE THEIR LATENT SPACE MODEL FOR EACH SUB-NETWORK.
# estimates     <- list()
# list.trait    <- list()
# obARD         <- list()
# G0            <- list()
# for (m in 1:M) {
#   #######################################################################################
#   #######                             SIMULATION STAGE                           ########
#   #######################################################################################
#   # ARD parameters
#   # Generate z (spherical coordinates)
#   genz    <- rvMF(N[m], rep(0,P))
#   # Generate nu  from a Normal(mu, sigma^2) (The gregariousness)
#   gennu   <- rnorm(N[m],mu,sigma)
#   # compute degrees
#   gend    <- N[m]*exp(gennu)*exp(mu+0.5*sigma^2)*exp(logCpvMF(P,0) - logCpvMF(P,genzeta))
#   # Link probabilities
#   distr   <- sim.dnetwork(gennu, gend, genzeta, genz)
#   # Adjacency matrix
#   G        <- sim.network(distr)
#   G0[[m]]  <- G
#   # Generate vk, the trait location
#   genv     <- rvMF(K, rep(0, P))
#   # set fixed some vk  distant
#   genv[1,] <- c(1, 0, 0)
#   genv[2,] <- c(0, 1, 0)
#   genv[3,] <- c(0, 0, 1)
#   # eta, the intensity parameter
#   geneta   <-abs(rnorm(K, 2, 1))
#   # Build traits matrix
#   densityatz  <- matrix(0, N[m], K)
#   for(k in 1:K){
#     densityatz[,k] <- dvMF(genz, genv[k,]*geneta[k])
#   }
#   trait       <- matrix(0, N[m], K)
#   NK          <- floor(runif(K, .8, .95)*colSums(densityatz)/apply(densityatz, 2, max))
#   for (k in 1:K) {
#     trait[,k] <- rbinom(N[m], 1, NK[k]*densityatz[,k]/sum(densityatz[,k]))
#   }
#   list.trait[[m]]  <- trait
#   # Build ADR
#   ARD         <- G %*% trait
#   # generate b
#   genb        <- numeric(K)
#   for(k in 1:K){
#     genb[k]  <- sum(G[,trait[,k]==1])/sum(G) + 1e-8
#   }
#   # sample with ARD
#   n          <- round(runif(1, .7, 1)*N[m])
#   # individual with ARD
#   iselect    <- sort(sample(1:N[m], n, replace = FALSE))
#   hasARD     <- (1:N[m]) %in% iselect
#   obARD[[m]] <- hasARD
#   ARDs       <- ARD[iselect,]
#   traits     <- trait[iselect,]
#   #######################################################################################
#   #######                             ESTIMATION STAGE                           ########
#   #######################################################################################
#   # initialization
#   d0         <- gend[iselect]; b0 <- exp(rnorm(K)); eta0 <- rep(1,K); zeta0 <- genzeta
#   z0         <- matrix(rvMF(n, rep(0,P)), n); v0 <- matrix(rvMF(K, rep(0, P)), K)
#   # We should fix some vk and bk
#   vfixcolumn     <- 1:5
#   bfixcolumn     <- c(1, 3, 5, 7, 9, 11)
#   b0[bfixcolumn] <- genb[bfixcolumn]; v0[vfixcolumn,] <- genv[vfixcolumn,]
#   start          <- list("z" = z0, "v" = v0, "d" = d0, "b" = b0, "eta" = eta0,
#                           "zeta" = zeta0)
#   estimates[[m]] <- mcmcARD(Y = ARDs, traitARD = traits, start = start, fixv = vfixcolumn,
#                             consb = bfixcolumn,  sim.d = FALSE, sim.zeta = FALSE,
#                             iteration = 5000, ctrl.mcmc = list(print = FALSE))
# }
# 
# # SIMULATE DATA FOR THE OUTCOME MODEL
# # individual effects
# beta          <- c(2,1,1.5)
# # contextual effects
# gamma         <- c(5,-3)
# # endogenous effects
# alpha         <- 0.4
# # std-dev errors
# se            <- 1
# # covariates
# X             <- cbind(rnorm(sum(N),0,5),rpois(sum(N),7))
# # Normalise G0
# G0norm        <- norm.network(G0)
# # GX
# GX            <- peer.avg(G0norm, X)
# # simulate dependent variable use an external package
# y             <- CDatanet::simsar(~ X + GX, Glist = G0norm,
#                                   theta = c(alpha, beta, gamma, se))
# y             <- y$y
# # dataset
# dataset       <- as.data.frame(cbind(y, X1 = X[,1], X2 = X[,2]))

## ----bayesard4s, eval=FALSE---------------------------------------------------
# mlinks       <- list(model = "latent space", estimates = estimates,
#                      mlinks.data = list.trait, obsARD = obARD)
# out.lspa2    <- mcmcSAR(formula = y ~ X1 + X2, contextual = TRUE, G0.obs = "none",
#                         data = dataset, mlinks = mlinks, iteration = 2e4)
# summary(out.lspa2)

## ----Bayesard4b, echo=FALSE---------------------------------------------------
summary(out.lspa2)

## ----Bayesard4c, fig.height = 3, fig.align = "center"-------------------------
plot(out.lspa2, plot.type = "sim", mar = c(3, 2.1, 1, 1))

## ----selb1, echo = TRUE, eval=FALSE-------------------------------------------
# rm(list = ls())
# library(PartialNetwork)
# library(dplyr)
# set.seed(123)
# # Number of groups
# M             <- 50
# # size of each group
# N             <- rep(30,M)
# # individual effects
# beta          <- c(2, 1, 1.5)
# # contextual effects
# gamma         <- c(5,-3)
# # endogenous effects
# alpha         <- 0.4
# # std-dev errors
# se            <- 2
# # parameters of the network formation model
# rho           <- c(-0.5, -.5, .4)
# # covariates
# X             <- cbind(rnorm(sum(N),0,5),rpois(sum(N),7))
# # compute distance between individuals
# tmp           <- c(0, cumsum(N))
# X1l           <- lapply(1:M, function(x) X[c(tmp[x] + 1):tmp[x+1],1])
# X2l           <- lapply(1:M, function(x) X[c(tmp[x] + 1):tmp[x+1],2])
# dist.net      <- function(x, y) abs(x - y)
# X1.mat        <- lapply(1:M, function(m) {
#   matrix(kronecker(X1l[[m]], X1l[[m]], FUN = dist.net), N[m])})
# X2.mat        <- lapply(1:M, function(m) {
#   matrix(kronecker(X2l[[m]], X2l[[m]], FUN = dist.net), N[m])})
# # true network
# Xnet          <- as.matrix(cbind("Const" = 1,
#                                  "dX1"   = mat.to.vec(X1.mat),
#                                  "dX2"   = mat.to.vec(X2.mat)))
# ynet          <- Xnet %*% rho
# ynet          <- c(1*((ynet + rlogis(length(ynet))) > 0))
# G0            <- vec.to.mat(ynet, N, normalise = FALSE)
# # number of friends
# nfriends      <- unlist(lapply(G0, function(x) rowSums(x)))
# # number of missing links
# nmislink      <- sapply(nfriends, function(x) sample(0:x, 1))

## ----selb2, echo = TRUE, eval=FALSE-------------------------------------------
# Gobs          <- list(M)  # The observed network
# G0.obs        <- list(M)  # Which information is true and doubtful
# for(x in 1:M){
#   Gx          <- G0[[x]]
#   G0.obsx     <- matrix(1, N[x], N[x]); diag(G0.obsx) <- 0
#   csum        <- cumsum(c(0, N))
#   nmis        <- nmislink[(csum[x] + 1):csum[x + 1]]
#   for (i in 1:N[x]) {
#     if(nmis[i] > 0){
#       tmp     <- which(c(Gx[i,]) == 1)
#       if(length(which(c(Gx[i,]) == 1)) > 1) {
#         tmp   <- sample(which(c(Gx[i,]) == 1), nmis[i])
#       }
#       Gx[i,tmp]   <- 0
#       G0.obsx[i,] <- 0
#     }
#   }
#   Gobs[[x]]   <- Gx
#   G0.obs[[x]] <- G0.obsx
# }
# G0norm        <- norm.network(G0)
# # GX
# GX            <- peer.avg(G0norm, X)
# # simulate dependent variable use an external package
# y             <- CDatanet::simsar(~ X + GX, Glist = G0norm,
#                                   theta = c(alpha, beta, gamma, se))
# y             <- y$y
# # data set
# dataset       <- as.data.frame(cbind(y, X1 = X[,1], X2 = X[,2]))

## ----selb3a, eval=FALSE, echo=TRUE--------------------------------------------
# mlinks        <- list(model = "logit", mlinks.formula = ~ dX1 + dX2,
#                       mlinks.data = as.data.frame(Xnet))
# out.selb1     <- mcmcSAR(formula = y ~ X1 + X2, contextual = TRUE, G0 = Gobs,
#                          G0.obs = G0.obs, data = dataset, mlinks = mlinks,
#                          iteration = 2e4)
# summary(out.selb1)

## ----selb3b, echo=FALSE-------------------------------------------------------
print(sout.selb1)

## ----selb3c, eval=FALSE, echo=TRUE--------------------------------------------
# G0.obsvec     <- as.logical(mat.to.vec(G0.obs))
# Gvec          <- mat.to.vec(Gobs, ceiled = TRUE)[G0.obsvec]
# W             <- unlist(data.frame(nfriends = nfriends, nmislink = nmislink) %>%
#                        group_by(nfriends) %>%
#                        summarise(w = length(nmislink)/sum(nmislink == 0)) %>%
#                        select(w))
# W             <- lapply(1:M, function(x){
#   matrix(rep(W[rowSums(G0[[x]]) + 1], each = N[x]), N[x], byrow = TRUE)})
# weights       <- mat.to.vec(W)[G0.obsvec]
# mlinks        <- list(model = "logit", mlinks.formula = ~ dX1 + dX2,
#                       mlinks.data = as.data.frame(Xnet), weights = weights)
# out.selb2     <- mcmcSAR(formula = y ~ X1 + X2, contextual = TRUE, G0 = Gobs,
#                          G0.obs = G0.obs, data = dataset, mlinks = mlinks,
#                          iteration = 2e4)
# summary(out.selb2)

## ----selb3d, echo=FALSE-------------------------------------------------------
print(sout.selb2)

