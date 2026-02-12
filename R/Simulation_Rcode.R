#########################################################
#### Packages 
#########################################################
library("splines")
library(SpatialTools)
library(doParallel)
library(eva)
library(extRemes)
library(MASS)
library(cgwtools)
library(mvtnorm)
library(fda)
library(rlang)
library(fields)
library(Iso)
library(qvalue)
library(sandwich)
library(lmtest)
library(geoR)
library(SpatialExtremes)
library(ismev)
library(mvtnorm)
library(cgwtools)
library(simcausal)
#########################################
### Import functions 
#########################################
source('/Functions/Nonparam_functions.R')
source('/Functions/GEVfunctions.R')
source('/Functions/Mirror_pi.R')
source('/Functions/Neigh_Function.R')

#### Import US map data #### 
load('/Data/NAM_Polygon_index.Rdata')

#### Set parameters #### 
n.obs <- 60 
n.site<- nrow(loc)
cov <- 'whitmat'
range <- 3
smooth <-0.2
t1 = 10
t2 = 20
t5 =50
t10 = 100
N = nrow(loc)

###### Create B-Spline basis function ########
sp1 <- create.bspline.basis(rangeval=c(-125,-66), nbasis=4, norder=3) ## What is coefficient?
sp2 <- create.bspline.basis(rangeval=c(25,50), nbasis=4, norder=3)

eval.sp1 <- eval.basis(loc[,1], sp1)
eval.sp2 <- eval.basis(loc[,2], sp2)

## create design matrix for least square esitmates of coeff ##
eval.sp <- matrix(NA, n.site, ncol(eval.sp1)*ncol(eval.sp2))
for (i in 1:n.site){
  eval.sp[i,] <- kronecker(eval.sp1[i,], eval.sp2[i,])
}

## Generate delta(s) / pi(s) 
beta_coeff = matrix(c(-34.2,20.4,4.7,-1.0,-0.2,-8.2
                      ,-0.6,2.6,14.4,-7.9,-3.8,-3.4,-40.3,2.4,10.1,-6.5), nrow=16)
pi_hat = plogis(as.numeric(eval.sp %*% beta_coeff))
delta.tr = sapply(pi_hat,function(x) rbern(1,x) )
table(delta.tr)/n.site

########################################################################
### Simulation 
### Generate Parameters : Fit Exp model for parameters
########################################################################
## load the actual precipitation of ERA5 and observation 
load('/Data/Winter_AnnualMax_Precip0.Rdata') 

#### Generate location parameter #### 
sigma2 = 155; phi = 12 # based on the 
covm.loc <- sigma2*exp(-distmat/phi)  ## compute covariance matrix
diag(covm.loc) =sigma2
Dchol = t(chol(covm.loc))
Dtmpx <- pmax(3, 12+Dchol %*%rnorm(N));
param.loc = Dtmpx
summary(param.loc)

#### Generate Scale Parameter #### 
sigma2 =30; phi = 10;
covm.sc <- sigma2*exp(-distmat/phi)  ## compute covariance matrix
diag(covm.sc) = sigma2
Dchol = t(chol(covm.sc))
Dtmpx <- pmax(3,7.5+Dchol %*%rnorm(N));
param.scale = Dtmpx
summary(param.scale)

#### Generate Shape Paramter #### 
sigma2 = 0.03; phi =30
covm.sh <- sigma2*exp(-distmat/phi)  ## compute covariance matrix
diag(covm.sh) = sigma2#+0.005
Dchol = t(chol(covm.sh))
Dtmpx <- pmax(-0.4999, 0.1+Dchol %*%rnorm(N));
Dtmpx = pmin(Dtmpx, 0.5); Dtmpx = pmax(Dtmpx, 0.01)
param.shape = Dtmpx
summary(param.shape)

### Bind into data and save
Param_SimX = data.frame(loc,scale = param.scale, loc = param.loc,
                        shape = param.shape, delta = delta.tr)

#save(Param_SimX, file="Simulated Data/Simulation_Parameters.Rdata")

############ Generate Signal (alternative parameters) ############
#### 1. Location signal ####
# basis function
load('Data/Simulation_Parameters.Rdata')
sp1 <- create.bspline.basis(rangeval=c(3, 53), nbasis=6, norder=6)
eval.loc <- eval.basis(param.loc0, sp1)

# Generate Signal
coeff.loc1 = c(0.96, 3.55, 1.55, 4.22, 4.61, 6.31)
k.sig = 1.5 # k.sig = 1.3, 1.5, 1.7
param.loc02 = Param_SimX$loc ## location parameter at alternative location
param.loc02 = Param_SimX$loc + (eval.loc %*% coeff.loc1)*k.sig*Param_SimX$delta  
Param_SimX = data.frame(loc,scale = param.scale, loc = param.loc,
                        shape = param.shape, delta = delta.tr, loc1 = param.loc02)

#save(Param_SimX, file="/Data/Simulated Data/SimData_Locdiff_K1.3_P35.Rdata")
#save(Param_SimX, file="/Data/Simulated Data/SimData_Locdiff_K1.5_P35.Rdata")
#save(Param_SimX, file="/Data/Simulated Data/SimData_Locdiff_K1.7_P35.Rdata")

#### 2. Scale signal ####
# basis function
sp1 <- create.bspline.basis(rangeval=c(3, 23.1), nbasis=6, norder=6)
eval.sc <- eval.basis(param.scale0, sp1)

# Generate Signal 
coeff.sc1 = c(0.53, 1.65, 1.85, 3.65, 4.58, 5.99)
k.sig = 2 # k.sig = 2, 2.5, 3
param.scale02 = Param_SimX$scale
param.scale02 = Param_SimX$scale + ((eval.sc %*% coeff.sc1))*k.sig*Param_SimX$delta
param.scale02 = pmax(1, param.scale02)
Param_SimX = data.frame(loc,scale = param.scale, loc = param.loc,
                        shape = param.shape, delta = delta.tr, scale1 = param.scale02)

#save(Param_SimX, file="/Data/Simulated Data/SimData_Scalediff_K2_P35.Rdata")
#save(Param_SimX, file="/Data/Simulated Data/SimData_Scalediff_K2.5_P35.Rdata")
#save(Param_SimX, file="/Data/Simulated Data/SimData_Scalediff_K3_P35.Rdata")

#### 3. Shape signal ####
# basis function
sp1 <- create.bspline.basis(rangeval=c(0.009, 0.421), nbasis=6, norder=6)
eval.sh <- eval.basis(param.shape0, sp1)

# Generate Signal 
coeff.sh = c(0.06, -0.15, 0.45, 0.15, 0.58, 0.56)
k.sig = 1.5 # k.sig = 0.8, 1, 1.2
param.shape02 = Param_SimX$shape
param.shape02 = Param_SimX$shape + (eval.sh %*% coeff.sh)*k.sig*Param_SimX$delta

## Set the alternative parameter
Param_SimX = data.frame(loc,scale = param.scale, loc = param.loc,
                        shape = param.shape, delta = delta.tr, shape1 = param.shape02)

#save(Param_SimX, file="/Data/Simulated Data/SimData_Shapediff_K1.2_P35_2.Rdata")# k.sig = 1.2
#save(Param_SimX, file="/Data/Simulated Data/SimData_Shapediff_K1_P35_2.Rdata")# k.sig = 1
#save(Param_SimX, file="/Data/Simulated Data/SimData_Shapediff_K0.8_P35_2.Rdata") # k.sig = 0.8
