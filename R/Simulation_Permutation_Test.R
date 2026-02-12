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
source('Functions/General_functions.R')
source('Functions/Nonparam_functions.R')
source('Functions/GEVfunctions.R')
source('Functions/Mirror_pi.R')
source('Functions/Neigh_Function.R')
source('Functions/Optimal_Basis_functions.R')
#############################################
# get args
#############################################

args = commandArgs(TRUE)
ii = as.integer(args[1])

print(ii)
seed = as.integer(8*ii)
print(class(seed))

cov <- 'whitmat'
r0 <- 3
s0 <- 0.2

######################################
### Simulate Max-Stable Data
#######################################
## load the actual precipitation of ERA5 and observation 
load('/Data/Winter_AnnualMax_Precip0.Rdata') 

loc = as.matrix(loc)
n.obs = ncol(annualmax_precip_obs)-2; 
N = n.site = nrow(loc)
cov <- 'whitmat'

## Generate max-stable data
library('SpatialExtremes')
set.seed(seed*1+1)
data1 <- rmaxstab(n.obs, loc, cov.mod = cov, nugget = 0,range =r0, smooth = s0)

set.seed(seed*2+2)
data2 <- rmaxstab(n.obs, loc, cov.mod = cov, nugget = 0,range =r0, smooth = s0)
########################################################################
### P-value Generation 
### Permutation test 
########################################################################
### Example : Moderate signal of scale parameter difference
### Load simulated data
load('/Data/Simulated Data/SimData_Scalediff_K2_P35.Rdata') 

delta.tr = Param_SimX$delta
param.loc.p0 = Param_SimX$loc; param.scale.p0 = Param_SimX$scale; param.shape.p0 = Param_SimX$shape
param.loc.p1 = Param_SimX$loc; param.scale.p1 = Param_SimX$scale1; param.shape.p1 = Param_SimX$shape

range <- 3
smooth <- 0.2
p=1
t0 = 5*p
t1 = 10*p
t2 = 20*p
t3 =50*p
t4 = 100*p

###############################################################################################
# Use two different maxstable model to estimate the param
###############################################################################################
### Generate 100 Max-stable data with r=3 / s =0.2 
# get args
args = commandArgs(TRUE)
ii = as.integer(args[1])

print(ii)
seed = as.integer(8*ii)
print(class(seed))

cov <- 'whitmat'
r0 <- 3
s0 <- 0.2
######################################
### Simulate Max-Stable Data
#######################################
loc = as.matrix(loc)
n.obs = ncol(annualmax_precip_obs)-2; N = n.site = nrow(loc)
cov <- 'whitmat'

### Generate max-stable data
library('SpatialExtremes')
set.seed(seed*1+1)
data1 <- rmaxstab(n.obs, loc, cov.mod = cov, nugget = 0,range =r0, smooth = s0)
s = Sys.time()
set.seed(seed*2+2)
data2 <- rmaxstab(n.obs, loc, cov.mod = cov, nugget = 0,range =r0, smooth = s0)
e = Sys.time(); e-s

### Frechet to GEV
dataxp = datayp = matrix(nrow=n.obs, ncol= n.site)
for (i in 1:n.site){
  dataxp[,i] <- frech2gev(data1[,i], param.loc.p0[i], param.scale.p0[i],param.shape.p0[i])
  datayp[,i] <- frech2gev(data2[,i], param.loc.p1[i], param.scale.p1[i],param.shape.p1[i])
  dataxp[,i] = pmin(dataxp[,i], 1000); dataxp[,i] = pmax(0.05, dataxp[,i])
  datayp[,i] = pmin(datayp[,i], 1000); datayp[,i] = pmax(0.05, datayp[,i])
}

#########--- FIT GEV at each location ----#############
rt.time = c(t0, t1, t2, t3, t4)
param_x = cbind(param.loc.p0, param.scale.p0, param.shape.p0);
param_y = cbind(param.loc.p1, param.scale.p1, param.shape.p1)

### Calculate the difference of parameters and return levels
GEV_Rt_Scalediff = Fit_GEV_RT(data1, data2, param_x, param_y,rt.time)
rt.diff.GEVP=GEV_Rt_Scalediff$Rt.diff
param.diff.GEVP=GEV_Rt_Scalediff$Param.diff

######################################################
# Set up parallel computing for permutation test 
######################################################
num_cores <- detectCores() - 1  # Use all but one core
registerDoParallel(cores = 12)

# Check if the parallel backend is properly registered
print(paste0("Running with ", getDoParWorkers(), " workers"))

comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}
##################################################################################################################
# Rt level and uncertainty calculation with aggregation
#set.seed(seed)
t0 = 5; t1 = 10; t2 = 20; t5= 50; t10 = 100
##  BOOT STRAP ##
s = Sys.time()

r <- foreach(ll=seq_len(1000),.combine='comb',.multicombine=TRUE,.init=list(list(), list())
             , .packages = c("extRemes", "eva","SpatialTools")) %dopar% {
               #source('~/scratch/GEVfunctions.R')
               set.seed(ll+seed*2)
               data.pool.p = cbind(t(dataxp), t(datayp)); len.p = ncol(data.pool.p)
               
               param.X.P = param.Y.P =c()
               
               for (j in 1:n.site){
                 
                 id.poolP = sample(1:len.p, replace=F)
                 
                 data.xbootP = c(data.pool.p[j,id.poolP[1:(length(id.poolP)/2)]])
                 data.ybootP = c(data.pool.p[j,id.poolP[(length(id.poolP)/2 +1):length(id.poolP)]])
                 # FIT GEV
                 modxP <-gev.fit(data.xbootP, show=FALSE)
                 modyP <-gev.fit(data.ybootP, show=FALSE)
                 
                 param.X.P = rbind(param.X.P,modxP$mle); param.Y.P = rbind(param.Y.P,modyP$mle)
                 
               }
               
               # GEV RL
               rl_5_GEVX <- gev.return_level(param.X.P[,1], param.X.P[,2], param.X.P[,3], n.site, time = t0)
               rl_10_GEVX <- gev.return_level(param.X.P[,1], param.X.P[,2], param.X.P[,3], n.site, time = t1)
               rl_20_GEVX <- gev.return_level(param.X.P[,1], param.X.P[,2], param.X.P[,3], n.site, time = t2)
               rl_50_GEVX <- gev.return_level(param.X.P[,1], param.X.P[,2], param.X.P[,3], n.site, time = t5)
               rl_100_GEVX <- gev.return_level(param.X.P[,1], param.X.P[,2], param.X.P[,3], n.site, time = t10)
               
               rl_5_GEVY <- gev.return_level(param.Y.P[,1], param.Y.P[,2], param.Y.P[,3], n.site, time = t0)
               rl_10_GEVY <- gev.return_level(param.Y.P[,1], param.Y.P[,2], param.Y.P[,3], n.site, time = t1)
               rl_20_GEVY <- gev.return_level(param.Y.P[,1], param.Y.P[,2], param.Y.P[,3], n.site, time = t2)
               rl_50_GEVY <- gev.return_level(param.Y.P[,1], param.Y.P[,2], param.Y.P[,3], n.site, time = t5)
               rl_100_GEVY <- gev.return_level(param.Y.P[,1], param.Y.P[,2], param.Y.P[,3], n.site, time = t10)
               
               rtX.boot.GEVP = cbind(rl_5_GEVX,rl_10_GEVX, rl_20_GEVX, rl_50_GEVX, rl_100_GEVX)
               rtY.boot.GEVP = cbind(rl_5_GEVY,rl_10_GEVY, rl_20_GEVY, rl_50_GEVY, rl_100_GEVY)
               rt.diff.bootP = rtY.boot.GEVP - rtX.boot.GEVP
               param.diff.P = param.Y.P - param.X.P
               
               list(rt.diff.bootP,  param.diff.P)
             }


e = Sys.time();
paste("Permuatation Test=", e-s)
##########################################################################################################################################
combine.listP =  list()
for ( m in 1:ncol(r[[1]][[1]])){
  combine.listP[[m]] = do.call(cbind, lapply(r[[1]], function(x) x[,m]))
}

combine.paramP = list()
for ( m in 1:ncol(r[[2]][[1]])){
  combine.paramP[[m]] = do.call(cbind, lapply(r[[2]], function(x) x[,m]))
}

#### P-values ####
sum.bootP =matrix(nrow= n.site, ncol=ncol(rt.diff.GEVP))
for ( m in 1:ncol(rt.diff.GEVP)) {
  for ( j in 1:nrow(combine.listP[[m]])) {
    sum.bootP[j,m] =  sum(ifelse(abs(combine.listP[[m]][j,])> abs(rt.diff.GEVP[j,m]), 1, 0))
  }
}

sum.bootParam =matrix(nrow= n.site, ncol=ncol(param.diff.GEVP))
for ( m in 1:ncol(param.diff.GEVP)) {
  for ( j in 1:nrow(combine.paramP[[m]])) {
    sum.bootParam[j,m] =  sum(ifelse(abs(combine.paramP[[m]][j,])> abs(param.diff.GEVP[j,m]), 1, 0))
  }
}

pvalP = sum.bootP/1000;
pvalP = apply(pvalP, 2, function(x) ifelse(x <= 1e-4, 1e-4,x))
pval.param.each = sum.bootParam/1000
pval.param.each = apply(pval.param.each, 2, function(x) ifelse(x <= 1e-4, 1e-4, x))

######################################################
## Search for the optimal number of basis function
######################################################
Param_opt_basis = OptimalBasis_BIC(pval.param.each[,2], 2, 3, 6)
Rt10_opt_basis = OptimalBasis_BIC(pvalP[,2], 2, 3, 6)
Rt50_opt_basis = OptimalBasis_BIC(pvalP[,4], 2, 3, 6)

basis_param =Param_opt_basis$basis_opt
basis_R10 = Rt10_opt_basis$basis_opt
basis_R50 = Rt50_opt_basis$basis_opt

### Return level Testing ### 
BiP10_2 <- try(BiModel_Cauchy_Test0(pvalP[,2],2,delta.tr, basis.tr[,1:4], basis_R10, index=TRUE))
if(class(BiP10_2)=="try-error"){BiP10_2 <- NA}
BiP10_4 <- try(BiModel_Cauchy_Test0(pvalP[,2],4,delta.tr, basis.tr[,1:4], basis_R10, index=TRUE))
if(class(BiP10_4)=="try-error"){BiP10_4 <- NA}
BiP10_12 <- try(BiModel_Cauchy_Test0(pvalP[,2],12,delta.tr, basis.tr[,1:4], basis_R10, index=TRUE))
if(class(BiP10_12)=="try-error"){BiP10_12 <- NA}

BiP50_2 <- try(BiModel_Cauchy_Test0(pvalP[,4],2,delta.tr, basis.tr[,1:4], basis_R50, index=TRUE))
if(class(BiP50_2)=="try-error"){BiP50_2 <- NA}
BiP50_4 <- try(BiModel_Cauchy_Test0(pvalP[,4],4,delta.tr, basis.tr[,1:4], basis_R50, index=TRUE))
if(class(BiP50_4)=="try-error"){BiP50_4 <- NA}
BiP50_12 <- try(BiModel_Cauchy_Test0(pvalP[,4],12,delta.tr, basis.tr[,1:4], basis_R50, index=TRUE))
if(class(BiP50_12)=="try-error"){BiP50_12 <- NA}

### Parameter Testing ###
BiParam_2 <- try(BiModel_Cauchy_Test0(pval.param.each[,2],2,delta.tr, basis.tr[,1:4], basis_param, index=TRUE))
if(class(BiParam_2)=="try-error"){BiParam_2 <- NA}
BiParam_4 <- try(BiModel_Cauchy_Test0(pval.param.each[,2],4,delta.tr, basis.tr[,1:4], basis_param, index=TRUE))
if(class(BiParam_4)=="try-error"){BiParam_4 <- NA}
BiParam_12 <- try(BiModel_Cauchy_Test0(pval.param.each[,2],12,delta.tr, basis.tr[,1:4], basis_param, index=TRUE))
if(class(BiParam_12)=="try-error"){BiParam_12 <- NA}

save(rt.diff.GEVP,param.diff.GEVP, pvalP, pval.param.each,BiParam_2,BiParam_4,BiParam_12
     , BiP10_2,BiP10_4,BiP10_12, BiP50_2, BiP50_4,BiP50_12
     ,file = paste0("/Data/Simulated Data/Permutation_Test_Result/Scale_Moderate/MLEP_ScaleRt_CauchyN4N12R3S0.2_", ii ,".Rdata"))
