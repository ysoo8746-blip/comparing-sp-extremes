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
setwd("/Users/syun/Baruch College Dropbox/Soo In Yun/Research/Comparing Extremes/JCGS/Rcode_Submit")
source('Functions/General_functions.R')
source('Functions/Nonparam_functions.R')
source('Functions/GEVfunctions.R')
source('Functions/Mirror_pi.R')
source('Functions/Neigh_Function.R')

#### Import US map data #### 
setwd("/Users/syun/Baruch College Dropbox/Soo In Yun/Research/Comparing Extremes/JCGS/Rcode_Submit")
load('Data/NAM_Polygon_index.Rdata')

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
setwd("/Users/syun/Baruch College Dropbox/Soo In Yun/Research/Comparing Extremes/JCGS/Rcode_Submit")
load('Data/Winter_AnnualMax_Precip0.Rdata') ## load the actual precipitation of ERA5 and observation 

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
setwd("/Users/syun/Baruch College Dropbox/Soo In Yun/Research/Comparing Extremes/JCGS/Rcode_Submit")
Param_SimX = data.frame(loc,scale = param.scale, loc = param.loc,
                        shape = param.shape, delta = delta.tr)

#save(Param_SimX, file="Data/Simulation_Parameters.Rdata")

############ Generate Signal (alternative parameters) ############
#### 1. Location signal ####
# basis function
setwd('/Users/syun/Baruch College Dropbox/Soo In Yun/Research/Comparing Extremes/JCGS/Rcode_Submit')
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

setwd('/Users/syun/Baruch College Dropbox/Soo In Yun/Research/Comparing Extremes/JCGS/Rcode_Submit')
#save(Param_SimX, file="SimData_Locdiff_K1.3_P35.Rdata")
#save(Param_SimX, file="SimData_Locdiff_K1.5_P35.Rdata")
#save(Param_SimX, file="SimData_Locdiff_K1.7_P35.Rdata")

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

setwd('/Users/syun/Baruch College Dropbox/Soo In Yun/Research/Comparing Extremes/JCGS/Rcode_Submit/Data')
#save(Param_SimX, file="SimData_Scalediff_K2_P35.Rdata")
#save(Param_SimX, file="SimData_Scalediff_K2.5_P35.Rdata")
#save(Param_SimX, file="SimData_Scalediff_K3_P35.Rdata")

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

setwd('/Users/syun/Baruch College Dropbox/Soo In Yun/Research/Comparing Extremes/JCGS/Rcode_Submit/Data')
#save(Param_SimX, file="SimData_Shapediff_K1.2_P35_2.Rdata")# k.sig = 1.2
#save(Param_SimX, file="SimData_Shapediff_K1_P35_2.Rdata")# k.sig = 1
#save(Param_SimX, file="SimData_Shapediff_K0.8_P35_2.Rdata") # k.sig = 0.8

########################################################################
### P-value Generation 
### Permutation test 
########################################################################
setwd('/Users/syun/Baruch College Dropbox/Soo In Yun/Research/Comparing Extremes/JCGS/Rcode_Submit')
load('Data/SimData_Scalediff_K2_P35.Rdata') ## Load simulated data

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
### Call Max-stable data with r=3 / s =0.2 
setwd('/Users/syun/Baruch College Dropbox/Soo In Yun/Research/Comparing Extremes/JCGS/Rcode_Submit')
data_name <- paste0("MaxStab_R3S0.2/MaxstabDataPrecip_R3_S0.2_", ii, ".Rdata") 
load(data_name)

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
s = Sys.time()
GEV_Rt_Scalediff = Fit_GEV_RT(data1, data2, param_x, param_y,rt.time)
rt.diff.GEVP=GEV_Rt_Scalediff$Rt.diff
param.diff.GEVP=GEV_Rt_Scalediff$Param.diff
e = Sys.time(); e-s

# Set up parallel computing
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
                 
                 id.poolP = sample(1:len.p, replace=F)#;id.poolT = sample(1:len.t, replace=F)
                 
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

pvalP = sum.bootP/1000; #pvalA = sum.bootA/500
pvalP = apply(pvalP, 2, function(x) ifelse(x <= 1e-4, 1e-4,x))
pval.param.each = sum.bootParam/1000
pval.param.each = apply(pval.param.each, 2, function(x) ifelse(x <= 1e-4, 1e-4, x))

######################################################
## Search for the optimal number of basis function
######################################################
setwd('/Users/syun/Baruch College Dropbox/Soo In Yun/Research/Comparing Extremes/JCGS/Rcode_Submit')
s = Sys.time()
source('/Functions/Optimal_Basis_functions.R')
Param_opt_basis = OptimalBasis_BIC(pval.param.each[,2], 2, 3, 6)
Rt10_opt_basis = OptimalBasis_BIC(pvalP[,2], 2, 3, 6)
Rt50_opt_basis = OptimalBasis_BIC(pvalP[,4], 2, 3, 6)

basis_param =Param_opt_basis$basis_opt
basis_R10 = Rt10_opt_basis$basis_opt
basis_R50 = Rt50_opt_basis$basis_opt
e = Sys.time();
paste("OptimalBasis Computation=", e-s)

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

# save(rt.diff.GEVP,param.diff.GEVP, pvalP, pval.param.each,BiParam_2,BiParam_4,BiParam_12
#      , BiP10_2,BiP10_4,BiP10_12, BiP50_2, BiP50_4,BiP50_12
#      ,file = paste0("Result/Shape_Rt_R3S0.2B9_Sig0.8/MLEP_ShapeRt_CauchyN4N12R3S0.2_", ii ,".Rdata"))

######################################################
## Plot Multiple test results
######################################################
source('/Functions/Optimal_Basis_functions.R')
m1= c('BH', 'JS', 'Mirror', 'BiCLfdr-1');m2= c('BH', 'JS', 'Mirror', 'BiCLfdr-2')
m4= c('BH', 'JS', 'Mirror', 'BiCLfdr-4');m12= c('BH', 'JS', 'Mirror', 'BiCLfdr-12')
m= c('BH', 'JS', 'Mirror', 'BiCLfdr')

# Location parameter
Loc_Multi_Strong = MyMulti_outputNeiLoc("Size2Diff_P33_R3S0.2B18_Sig1.7","Loc",Method=m2, Nei=2)
Loc_Multi_Mod = MyMulti_outputNeiLoc("Size2Diff_P33_R3S0.2B18_Sig1.5","Loc",Method=m2, Nei=2)
Loc_Multi_Weak= MyMulti_outputNeiLoc("Size2Diff_P33_R3S0.2B18_Sig1.3","Loc",Method=m2, Nei=2)
Loc_ResOpt_4 = list(Weak = Loc_Multi_Weak, Moderate= Loc_Multi_Mod, Strong=Loc_Multi_Strong)
#Loc_ResOpt_12 = list(Weak = Loc_Multi_Weak, Moderate= Loc_Multi_Mod, Strong=Loc_Multi_Strong)

# Scale Parameter
Scale_Multi_Strong = MyMulti_outputNeiN(Param="Scale",N1=2,Method=m2, name="Sig3", Nei=2)
Scale_Multi_Mod = MyMulti_outputNeiN(Param="Scale",N1 = 2,Method=m2, name="Sig2.5", Nei=2)
Scale_Multi_Weak = MyMulti_outputNeiN(Param="Scale",N1 = 2,Method=m2, name="Sig2", Nei=2)
Scale_ResOpt_4 = list(Weak = Scale_Multi_Weak, Moderate= Scale_Multi_Mod, Strong=Scale_Multi_Strong)
#Scale_ResOpt_12 = list(Weak = Scale_Multi_Weak, Moderate= Scale_Multi_Mod, Strong=Scale_Multi_Strong)

# Shape Parameter
Shape_Multi_Strong = MyMulti_outputNeiN(Param="Shape",N1=2,Method=m2, name="Sig1.2", Nei=2)
Shape_Multi_Mod = MyMulti_outputNeiN(Param="Shape",N1 = 2,Method=m2, name="Sig1", Nei=2)
Shape_Multi_Weak = MyMulti_outputNeiN(Param="Shape",N1 = 2,Method=m2, name="Sig0.8", Nei=2)
 Shape_ResOpt_4 = list(Weak = Shape_Multi_Weak, Moderate= Shape_Multi_Mod, Strong=Shape_Multi_Strong)
# Shape_ResOpt_12 = list(Weak = Shape_Multi_Weak, Moderate= Shape_Multi_Mod, Strong=Shape_Multi_Strong)

### PLOT ### 
# Figure 1
G.Loc_N4_12 = My_GG_MultiTestLoc(Loc_ResOpt_4, Loc_ResOpt_12, "Loc", Method=m_n0, Sig=s)
# Figure 2
G.Scale_N4_12 = My_GG_MultiTest(Scale_ResOpt_4,  Scale_ResOpt_12, "Scale", Method=m_n0, Sig=s)
# Figure 4
G.Shape_N4_12 = My_GG_MultiTest(Shape_ResOpt_4,  Shape_ResOpt_12, "Shape", Method=m_n0, Sig=s)


# Parameter and return level estimation
##### ---- Figure 3 : Scale Difference ----- ####
# Weak Signal
setwd('Result/Scale_Rt_Size2Diff_P33_R3S0.2B9_Sig2')
datafile_names0 = list.files();load(datafile_names0[7])
R10_Scale_diff_Weak = rt.diff.GEVP[,2];R50_Scale_diff_Weak = rt.diff.GEVP[,4] ## j = 1/ j = 7
Scale_diff_Weak = param.diff.GEVP[,2] # j = 8

# Moderate Signal 
setwd('Result/Scale_Rt_Size2Diff_P33_R3S0.2B9_Sig2.5')
datafile_names0 = list.files();load(datafile_names0[7])
R10_Scale_diff_Moderate = rt.diff.GEVP[,2];R50_Scale_diff_Moderate = rt.diff.GEVP[,4] ## j = 1/ j = 7
Scale_diff_Moderate = param.diff.GEVP[,2] # j = 8

# Strong Signal 
setwd('Result/Scale_Rt_Size2Diff_P33_R3S0.2B9_Sig3')
datafile_names0 = list.files();load(datafile_names0[7])
R10_Scale_diff_Strong = rt.diff.GEVP[,2];R50_Scale_diff_Strong = rt.diff.GEVP[,4] ## j = 1/ j = 7
Scale_diff_Strong = param.diff.GEVP[,2] # j = 8
# Apply the function to each case
R10_df_weak     <- reshape_loc_diff(R10_Scale_diff_Weak, "Weak", delta.tr)
R10_df_moderate <- reshape_loc_diff(R10_Scale_diff_Moderate, "Moderate", delta.tr)
R10_df_strong   <- reshape_loc_diff(R10_Scale_diff_Strong, "Strong", delta.tr)

R50_df_weak     <- reshape_loc_diff(R50_Scale_diff_Weak, "Weak", delta.tr)
R50_df_moderate <- reshape_loc_diff(R50_Scale_diff_Moderate, "Moderate", delta.tr)
R50_df_strong   <- reshape_loc_diff(R50_Scale_diff_Strong, "Strong", delta.tr)

# Build one combined data frame
plot_df <- bind_rows(
  make_df(R10_Scale_diff_Weak,     R50_Scale_diff_Weak,     "Weak",     delta.tr),
  make_df(R10_Scale_diff_Moderate, R50_Scale_diff_Moderate, "Moderate", delta.tr),
  make_df(R10_Scale_diff_Strong,   R50_Scale_diff_Strong,   "Strong",   delta.tr)
)

# Plot: 2 (rows: 10/50-year) × 3 (cols: Weak/Moderate/Strong)
# Note: `scales = "free_y"` lets the 10-year and 50-year rows use different y-axis ranges.
plot_df$delta <- factor(plot_df$delta, levels = c("Null", "Alternative"))
RT_scalediff = ggplot(plot_df, aes(x = loc_id, y = diff)) +
  
  # Null locations (background)
  geom_point(
    data = subset(plot_df, delta == "Null"),
    aes(shape = delta),
    color = "grey60",
    size  = 1.0,
    alpha = 0.7
  ) +
  
  # Alternative locations (foreground)
  geom_point(
    data = subset(plot_df, delta == "Alternative"),
    aes(shape = delta),
    color = "black",
    size  = 1.6,
    alpha = 0.9
  ) +
  
  scale_shape_manual(
    values = c("Null" = 1,        # open circle
               "Alternative" = 16 # filled circle
    ),
    name = "Location Type"
  ) +
  
  facet_grid(horizon ~ strength, scales = "free_y") +
  
  labs(
    title = NULL,
    x = "Location Index",
    y = "Return-Level Difference"
  ) +
  
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(size = 20, face = "bold"),
    strip.background = element_rect(fill = "grey85", color = "black"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 20),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text.x  = element_text(size = 14),
    axis.text.y  = element_text(size = 14),
    legend.position = "bottom"
  ) +
  guides(shape = guide_legend(override.aes = list(size = 4)))

##### ---- Figure 5 : Shape Difference ----- ####
# Weak signal
setwd('Result/Shape_Rt_Size2Diff_P33_R3S0.2B9_Sig0.8')
datafile_names0 = list.files();load(datafile_names0[4])
R10_Shape_diff_Weak = rt.diff.GEVP[,2];R50_Shape_diff_Weak = rt.diff.GEVP[,4] # j =3/j= 4
Shape_diff_Weak = param.diff.GEVP[,3] # j = 3

# Moderate Signal
setwd('Result/Shape_Rt_Size2Diff_P33_R3S0.2B9_Sig1')
datafile_names0 = list.files();load(datafile_names0[4])
R10_Shape_diff_Moderate = rt.diff.GEVP[,2];R50_Shape_diff_Moderate = rt.diff.GEVP[,4] # j =3/j= 4
Shape_diff_Moderate = param.diff.GEVP[,3] # j = 3

# Strong Signal 
setwd('Result/Shape_Rt_Size2Diff_P33_R3S0.2B9_Sig1.2')
datafile_names0 = list.files();load(datafile_names0[4])
R10_Shape_diff_Strong = rt.diff.GEVP[,2];R50_Shape_diff_Strong = rt.diff.GEVP[,4] # j =3/j= 4
Shape_diff_Strong = param.diff.GEVP[,3] # j = 3

# Apply the function to each case
R10_df_weak     <- reshape_loc_diff(R10_Shape_diff_Weak, "Weak", delta.tr)
R10_df_moderate <- reshape_loc_diff(R10_Shape_diff_Moderate, "Moderate", delta.tr)
R10_df_strong   <- reshape_loc_diff(R10_Shape_diff_Strong, "Strong", delta.tr)

R50_df_weak     <- reshape_loc_diff(R50_Shape_diff_Weak, "Weak", delta.tr)
R50_df_moderate <- reshape_loc_diff(R50_Shape_diff_Moderate, "Moderate", delta.tr)
R50_df_strong   <- reshape_loc_diff(R50_Shape_diff_Strong, "Strong", delta.tr)

# Build one combined data frame
plot_df <- bind_rows(
  make_df(R10_Shape_diff_Weak,     R50_Shape_diff_Weak,     "Weak",     delta.tr),
  make_df(R10_Shape_diff_Moderate, R50_Shape_diff_Moderate, "Moderate", delta.tr),
  make_df(R10_Shape_diff_Strong,   R50_Shape_diff_Strong,   "Strong",   delta.tr)
)

# Plot: 2 (rows: 10/50-year) × 3 (cols: Weak/Moderate/Strong)
# Note: `scales = "free_y"` lets the 10-year and 50-year rows use different y-axis ranges.
plot_df$delta <- factor(plot_df$delta, levels = c("Null", "Alternative"))
RT_shapediff = ggplot(plot_df, aes(x = loc_id, y = diff)) +
  
  # Null locations (background)
  geom_point(
    data = subset(plot_df, delta == "Null"),
    aes(shape = delta),
    color = "grey60",
    size  = 1.0,
    alpha = 0.7
  ) +
  
  # Alternative locations (foreground)
  geom_point(
    data = subset(plot_df, delta == "Alternative"),
    aes(shape = delta),
    color = "black",
    size  = 1.6,
    alpha = 0.9
  ) +
  
  scale_shape_manual(
    values = c("Null" = 1,        # open circle
               "Alternative" = 16 # filled circle
    ),
    name = "Location Type"
  ) +
  
  facet_grid(horizon ~ strength, scales = "free_y") +
  
  labs(
    title = NULL,
    x = "Location Index",
    y = "Return-Level Difference"
  ) +
  
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(size = 20, face = "bold"),
    strip.background = element_rect(fill = "grey85", color = "black"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 20),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text.x  = element_text(size = 14),
    axis.text.y  = element_text(size = 14),
    legend.position = "bottom"
  ) +
  guides(shape = guide_legend(override.aes = list(size = 4)))

