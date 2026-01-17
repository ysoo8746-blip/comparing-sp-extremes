##### Load Functions and data 
setwd('/Users/syun/Baruch College Dropbox/Soo In Yun/Research/Comparing Extremes/JCGS/Rcode_Submit/Functions')
source('General_functions.R')
source('Nonparam_functions.R')
source('GEVfunctions.R')
source('Mirror_pi.R')
source('Neigh_Function.R')
setwd('/Users/syun/Baruch College Dropbox/Soo In Yun/Research/Comparing Extremes/JCGS/Rcode_Submit/Data')
load('Winter_AnnualMax_Precip0.Rdata')

n.obs <- ncol(annualmax_precip_obs)-2
n.site<- nrow(loc)
loc = as.matrix(loc)
cov <- 'whitmat'
range <- 3
smooth <-0.2
p=1
t0 = 5*p
t1 = 10*p
t2 = 20*p
t5 =50*p
t10 = 100*p

N = nrow(loc)

###############################################################################################
# Use two different maxstable model to estimate the param
###############################################################################################
loc0 = loc

#### GEV function #####
df_dg0 = function(param, return){
  pm = param
  rt = return
  yRes = c(1, -pm[3]^(-1)*(1-rt^(-pm[3])), pm[2]*pm[3]^(-2)*(1-rt^(-pm[3]))-pm[2]*pm[3]^(-1)*rt^(-pm[3])*log(rt))
  
  return(yRes)
}
############################################################################################
############################################################################################################
# combine function
#set.seed(seed)
library(doParallel)
registerDoParallel(cores=12)
detectCores()

comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

datP1 = annualmax_precip_obs[,-c(1:2)]
datP2 = annualmax_precip_era5[,-c(1:2)]

loc = as.matrix(loc)

paramx.mle.p = paramy.mle.p =c()
for ( i in 1:n.site){
  datX.P = as.numeric(datP1[i,])
  datY.P = as.numeric(datP2[i,])

  modxp <-gev.fit(datX.P, show=FALSE)
  modyp <-gev.fit(datY.P, show=FALSE)

  paramx.mle.p = rbind(paramx.mle.p, modxp$mle)
  paramy.mle.p = rbind(paramy.mle.p, modyp$mle)

  if (i %% 100 ==0) print(i)
}

# GEV RL for precip
rl_5_GEVX <- gev.return_level(paramx.mle.p[,1], paramx.mle.p[,2], paramx.mle.p[,3], n.site, time = t0)
rl_10_GEVX <- gev.return_level(paramx.mle.p[,1], paramx.mle.p[,2], paramx.mle.p[,3], n.site, time =  t1)
rl_20_GEVX <- gev.return_level(paramx.mle.p[,1], paramx.mle.p[,2], paramx.mle.p[,3], n.site, time =  t2)
rl_50_GEVX <- gev.return_level(paramx.mle.p[,1], paramx.mle.p[,2], paramx.mle.p[,3], n.site, time =  t5)
rl_100_GEVX <- gev.return_level(paramx.mle.p[,1], paramx.mle.p[,2], paramx.mle.p[,3], n.site, time =  t10)
rtX.GEV.P = cbind(rl_5_GEVX, rl_10_GEVX, rl_20_GEVX, rl_50_GEVX, rl_100_GEVX)
colnames(rtX.GEV.P) = c('Rt5','Rt10', 'Rt20', 'Rt50', 'Rt100')
summary(rtX.GEV.P) ## ESTIMATED RETURN LEVEL OF X
rl_5_GEVY <- gev.return_level(paramy.mle.p[,1], paramy.mle.p[,2], paramy.mle.p[,3], n.site, time = t0)
rl_10_GEVY <- gev.return_level(paramy.mle.p[,1], paramy.mle.p[,2], paramy.mle.p[,3], n.site, time = t1)
rl_20_GEVY <- gev.return_level(paramy.mle.p[,1], paramy.mle.p[,2], paramy.mle.p[,3], n.site, time = t2)
rl_50_GEVY <- gev.return_level(paramy.mle.p[,1], paramy.mle.p[,2], paramy.mle.p[,3], n.site, time = t5)
rl_100_GEVY <- gev.return_level(paramy.mle.p[,1], paramy.mle.p[,2], paramy.mle.p[,3], n.site, time = t10)
rtY.GEV.P = cbind(rl_5_GEVY,rl_10_GEVY, rl_20_GEVY, rl_50_GEVY, rl_100_GEVY)
colnames(rtY.GEV.P) = c('Rt5','Rt10', 'Rt20', 'Rt50', 'Rt100')
summary(rtY.GEV.P) ## ESTIMATED RETURN LEVEL OF Y

param.diffpercip = paramy.mle.p- paramx.mle.p
rt.diff.GEVP = rtY.GEV.P -rtX.GEV.P
colnames(rt.diff.GEVP) = c('Rt5_diff','Rt10_diff', 'Rt20_diff', 'Rt50_diff', 'Rt100_diff')
summary(rt.diff.GEVP) ## ESTIMATED RETURN LEVEL DIFF
trend.diff = trend.y - trend.x

##################################################################################################################
# Permutation test
#set.seed(seed)
s = Sys.time()
ptime <- system.time({
  r <- foreach(ll=seq_len(1000),.combine='comb',.multicombine=TRUE,.init=list(list(), list(),list())
               , .packages = c("extRemes", "eva","SpatialTools")) %dopar% {
                 
                 set.seed(ll+seed*2)
                 data.pool.p = cbind(datP1, datP2); len.p = ncol(data.pool.p)
                 
                 param.X.P0 = param.Y.P0 =c()
                 
                 for (j in 1:n.site){
                   
                   set.seed(ll+seed*2+j)
                   id.poolP = sample(1:len.p, replace=F)#;id.poolT = sample(1:len.t, replace=F)
                   
                   data.xbootP = as.numeric(c(data.pool.p[j,id.poolP[1:n.obs]]))
                   data.ybootP = as.numeric(c(data.pool.p[j,id.poolP[(n.obs+1):(2*n.obs)]]))
                   # FIT GEV
                   modxP <-gev.fit(data.xbootP, ydat = cov, mul =1,show=FALSE)
                   modyP <-gev.fit(data.ybootP, ydat = cov, mul =1,show=FALSE)
                   
                   param.X.P0 = rbind(param.X.P0,modxP$mle); param.Y.P0 = rbind(param.Y.P0,modyP$mle)
                   if (j %% 100==0) print(j)
                 }
                 newlocx = param.X.P0[,1] + param.X.P0[,2]*cov[74]
                 newlocy = param.Y.P0[,1] + param.Y.P0[,2]*cov[74]
                 param.X.P = cbind(newlocx, param.X.P0[,-c(1:2)])
                 param.Y.P = cbind(newlocy, param.Y.P0[,-c(1:2)])
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
                 trend.diff.P = param.Y.P0[,2]- param.X.P0[,2]
                 
                 list(rt.diff.bootP,  param.diff.P,trend.diff.P)
                 #list(rt.diff.boot, rt.diff.bootA)
               }
  
})
e = Sys.time(); e-s

#paramy.mle
##########################################################################################################################################
param.diffpercip = paramy.mle.p- paramx.mle.p
rt.diff.GEVP = rtY.GEV.P -rtX.GEV.P

combine.listP =  list()
for ( m in 1:ncol(r[[1]][[1]])){
  combine.listP[[m]] = do.call(cbind, lapply(r[[1]], function(x) x[,m]))
}

combine.paramP = list()
for ( m in 1:ncol(r[[2]][[1]])){
  combine.paramP[[m]] = do.call(cbind, lapply(r[[2]], function(x) x[,m]))
}

combine.paramT = do.call(cbind, r[[3]])

## P-values
sum.bootP =matrix(nrow= n.site, ncol=ncol(rt.diff.GEVP))
for ( m in 1:ncol(rt.diff.GEVP)) {
  for ( j in 1:nrow(combine.listP[[m]])) {
    sum.bootP[j,m] =  sum(ifelse(abs(combine.listP[[m]][j,])> abs(rt.diff.GEVP[j,m]), 1, 0))
  }
}

sum.bootParam =matrix(nrow= n.site, ncol=ncol(param.diffpercip))
for ( m in 1:ncol(param.diffpercip)) {
  for ( j in 1:nrow(combine.paramP[[m]])) {
    sum.bootParam[j,m] =  sum(ifelse(abs(combine.paramP[[m]][j,])> abs(param.diffpercip[j,m]), 1, 0))
  }
}

sum.bootTrend =c()
for ( j in 1:nrow(combine.paramT)) {
  sum.bootTrend[j] =  sum(ifelse(abs(combine.paramT[j,])> abs(trend.diff[j]), 1, 0))
}


pvalP = sum.bootP/1000;#pvalA = sum.bootA/500
pvalP = apply(pvalP, 2, function(x) ifelse(x <= 1e-4, 1e-4,x))
pval.param.each = sum.bootParam/1000
pval.param.each = apply(pval.param.each, 2, function(x) ifelse(x <= 1e-4, 1e-4, x))
pval.trend = sum.bootTrend/1000
pval.trend = ifelse(pval.trend <= 1e-4, 1e-4, pval.trend)

##############################################
### Multiple testing #########################
##############################################
### Basis functions

delta.tr0 = sample(c(0,1), n.site,replace=TRUE, prob=c(2/3, 1/3))
source('Functions/Optimal_Basis_functions.R')
Trend_opt_basis = OptimalBasis_BIC(pval.trend, 4, 3, 6)
Loc_opt_basis = OptimalBasis_BIC(pval.param.each[,1], 4, 3, 6)
Scale_opt_basis = OptimalBasis_BIC(pval.param.each[,2], 4, 3, 6)
Shape_opt_basis = OptimalBasis_BIC(pval.param.each[,3], 4, 3, 6)
Rt10_opt_basis = OptimalBasis_BIC(pvalP[,2], 4, 3, 6)
Rt50_opt_basis = OptimalBasis_BIC(pvalP[,4], 4, 3, 6)

basis_trend = Trend_opt_basis$basis_opt
basis_loc =Loc_opt_basis$basis_opt;basis_sc =Scale_opt_basis$basis_opt
basis_sh =Shape_opt_basis$basis_opt
basis_R10 = Rt10_opt_basis$basis_opt
basis_R50 = Rt50_opt_basis$basis_opt
basis.tr =Create_Basis(c(-125,-66),c(25, 50), 3, 3, loc )

##### Multiple testing of trend, parameters, and RT ######
trend_4 =try(BiModel_Cauchy_Real(pval.trend, 4, basis.tr[, 1:6], basis_trend))
if (class(trend_4)=="try-error"){trend_2 <- NA}

loc_4 =try(BiModel_Cauchy_Real(pval.param.each[,1], 4, basis.tr[, 1:6], basis_loc))
if (class(loc_4)=="try-error"){loc_2 <- NA}

Scale_4 =try(BiModel_Cauchy_Real(pval.param.each[,2], 4, basis.tr[, 1:6], basis_sc))
if (class(Scale_4)=="try-error"){Scale_2 <- NA}

Shape_4 =try(BiModel_Cauchy_Real(pval.param.each[,3], 4, basis.tr[, 1:6], basis_sh))
if (class(Shape_4)=="try-error"){Shape_2 <- NA}

test10_4 =try(BiModel_Cauchy_Real(pvalP[,2], 4, basis.tr[, 1:6], basis_R10))
if (class(test10_4)=="try-error"){test10_4 <- NA}
test50_4 =try(BiModel_Cauchy_Real(pvalP[,4], 2, basis.tr[, 1:6], basis_R50))
if (class(test50_4)=="try-error"){test50_4 <- NA}

