########Function #############
## neigh 
neigh = function(x, nei){
  idres = c()
  for ( j in 1:n.site){
    a = round((dist2(matrix(loc[j,], ncol=2), loc))[1,],2)
    id = sort(a, index=TRUE)$ix[2:(nei+1)]
    idres = rbind(idres, id)
  }
  
  xres= c()
  for ( j in 1:n.site) { 
    xres = rbind(xres, c(x[j],x[idres[j,]] ))
  }
  return(xres)
}
## pvalue to stat
pvalue_to_teststat = function(x){
  
  out = sort(x, index.return=TRUE)
  pval.s = out$x
  ind = out$ix
  
  ecdf.p1 = ecdf(pval.s);ecdf.p1 = ecdf.p1(pval.s);
  cdf.p1 = ifelse(ecdf.p1 > 1-1e-5, 1-1e-5, ecdf.p1)
  tstat = rep(NA, n.site);tstat[ind] = qnorm(cdf.p1);
  return(tstat)
}

## Marginal 
marg_ftn= function(nei, null, nullalt, altnull, alt, pi){
  pi_0 = pi[,1]; pi_n = pi[,-1]
  marg = nullmarg = 0
  for ( i in 1:nei) {
    margtmp = pi_0*pi_n[,i]*null[,i] + pi_0*(1-pi_n[,i])*nullalt[,i] +(1-pi_0)*pi_n[,i]*altnull[,i] +(1-pi_0)*(1-pi_n[,i])*alt[,i] 
    marg = marg + margtmp
    nullmargtmp = pi_0*pi_n[,i]*null[,i] + pi_0*(1-pi_n[,i])*nullalt[,i]
    nullmarg = nullmarg +nullmargtmp
  }
  
  return (list(Marg =marg, NMarg = nullmarg/marg))
}

marg_ftn_each= function(nei, null, nullalt, altnull, alt, pi){
  pi_0 = pi[,1]; pi_n = pi[,-1]
  marg = nullmarg = c()
  for ( i in 1:nei) {
    margtmp = pi_0*pi_n[,i]*null[,i] + pi_0*(1-pi_n[,i])*nullalt[,i] +(1-pi_0)*pi_n[,i]*altnull[,i] +(1-pi_0)*(1-pi_n[,i])*alt[,i] 
    marg =  cbind(marg ,margtmp)
    nullmargtmp = pi_0*pi_n[,i]*null[,i] + pi_0*(1-pi_n[,i])*nullalt[,i]
    nullmarg = cbind(nullmarg, nullmargtmp)
  }
  
  return (list(Marg =marg, NMarg = nullmarg/marg))
}
## Null alt joint dist
null_joint = function(x){ 
  x0 = x[,1]
  xx = x[,-1]
  f00_res = c()
  if (ncol(x) > 2){
  for (i in 1:ncol(xx)){ 
    f00_res = cbind(f00_res,dmvnorm(cbind(x0, xx[,i]), cbind(0,0), diag(1,2)))
  }
  } else if (ncol(x)<=2){
    f00_res = dmvnorm(cbind(x0, xx), cbind(0,0), diag(1,2))
  }
  return(f00_res)
}

nullalt_joint = function(x, basis, delta) {
  x0 = x[,1]; xx =x[,-1]
  A = Mirror_Piest(x[,1],basis, delta )
  piest = A$piest
  if (ncol(x) > 2) {
    f10_res = matrix(rep(A$f1est, ncol(xx)), ncol= ncol(xx))
  } else if (ncol(x) <= 2) {
    f10_res = matrix(A$f1est, ncol=1)
  }
  f01_res = c()
  if (ncol(x) > 2){
  for ( i in 1:ncol(xx)){ 
    
    R =Mirror_Piest(xx[,i], basis, delta)
    f01_res = cbind(f01_res,R$f1est)
  }
  } else if (ncol(x)<=2){
    R =Mirror_Piest(xx, basis, delta)
    f01_res = R$f1est
  }
  f11_res = f10_res*f01_res
  return(list(f10est = f10_res,f01est = f01_res, f11est= f11_res, piest =piest, M.res = A$MirrRes))
}

## GMM
GMM_mean_cov = function(M, null, nullalt, altnull, alt,pi,nei, stat){
  
  pi_0 = pi[,1]; pi_n = pi[,-1]
  w00 = w01 = w10 = w11 = c()
  for (i in 1:nei){
    w00 = cbind(w00, as.numeric(pi_0*pi_n[,i]*null[,i]/M))
    w01 = cbind(w01, as.numeric(pi_0*(1-pi_n[,i])*nullalt[,i]/M))
    w10 = cbind(w10, as.numeric((1-pi_0)*pi_n[,i]*altnull[,i]/M))
    w11 = cbind(w11, as.numeric((1-pi_0)*(1-pi_n[,i])*alt[,i]/M))
  }
  N00 = apply(w00, 2, sum);N01 = apply(w01, 2, sum);N10 = apply(w10, 2, sum);N11 = apply(w11, 2, sum)
  
  # Mean est * Cov est
  mu_00 = mu_01 = mu_10 = mu_11 = c()
  f00_res = f10_res = f01_res = f11_res = c()
  for (i in 1:nei){
    mu_00 = c(0,0)
    mu_01tmp = (1/N01[i])*apply(w01[,i]*stat[,c(1,i+1)], 2, sum); mu_01tmp = c(0,mu_01tmp[2])
    mu_10tmp = (1/N10[i])*apply(w10[,i]*stat[,c(1,i+1)], 2, sum); mu_10tmp = c(mu_10tmp[1],0)
    mu_11tmp =  (1/N11[i])*apply(w11[,i]*stat[,c(1,i+1)], 2, sum)
    mu_01 = cbind(mu_01,mu_01tmp);mu_10 = cbind(mu_10,mu_10tmp); mu_11 = cbind(mu_11, mu_11tmp)
    
    #Sig_00 = matrix(0,ncol=2, nrow=2);
    Sig_00 = matrix((1/N00[i])*apply(diag(w00[,i]) %*% t(apply(stat[,c(1,i+1)], 1, function(x) (x) %*% t(x))), 2,sum), ncol=2)
    Sig_01 = matrix((1/N01[i])*apply(diag(w01[,i]) %*% t(apply(stat[,c(1,i+1)]-mu_01tmp, 1, function(x) (x) %*% t(x))), 2,sum), ncol=2)
    Sig_10 = matrix((1/N10[i])*apply(diag(w10[,i]) %*% t(apply(stat[,c(1,i+1)]-mu_10tmp, 1, function(x) (x) %*% t(x))), 2,sum),ncol=2)
    Sig_11 = matrix((1/N11[i])*apply(diag(w11[,i]) %*% t(apply(stat[,c(1,i+1)]-mu_11tmp, 1, function(x) (x) %*% t(x))), 2,sum),ncol=2)
    
    Sig_00[1,1] =Sig_00[2,2] = 1 ; Sig_01[1,1] = 1; Sig_10[2,2]=1
    
    f00_res = cbind(f00_res, dmvnorm(stat[,c(1,i+1)], c(0,0), Sig_00))
    f01_res = cbind(f01_res, dmvnorm(stat[,c(1,i+1)], mu_01tmp, Sig_01))
    f10_res = cbind(f10_res, dmvnorm(stat[,c(1,i+1)], mu_10tmp, Sig_10))
    f11_res = cbind(f11_res, dmvnorm(stat[,c(1,i+1)], mu_11tmp, Sig_11))
  }
  return(list(f00 = f00_res, f01 = f01_res, f10 = f10_res, f11 = f11_res))
  
}