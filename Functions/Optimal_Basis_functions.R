Mirror_Piest_noMT = function(p_value, basis.tr0, fdr.c=0.1, JS=TRUE){
  N = length(p_value); Nm = dim(basis.tr0)[2]
  pval = p_value
  out = sort(pval, index.return=TRUE)
  pval = out$x
  ind = out$ix
  basis = basis.tr0[ind,]
 # delta = delta.tr0[ind]
  df = c(pval[1],diff(pval))       # z(i)-z(i-1)
  f1.int = f1 = dexp(pval,100)     # initial yi=f1(z)
  pi0  = runif(N,0,1)          # initial pi(s)
  f0 = rep(1,N)                # known   f0(zi)
  
  ### ESTIMATION METHOD
  q = pi0*f0/(pi0*f0+(1-pi0)*f1);d = ncol(basis)
  beta.int = beta.cur = rep(1, d)  # initial value
  K = 100
  lik = -1e5
  iter = 1
  epsilon = 1
  while(iter <= K & epsilon > 1e-6)
  {
    # E-step:
    q = pi0*f0/(pi0*f0+(1-pi0)*f1)
    
    # M-step (Newton method):
    Q = function(beta) {
      sum(-q*(basis %*% beta)-log(1/(exp(basis %*% beta)+1)))
    }
    grr = function(beta) {
      A = matrix(-q, N, Nm)*basis
      C = exp(basis %*% beta)
      C[C==Inf] = 1e+30
      B = basis*matrix(C/(C+1), N, Nm)
      return(apply(A+B,2,sum))
    }
    fit.optim = nlminb(beta.int, Q, gradient=grr)   ### Initial value !!!
    beta.cur = fit.optim$par
    pi0 = plogis(basis %*% beta.cur)
    pi0 = ifelse(pi0 < 1e-5, 1e-5, pi0); pi0 = ifelse(pi0 > 1-1e-5, 1-1e-5, pi0)
    # PAVA
    q1 = 1-q
    qsum = sum(q1)
    y = -df*qsum/q1
    
    ### modify the code slightly
    y[is.na(y)] = -Inf
    y[y==0] = max(y[y!=0])
    y[y==-Inf] = min(y[y!=-Inf])
    
    out = pava(y, q1, decreasing=TRUE, long.out=FALSE, stepfun=FALSE)
    f1 = -1/out
    f1 = f1/sum( f1*df )
    
    iter = iter + 1
    lik.cur = sum(log(pi0*f0+(1-pi0)*f1))
    epsilon = (lik.cur - lik[length(lik)])/abs(lik[length(lik)])
    lik = c(lik, lik.cur)
    if(iter%%200==0) print(iter)
  }
  
  f0.marg = f1.marg = pi0.marg = p.marg = rep(NA, N)
  f0.marg[ind] = f0;f1.marg[ind] = f1;pi0.marg[ind] = pi0
  p.marg = pi0.marg*f0.marg+(1-pi0.marg)*f1.marg
  p.cond = pi0.marg*f0.marg/p.marg
  
  return(list(piest=pi0.marg, lfdr= p.cond, f0est = f0.marg, f1est = f1.marg))
}
CCT <- function(pvals, weights=NULL){
  #### check if there is NA
  if(sum(is.na(pvals)) > 0){
    stop("Cannot have NAs in the p-values!")
  }
  
  #### check if all p-values are between 0 and 1
  if((sum(pvals<0) + sum(pvals>1)) > 0){
    stop("All p-values must be between 0 and 1!")
  }
  
  #### check if there are p-values that are either exactly 0 or 1.
  is.zero <- (sum(pvals==0)>=1)
  is.one <- (sum(pvals==1)>=1)
  if(is.zero && is.one){
    stop("Cannot have both 0 and 1 p-values!")
  }
  if(is.zero){
    return(0)
  }
  if(is.one){
    warning("There are p-values that are exactly 1!")
    return(1)
  }
  
  #### check the validity of weights (default: equal weights) and standardize them.
  if(is.null(weights)){
    weights <- rep(1/length(pvals),length(pvals))
  }else if(length(weights)!=length(pvals)){
    stop("The length of weights should be the same as that of the p-values!")
  }else if(sum(weights < 0) > 0){
    stop("All the weights must be positive!")
  }else{
    weights <- weights/sum(weights)
  }
  
  #### check if there are very small non-zero p-values
  is.small <- (pvals < 1e-16)
  if (sum(is.small) == 0){
    cct.stat <- sum(weights*tan((0.5-pvals)*pi))
  }else{
    cct.stat <- sum((weights[is.small]/pvals[is.small])/pi)
    cct.stat <- cct.stat + sum(weights[!is.small]*tan((0.5-pvals[!is.small])*pi))
  }
  
  #### check if the test statistic is very large.
  if(cct.stat > 1e+15){
    pval <- 2*(1/cct.stat)/pi
  }else{
    pval <- 2*(1-pcauchy(abs(cct.stat)))
  }
  return(pval)
}
Create_Basis = function(range1, range2, n1, n2, loc){
  ### B-spline basis function
  sp1 <- create.bspline.basis(rangeval=range1, nbasis=n1, norder=n1) ## What is coefficient?
  sp2 <- create.bspline.basis(rangeval=range2, nbasis=n2, norder=n2)
  
  eval.sp1 <- eval.basis(loc[,1], sp1)
  eval.sp2 <- eval.basis(loc[,2], sp2)
  
  ## creat design matrix for least square esitmates of coeff ##
  eval.sp <- matrix(NA, n.site, ncol(eval.sp1)*ncol(eval.sp2))
  for (i in 1:n.site){
    eval.sp[i,] <- kronecker(eval.sp1[i,], eval.sp2[i,])
  }
  basis.tr = eval.sp
  return(basis.tr)
}

OptimalBasis_BIC = function(pval, nei, n1,n2){
  
  eval.sp = Create_Basis(c(-125,-66),c(25, 50), 3, 3, loc )
  basis.tr00 = eval.sp[,1:4]
  
  ### different numb of basis functions
  numb = expand.grid(n1:n2, n1:(n2-1))
  numb = cbind(numb, numb[,1]*numb[,2]) 
  numb = data.frame(numb)
  names(numb) = c('n1', 'n2', 'prod')
  numb = numb %>% arrange(prod, n1, n2)
  bicvals = c()
  for ( ii in 1:nrow(numb)){
    
    basis.tr = Create_Basis(c(-125,-66),c(25, 50), numb[ii,1], numb[ii,2], loc )
    idres = distres =c()
    pnei = neigh(pval,nei)
    pnei = apply(pnei, 2, function(x) ifelse(x < 1e-5, 1e-5, ifelse(x > 1-1e-5, 1-1e-5,x)))
    
    for ( j in 1:n.site){
      a = round((dist2(matrix(loc[j,], ncol=2), loc))[1,],2)
      id = sort(a, index=TRUE)$ix[2:(nei+1)]
      idres = rbind(idres, id)
      distres = rbind(distres, sort(a, index=TRUE)$x[2:(nei+1)])
    }
    pnei.comb = apply(pnei, 1, CCT)
    pnei = cbind(pval,pnei.comb)
    pnei = apply(pnei, 2, function(x) ifelse(x > 1-1e-5, 1-1e-5, x))
    statnei = apply(pnei, 2, pvalue_to_teststat)
    stat = statnei
    
    N = nrow(stat); Nm = dim(basis.tr)[2]
    basis = basis.tr
    #delta = delta.tr0
    
    f00 = null_joint(stat)
    AA = Mirror_Piest_noMT(pval, basis.tr00)
    
    if (nei!=1) {
      pi.est = neigh(AA$piest, nei)
      pi.est = cbind(pi.est[,1], apply(pi.est[,-1], 1, prod))
    } else if (nei==1){
      pi.est = AA$piest
      pi.est =neigh(pi.est, 1)
    }
    
    Sig10est = matrix()
    f10 = dnorm(stat[,1], -1, 0.5)*f00
    f01 = dnorm(stat[,2], -1, 0.5)*f00
    f11 = dnorm(stat[,1], -1, 0.5)*dnorm(stat[,2], -1, 0.5)
    
    # initial values
    pi0.1 = pi.est[,1]; pi0.2 = pi.est[,2]
    ### ESTIMATION METHOD
    fmarg = pi0.1*pi0.2*f00 + pi0.1*(1-pi0.2)*f01 + (1-pi0.1)*pi0.2*f10 + (1-pi0.1)*(1-pi0.2)*f11
    q1 = (pi0.1*pi0.2*f00 + pi0.1*(1-pi0.2)*f01)/fmarg
    q2 = (pi0.1*pi0.2*f00 + (1-pi0.1)*pi0.2*f10)/fmarg
    
    gamma.cur =  rep(1, 2*ncol(basis)) # initial value
    param.cur = c(rep(-1,4), 0.5, 1, 0.5,1, 1, 1, 0.5, 0.5)
    K = 100
    lik = -1e5
    iter = 1
    epsilon = 1
    d = Nm = ncol(basis)
    while(iter <= K & epsilon > 1e-6)
    {
      ######## EM 1: Incorporate spatial variability to estimate \pi(s) through basis ftns ########
      ###--- E-step: probability of being null at s_0 given other ftns ---- ###
      q1 = (pi0.1*pi0.2*f00 + pi0.1*(1-pi0.2)*f01)/fmarg
      q2 = (pi0.1*pi0.2*f00 + (1-pi0.1)*pi0.2*f10)/fmarg
      
      ###--- M-step (Newton method): estimate \pi through gamma --- ###
      Q_g = function(gamma) {
        #gamma1 = gamma[1:d]; gamma2 = gamma[-c(1:d)]
        sum(-q1*(basis %*% gamma[1:d])-q2*(basis %*% gamma[-c(1:d)])-log(1/(exp(basis %*% gamma[(1:d)])+1))-log(1/(exp(basis %*% gamma[-c(1:d)])+1)))
      }
      grr_g = function(gamma) {
        gamma1 = gamma[1:d]; gamma2 = gamma[-c(1:d)]
        A1 = matrix(-q1, N, Nm)*basis
        A2 = matrix(-q2, N, Nm)*basis
        A = cbind(A1, A2)
        C1 = exp(basis %*% gamma1);C2 = exp(basis %*% gamma2)
        C1[C1==Inf] = 1e+30;C2[C2==Inf] = 1e+30
        B1 =  basis*matrix(C1/(C1+1), N, Nm);B2 =  basis*matrix(C2/(C2+1), N, Nm)
        B = cbind(B1, B2)
        return(apply(A+B,2,sum))
      }
      
      fit.optim.g = nlminb(gamma.cur, Q_g, gradient=grr_g)   ### Initial value !!!
      gamma.cur = fit.optim.g$par
      pi0.1 = plogis(basis %*% gamma.cur[1:d])
      pi0.2 = plogis(basis %*% gamma.cur[-c(1:d)])
      
      fmarg = pi0.1*pi0.2*f00 + pi0.1*(1-pi0.2)*f01 + (1-pi0.1)*pi0.2*f10 + (1-pi0.1)*(1-pi0.2)*f11
      q1 = (pi0.1*pi0.2*f00 + pi0.1*(1-pi0.2)*f01)/fmarg
      q2 = (pi0.1*pi0.2*f00 + (1-pi0.1)*pi0.2*f10)/fmarg
      
      ######## EM 2: Estimate joint density through GMM ########
      logBivNorm = function(param,x){
        mu_00 = mu_01 = mu_10 = mu_11 = c(0,0)
        Sig_00 = Sig_01 = Sig_10 = Sig_11 = matrix(1, 2,2)
        mu_01[2] = param[1]; mu_10[1]= param[2]; mu_11[1] = param[3]; mu_11[2] = param[4]
        Sig_01[1,2] = Sig_01[2,1] = param[5]; Sig_01[2,2] = param[6]
        Sig_10[1,2] = Sig_10[2,1] = param[7]; Sig_10[1,1] = param[8]
        Sig_11[1,1] = param[9]; Sig_11[2,2] =param[10]; Sig_11[1,2]= Sig_11[2,1] =param[11]
        Sig_00[1,2] = Sig_00[2,1] = param[12]
        df00 = dmvnorm(x,mu_00, Sig_00)
        df01 = dmvnorm(x,mu_01, Sig_01)
        df10 = dmvnorm(x,mu_10, Sig_10)
        df11 = dmvnorm(x,mu_11, Sig_11)
        -sum(log(pi0.1*pi0.2*df00+pi0.1*(1-pi0.2)*df01+(1-pi0.1)*pi0.2*df10+(1-pi0.1)*(1-pi0.2)*df11))
      }
      logBivNormX <- function(param) {logBivNorm(param,stat)}
      #optim.nlm = nlm(f=logBivNorm, p=c(rep(-1,4), 0.5, 1, 0.5,1, 1, 1, 0.5, 0.5), x=stat)
      optim.op = optim(par=param.cur, fn=logBivNormX)
      param.cur = optim.op$par
      ## Mean estimation : mean of the marginal null is 0
      mu_00 = c(0,0); mu_01=c(0,param.cur[1]); mu_10 = c(param.cur[2],0)
      mu_11 = c(param.cur[3], param.cur[4])
      Sig_00 = Sig_01 = Sig_10 = Sig_11 = matrix(1, 2,2)
      Sig_01[1,2] = Sig_01[2,1] = param.cur[5]; Sig_01[2,2] = param.cur[6]
      Sig_10[1,2] = Sig_10[2,1] = param.cur[7]; Sig_10[1,1] = param.cur[8]
      Sig_11[1,1] = param.cur[9]; Sig_11[2,2] =param.cur[10]; Sig_11[1,2]= Sig_11[2,1] =param.cur[11]
      Sig_00[1,2] = Sig_00[2,1] = param.cur[12]
      
      ### update joint density ###
      f00 = dmvnorm(stat, mu_00, Sig_00)
      f01 = dmvnorm(stat, mu_01, Sig_01)
      f10 = dmvnorm(stat, mu_10, Sig_10)
      f11 = dmvnorm(stat, mu_11, Sig_11)
      
      iter = iter + 1
      fmarg = pi0.1*pi0.2*f00 + pi0.1*(1-pi0.2)*f01 + (1-pi0.1)*pi0.2*f10 + (1-pi0.1)*(1-pi0.2)*f11
      lik.cur = sum(log(fmarg))
      epsilon = abs(lik.cur - lik[length(lik)])/abs(lik[length(lik)])
      lik = c(lik, lik.cur)
      if(iter%%10==0) print(iter)
    }
    
    bic = -2*lik.cur+(9+numb[ii,3])*log(n.site)
    print(ii)
    bicvals = c(bicvals, bic)
  }
  numb_basis_opt = numb[which.min(bicvals),]
  basis_opt = Create_Basis(c(-125,-66),c(25, 50), as.numeric(numb_basis_opt[1])
                           , as.numeric(numb_basis_opt[2]), loc )
  return(list(res=numb_basis_opt, bic = bicvals, basis_opt = basis_opt))
}

OptimalBasis_Nei_BIC = function(pval, nei_vec, n1,n2, fixed=FALSE){
  
  eval.sp = Create_Basis(c(-125,-66),c(25, 50), 3, 3, loc )
  basis.tr00 = eval.sp[,1:4]
  
  ### different numb of basis functions
  if (fixed==FALSE){
  numb = expand.grid(n1:n2, n1:(n2-1), nei_vec)
  numb = cbind(numb, numb[,1]*numb[,2]) 
  numb = data.frame(numb)
  names(numb) = c('n1', 'n2','nei', 'prod')
  numb = numb %>% arrange(prod, n1, n2)
  } else if (fixed==TRUE){
    numb = cbind(rep(n1, length(nei_vec)), rep(n2, length(nei_vec)), nei_vec)
    numb = cbind(numb, numb[,1]*numb[,2]) 
    numb = data.frame(numb)
    names(numb) = c('n1', 'n2','nei', 'prod')
    numb = numb %>% arrange(prod, n1, n2)
  }
  bicvals = c()
  for ( ii in 1:nrow(numb)){
    
    basis.tr = Create_Basis(c(-125,-66),c(25, 50), numb[ii,1], numb[ii,2], loc )
    nei = numb[ii,3]
    idres = distres =c()
    pnei = neigh(pval,nei)
    pnei = apply(pnei, 2, function(x) ifelse(x < 1e-5, 1e-5, ifelse(x > 1-1e-5, 1-1e-5,x)))
    
    for ( j in 1:n.site){
      a = round((dist2(matrix(loc[j,], ncol=2), loc))[1,],2)
      id = sort(a, index=TRUE)$ix[2:(nei+1)]
      idres = rbind(idres, id)
      distres = rbind(distres, sort(a, index=TRUE)$x[2:(nei+1)])
    }
    pnei.comb = apply(pnei, 1, CCT)
    pnei = cbind(pval,pnei.comb)
    pnei = apply(pnei, 2, function(x) ifelse(x > 1-1e-5, 1-1e-5, x))
    statnei = apply(pnei, 2, pvalue_to_teststat)
    stat = statnei
    
    N = nrow(stat); Nm = dim(basis.tr)[2]
    basis = basis.tr
    #delta = delta.tr0
    
    f00 = null_joint(stat)
    AA = Mirror_Piest_noMT(pval, basis.tr00)
    
    if (nei!=1) {
      pi.est = neigh(AA$piest, nei)
      pi.est = cbind(pi.est[,1], apply(pi.est[,-1], 1, prod))
    } else if (nei==1){
      pi.est = AA$piest
      pi.est =neigh(pi.est, 1)
    }
    
    Sig10est = matrix()
    f10 = dnorm(stat[,1], -1, 0.5)*f00
    f01 = dnorm(stat[,2], -1, 0.5)*f00
    f11 = dnorm(stat[,1], -1, 0.5)*dnorm(stat[,2], -1, 0.5)
    
    # initial values
    pi0.1 = pi.est[,1]; pi0.2 = pi.est[,2]
    ### ESTIMATION METHOD
    fmarg = pi0.1*pi0.2*f00 + pi0.1*(1-pi0.2)*f01 + (1-pi0.1)*pi0.2*f10 + (1-pi0.1)*(1-pi0.2)*f11
    q1 = (pi0.1*pi0.2*f00 + pi0.1*(1-pi0.2)*f01)/fmarg
    q2 = (pi0.1*pi0.2*f00 + (1-pi0.1)*pi0.2*f10)/fmarg
    
    gamma.cur =  rep(1, 2*ncol(basis)) # initial value
    param.cur = c(rep(-1,4), 0.5, 1, 0.5,1, 1, 1, 0.5, 0.5)
    K = 100
    lik = -1e5
    iter = 1
    epsilon = 1
    d = Nm = ncol(basis)
    while(iter <= K & epsilon > 1e-6)
    {
      ######## EM 1: Incorporate spatial variability to estimate \pi(s) through basis ftns ########
      ###--- E-step: probability of being null at s_0 given other ftns ---- ###
      q1 = (pi0.1*pi0.2*f00 + pi0.1*(1-pi0.2)*f01)/fmarg
      q2 = (pi0.1*pi0.2*f00 + (1-pi0.1)*pi0.2*f10)/fmarg
      
      ###--- M-step (Newton method): estimate \pi through gamma --- ###
      Q_g = function(gamma) {
        #gamma1 = gamma[1:d]; gamma2 = gamma[-c(1:d)]
        sum(-q1*(basis %*% gamma[1:d])-q2*(basis %*% gamma[-c(1:d)])-log(1/(exp(basis %*% gamma[(1:d)])+1))-log(1/(exp(basis %*% gamma[-c(1:d)])+1)))
      }
      grr_g = function(gamma) {
        gamma1 = gamma[1:d]; gamma2 = gamma[-c(1:d)]
        A1 = matrix(-q1, N, Nm)*basis
        A2 = matrix(-q2, N, Nm)*basis
        A = cbind(A1, A2)
        C1 = exp(basis %*% gamma1);C2 = exp(basis %*% gamma2)
        C1[C1==Inf] = 1e+30;C2[C2==Inf] = 1e+30
        B1 =  basis*matrix(C1/(C1+1), N, Nm);B2 =  basis*matrix(C2/(C2+1), N, Nm)
        B = cbind(B1, B2)
        return(apply(A+B,2,sum))
      }
      
      fit.optim.g = nlminb(gamma.cur, Q_g, gradient=grr_g)   ### Initial value !!!
      gamma.cur = fit.optim.g$par
      pi0.1 = plogis(basis %*% gamma.cur[1:d])
      pi0.2 = plogis(basis %*% gamma.cur[-c(1:d)])
      
      fmarg = pi0.1*pi0.2*f00 + pi0.1*(1-pi0.2)*f01 + (1-pi0.1)*pi0.2*f10 + (1-pi0.1)*(1-pi0.2)*f11
      q1 = (pi0.1*pi0.2*f00 + pi0.1*(1-pi0.2)*f01)/fmarg
      q2 = (pi0.1*pi0.2*f00 + (1-pi0.1)*pi0.2*f10)/fmarg
      
      ######## EM 2: Estimate joint density through GMM ########
      logBivNorm = function(param,x){
        mu_00 = mu_01 = mu_10 = mu_11 = c(0,0)
        Sig_00 = Sig_01 = Sig_10 = Sig_11 = matrix(1, 2,2)
        mu_01[2] = param[1]; mu_10[1]= param[2]; mu_11[1] = param[3]; mu_11[2] = param[4]
        Sig_01[1,2] = Sig_01[2,1] = param[5]; Sig_01[2,2] = param[6]
        Sig_10[1,2] = Sig_10[2,1] = param[7]; Sig_10[1,1] = param[8]
        Sig_11[1,1] = param[9]; Sig_11[2,2] =param[10]; Sig_11[1,2]= Sig_11[2,1] =param[11]
        Sig_00[1,2] = Sig_00[2,1] = param[12]
        df00 = dmvnorm(x,mu_00, Sig_00)
        df01 = dmvnorm(x,mu_01, Sig_01)
        df10 = dmvnorm(x,mu_10, Sig_10)
        df11 = dmvnorm(x,mu_11, Sig_11)
        -sum(log(pi0.1*pi0.2*df00+pi0.1*(1-pi0.2)*df01+(1-pi0.1)*pi0.2*df10+(1-pi0.1)*(1-pi0.2)*df11))
      }
      logBivNormX <- function(param) {logBivNorm(param,stat)}
      #optim.nlm = nlm(f=logBivNorm, p=c(rep(-1,4), 0.5, 1, 0.5,1, 1, 1, 0.5, 0.5), x=stat)
      optim.op = optim(par=param.cur, fn=logBivNormX)
      param.cur = optim.op$par
      ## Mean estimation : mean of the marginal null is 0
      mu_00 = c(0,0); mu_01=c(0,param.cur[1]); mu_10 = c(param.cur[2],0)
      mu_11 = c(param.cur[3], param.cur[4])
      Sig_00 = Sig_01 = Sig_10 = Sig_11 = matrix(1, 2,2)
      Sig_01[1,2] = Sig_01[2,1] = param.cur[5]; Sig_01[2,2] = param.cur[6]
      Sig_10[1,2] = Sig_10[2,1] = param.cur[7]; Sig_10[1,1] = param.cur[8]
      Sig_11[1,1] = param.cur[9]; Sig_11[2,2] =param.cur[10]; Sig_11[1,2]= Sig_11[2,1] =param.cur[11]
      Sig_00[1,2] = Sig_00[2,1] = param.cur[12]
      
      ### update joint density ###
      f00 = dmvnorm(stat, mu_00, Sig_00)
      f01 = dmvnorm(stat, mu_01, Sig_01)
      f10 = dmvnorm(stat, mu_10, Sig_10)
      f11 = dmvnorm(stat, mu_11, Sig_11)
      
      iter = iter + 1
      fmarg = pi0.1*pi0.2*f00 + pi0.1*(1-pi0.2)*f01 + (1-pi0.1)*pi0.2*f10 + (1-pi0.1)*(1-pi0.2)*f11
      lik.cur = sum(log(fmarg))
      epsilon = abs(lik.cur - lik[length(lik)])/abs(lik[length(lik)])
      lik = c(lik, lik.cur)
      if(iter%%10==0) print(iter)
    }
    
    bic = -2*lik.cur+(12+numb[ii,3])*log(n.site)
    print(ii)
    bicvals = c(bicvals, bic)
  }
  numb_basis_opt = numb[which.min(bicvals),]
  basis_opt = Create_Basis(c(-125,-66),c(25, 50), as.numeric(numb_basis_opt[1])
                           , as.numeric(numb_basis_opt[2]), loc )
  numb_nei = numb_basis_opt[3]
  return(list(res=numb_basis_opt, bic = bicvals, basis_opt = basis_opt, nnei = numb_nei))
}

OptimalNei_BIC = function(pval, nei_vec, n1, n2){
  
  eval.sp = Create_Basis(c(-125,-66),c(25, 50), n1, n2, loc )
  basis.tr00 = eval.sp[,1:4]
  basis.tr = eval.sp
  
  ### different numb of nei
  #numb = seq(4, 24, by=4)
  bicvals = c()
  for ( ii in 1:length(nei_vec)){
    
    idres = distres =c()
    nn = nei = nei_vec[ii]
    pnei = neigh(pval,nn)
    pnei = apply(pnei, 2, function(x) ifelse(x < 1e-5, 1e-5, ifelse(x > 1-1e-5, 1-1e-5,x)))
    
    for ( j in 1:n.site){
      a = round((dist2(matrix(loc[j,], ncol=2), loc))[1,],2)
      id = sort(a, index=TRUE)$ix[2:(nn+1)]
      idres = rbind(idres, id)
      distres = rbind(distres, sort(a, index=TRUE)$x[2:(nn+1)])
    }
    pnei.comb = apply(pnei, 1, CCT)
    pnei = cbind(pval,pnei.comb)
    pnei = apply(pnei, 2, function(x) ifelse(x > 1-1e-5, 1-1e-5, x))
    statnei = apply(pnei, 2, pvalue_to_teststat)
    stat = statnei
    
    N = nrow(stat); Nm = dim(basis.tr)[2]
    basis = basis.tr
    #delta = delta.tr0
    
    f00 = null_joint(stat)
    AA = Mirror_Piest_noMT(pval, basis.tr00)
    
    if (nei!=1) {
      pi.est = neigh(AA$piest, nn)
      pi.est = cbind(pi.est[,1], apply(pi.est[,-1], 1, prod))
    } else if (nei==1){
      pi.est = AA$piest
      pi.est =neigh(pi.est, 1)
    }
    
    Sig10est = matrix()
    f10 = dnorm(stat[,1], -1, 0.5)*f00
    f01 = dnorm(stat[,2], -1, 0.5)*f00
    f11 = dnorm(stat[,1], -1, 0.5)*dnorm(stat[,2], -1, 0.5)
    
    # initial values
    pi0.1 = pi.est[,1]; pi0.2 = pi.est[,2]
    ### ESTIMATION METHOD
    fmarg = pi0.1*pi0.2*f00 + pi0.1*(1-pi0.2)*f01 + (1-pi0.1)*pi0.2*f10 + (1-pi0.1)*(1-pi0.2)*f11
    q1 = (pi0.1*pi0.2*f00 + pi0.1*(1-pi0.2)*f01)/fmarg
    q2 = (pi0.1*pi0.2*f00 + (1-pi0.1)*pi0.2*f10)/fmarg
    
    gamma.cur =  rep(1, 2*ncol(basis)) # initial value
    param.cur = c(rep(-1,4), 0.5, 1, 0.5,1, 1, 1, 0.5, 0.5)
    K = 100
    lik = -1e5
    iter = 1
    epsilon = 1
    d = Nm = ncol(basis)
    while(iter <= K & epsilon > 1e-6)
    {
      ######## EM 1: Incorporate spatial variability to estimate \pi(s) through basis ftns ########
      ###--- E-step: probability of being null at s_0 given other ftns ---- ###
      q1 = (pi0.1*pi0.2*f00 + pi0.1*(1-pi0.2)*f01)/fmarg
      q2 = (pi0.1*pi0.2*f00 + (1-pi0.1)*pi0.2*f10)/fmarg
      
      ###--- M-step (Newton method): estimate \pi through gamma --- ###
      Q_g = function(gamma) {
        #gamma1 = gamma[1:d]; gamma2 = gamma[-c(1:d)]
        sum(-q1*(basis %*% gamma[1:d])-q2*(basis %*% gamma[-c(1:d)])-log(1/(exp(basis %*% gamma[(1:d)])+1))-log(1/(exp(basis %*% gamma[-c(1:d)])+1)))
      }
      grr_g = function(gamma) {
        gamma1 = gamma[1:d]; gamma2 = gamma[-c(1:d)]
        A1 = matrix(-q1, N, Nm)*basis
        A2 = matrix(-q2, N, Nm)*basis
        A = cbind(A1, A2)
        C1 = exp(basis %*% gamma1);C2 = exp(basis %*% gamma2)
        C1[C1==Inf] = 1e+30;C2[C2==Inf] = 1e+30
        B1 =  basis*matrix(C1/(C1+1), N, Nm);B2 =  basis*matrix(C2/(C2+1), N, Nm)
        B = cbind(B1, B2)
        return(apply(A+B,2,sum))
      }
      
      fit.optim.g = nlminb(gamma.cur, Q_g, gradient=grr_g)   ### Initial value !!!
      gamma.cur = fit.optim.g$par
      pi0.1 = plogis(basis %*% gamma.cur[1:d])
      pi0.2 = plogis(basis %*% gamma.cur[-c(1:d)])
      
      fmarg = pi0.1*pi0.2*f00 + pi0.1*(1-pi0.2)*f01 + (1-pi0.1)*pi0.2*f10 + (1-pi0.1)*(1-pi0.2)*f11
      q1 = (pi0.1*pi0.2*f00 + pi0.1*(1-pi0.2)*f01)/fmarg
      q2 = (pi0.1*pi0.2*f00 + (1-pi0.1)*pi0.2*f10)/fmarg
      
      ######## EM 2: Estimate joint density through GMM ########
      logBivNorm = function(param,x){
        mu_00 = mu_01 = mu_10 = mu_11 = c(0,0)
        Sig_00 = Sig_01 = Sig_10 = Sig_11 = matrix(1, 2,2)
        mu_01[2] = param[1]; mu_10[1]= param[2]; mu_11[1] = param[3]; mu_11[2] = param[4]
        Sig_01[1,2] = Sig_01[2,1] = param[5]; Sig_01[2,2] = param[6]
        Sig_10[1,2] = Sig_10[2,1] = param[7]; Sig_10[1,1] = param[8]
        Sig_11[1,1] = param[9]; Sig_11[2,2] =param[10]; Sig_11[1,2]= Sig_11[2,1] =param[11]
        Sig_00[1,2] = Sig_00[2,1] = param[12]
        df00 = dmvnorm(x,mu_00, Sig_00)
        df01 = dmvnorm(x,mu_01, Sig_01)
        df10 = dmvnorm(x,mu_10, Sig_10)
        df11 = dmvnorm(x,mu_11, Sig_11)
        -sum(log(pi0.1*pi0.2*df00+pi0.1*(1-pi0.2)*df01+(1-pi0.1)*pi0.2*df10+(1-pi0.1)*(1-pi0.2)*df11))
      }
      logBivNormX <- function(param) {logBivNorm(param,stat)}
      #optim.nlm = nlm(f=logBivNorm, p=c(rep(-1,4), 0.5, 1, 0.5,1, 1, 1, 0.5, 0.5), x=stat)
      optim.op = optim(par=param.cur, fn=logBivNormX)
      param.cur = optim.op$par
      ## Mean estimation : mean of the marginal null is 0
      mu_00 = c(0,0); mu_01=c(0,param.cur[1]); mu_10 = c(param.cur[2],0)
      mu_11 = c(param.cur[3], param.cur[4])
      Sig_00 = Sig_01 = Sig_10 = Sig_11 = matrix(1, 2,2)
      Sig_01[1,2] = Sig_01[2,1] = param.cur[5]; Sig_01[2,2] = param.cur[6]
      Sig_10[1,2] = Sig_10[2,1] = param.cur[7]; Sig_10[1,1] = param.cur[8]
      Sig_11[1,1] = param.cur[9]; Sig_11[2,2] =param.cur[10]; Sig_11[1,2]= Sig_11[2,1] =param.cur[11]
      Sig_00[1,2] = Sig_00[2,1] = param.cur[12]
      
      ### update joint density ###
      f00 = dmvnorm(stat, mu_00, Sig_00)
      f01 = dmvnorm(stat, mu_01, Sig_01)
      f10 = dmvnorm(stat, mu_10, Sig_10)
      f11 = dmvnorm(stat, mu_11, Sig_11)
      
      iter = iter + 1
      fmarg = pi0.1*pi0.2*f00 + pi0.1*(1-pi0.2)*f01 + (1-pi0.1)*pi0.2*f10 + (1-pi0.1)*(1-pi0.2)*f11
      lik.cur = sum(log(fmarg))
      epsilon = abs(lik.cur - lik[length(lik)])/abs(lik[length(lik)])
      lik = c(lik, lik.cur)
      if(iter%%10==0) print(iter)
    }
    
    bic = -2*lik.cur+(12+n1*n2)*log(n.site)
    print(ii)
    bicvals = c(bicvals, bic)
  }
  numb_nei_opt = nei_vec[which.min(bicvals)]
  
  return(list(res=numb_nei_opt, bic = bicvals))
}

#### function
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

    Sig_00 = matrix((1/N00[i])*apply(diag(w00[,i]) %*% t(apply(stat[,c(1,i+1)], 1, function(x) (x) %*% t(x))), 2,sum), ncol=2)
    Sig_01 = matrix((1/N01[i])*apply(diag(w01[,i]) %*% t(apply(stat[,c(1,i+1)]-mu_01tmp, 1, function(x) (x) %*% t(x))), 2,sum), ncol=2)
    Sig_10 = matrix((1/N10[i])*apply(diag(w10[,i]) %*% t(apply(stat[,c(1,i+1)]-mu_10tmp, 1, function(x) (x) %*% t(x))), 2,sum),ncol=2)
    Sig_11 = matrix((1/N11[i])*apply(diag(w11[,i]) %*% t(apply(stat[,c(1,i+1)]-mu_11tmp, 1, function(x) (x) %*% t(x))), 2,sum),ncol=2)
    Sig_00[1,1] =Sig_00[2,2] = 1 ;
    Sig_01[1,1] = 1; Sig_10[2,2]=1
    
    f00_res = cbind(f00_res, dmvnorm(stat[,c(1,i+1)], c(0,0), Sig_00))
    f01_res = cbind(f01_res, dmvnorm(stat[,c(1,i+1)], mu_01tmp, Sig_01))
    f10_res = cbind(f10_res, dmvnorm(stat[,c(1,i+1)], mu_10tmp, Sig_10))
    f11_res = cbind(f11_res, dmvnorm(stat[,c(1,i+1)], mu_11tmp, Sig_11))
  }
  return(list(f00 = f00_res, f01 = f01_res, f10 = f10_res, f11 = f11_res))
  
}

######--- BiCLfdr_Cauhcy -----#######
BiModel_Cauchy_Test0 = function(pval, nei, delta.tr0,basis.tr00, basis.tr, fdr.c=0.1, index=FALSE, den=FALSE) {
  # FDR control using nonparametric estimate for f1
  idres = distres =c()
  pnei = neigh(pval,nei)
  pnei = apply(pnei, 2, function(x) ifelse(x < 1e-5, 1e-5, ifelse(x > 1-1e-5, 1-1e-5,x)))
  
  for ( j in 1:n.site){
    a = round((dist2(matrix(loc[j,], ncol=2), loc))[1,],2)
    id = sort(a, index=TRUE)$ix[2:(nei+1)]
    idres = rbind(idres, id)
    distres = rbind(distres, sort(a, index=TRUE)$x[2:(nei+1)])
  }
  pnei.comb = apply(pnei, 1, CCT)
  pnei = cbind(pval,pnei.comb)
  pnei = apply(pnei, 2, function(x) ifelse(x > 1-1e-5, 1-1e-5, x))
  
  statnei = apply(pnei, 2, pvalue_to_teststat)
  stat = statnei
  
  N = nrow(stat); Nm = dim(basis.tr)[2]
  basis = basis.tr
  delta = delta.tr0
  
  f00 = null_joint(stat)
  AA = Mirror_Piest(pval, basis.tr00, delta.tr0)
  
  f10 = dnorm(stat[,1], -1, 0.5)*f00
  f01 = dnorm(stat[,2], -1, 0.5)*f00
  f11 = dnorm(stat[,1], -1, 0.5)*dnorm(stat[,2], -1, 0.5)
  
  # initial values
  pi0.1 = AA$piest;
  if (nei!=1) {
    pitmp = neigh(pi0.1, nei)
    pi0.2 = apply(pitmp[,-1], 1, prod)
  } else if (nei==1){
    pitmp = neigh(pi0.1, 1)
    pi0.2 = pitmp[,2]
  }
  ### ESTIMATION METHOD
  fmarg = pi0.1*pi0.2*f00 + pi0.1*(1-pi0.2)*f01 + (1-pi0.1)*pi0.2*f10 + (1-pi0.1)*(1-pi0.2)*f11
  q1 = (pi0.1*pi0.2*f00 + pi0.1*(1-pi0.2)*f01)/fmarg
  q2 = (pi0.1*pi0.2*f00 + (1-pi0.1)*pi0.2*f10)/fmarg
  
  gamma.cur =  rep(1, ncol(basis)) # initial value
  K = 100
  lik = -1e5
  iter = 1
  epsilon = 1
  d = Nm = ncol(basis)
  while(iter <= K & epsilon > 1e-6)
  {
    # E-step:
    q1 = (pi0.1*pi0.2*f00 + pi0.1*(1-pi0.2)*f01)/fmarg
    q2 = (pi0.1*pi0.2*f00 + (1-pi0.1)*pi0.2*f10)/fmarg
    
    # M-step (Newton method):
    Q_g = function(gamma) {
      sum(-q1*(basis %*% gamma)-log(1/(exp(basis %*% gamma)+1)))
      #sum(-q1*(basis %*% beta.cur)-q2*(basis %*% gamma)-log(1/((exp(basis %*% beta.cur)+1))*(exp(basis %*% gamma)+1)))
    }
    grr_g = function(gamma) {
      A = matrix(-q1, N, Nm)*basis
      C = exp(basis %*% gamma)
      C[C==Inf] = 1e+30#;C2[C2==Inf] = 1e+30
      B =  basis*matrix(C/(C+1), N, Nm)
      return(apply(A+B,2,sum))
    }
    
    fit.optim.g = nlminb(gamma.cur, Q_g, gradient=grr_g)   ### Initial value !!!
    gamma.cur = fit.optim.g$par
    pi0.1 = plogis(basis %*% gamma.cur)
    
    if (nei!=1) {
      pitmp = neigh(pi0.1, nei)
      pi0.2 = apply(pitmp[,-1], 1, prod)
      #pi0.2 = apply(pitmp[,-1], 1, median)
    } else if (nei==1){
      pitmp = neigh(pi0.1, 1)
      pi0.2 = pitmp[,2]
    }
    
    #fmarg = pi0.1*pi0.2*f00 + pi0.1*(1-pi0.2)*f01 + (1-pi0.1)*pi0.2*f10 + (1-pi0.1)*(1-pi0.2)*f11
    w00 = as.numeric(q1*q2); w11 = as.numeric((1-q1)*(1-q2))
    w10 = as.numeric((1-q1)*q2); w01 = as.numeric(q1*(1-q2))
    
    #w00 = as.numeric(pi0.1*pi0.2*f00/fmarg); w11 = as.numeric((1-pi0.1)*(1-pi0.2)*f11/fmarg)
    #w10 = as.numeric((1-pi0.1)*pi0.2*f10/fmarg); w01 = as.numeric(pi0.1*(1-pi0.2)*f01/fmarg)
    N00 = sum(w00);N10 = sum(w10);N01 = sum(w01);N11 = sum(w11)
    mu_00 = (1/N00)*apply(w00*stat, 2, sum); mu_00 = c(0,0)
    mu_11 = (1/N11)*apply(w11*stat, 2, sum)
    mu_01 = (1/N01)*apply(w01*stat, 2, sum);
    mu_01 = c(0,mu_01[2])
    mu_10 = (1/N10)*apply(w10*stat, 2, sum);
    mu_10 = c(mu_10[1],0)
    
    Sig_00 = matrix((1/N00)*apply(diag(w00) %*% t(apply(stat-mu_00, 1, function(x) (x) %*% t(x))), 2,sum), ncol=2)
    #Sig_00 = matrix(0,ncol=2, nrow=2)
    #Sig_00[1,1] = Sig_00[2,2] = 1#;
    #Sig_00[1,2] = Sig_00[2,1] = 0
    Sig_01 = matrix((1/N01)*apply(diag(w01) %*% t(apply(stat-mu_01, 1, function(x) (x) %*% t(x))), 2,sum), ncol=2)
    Sig_10 = matrix((1/N10)*apply(diag(w10) %*% t(apply(stat-mu_10, 1, function(x) (x) %*% t(x))), 2,sum),ncol=2)
    Sig_11 = matrix((1/N11)*apply(diag(w11) %*% t(apply(stat-mu_11, 1, function(x) (x) %*% t(x))), 2,sum),ncol=2)
    Sig_00[1,1] =Sig_00[2,2] = 1 ;
    Sig_01[1,1] = 1; Sig_10[2,2]=1
    
    f00 = dmvnorm(stat, mu_00, Sig_00)
    f01 = dmvnorm(stat, mu_01, Sig_01)
    f10 = dmvnorm(stat, mu_10, Sig_10)
    f11 = dmvnorm(stat, mu_11, Sig_11)
    
    iter = iter + 1
    # pi0.2 = pi0.1
    fmarg = pi0.1*pi0.2*f00 + pi0.1*(1-pi0.2)*f01 + (1-pi0.1)*pi0.2*f10 + (1-pi0.1)*(1-pi0.2)*f11
    lik.cur = sum(log(fmarg))
    epsilon = abs(lik.cur - lik[length(lik)])/abs(lik[length(lik)])
    lik = c(lik, lik.cur)
    #
    if(iter%%200==0) print(iter)
  }
  CLFDR =  (pi0.1*pi0.2*f00 + pi0.1*(1-pi0.2)*f01)/fmarg
  bic =  -2*lik.cur+(ncol(basis))*log(N)
  # Local approach
  k <- sum(cumsum(sort(CLFDR))/(1:(N)) < fdr.c )               # Nominal level 10%
  ind1 <- sort(CLFDR, index.return=TRUE)$ix[1:k]
  FDP.1 <- mean(delta[ind1] == 0)             # FDP (should be close to the nominal level)
  Power.1 <- sum(delta[ind1]) / (sum(delta==1))            # Power (the larger the better)
  
  if (index==FALSE & den==FALSE){
    Res_B = c(AA$MirrRes,AA$JSRes, AA$BHRes, FDP.1, Power.1)
    return(Res_B)
  }else if (index==TRUE & den==FALSE){
    Res_B = c(AA$MirrRes,AA$JSRes, AA$BHRes, FDP.1, Power.1)
    return(list(Res = Res_B, RejM=AA$Rej, RejBi = ind1, BLFDR = CLFDR, BIC = bic))
  } else if (index==FALSE & den==TRUE){
    Res_B = c(AA$MirrRes,AA$JSRes, AA$BHRes, FDP.1, Power.1)
    Den = cbind(f00, f01, f10, f11)
    return(list(Res=Res_B, density=Den))
  }
}

BiModel_Cauchy_Real = function(pval, nei,  basis.tr00, basis.tr, fdr.c=0.1, index=F, den=F) {
  
  pnei = neigh(pval,nei)
  pnei = apply(pnei, 2, function(x) ifelse(x < 1e-5, 1e-5, ifelse(x > 1-1e-5, 1-1e-5,x)))
  
  pnei.comb = apply(pnei, 1, CCT)
  pnei = cbind(pval,pnei.comb)
  pnei = apply(pnei, 2, function(x) ifelse(x > 1-1e-5, 1-1e-5, x))
  
  statnei = apply(pnei, 2, pvalue_to_teststat)
  stat = statnei
  
  N = nrow(stat); Nm = dim(basis.tr)[2]#; Nm2 = dim(basis.tr2)[2]
  basis = basis.tr#; basis2 = basis.tr2
  delta = delta.tr0
  
  f00 = null_joint(stat)
  AA = Mirror_Piest(pval, basis.tr00, delta.tr0)
  
  f10 = dnorm(stat[,1], -1, 0.5)*f00
  f01 = dnorm(stat[,2], -1, 0.5)*f00
  f11 = dnorm(stat[,1], -1, 0.5)*dnorm(stat[,2], -1, 0.5)
  
  # initial values
  pi0.1 = AA$piest;
  if (nei!=1) {
    pitmp = neigh(pi0.1, nei)
    pi0.2 = apply(pitmp[,-1], 1, prod)
  } else if (nei==1){
    pitmp = neigh(pi0.1, 1)
    pi0.2 = pitmp[,2]
  }
  ### ESTIMATION METHOD
  fmarg = pi0.1*pi0.2*f00 + pi0.1*(1-pi0.2)*f01 + (1-pi0.1)*pi0.2*f10 + (1-pi0.1)*(1-pi0.2)*f11
  q1 = (pi0.1*pi0.2*f00 + pi0.1*(1-pi0.2)*f01)/fmarg
  q2 = (pi0.1*pi0.2*f00 + (1-pi0.1)*pi0.2*f10)/fmarg
  
  gamma.cur =  rep(1, ncol(basis)) # initial value
  K = 500
  lik = -1e5
  iter = 1
  epsilon = 1
  d = Nm = ncol(basis)
  while(iter <= K & epsilon > 1e-6)
  {
    # E-step:
    q1 = (pi0.1*pi0.2*f00 + pi0.1*(1-pi0.2)*f01)/fmarg
    q2 = (pi0.1*pi0.2*f00 + (1-pi0.1)*pi0.2*f10)/fmarg
    
    # M-step (Newton method):
    Q_g = function(gamma) {
      sum(-q1*(basis %*% gamma)-log(1/(exp(basis %*% gamma)+1)))
          }
    grr_g = function(gamma) {
      A = matrix(-q1, N, Nm)*basis
      C = exp(basis %*% gamma)
      C[C==Inf] = 1e+30
      B =  basis*matrix(C/(C+1), N, Nm)
      return(apply(A+B,2,sum))
    }
    
    fit.optim.g = nlminb(gamma.cur, Q_g, gradient=grr_g)   ### Initial value !!!
    gamma.cur = fit.optim.g$par
    pi0.1 = plogis(basis %*% gamma.cur)
    
    if (nei!=1) {
      pitmp = neigh(pi0.1, nei)
      pi0.2 = apply(pitmp[,-1], 1, prod)
    } else if (nei==1){
      pitmp = neigh(pi0.1, 1)
      pi0.2 = pitmp[,2]
    }
    
    fmarg = pi0.1*pi0.2*f00 + pi0.1*(1-pi0.2)*f01 + (1-pi0.1)*pi0.2*f10 + (1-pi0.1)*(1-pi0.2)*f11
    
    w00 = as.numeric(pi0.1*pi0.2*f00/fmarg); w11 = as.numeric((1-pi0.1)*(1-pi0.2)*f11/fmarg)
    w10 = as.numeric((1-pi0.1)*pi0.2*f10/fmarg); w01 = as.numeric(pi0.1*(1-pi0.2)*f01/fmarg)
    N00 = sum(w00);N10 = sum(w10);N01 = sum(w01);N11 = sum(w11)
    mu_00 = (1/N00)*apply(w00*stat, 2, sum); mu_00 = c(0,0)
    mu_11 = (1/N11)*apply(w11*stat, 2, sum)
    mu_01 = (1/N01)*apply(w01*stat, 2, sum);
    mu_01 = c(0,mu_01[2])
    mu_10 = (1/N10)*apply(w10*stat, 2, sum);
    mu_10 = c(mu_10[1],0)
    
    Sig_00 = matrix((1/N00)*apply(diag(w00) %*% t(apply(stat-mu_00, 1, function(x) (x) %*% t(x))), 2,sum), ncol=2)
    Sig_01 = matrix((1/N01)*apply(diag(w01) %*% t(apply(stat-mu_01, 1, function(x) (x) %*% t(x))), 2,sum), ncol=2)
    Sig_10 = matrix((1/N10)*apply(diag(w10) %*% t(apply(stat-mu_10, 1, function(x) (x) %*% t(x))), 2,sum),ncol=2)
    Sig_11 = matrix((1/N11)*apply(diag(w11) %*% t(apply(stat-mu_11, 1, function(x) (x) %*% t(x))), 2,sum),ncol=2)
    Sig_00[1,1] =Sig_00[2,2] = 1 ;
    Sig_01[1,1] = 1; Sig_10[2,2]=1
    
    f00 = dmvnorm(stat, mu_00, Sig_00)
    f01 = dmvnorm(stat, mu_01, Sig_01)
    f10 = dmvnorm(stat, mu_10, Sig_10)
    f11 = dmvnorm(stat, mu_11, Sig_11)
    
    iter = iter + 1
    # pi0.2 = pi0.1
    fmarg = pi0.1*pi0.2*f00 + pi0.1*(1-pi0.2)*f01 + (1-pi0.1)*pi0.2*f10 + (1-pi0.1)*(1-pi0.2)*f11
    lik.cur = sum(log(fmarg))
    epsilon = abs(lik.cur - lik[length(lik)])/abs(lik[length(lik)])
    lik = c(lik, lik.cur)
    if(iter%%10==0) print(iter)
  }
  CLFDR =  (pi0.1*pi0.2*f00 + pi0.1*(1-pi0.2)*f01)/fmarg
  
  # Local approach
  k <- sum(cumsum(sort(CLFDR))/(1:(N)) < fdr.c )               # Nominal level 10%
  ind1 <- sort(CLFDR, index.return=TRUE)$ix[1:k]
  return(list(RejM=AA$Rej, RejBi = ind1))
  
}
