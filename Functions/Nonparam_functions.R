## Mirror and local FDR procedure
NonparamEst0 =function(p_value, delta.tr0, delta.tr1, basis.tr0, fdr.c=0.1
                       , index=FALSE, JS=TRUE, local =FALSE) { 
  # FDR control using nonparametric estimate for f1
  N = length(p_value); Nm = dim(basis.tr0)[2]
  pval = p_value
  out = sort(pval, index.return=TRUE)
  pval = out$x
  ind = out$ix
  basis = basis.tr0[ind,]
  delta = delta.tr0[ind]		
  df = c(pval[1],diff(pval))       # z(i)-z(i-1)
  f1.int = f1 = dexp(pval,100)     # initial yi=f1(z)
  pi0  = runif(N,0,1)		  # initial pi(s)
  f0 = rep(1,N)		        # known   f0(zi)
  
  ### ESTIMATION METHOD
  q = pi0*f0/(pi0*f0+(1-pi0)*f1);d = ncol(basis)
  beta.int = beta.cur = rep(1, d)  # initial value
  K = 100
  lik = -1e5
  iter = 1
  epsilon = 1
  while(iter <= K & epsilon > 1e-9)
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
  
  # Local approach
  #fdr.c = 0.05
  lfdr <- pi0*f0/(pi0*f0+(1-pi0)*f1)
  k <- sum(cumsum(sort(lfdr))/(1:(N)) < fdr.c )               # Nominal level 10%
  ind1 <- sort(lfdr, index.return=TRUE)$ix[1:k]
  sum(delta == 0)/length(delta)
  FDP.1 <- mean(delta[ind1] == 0)             # FDP (should be close to the nominal level)
  Power.1 <- sum(delta[ind1]) / (sum(delta.tr1==1))            # Power (the larger the better) 
  
  # Shape-adaptive approach
  f.step = stepfun(pval, c(f1[1],f1), right = FALSE) 
  fp = f.step(pval)
  f1p = f.step(1-pval)
  nu = (1-pi0)/pi0*f1p
  de = (1-pi0)/pi0*fp
  TT = sort(c(nu,de))
  i = 1
  hit = FALSE
  while(i <= length(TT) & hit == FALSE)
  {
    if ( sum(nu > TT[i])/max(1,sum(de > TT[i])) < fdr.c & sum(de > TT[i])>0) hit=TRUE
    i = i+1
  }
  if(hit == TRUE) {
    ind2 = (1:N)[de > TT[i]]
    FDP.2 = mean( delta[(1:N)[de > TT[i]]] ==0 )
    Power.2 = sum( delta[(1:N)[de > TT[i]]] ==1 )/(sum(delta.tr1==1))  
  } else
  {
    FDP.2 = 0;
    Power.2 = 0;
  }  
  
  # John Storey approach
  if (JS==TRUE) p.adj = qvalue(pval, fdr.level=fdr.c)
  else if(JS==FALSE) p.adj = qvalue(pval, fdr.level=fdr.c, pi0=1)
  indQ = which(p.adj$significant == TRUE) 	
  FDP.J = mean(delta[indQ] == 0)                   
  Power.J = sum(delta[indQ]) / (sum(delta.tr1 == 1))
  # BH approach
  
  p1.adj = p.adjust(pval, method='BH', n =length(pval))
  indBH = which(p1.adj <fdr.c)
  
  FDP.B <- mean(delta[indBH] == 0)                    # FDP (should be close to the nominal level)
  Power.B <- sum(delta[indBH]) / (sum(delta.tr1==1));
  
  if(index==FALSE & local==FALSE) {
    
    return(list(FDP.L = FDP.1, FDP.A = FDP.2, Power.L = Power.1, Power.A = Power.2
                , FDP.J = FDP.J, FDP.B = FDP.B, Power.B = Power.B, Power.J = Power.J))
  } else if (index == TRUE & local==FALSE){
    return(list(FDP.L = FDP.1, FDP.A = FDP.2, Power.L = Power.1, Power.A = Power.2
                , FDP.J = FDP.J, FDP.B = FDP.B, Power.B = Power.B, Power.J = Power.J
                , Rej.L = ind1, Rej.M =ind2 , Rej.B = indBH, Rej.J = indQ, LFDR = lfdr, Idx = ind))
  } else if (index == TRUE & local==TRUE){
    return(list(FDP.L = FDP.1, FDP.A = FDP.2, Power.L = Power.1, Power.A = Power.2
                , FDP.J = FDP.J, FDP.B = FDP.B, Power.B = Power.B, Power.J = Power.J
                , Rej.L = ind1, Rej.B = indBH, Rej.J = indQ, LFDR = lfdr, Idx = ind))
  }
  #  , Rej.L = ind1, Rej.M =ind2 , Rej.B = indBH, Rej.J = indQ, LFDR = lfdr, Idx = ind))
}

NonparamEst0_cov =function(p_value, delta.tr0, delta.tr1, basis.tr0, fdr.c=0.1
                       , index=FALSE, JS=TRUE, local =FALSE) { 
  
  #median pvalue
  nn <- length(p_value)
  pstar = pmin  =  pmax = numeric(length = nn)
  
  for( ii in 1:nn ) {
    # neighborhood
    a = (dist2(matrix(loc0[ii,], ncol=2), loc0))[1,]
    short.dist = which(round(a,1) <= wi)
    
    #print(length(short.dist))
    p0 = ifelse(p_value[short.dist] < 1e-5, 1e-5, p_value[short.dist])
    pstar[ii] <- -2*sum(log(p0))
    pmin[ii] = median(p0)#; pmax[ii] = max(p0)
  }
  pstar = cbind(1, pstar, loc0)
  
  # FDR control using nonparametric estimate for f1
  N = length(p_value);
  pval = p_value
  out = sort(pval, index.return=TRUE)
  pval = out$x
  ind = out$ix
  pstar0 = pstar[ind,]; Nm = dim(pstar0)[2]
  basis = basis.tr0[ind,]
  delta = delta.tr0[ind]		
  df = c(pval[1],diff(pval))       # z(i)-z(i-1)
  f1.int = f1 = dexp(pval,100)     # initial yi=f1(z)
  pi0  = runif(N,0,1)		  # initial pi(s)
  f0 = rep(1,N)		        # known   f0(zi)
  
  ### ESTIMATION METHOD
  q = pi0*f0/(pi0*f0+(1-pi0)*f1);d = ncol(pstar0)
  beta.int = beta.cur = rep(1, d)  # initial value
  K = 100
  lik = -1e5
  iter = 1
  epsilon = 1
  while(iter <= K & epsilon > 1e-9)
  {
    # E-step: 
    q = pi0*f0/(pi0*f0+(1-pi0)*f1)
    
    # M-step (Newton method):
    Q = function(ps) {
      sum(-q*(pstar0 %*% ps)-log(1/(exp(pstar0 %*% ps)+1)))
    }
    grr = function(ps) {
      A = matrix(-q, N, Nm)*pstar0 
      C = exp(pstar0 %*% ps)
      C[C==Inf] = 1e+30
      B = pstar0*matrix(C/(C+1), N, Nm)
      return(apply(A+B,2,sum))
    }
    fit.optim = nlminb(beta.int, Q, gradient=grr)   ### Initial value !!!
    beta.cur = fit.optim$par
    pi0 = plogis(pstar0 %*% beta.cur)
    pi0 = ifelse(pi0<= 1e-5, 1e-5, pi0); pi0 = ifelse(pi0 >= 1-1e-5, 1-1e-5, pi0)
    
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
  
  # Local approach
  #fdr.c = 0.05
  lfdr <- pi0*f0/(pi0*f0+(1-pi0)*f1)
  k <- sum(cumsum(sort(lfdr))/(1:(N)) < fdr.c )               # Nominal level 10%
  ind1 <- sort(lfdr, index.return=TRUE)$ix[1:k]
  sum(delta == 0)/length(delta)
  FDP.1 <- mean(delta[ind1] == 0)             # FDP (should be close to the nominal level)
  Power.1 <- sum(delta[ind1]) / (sum(delta.tr1==1))            # Power (the larger the better) 
  
  # Shape-adaptive approach
  f.step = stepfun(pval, c(f1[1],f1), right = FALSE) 
  fp = f.step(pval)
  f1p = f.step(1-pval)
  nu = (1-pi0)/pi0*f1p
  de = (1-pi0)/pi0*fp
  TT = sort(c(nu,de))
  i = 1
  hit = FALSE
  while(i <= length(TT) & hit == FALSE)
  {
    if ( sum(nu > TT[i])/max(1,sum(de > TT[i])) < fdr.c & sum(de > TT[i])>0) hit=TRUE
    i = i+1
  }
  if(hit == TRUE) {
    ind2 = (1:N)[de > TT[i]]
    FDP.2 = mean( delta[(1:N)[de > TT[i]]] ==0 )
    Power.2 = sum( delta[(1:N)[de > TT[i]]] ==1 )/(sum(delta.tr1==1))  
  } else
  {
    FDP.2 = 0;
    Power.2 = 0;
  }  
  
  # John Storey approach
  if (JS==TRUE) p.adj = qvalue(pval, fdr.level=fdr.c)
  else if(JS==FALSE) p.adj = qvalue(pval, fdr.level=fdr.c, pi0=1)
  indQ = which(p.adj$significant == TRUE) 	
  FDP.J = mean(delta[indQ] == 0)                   
  Power.J = sum(delta[indQ]) / (sum(delta.tr1 == 1))
  # BH approach
  
  p1.adj = p.adjust(pval, method='BH', n =length(pval))
  indBH = which(p1.adj <fdr.c)
  
  FDP.B <- mean(delta[indBH] == 0)                    # FDP (should be close to the nominal level)
  Power.B <- sum(delta[indBH]) / (sum(delta.tr1==1));
  
  if(index==FALSE & local==FALSE) {
    
    return(list(FDP.L = FDP.1, FDP.A = FDP.2, Power.L = Power.1, Power.A = Power.2
                , FDP.J = FDP.J, FDP.B = FDP.B, Power.B = Power.B, Power.J = Power.J))
  } else if (index == TRUE & local==FALSE){
    return(list(FDP.L = FDP.1, FDP.A = FDP.2, Power.L = Power.1, Power.A = Power.2
                , FDP.J = FDP.J, FDP.B = FDP.B, Power.B = Power.B, Power.J = Power.J
                , Rej.L = ind1, Rej.M =ind2 , Rej.B = indBH, Rej.J = indQ, LFDR = lfdr, Idx = ind))
  } else if (index == TRUE & local==TRUE){
    return(list(FDP.L = FDP.1, FDP.A = FDP.2, Power.L = Power.1, Power.A = Power.2
                , FDP.J = FDP.J, FDP.B = FDP.B, Power.B = Power.B, Power.J = Power.J
                , Rej.L = ind1, Rej.B = indBH, Rej.J = indQ, LFDR = lfdr, Idx = ind))
  }
  #  , Rej.L = ind1, Rej.M =ind2 , Rej.B = indBH, Rej.J = indQ, LFDR = lfdr, Idx = ind))
}

## covariates - loc, pvalA
NonparamEst0_cov2=function(p_value, delta.tr0, delta.tr1, basis.tr0, fdr.c=0.1
                           , index=FALSE, JS=TRUE, local =FALSE, method='nonparam') { 
  
  #median pvalue
  nn <- length(p_value)
  pstar = pmin  =  pmax = numeric(length = nn)
  
  pstar = cbind(1, p_value_neigh, sqrt(loc0), sqrt(loc0[,1]*loc0[,2]), loc0)
  # FDR control using nonparametric estimate for f1
  N = length(p_value);
  pval = p_value
  out = sort(pval, index.return=TRUE)
  pval = out$x
  ind = out$ix
  pstar0 = pstar[ind,]; Nm = dim(pstar0)[2]
  basis = basis.tr0[ind,]
  delta = delta.tr0[ind]		
  df = c(pval[1],diff(pval))       # z(i)-z(i-1)
  f1.int = f1 = dexp(pval,100)     # initial yi=f1(z)
  pi0  = (1-pis.hat)[ind]	  # initial pi(s)
  f0 = rep(1,N)		        # known   f0(zi)
  pval = ifelse(pval < 1e-5, 1e-5, pval); pval = ifelse(pval > 1-1e-5, 1-1e-5, pval)
  ### ESTIMATION METHOD
  theta0 = seq(2, 100, length =N); 	alpha.cur = matrix(lm((log(theta0))~ -1+pstar0)$coef)
  q = pi0*f0/(pi0*f0+(1-pi0)*f1);d = ncol(pstar0)
  beta.int = beta.cur = rep(1, d)  # initial value
  K = 100
  lik = -1e5
  iter = 1
  epsilon = 1
  while(iter <= K & epsilon > 1e-9)
  {
    # E-step: 
    q = pi0*f0/(pi0*f0+(1-pi0)*f1)
    
    # M-step (Newton method):
    Q = function(ps) {
      sum(-q*(pstar0 %*% ps)-log(1/(exp(pstar0 %*% ps)+1)))
    }
    grr = function(ps) {
      A = matrix(-q, N, Nm)*pstar0 
      C = exp(pstar0 %*% ps)
      C[C==Inf] = 1e+30
      B = pstar0*matrix(C/(C+1), N, Nm)
      return(apply(A+B,2,sum))
    }
    fit.optim = nlminb(beta.int, Q, gradient=grr)   ### Initial value !!!
    beta.cur = fit.optim$par
    pi0 = plogis(pstar0 %*% beta.cur)
    pi0 = ifelse(pi0<= 1e-5, 1e-5, pi0); pi0 = ifelse(pi0 >= 1-1e-5, 1-1e-5, pi0)
    
    if (method=='param'){
      # F1 parametric
      # f1 update
      Q.f1 = function(alpha){  
        sum(-(1-q)*(pstar0%*%alpha+(exp(pstar0%*%alpha)-1)*log(1-pval)))
      }
      fit.optim.theta = optim(alpha.cur, Q.f1, method ='Nelder-Mead')
      alpha.cur = fit.optim.theta$par;
      theta0    = ifelse(exp(pstar0%*%alpha.cur)<1, 1, exp(pstar0%*%alpha.cur))
      
      # Estimating f1
      f1 = dbeta(pval, 1, theta0)
    } else if (method=='nonparam'){
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
    }
    iter = iter + 1
    lik.cur = sum(log(pi0*f0+(1-pi0)*f1))
    epsilon = (lik.cur - lik[length(lik)])/abs(lik[length(lik)])
    lik = c(lik, lik.cur)
    if(iter%%200==0) print(iter)
  }
  
  # LAW approach 
  law.dd.res<-law.func(pvs=pval, pi0, 0.1)
  law.dd.de<-law.dd.res$de
  ind.law = which(law.dd.de==1)
  law.fdr0 = sum((1-delta)*law.dd.de)/max(sum(law.dd.de), 1)
  law.power0 = sum(delta*law.dd.de)/sum(delta)
  
  # Local approach
  #fdr.c = 0.05
  lfdr <- pi0*f0/(pi0*f0+(1-pi0)*f1)
  k <- sum(cumsum(sort(lfdr))/(1:(N)) < fdr.c )               # Nominal level 10%
  ind1 <- sort(lfdr, index.return=TRUE)$ix[1:k]
  sum(delta == 0)/length(delta)
  FDP.1 <- mean(delta[ind1] == 0)             # FDP (should be close to the nominal level)
  Power.1 <- sum(delta[ind1]) / (sum(delta.tr1==1))            # Power (the larger the better) 
  
  # Shape-adaptive approach
  f.step = stepfun(pval, c(f1[1],f1), right = FALSE) 
  fp = f.step(pval)
  f1p = f.step(1-pval)
  nu = (1-pi0)/pi0*f1p
  de = (1-pi0)/pi0*fp
  TT = sort(c(nu,de))
  i = 1
  hit = FALSE
  while(i <= length(TT) & hit == FALSE)
  {
    if ( sum(nu > TT[i])/max(1,sum(de > TT[i])) < fdr.c & sum(de > TT[i])>0) hit=TRUE
    i = i+1
  }
  if(hit == TRUE) {
    ind2 = (1:N)[de > TT[i]]
    FDP.2 = mean( delta[(1:N)[de > TT[i]]] ==0 )
    Power.2 = sum( delta[(1:N)[de > TT[i]]] ==1 )/(sum(delta.tr1==1))  
  } else
  {
    FDP.2 = 0;
    Power.2 = 0;
  }  
  
  # # John Storey approach
  # if (JS==TRUE) p.adj = qvalue(pval, fdr.level=fdr.c)
  # else if(JS==FALSE) p.adj = qvalue(pval, fdr.level=fdr.c, pi0=1)
  # indQ = which(p.adj$significant == TRUE) 	
  # FDP.J = mean(delta[indQ] == 0)                   
  # Power.J = sum(delta[indQ]) / (sum(delta.tr1 == 1))
  # # BH approach
  # 
  # p1.adj = p.adjust(pval, method='BH', n =length(pval))
  # indBH = which(p1.adj <fdr.c)
  # 
  # FDP.B <- mean(delta[indBH] == 0)                    # FDP (should be close to the nominal level)
  # Power.B <- sum(delta[indBH]) / (sum(delta.tr1==1));
  
  if(index==FALSE & local==FALSE) {
    
    return(list(FDP.L = FDP.1, FDP.A = FDP.2, Power.L = Power.1, Power.A = Power.2
                , FDP.LA = law.fdr0, Power.LA = law.power0))
  } else if (index == TRUE & local==FALSE){
    return(list(FDP.L = FDP.1, FDP.A = FDP.2, Power.L = Power.1, Power.A = Power.2
                , FDP.LA = law.fdr0, Power.LA = law.power0
                , Rej.L = ind1, Rej.M =ind2 , Rej.LA =ind.law , LFDR = lfdr, Idx = ind))
  } else if (index == TRUE & local==TRUE){
    return(list(FDP.L = FDP.1, FDP.A = FDP.2, Power.L = Power.1, Power.A = Power.2
                , FDP.LA = law.fdr0, Power.LA = law.power0
                , Rej.L = ind1, Rej.M =ind2 , Rej.LA =ind.law , LFDR = lfdr, Idx = ind))
  }
  #  , Rej.L = ind1, Rej.M =ind2 , Rej.B = indBH, Rej.J = indQ, LFDR = lfdr, Idx = ind))
}

NonparamEst0_LAWS =function(p_value, delta.tr0, delta.tr1, basis.tr0, fdr.c=0.1
                               , index=FALSE, JS=TRUE, local =FALSE) { 
  
  #median pvalue
  nn <- length(p_value)
  pstar = pmin  =  pmax = numeric(length = nn)
  
  # for( ii in 1:nn ) {
  #   # neighborhood
  #   a = (dist2(matrix(loc0[ii,], ncol=2), loc0))[1,]
  #   short.dist = which(round(a,1) <= wi)
  # 
  #   #print(length(short.dist))
  #   p0 = ifelse(p_value[short.dist] < 1e-5, 1e-5, p_value[short.dist])
  #   pstar[ii] <- -2*sum(log(p0))
  #   pmin[ii] = median(p0)#; pmax[ii] = max(p0)
  # }
  pstar = cbind(1, pvalA, loc0)
  
  # FDR control using nonparametric estimate for f1
  N = length(p_value);
  pval = p_value
  out = sort(pval, index.return=TRUE)
  pval = out$x
  ind = out$ix
  pstar0 = pstar[ind,]; Nm = dim(pstar0)[2]
  basis = basis.tr0[ind,]
  delta = delta.tr0[ind]		
  df = c(pval[1],diff(pval))       # z(i)-z(i-1)
  f1.int = f1 = dexp(pval,100)     # initial yi=f1(z)
  pi0  = runif(N,0,1)		  # initial pi(s)
  f0 = rep(1,N)		        # known   f0(zi)
  
  ### ESTIMATION METHOD
  q = pi0*f0/(pi0*f0+(1-pi0)*f1);d = ncol(pstar0)
  beta.int = beta.cur = rep(1, d)  # initial value
  K = 100
  lik = -1e5
  iter = 1
  epsilon = 1
  while(iter <= K & epsilon > 1e-9)
  {
    # E-step: 
    q = pi0*f0/(pi0*f0+(1-pi0)*f1)
    
    # M-step (Newton method):
    Q = function(ps) {
      sum(-q*(pstar0 %*% ps)-log(1/(exp(pstar0 %*% ps)+1)))
    }
    grr = function(ps) {
      A = matrix(-q, N, Nm)*pstar0 
      C = exp(pstar0 %*% ps)
      C[C==Inf] = 1e+30
      B = pstar0*matrix(C/(C+1), N, Nm)
      return(apply(A+B,2,sum))
    }
    fit.optim = nlminb(beta.int, Q, gradient=grr)   ### Initial value !!!
    beta.cur = fit.optim$par
    pi0 = plogis(pstar0 %*% beta.cur)
    pi0 = ifelse(pi0<= 1e-5, 1e-5, pi0); pi0 = ifelse(pi0 >= 1-1e-5, 1-1e-5, pi0)
    
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
  
  # LAW approach 
  law.dd.res<-law.func(pvs=pval, pi0, 0.1)
  law.dd.de<-law.dd.res$de
  ind.law = which(law.dd.de==1)
  law.fdr0 = sum((1-delta)*law.dd.de)/max(sum(law.dd.de), 1)
  law.power0 = sum(delta*law.dd.de)/sum(delta)
  
  # Local approach
  #fdr.c = 0.05
  lfdr <- pi0*f0/(pi0*f0+(1-pi0)*f1)
  k <- sum(cumsum(sort(lfdr))/(1:(N)) < fdr.c )               # Nominal level 10%
  ind1 <- sort(lfdr, index.return=TRUE)$ix[1:k]
  sum(delta == 0)/length(delta)
  FDP.1 <- mean(delta[ind1] == 0)             # FDP (should be close to the nominal level)
  Power.1 <- sum(delta[ind1]) / (sum(delta.tr1==1))            # Power (the larger the better) 
  
  # Shape-adaptive approach
  f.step = stepfun(pval, c(f1[1],f1), right = FALSE) 
  fp = f.step(pval)
  f1p = f.step(1-pval)
  nu = (1-pi0)/pi0*f1p
  de = (1-pi0)/pi0*fp
  TT = sort(c(nu,de))
  i = 1
  hit = FALSE
  while(i <= length(TT) & hit == FALSE)
  {
    if ( sum(nu > TT[i])/max(1,sum(de > TT[i])) < fdr.c & sum(de > TT[i])>0) hit=TRUE
    i = i+1
  }
  if(hit == TRUE) {
    ind2 = (1:N)[de > TT[i]]
    FDP.2 = mean( delta[(1:N)[de > TT[i]]] ==0 )
    Power.2 = sum( delta[(1:N)[de > TT[i]]] ==1 )/(sum(delta.tr1==1))  
  } else
  {
    FDP.2 = 0;
    Power.2 = 0;
  }  
  
  # # John Storey approach
  # if (JS==TRUE) p.adj = qvalue(pval, fdr.level=fdr.c)
  # else if(JS==FALSE) p.adj = qvalue(pval, fdr.level=fdr.c, pi0=1)
  # indQ = which(p.adj$significant == TRUE) 	
  # FDP.J = mean(delta[indQ] == 0)                   
  # Power.J = sum(delta[indQ]) / (sum(delta.tr1 == 1))
  # # BH approach
  # 
  # p1.adj = p.adjust(pval, method='BH', n =length(pval))
  # indBH = which(p1.adj <fdr.c)
  # 
  # FDP.B <- mean(delta[indBH] == 0)                    # FDP (should be close to the nominal level)
  # Power.B <- sum(delta[indBH]) / (sum(delta.tr1==1));
  
  if(index==FALSE & local==FALSE) {
    
    return(list(FDP.L = FDP.1, FDP.A = FDP.2, Power.L = Power.1, Power.A = Power.2
                , FDP.LA = law.fdr0, Power.LA = law.power0))
  } else if (index == TRUE & local==FALSE){
    return(list(FDP.L = FDP.1, FDP.A = FDP.2, Power.L = Power.1, Power.A = Power.2
                , FDP.LA = law.fdr0, Power.LA = law.power0
                , Rej.L = ind1, Rej.M =ind2 , Rej.LA =ind.law , LFDR = lfdr, Idx = ind))
  } else if (index == TRUE & local==TRUE){
    return(list(FDP.L = FDP.1, FDP.A = FDP.2, Power.L = Power.1, Power.A = Power.2
                , FDP.LA = law.fdr0, Power.LA = law.power0
                , Rej.L = ind1, Rej.M =ind2 , Rej.LA =ind.law , LFDR = lfdr, Idx = ind))
  }
  #  , Rej.L = ind1, Rej.M =ind2 , Rej.B = indBH, Rej.J = indQ, LFDR = lfdr, Idx = ind))
}

## Estimate f0 using FDRL method (beta distribution) on median p-values
NonparamEst0_f0beta =function(p_value, delta.tr0, delta.tr1, basis.tr0, fdr.c=0.1
                              , index=FALSE, JS=TRUE, local =FALSE, beta=7) { 
  # FDR control using nonparametric estimate for f1
  N = length(p_value); Nm = dim(basis.tr0)[2]
  
  #median pvalue
  nn <- length(p_value)
  
  # pcov
  #pcov = cbind(1, pvalA, loc)
  pstar <- numeric(length = nn)
  
  for( ii in 1:nn ) {
    # neighborhood
    a = (dist2(matrix(loc0[ii,], ncol=2), loc0))[1,]
    short.dist = which(round(a,1) <= wi)
    
    #print(length(short.dist))
    pstar[ii] <- median(p_value[short.dist])
  }
  
  pval = pstar
  out = sort(pval, index.return=TRUE)
  pval = out$x
  ind = out$ix
  #pcov0 = pcov[ind,]; Nm = d =dim(pcov0)[2]
  basis = basis.tr0[ind,]
  delta = delta.tr0[ind]		
  df = c(pval[1],diff(pval))       # z(i)-z(i-1)
  f1.int = f1 = dexp(pval,100)     # initial yi=f1(z)
  pi0  = runif(N,0,1)		  # initial pi(s)W
  #f0 = rep(1,N)		        # known   f0(zi)
  f0 = dbeta(pval, (beta+1)/2,(beta+1)/2)
  #f0 = dnorm(pstar,0.5, 0.1)#/sum(dnorm(pstar,0.5, 0.1))
  
  ### ESTIMATION METHOD
  q = pi0*f0/(pi0*f0+(1-pi0)*f1);d = ncol(basis)
  beta.int = beta.cur = rep(1, d)  # initial value
  K = 100
  lik = -1e5
  iter = 1
  epsilon = 1
  while(iter <= K & epsilon > 1e-9)
  {
    # E-step: 
    q = pi0*f0/(pi0*f0+(1-pi0)*f1)
    
    # M-step (Newton method):
    # Q = function(ps) {
    #   sum(-q*(pcov0 %*% ps)-log(1/(exp(pcov0 %*% ps)+1)))
    # }
    # grr = function(ps) {
    #   A = matrix(-q, N, Nm)*pcov0 
    #   C = exp(pcov0 %*% ps)
    #   C[C==Inf] = 1e+30
    #   B = pcov0*matrix(C/(C+1), N, Nm)
    #   return(apply(A+B,2,sum))
    # }
    # fit.optim = nlminb(beta.int, Q, gradient=grr)   ### Initial value !!!
    # beta.cur = fit.optim$par
    # pi0 = plogis(pcov0 %*% beta.cur)
    # pi0 = ifelse(pi0<= 1e-5, 1e-5, pi0); pi0 = ifelse(pi0 >= 1-1e-5, 1-1e-5, pi0)
    
    
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
  
  # Local approach
  #fdr.c = 0.05
  lfdr <- pi0*f0/(pi0*f0+(1-pi0)*f1)
  k <- sum(cumsum(sort(lfdr))/(1:(N)) < fdr.c )               # Nominal level 10%
  ind1 <- sort(lfdr, index.return=TRUE)$ix[1:k]
  sum(delta == 0)/length(delta)
  FDP.1 <- mean(delta[ind1] == 0)             # FDP (should be close to the nominal level)
  Power.1 <- sum(delta[ind1]) / (sum(delta.tr1==1))            # Power (the larger the better) 
  
  # Shape-adaptive approach
  f.step = stepfun(pval, c(f1[1],f1), right = FALSE) 
  fp = f.step(pval)
  f1p = f.step(1-pval)
  nu = (1-pi0)/pi0*f1p
  de = (1-pi0)/pi0*fp
  TT = sort(c(nu,de))
  i = 1
  hit = FALSE
  while(i <= length(TT) & hit == FALSE)
  {
    if ( sum(nu > TT[i])/max(1,sum(de > TT[i])) < fdr.c & sum(de > TT[i])>0) hit=TRUE
    i = i+1
  }
  if(hit == TRUE) {
    ind2 = (1:N)[de > TT[i]]
    FDP.2 = mean( delta[(1:N)[de > TT[i]]] ==0 )
    Power.2 = sum( delta[(1:N)[de > TT[i]]] ==1 )/(sum(delta.tr1==1))  
  } else
  {
    FDP.2 = 0;
    Power.2 = 0;
  }  
  
  # John Storey approach
  if (JS==TRUE) p.adj = qvalue(pval, fdr.level=fdr.c)
  else if(JS==FALSE) p.adj = qvalue(pval, fdr.level=fdr.c, pi0=1)
  indQ = which(p.adj$significant == TRUE) 	
  FDP.J = mean(delta[indQ] == 0)                   
  Power.J = sum(delta[indQ]) / (sum(delta.tr1 == 1))
  # BH approach
  
  p1.adj = p.adjust(pval, method='BH', n =length(pval))
  indBH = which(p1.adj <fdr.c)
  
  FDP.B <- mean(delta[indBH] == 0)                    # FDP (should be close to the nominal level)
  Power.B <- sum(delta[indBH]) / (sum(delta.tr1==1));
  
  if(index==FALSE & local==FALSE) {
    
    return(list(FDP.L = FDP.1, FDP.A = FDP.2, Power.L = Power.1, Power.A = Power.2
                , FDP.J = FDP.J, FDP.B = FDP.B, Power.B = Power.B, Power.J = Power.J))
  } else if (index == TRUE & local==FALSE){
    return(list(FDP.L = FDP.1, FDP.A = FDP.2, Power.L = Power.1, Power.A = Power.2
                , FDP.J = FDP.J, FDP.B = FDP.B, Power.B = Power.B, Power.J = Power.J
                , Rej.L = ind1, Rej.M =ind2 , Rej.B = indBH, Rej.J = indQ, LFDR = lfdr, Idx = ind))
  } else if (index == TRUE & local==TRUE){
    return(list(FDP.L = FDP.1, FDP.A = FDP.2, Power.L = Power.1, Power.A = Power.2
                , FDP.J = FDP.J, FDP.B = FDP.B, Power.B = Power.B, Power.J = Power.J
                , Rej.L = ind1, Rej.B = indBH, Rej.J = indQ, LFDR = lfdr, Idx = ind))
  }
  #  , Rej.L = ind1, Rej.M =ind2 , Rej.B = indBH, Rej.J = indQ, LFDR = lfdr, Idx = ind))
}

NonparamEst0_f0beta_2 =function(p_value, delta.tr0, delta.tr1, basis.tr0, fdr.c=0.1
                              , index=FALSE, JS=TRUE, local =FALSE, beta=7) { 
  # FDR control using nonparametric estimate for f1
  N = length(p_value); Nm = dim(basis.tr0)[2]
  
  #median pvalue
  nn <- length(p_value)
  
  # pcov
  #pcov = cbind(1, pvalA, loc)
  pstar <- p_value
  
  # for( ii in 1:nn ) {
  #   # neighborhood
  #   a = (dist2(matrix(loc0[ii,], ncol=2), loc0))[1,]
  #   short.dist = which(round(a,1) <= wi)
  #   
  #   #print(length(short.dist))
  #   pstar[ii] <- median(p_value[short.dist])
  # }
  
  pval = pstar
  out = sort(pval, index.return=TRUE)
  pval = out$x
  ind = out$ix
  #pcov0 = pcov[ind,]; Nm = d =dim(pcov0)[2]
  basis = basis.tr0[ind,]
  delta = delta.tr0[ind]		
  df = c(pval[1],diff(pval))       # z(i)-z(i-1)
  f1.int = f1 = dexp(pval,100)     # initial yi=f1(z)
  pi0  = runif(N,0,1)		  # initial pi(s)W
  #f0 = rep(1,N)		        # known   f0(zi)
  f0 = dbeta(pval, (beta+1)/2,(beta+1)/2)
  #f0 = dnorm(pstar,0.5, 0.1)#/sum(dnorm(pstar,0.5, 0.1))
  
  ### ESTIMATION METHOD
  q = pi0*f0/(pi0*f0+(1-pi0)*f1);d = ncol(basis)
  beta.int = beta.cur = rep(1, d)  # initial value
  K = 100
  lik = -1e5
  iter = 1
  epsilon = 1
  while(iter <= K & epsilon > 1e-9)
  {
    # E-step: 
    q = pi0*f0/(pi0*f0+(1-pi0)*f1)
    
    # M-step (Newton method):
    # Q = function(ps) {
    #   sum(-q*(pcov0 %*% ps)-log(1/(exp(pcov0 %*% ps)+1)))
    # }
    # grr = function(ps) {
    #   A = matrix(-q, N, Nm)*pcov0 
    #   C = exp(pcov0 %*% ps)
    #   C[C==Inf] = 1e+30
    #   B = pcov0*matrix(C/(C+1), N, Nm)
    #   return(apply(A+B,2,sum))
    # }
    # fit.optim = nlminb(beta.int, Q, gradient=grr)   ### Initial value !!!
    # beta.cur = fit.optim$par
    # pi0 = plogis(pcov0 %*% beta.cur)
    # pi0 = ifelse(pi0<= 1e-5, 1e-5, pi0); pi0 = ifelse(pi0 >= 1-1e-5, 1-1e-5, pi0)
    
    
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
  
  # Local approach
  #fdr.c = 0.05
  lfdr <- pi0*f0/(pi0*f0+(1-pi0)*f1)
  k <- sum(cumsum(sort(lfdr))/(1:(N)) < fdr.c )               # Nominal level 10%
  ind1 <- sort(lfdr, index.return=TRUE)$ix[1:k]
  sum(delta == 0)/length(delta)
  FDP.1 <- mean(delta[ind1] == 0)             # FDP (should be close to the nominal level)
  Power.1 <- sum(delta[ind1]) / (sum(delta.tr1==1))            # Power (the larger the better) 
  
  # Shape-adaptive approach
  f.step = stepfun(pval, c(f1[1],f1), right = FALSE) 
  fp = f.step(pval)
  f1p = f.step(1-pval)
  nu = (1-pi0)/pi0*f1p
  de = (1-pi0)/pi0*fp
  TT = sort(c(nu,de))
  i = 1
  hit = FALSE
  while(i <= length(TT) & hit == FALSE)
  {
    if ( sum(nu > TT[i])/max(1,sum(de > TT[i])) < fdr.c & sum(de > TT[i])>0) hit=TRUE
    i = i+1
  }
  if(hit == TRUE) {
    ind2 = (1:N)[de > TT[i]]
    FDP.2 = mean( delta[(1:N)[de > TT[i]]] ==0 )
    Power.2 = sum( delta[(1:N)[de > TT[i]]] ==1 )/(sum(delta.tr1==1))  
  } else
  {
    FDP.2 = 0;
    Power.2 = 0;
  }  
  
  # John Storey approach
  if (JS==TRUE) p.adj = qvalue(pval, fdr.level=fdr.c)
  else if(JS==FALSE) p.adj = qvalue(pval, fdr.level=fdr.c, pi0=1)
  indQ = which(p.adj$significant == TRUE) 	
  FDP.J = mean(delta[indQ] == 0)                   
  Power.J = sum(delta[indQ]) / (sum(delta.tr1 == 1))
  # BH approach
  
  p1.adj = p.adjust(pval, method='BH', n =length(pval))
  indBH = which(p1.adj <fdr.c)
  
  FDP.B <- mean(delta[indBH] == 0)                    # FDP (should be close to the nominal level)
  Power.B <- sum(delta[indBH]) / (sum(delta.tr1==1));
  
  if(index==FALSE & local==FALSE) {
    
    return(list(FDP.L = FDP.1, FDP.A = FDP.2, Power.L = Power.1, Power.A = Power.2
                , FDP.J = FDP.J, FDP.B = FDP.B, Power.B = Power.B, Power.J = Power.J))
  } else if (index == TRUE & local==FALSE){
    return(list(FDP.L = FDP.1, FDP.A = FDP.2, Power.L = Power.1, Power.A = Power.2
                , FDP.J = FDP.J, FDP.B = FDP.B, Power.B = Power.B, Power.J = Power.J
                , Rej.L = ind1, Rej.M =ind2 , Rej.B = indBH, Rej.J = indQ, LFDR = lfdr, Idx = ind))
  } else if (index == TRUE & local==TRUE){
    return(list(FDP.L = FDP.1, FDP.A = FDP.2, Power.L = Power.1, Power.A = Power.2
                , FDP.J = FDP.J, FDP.B = FDP.B, Power.B = Power.B, Power.J = Power.J
                , Rej.L = ind1, Rej.B = indBH, Rej.J = indQ, LFDR = lfdr, Idx = ind))
  }
  #  , Rej.L = ind1, Rej.M =ind2 , Rej.B = indBH, Rej.J = indQ, LFDR = lfdr, Idx = ind))
}

## 1. Estimate f1 / After trim / 2. estimate f0 
NonparamEst0_f0 =function(p_value0, p_valueA,trim, delta.tr0, basis.tr0, fdr.c=0.1
                          , index=FALSE, JS=TRUE, local =FALSE, Agg=FALSE) { 
  # FDR control using nonparametric estimate for f1
  
  # p-value before trimming 
  if (Agg==T) { 
    N = length(p_valueA); Nm = dim(basis.tr0)[2]
    pval = p_valueA
  } else if (Agg==F) { 
    N = length(p_value0); Nm = dim(basis.tr0)[2]
    pval = p_value0
  }
  out = sort(pval, index.return=TRUE)
  pval = out$x
  ind = out$ix
  basis = basis.tr0[ind,]
  delta = delta.tr0[ind]		
  df = c(pval[1],diff(pval))       # z(i)-z(i-1)
  f1.int = f1 = dbeta(pval, 1, 50)     # initial yi=f1(z)
  #f0.int = f0 = dexp(pval, 100)
  pi0  = runif(N,0,1)		  # initial pi(s)
  f0 = rep(1,N)		        # known   f0(zi)
  trim.ind = match(trim, ind)
  
  ### ESTIMATION METHOD
  q = pi0*f0/(pi0*f0+(1-pi0)*f1);d = ncol(basis)
  beta.int = beta.cur = rep(1, d)  # initial value
  K = 200
  lik = -1e5
  iter = 1
  epsilon = 1
  while(iter <= K & epsilon > 1e-9)
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
    
    # PAVA - f1
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
  
  f.step = stepfun(pval, c(f1[1],f1), right = FALSE) 
  #pi0 = pi0[trim.ind]
  
  # Estimate the null after trimming 
  #p_value = SimGEV.boot[[7]][trim10, 1]
  delta.tr00 = delta.tr0[trim]
  basis.tr0 = basis.tr0[trim,]
  delta.tr1 = delta.tr0
  N = length(trim10); Nm = dim(basis.tr0)[2]
  pval = p_value0[trim]
  
  out = sort(pval, index.return=TRUE)
  pval = out$x
  
  f1.c = f.step(pval)
  
  ind = out$ix
  basis = basis.tr0[ind,]
  delta = delta.tr00[ind]		
  df = c(pval[1],diff(pval))       # z(i)-z(i-1)
  #f1.c = f1.c[ind]
  f0.int = f0 = dexp(pval, 100)
  pi0  = runif(N,0,1)		  # initial pi(s)
  
  ### ESTIMATION METHOD
  q = pi0*f0/(pi0*f0+(1-pi0)*f1.c);d = ncol(basis)
  beta.int = beta.cur = rep(1, d)  # initial value
  K = 200
  lik = -1e5
  iter = 1
  epsilon = 1
  while(iter <= K & epsilon > 1e-9)
  {
    # E-step: 
    q = pi0*f0/(pi0*f0+(1-pi0)*f1.c)
    
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
    
    # PAVA - f0
    q0 = q
    qsum0 = sum(q0)
    y0 = -df*qsum0/q0
    
    ### modify the code slightly
    y0[is.na(y0)] = -Inf
    y0[y0==0] = max(y0[y0!=0])
    y0[y0==-Inf] = min(y0[y0!=-Inf])
    
    out = pava(y0, q0, decreasing=TRUE, long.out=FALSE, stepfun=FALSE)
    f0 = -1/out
    f0 = f0/sum( f0*df )
    
    iter = iter + 1
    lik.cur = sum(log(pi0*f0+(1-pi0)*f1.c))
    epsilon = (lik.cur - lik[length(lik)])/abs(lik[length(lik)])
    lik = c(lik, lik.cur)
    if(iter%%200==0) print(iter)
  }
  par(mfrow=c(1,2))
  plot(pval, f1.c, main = 'Estimated f1')
  plot(pval, f0, main = 'Estimated f0 after trimming')
  # Local approach
  delta.tr1 = delta.tr0
  lfdr <- pi0*f0/(pi0*f0+(1-pi0)*f1.c)
  k <- sum(cumsum(sort(lfdr))/(1:(N)) < fdr.c )               # Nominal level 10%
  ind1 <- sort(lfdr, index.return=TRUE)$ix[1:k]
  sum(delta == 0)/length(delta)
  FDP.1 <- mean(delta[ind1] == 0)             # FDP (should be close to the nominal level)
  Power.1 <- sum(delta[ind1]) / (sum(delta.tr1==1))            # Power (the larger the better) 
  
  # Shape-adaptive approach
  f.step = stepfun(pval, c(f1.c[1],f1.c), right = FALSE) 
  fp = f.step(pval)
  f1p = f.step(1-pval)
  nu = (1-pi0)/pi0*f1p
  de = (1-pi0)/pi0*fp
  TT = sort(c(nu,de))
  i = 1
  hit = FALSE
  
  while(i <= length(TT) & hit == FALSE)
  {
    if ( sum(nu > TT[i])/max(1,sum(de > TT[i])) < fdr.c & sum(de > TT[i])>0) hit=TRUE
    i = i+1
  }
  if(hit == TRUE) {
    ind2 = (1:N)[de > TT[i]]
    FDP.2 = mean( delta[(1:N)[de > TT[i]]] ==0 )
    Power.2 = sum( delta[(1:N)[de > TT[i]]] ==1 )/(sum(delta.tr1==1))  
  } else
  {
    FDP.2 = 0;
    Power.2 = 0;
  }  
  
  # John Storey approach
  if (JS==TRUE) p.adj = qvalue(pval, fdr.level=fdr.c)
  else if(JS==FALSE) p.adj = qvalue(pval, fdr.level=fdr.c, pi0=1)
  indQ = which(p.adj$significant == TRUE) 	
  FDP.J = mean(delta[indQ] == 0)                   
  Power.J = sum(delta[indQ]) / (sum(delta.tr1 == 1))
  # BH approach
  
  p1.adj = p.adjust(pval, method='BH', n =length(pval))
  indBH = which(p1.adj <fdr.c)
  
  FDP.B <- mean(delta[indBH] == 0)                    # FDP (should be close to the nominal level)
  Power.B <- sum(delta[indBH]) / (sum(delta.tr1==1));
  
  if(index==FALSE) {
    
    return(list(FDP.L = FDP.1, FDP.A = FDP.2, Power.L = Power.1, Power.A = Power.2
                , FDP.J = FDP.J, FDP.B = FDP.B, Power.B = Power.B, Power.J = Power.J))
  } else if (index == TRUE & local==FALSE){
    return(list(FDP.L = FDP.1, FDP.A = FDP.2, Power.L = Power.1, Power.A = Power.2
                , FDP.J = FDP.J, FDP.B = FDP.B, Power.B = Power.B, Power.J = Power.J
                , Rej.L = ind1, Rej.M =ind2 , Rej.B = indBH, Rej.J = indQ, LFDR = lfdr, Idx = ind))
  } else if (index == TRUE & local==TRUE){
    return(list(FDP.L = FDP.1, FDP.A = FDP.2, Power.L = Power.1, Power.A = Power.2
                , FDP.J = FDP.J, FDP.B = FDP.B, Power.B = Power.B, Power.J = Power.J
                , Rej.L = ind1, Rej.B = indBH, Rej.J = indQ, LFDR = lfdr, Idx = ind))
  }
  #  , Rej.L = ind1, Rej.M =ind2 , Rej.B = indBH, Rej.J = indQ, LFDR = lfdr, Idx = ind))
}
gev.fitF = function (xdat, ydat = NULL, mul = NULL, sigl = NULL, shl = NULL, 
                     mulink = identity, siglink = identity, shlink = identity, 
                     muinit = NULL, siginit = NULL, shinit = NULL, show = TRUE, 
                     method = "Nelder-Mead", maxit = 10000, ...) 
{
  z <- list()
  npmu <- length(mul) + 1
  npsc <- length(sigl) + 1
  npsh <- length(shl) + 1
  z$trans <- FALSE
  in2 <- sqrt(6 * var(xdat))/pi
  in1 <- mean(xdat) - 0.57722 * in2
  if (is.null(mul)) {
    mumat <- as.matrix(rep(1, length(xdat)))
    if (is.null(muinit)) 
      muinit <- in1
  }
  else {
    z$trans <- TRUE
    mumat <- cbind(rep(1, length(xdat)), ydat[, mul])
    if (is.null(muinit)) 
      muinit <- c(in1, rep(0, length(mul)))
  }
  if (is.null(sigl)) {
    sigmat <- as.matrix(rep(1, length(xdat)))
    if (is.null(siginit)) 
      siginit <- in2
  }
  else {
    z$trans <- TRUE
    sigmat <- cbind(rep(1, length(xdat)), ydat[, sigl])
    if (is.null(siginit)) 
      siginit <- c(in2, rep(0, length(sigl)))
  }
  if (is.null(shl)) {
    shmat <- as.matrix(rep(1, length(xdat)))
    if (is.null(shinit)) 
      shinit <- 0.1
  }
  else {
    z$trans <- TRUE
    shmat <- cbind(rep(1, length(xdat)), ydat[, shl])
    if (is.null(shinit)) 
      shinit <- c(0.1, rep(0, length(shl)))
  }
  z$model <- list(mul, sigl, shl)
  z$link <- deparse(substitute(c(mulink, siglink, shlink)))
  init <- c(muinit, siginit, shinit)
  gev.lik.Fsh <- function(a) {
    mu <- mulink(mumat %*% (a[1:npmu]))
    sc <- siglink(sigmat %*% (a[seq(npmu + 1, length = npsc)]))
    xi <- 0.24#shlink(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)]))
    y <- (xdat - mu)/sc
    y <- 1 + xi * y
    if (any(y <= 0) || any(sc <= 0))
      return(10^6)
    sum(log(sc)) + sum(y^(-1/xi)) + sum(log(y) * (1/xi + 
                                                    1))
  }
  x <- optim(init[-3], gev.lik.Fsh, hessian = TRUE, method = method, 
             control = list(maxit = maxit, ...))
  z$conv <- x$convergence
  mu <- mulink(mumat %*% (x$par[1:npmu]))
  sc <- siglink(sigmat %*% (x$par[seq(npmu + 1, length = npsc)]))
  xi <- 0.24 #shlink(shmat %*% (x$par[seq(npmu + npsc + 1, length = npsh)]))
  z$nllh <- x$value
  z$data <- xdat
  if (z$trans) {
    z$data <- -log(as.vector((1 + (xi * (xdat - mu))/sc)^(-1/xi)))
  }
  z$mle <- x$par
  z$cov <- solve(x$hessian)
  z$se <- sqrt(diag(z$cov))
  z$vals <- cbind(mu, sc, xi)
  if (show) {
    if (z$trans) 
      print(z[c(2, 3, 4)])
    else print(z[4])
    if (!z$conv) 
      print(z[c(5, 7, 9)])
  }
  class(z) <- "gev.fit"
  invisible(z)
}


## Return level test 
Rt_test = function(pval , nn, trim=FALSE, JS=TRUE, pval2) { 
  
  if ( trim==FALSE) {
    test.res = c()
    for ( i in 1:nn) { 
      if (JS==TRUE) test.tmp = NonparamEst0(pval[,i], delta.tr0, delta.tr0, basis.tr)
      else if (JS==FALSE) test.tmp = NonparamEst0(pval[,i], delta.tr0, delta.tr0, basis.tr, JS=F)
      test.tmp = unlist(test.tmp)
      test.res = rbind(test.res, test.tmp)
    }
    FDR = test.res[, c(1,2,5,6)]; Power = test.res[, -c(1,2,5,6)]
  } else if (trim==TRUE){
    # trim.tmp = NonparamEst0(pval2[,1],  delta.tr0, delta.tr0, basis.tr[ab, ], ind=T, local=F)
    # trim.rej = trim.tmp$Idx[trim.tmp$Rej.M]
    test.res = c()
    for ( i in 1:nn) { 
      if (JS==TRUE) trim.tmp = NonparamEst0(pval2[,i],  delta.tr0, delta.tr0, basis.tr, ind=T, local=F)
      else if (JS==FALSE) trim.tmp = NonparamEst0(pval2[,i],  delta.tr0, delta.tr0, basis.tr, ind=T, local=F, JS=F)
      trim.rej = trim.tmp$Idx[trim.tmp$Rej.M]
      
      test.tmp = NonparamEst0(pval[trim.rej,i], delta.tr0[trim.rej], delta.tr0
                              , basis.tr[ab[trim.rej], ], JS = JS)
      test.tmp = unlist(test.tmp)
      test.res = rbind(test.res, test.tmp)
    }
    FDR = test.res[, c(1,2,5,6)]; Power = test.res[, -c(1,2,5,6)]
  }
  return(list(FDR = FDR, Power=Power))
  
  
  
}
