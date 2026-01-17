Mirror_Piest = function(p_value, basis.tr0, delta.tr0,fdr.c=0.1, JS=TRUE){
  N = length(p_value); Nm = dim(basis.tr0)[2]
  pval = p_value
  out = sort(pval, index.return=TRUE)
  pval = out$x
  ind = out$ix
  basis = basis.tr0[ind,]
  delta = delta.tr0[ind]
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
    Power.2 = sum( delta[(1:N)[de > TT[i]]] ==1 )/(sum(delta==1))
  } else
  {
    FDP.2 = 0;
    Power.2 = 0;
  }
  
  
  p1.adj = p.adjust(pval, method='BH', n =length(pval))
  indBH = which(p1.adj < fdr.c)
  #indBH2 = which(p1.adj < fdr.c2)
  
  FDP.B <- mean(delta[indBH] == 0)                    # FDP (should be close to the nominal level)
  Power.B <- sum(delta[indBH]) / (sum(delta==1));
  #FDP.B2 <- mean(delta[indBH2] == 0)                    # FDP (should be close to the nominal level)
  #Power.B2 <- sum(delta[indBH2]) / (sum(delta==1));
  # John Storey approach
  if(JS==TRUE) {
    p.adj = qvalue(pval, fdr.level=fdr.c)
    #p.adj = qvalue(pval, fdr.level=fdr.c, pi0=1)
    indQ = which(p.adj$significant == TRUE)
    FDP.J = mean(delta[indQ] == 0)
    Power.J = sum(delta[indQ]) / (sum(delta == 1))
    # BH approach
  } else
  {
    FDP.J = FDP.B
    Power.J = Power.B
  }
  
  f0.marg = f1.marg = pi0.marg = p.marg = rep(NA, N)
  f0.marg[ind] = f0;f1.marg[ind] = f1;pi0.marg[ind] = pi0
  p.marg = pi0.marg*f0.marg+(1-pi0.marg)*f1.marg
  p.cond = pi0.marg*f0.marg/p.marg
  
  Rej.M = ind[ind2]; Rej.B = ind[indBH]#;Rej.B0 = ind[indBH]
  MirrRes = c(FDP.2, Power.2);  BHRes = c(FDP.B, Power.B);  #BHRes2 = c(FDP.B2, Power.B2)
  JSRes = c(FDP.J, Power.J);#BHRes2 = c(FDP.B2, Power.B2)
  return(list(piest=pi0.marg, lfdr= p.cond, f0est = f0.marg, f1est = f1.marg,Rej = Rej.M, RejB = Rej.B,MirrRes = MirrRes, JSRes = JSRes,BHRes = BHRes))
}

#### GEV function #####
df_dg0 = function(param, return){
  pm = param
  rt = return
  yRes = c(1, -pm[3]^(-1)*(1-rt^(-pm[3])), pm[2]*pm[3]^(-2)*(1-rt^(-pm[3]))-pm[2]*pm[3]^(-1)*rt^(-pm[3])*log(rt))
  
  return(yRes)
}

gev.fit = function (xdat, ydat = NULL, mul = NULL, sigl = NULL, shl = NULL,
                    mulink = identity, siglink = identity, shlink = identity,
                    muinit = NULL, siginit = NULL, shinit = NULL, show = TRUE,
                    method = "L-BFGS-B", maxit = 10000, ...)
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
  gev.lik <- function(a) {
    mu <- mulink(mumat %*% (a[1:npmu]))
    sc <- siglink(sigmat %*% (a[seq(npmu + 1, length = npsc)]))
    xi <- shlink(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)]))
    y <- (xdat - mu)/sc
    y <- 1 + xi * y
    if (any(y <= 0) || any(sc <= 0))
      return(10^6)
    sum(log(sc)) + sum(y^(-1/xi)) + sum(log(y) * (1/xi +
                                                    1))
  }
  x <- optim(init, gev.lik, hessian = TRUE, method = method,
             control = list(maxit = maxit, ...), lower = c(1, 1.5, 0.01), upper=c(60, 60, 0.999))
  z$conv <- x$convergence
  mu <- mulink(mumat %*% (x$par[1:npmu]))
  sc <- siglink(sigmat %*% (x$par[seq(npmu + 1, length = npsc)]))
  xi <- shlink(shmat %*% (x$par[seq(npmu + npsc + 1, length = npsh)]))
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
    
    #Sig_00 = matrix(0,ncol=2, nrow=2);
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
###########################################################################################
BiModel_Test = function(pval, nei, delta.tr0,delta.tr00, basis.tr0, fdr.c=0.1, index=FALSE, den=FALSE) {
  # FDR control using nonparametric estimate for f1
  stat = pvalue_to_teststat(pval)
  pnei = neigh(pval,nei)
  #pdiff = abs(pnei[,2] - pnei[,1])
  #idp = which(pdiff < 0.05)
  statnei = neigh(stat, nei)
  
   if (nei!=1) {
   stat.neigh = apply(statnei[,-1], 1, median)
   stat = cbind(stat, stat.neigh)
   } else if (nei==1){
  stat = statnei
  }
  #stat.neigh = pvalue_to_teststat(pnei.mean)
  
  AA = Mirror_Piest(pval, basis.tr0, delta.tr0)
  #stat = neigh(stat, nei);
  
  N = nrow(stat); Nm = dim(basis.tr0)[2]
  
  basis = basis.tr0
  delta = delta.tr0
  
  f00 = null_joint(stat)
  #nullalt.est = nullalt_joint(pnei, basis, delta)
  
   if (nei!=1) {
     pi.est = neigh(AA$piest, nei)
     pi.est = cbind(pi.est[,1], apply(pi.est[,-1], 1, median))
   } else if (nei==1){
  #pi.est = neigh0(A$piest, nei, idp)
  pi.est = neigh(AA$piest, nei)
  pi.est =pi.est
   }
  
  Sig10est = matrix()
  f10 = dnorm(stat[,1], -1, 0.5)*f00
  f01 = dnorm(stat[,2], -1, 0.5)*f00
  f11 = dnorm(stat[,1], -1, 0.5)*dnorm(stat[,2], -1, 0.5)
  #f10 = nullalt.est$f10est; f01 = nullalt.est$f01est; f11 = nullalt.est$f11est
  
  # initial values
  pi0.1 = pi.est[,1]; pi0.2 = pi.est[,2]
  ### ESTIMATION METHOD
  fmarg = pi0.1*pi0.2*f00 + pi0.1*(1-pi0.2)*f01 + (1-pi0.1)*pi0.2*f10 + (1-pi0.1)*(1-pi0.2)*f11
  q1 = (pi0.1*pi0.2*f00 + pi0.1*(1-pi0.2)*f01)/fmarg
  q2 = (pi0.1*pi0.2*f00 + (1-pi0.1)*pi0.2*f10)/fmarg
  
  beta.int = beta.cur = gamma.cur = gamma.int = rep(1, ncol(basis)) # initial value
  K = 100
  lik = -1e5
  iter = 1
  epsilon = 1
  d = Nm = ncol(basis)
  while(iter <= K & epsilon > 1e-6)
  {
    # E-step:
    q1 = (pi0.1*pi0.2*f00 + pi0.1*(1-pi0.2)*f01)/fmarg
    
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
    
    # if (nei!=1) {
    #   pitmp = neighP(pi0.1, nei)
    #   pi0.2 = apply(pitmp[,-1], 1, median)
    # } else if (nei==1){
    pitmp = neigh(pi0.1, nei)
    pi0.2 = pitmp[,2]
    #}
    
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
    #Sig_00 = matrix(0,ncol=2, nrow=2)
    Sig_00[1,1] = Sig_00[2,2] = 1#;
    #Sig_00[1,2] = Sig_00[2,1] = 0
    Sig_01 = matrix((1/N01)*apply(diag(w01) %*% t(apply(stat-mu_01, 1, function(x) (x) %*% t(x))), 2,sum), ncol=2)
    Sig_10 = matrix((1/N10)*apply(diag(w10) %*% t(apply(stat-mu_10, 1, function(x) (x) %*% t(x))), 2,sum),ncol=2)
    Sig_11 = matrix((1/N11)*apply(diag(w11) %*% t(apply(stat-mu_11, 1, function(x) (x) %*% t(x))), 2,sum),ncol=2)
    Sig_01[1,1] = Sig_10[2,2]=1
    #Sig_01[2,2] = Sig_11[2,2];Sig_10[1,1] = Sig_11[1,1]
    
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
    if(iter%%200==0) print(iter)
  }
  CLFDR =  (pi0.1*pi0.2*f00 + pi0.1*(1-pi0.2)*f01)/fmarg
  
  #CLFDR =  (pi0.1*pi0.2*f00)/fmarg
  # Local approach
  k <- sum(cumsum(sort(CLFDR))/(1:(N)) < fdr.c )               # Nominal level 10%
  ind1 <- sort(CLFDR, index.return=TRUE)$ix[1:k]
  FDP.1 <- mean(delta[ind1] == 0)             # FDP (should be close to the nominal level)
  Power.1 <- sum(delta[ind1]) / (sum(delta==1))            # Power (the larger the better)
  
  if (index==FALSE & den==FALSE){
    Res_B = c(AA$MirrRes,AA$JSRes, AA$BHRes, FDP.1, Power.1)
    return(Res_B)
  }else if (index==TRUE& den==FALSE){
    Res_B = c(AA$MirrRes,AA$JSRes, AA$BHRes, FDP.1, Power.1)
    return(list(Res = Res_B, RejM=AA$Rej, RejBi = ind1, BLFDR = CLFDR))
  } else if (index==FALSE & den==TRUE){
    Res_B = c(AA$MirrRes,AA$JSRes, AA$BHRes, FDP.1, Power.1)
    Den = cbind(f00, f01, f10, f11)
    return(list(Res=Res_B, density=Den))
  }
}
