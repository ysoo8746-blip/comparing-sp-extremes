gev.fit <- function(xdat, ydat = NULL,
                    mul = NULL, sigl = NULL, shl = NULL,
                    mulink = identity, siglink = identity, shlink = identity,
                    muinit = NULL, siginit = NULL, shinit = NULL,
                    show = TRUE, method = "L-BFGS-B", maxit = 10000, ...) {
  
  z <- list()
  npmu <- length(mul) + 1
  npsc <- length(sigl) + 1
  npsh <- length(shl) + 1
  z$trans <- FALSE
  
  in2 <- sqrt(6 * var(xdat, na.rm = TRUE)) / pi
  in1 <- mean(xdat, na.rm = TRUE) - 0.57722 * in2
  
  if (is.null(mul)) {
    mumat <- matrix(1, length(xdat), 1)
    if (is.null(muinit)) muinit <- in1
  } else {
    z$trans <- TRUE
    mumat <- cbind(1, ydat[, mul])
    if (is.null(muinit)) muinit <- c(in1, rep(0, length(mul)))
  }
  
  if (is.null(sigl)) {
    sigmat <- matrix(1, length(xdat), 1)
    if (is.null(siginit)) siginit <- in2
  } else {
    z$trans <- TRUE
    sigmat <- cbind(1, ydat[, sigl])
    if (is.null(siginit)) siginit <- c(in2, rep(0, length(sigl)))
  }
  
  if (is.null(shl)) {
    shmat <- matrix(1, length(xdat), 1)
    if (is.null(shinit)) shinit <- 0.1
  } else {
    z$trans <- TRUE
    shmat <- cbind(1, ydat[, shl])
    if (is.null(shinit)) shinit <- c(0.1, rep(0, length(shl)))
  }
  
  z$model <- list(mul, sigl, shl)
  z$link <- deparse(substitute(c(mulink, siglink, shlink)))
  
  init <- c(muinit, siginit, shinit)
  
  ## ---- Robust likelihood ---- ##
  gev.lik <- function(a) {
    mu <- as.vector(mulink(mumat %*% a[1:npmu]))
    sc <- as.vector(siglink(sigmat %*% a[seq(npmu + 1, length = npsc)]))
    xi <- as.vector(shlink(shmat %*% a[seq(npmu + npsc + 1, length = npsh)]))
    
    if (any(!is.finite(mu)) || any(!is.finite(sc)) || any(!is.finite(xi)))
      return(1e6)
    
    if (any(sc <= 0))
      return(1e6)
    
    y <- (xdat - mu) / sc
    z <- 1 + xi * y
    
    if (any(z <= 0) || any(!is.finite(z)))
      return(1e6)
    
    loglik <- sum(log(sc)) + sum(z^(-1/xi)) + sum(log(z) * (1/xi + 1))
    if (!is.finite(loglik)) return(1e6)
    return(loglik)
  }
  
  ## ---- Optimization ---- ##
  if (!is.null(ydat)) {
    lower <- c(rep(-50, npmu), rep(0.001, npsc), rep(-0.1, npsh))
    upper <- c(rep(50, npmu),  rep(50, npsc),    rep(0.3, npsh))
  } else {
    #lower <- c(rep(-200, npmu), rep(0.0001, npsc), rep(-1.5, npsh))
    #upper <- c(rep(200, npmu),  rep(200, npsc),   rep(1.5, npsh))
    lower <- c(rep(1, npmu), rep(1, npsc), rep(-0.05, npsh))
    upper <- c(rep(100, npmu),  rep(50, npsc),   rep(0.5, npsh))
  }
  
  x <- tryCatch({
    optim(
      par = init,
      fn = gev.lik,
      hessian = TRUE,
      method = method,
      control = list(maxit = maxit, ...),
      lower = lower,
      upper = upper
    )
  }, error = function(e) {
    warning("Optimization failed: ", conditionMessage(e))
    return(NULL)
  })
  
  if (is.null(x)) return(NULL)
  
  z$conv <- x$convergence
  
  mu <- as.vector(mulink(mumat %*% x$par[1:npmu]))
  sc <- as.vector(siglink(sigmat %*% x$par[seq(npmu + 1, length = npsc)]))
  xi <- as.vector(shlink(shmat %*% x$par[seq(npmu + npsc + 1, length = npsh)]))
  
  z$nllh <- x$value
  z$data <- xdat
  
  if (z$trans) {
    z$data <- -log((1 + (xi * (xdat - mu)) / sc)^(-1 / xi))
  }
  
  z$mle <- x$par
  z$cov <- tryCatch(solve(x$hessian), error = function(e) matrix(NA, length(x$par), length(x$par)))
  z$se <- sqrt(diag(z$cov))
  z$vals <- cbind(mu, sc, xi)
  
  if (show) {
    if (z$trans) print(z[c("mle", "se", "vals")])
    else print(z["vals"])
    if (!z$conv) print(z[c("mle", "se", "nllh")])
  }
  
  class(z) <- "gev.fit"
  invisible(z)
}


gev.return_level = function(loc,scale,shape,n.site,time){
  if (length(time)==1){
    p=1-1/time
    y<- rep(-log(p), n.site)
    tr.return <- loc - ifelse(abs(shape)>1e-5, scale*(1-y^(-shape))/shape, scale*log(y))
  } else if (length(time)!=0){
    tr.return=c()
    for (i in 1:length(time)){
      p=1-1/time[i]
      y<- -log(p)#rep(-log(p), n.site)
      tr.return <- cbind(tr.return, loc - ifelse(abs(shape)>1e-5, scale*(1-y^(-shape))/shape, scale*log(y)))
    }
    
  }
  
  return(tr.return)
}
##################################################################
# #########--- FIT GEV at each location ----#############
##################################################################
Fit_GEV_RT = function(SimX, SimY,paramx, paramy,T) {
  
  dataxp = datayp = matrix(nrow=n.obs, ncol=n.site )
  for (i in 1:n.site){
    dataxp[,i] <- frech2gev(SimX[,i], paramx[i,1], paramx[i,2],paramx[i,3])
    datayp[,i] <- frech2gev(SimY[,i], paramy[i,1], paramy[i,2],paramy[i,3])
    dataxp[,i] = pmin(dataxp[,i], 1000); dataxp[,i] = pmax(0.05, dataxp[,i])
    datayp[,i] = pmin(datayp[,i], 1000); datayp[,i] = pmax(0.05, datayp[,i])
  }
  
  paramx.mle.p = paramy.mle.p =c()
  for ( i in 1:n.site){
    datX.P = dataxp[,i];#datX.T = dataxt[,i]
    datY.P = datayp[,i];#datY.T = datayt[,i]
    
    modxp <-gev.fit(datX.P, show=FALSE);#modxt <-gev.fit(datX.T, show=FALSE)
    modyp <-gev.fit(datY.P, show=FALSE);#modyt <-gev.fit(datY.T, show=FALSE)
    
    paramx.mle.p = rbind(paramx.mle.p, modxp$mle)
    paramy.mle.p = rbind(paramy.mle.p, modyp$mle)
    
    if (i %% 100 ==0) print(i)
  }
  
  # GEV RL for precip
  rl_5_GEVX <- gev.return_level(paramx.mle.p[,1], paramx.mle.p[,2], paramx.mle.p[,3], n.site, time = T[1])
  rl_10_GEVX <- gev.return_level(paramx.mle.p[,1], paramx.mle.p[,2], paramx.mle.p[,3], n.site, time =  T[2])
  rl_20_GEVX <- gev.return_level(paramx.mle.p[,1], paramx.mle.p[,2], paramx.mle.p[,3], n.site, time =  T[3])
  rl_50_GEVX <- gev.return_level(paramx.mle.p[,1], paramx.mle.p[,2], paramx.mle.p[,3], n.site, time =  T[4])
  rl_100_GEVX <- gev.return_level(paramx.mle.p[,1], paramx.mle.p[,2], paramx.mle.p[,3], n.site, time =  T[5])
  rtX.GEV.P = cbind(rl_5_GEVX, rl_10_GEVX, rl_20_GEVX, rl_50_GEVX, rl_100_GEVX)
  colnames(rtX.GEV.P) = c('Rt5','Rt10', 'Rt20', 'Rt50', 'Rt100')
  summary(rtX.GEV.P) ## ESTIMATED RETURN LEVEL OF X
  rl_5_GEVY <- gev.return_level(paramy.mle.p[,1], paramy.mle.p[,2], paramy.mle.p[,3], n.site, time = T[1])
  rl_10_GEVY <- gev.return_level(paramy.mle.p[,1], paramy.mle.p[,2], paramy.mle.p[,3], n.site, time = T[2])
  rl_20_GEVY <- gev.return_level(paramy.mle.p[,1], paramy.mle.p[,2], paramy.mle.p[,3], n.site, time = T[3])
  rl_50_GEVY <- gev.return_level(paramy.mle.p[,1], paramy.mle.p[,2], paramy.mle.p[,3], n.site, time = T[4])
  rl_100_GEVY <- gev.return_level(paramy.mle.p[,1], paramy.mle.p[,2], paramy.mle.p[,3], n.site, time = T[5])
  rtY.GEV.P = cbind(rl_5_GEVY,rl_10_GEVY, rl_20_GEVY, rl_50_GEVY, rl_100_GEVY)
  colnames(rtY.GEV.P) = c('Rt5','Rt10', 'Rt20', 'Rt50', 'Rt100')
  summary(rtY.GEV.P) ## ESTIMATED RETURN LEVEL OF Y
  
  param.diff.GEVP = paramy.mle.p- paramx.mle.p
  rt.diff.GEVP = rtY.GEV.P -rtX.GEV.P
  colnames(rt.diff.GEVP) = c('Rt5_diff','Rt10_diff', 'Rt20_diff', 'Rt50_diff', 'Rt100_diff')
  summary(rt.diff.GEVP) ## ESTIMATED RETURN LEVEL DIFF
  
  return(list(Rt.X =rtX.GEV.P, Rt.Y = rtY.GEV.P, Param.X = paramx.mle.p, Param.Y = paramy.mle.p
              , Rt.diff = rt.diff.GEVP, Param.diff =param.diff.GEVP ))
}

######################
### Mirror Piest
######################
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
