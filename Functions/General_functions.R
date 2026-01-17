###################################################################################
# Functions
##################################################################################

#covariance for multivariate gaussian parameters
gp.cov <- function(a,b,c){
  gp.cov = a*exp(-(loc.dist/b)^c)
  return(gp.cov)
}


## Calculate return level using GEV
gev.return_level<-function(loc,scale,shape,n.site,time){
  p=1-1/time
  y<- rep(-log(p), n.site)
  ture.return <- loc - ifelse(abs(shape)>1e-5, scale*(1-y^(-shape))/shape, scale*log(y))
  return(ture.return)
}


## Negative log-likelihood function for spatial GPD 
gpd.neglog.l <- function(params){
  gpd.neglog.l <- 0
  shape <- params[4]
  scale <- params[1] + params[2]*locations[,1] + params[3]*locations[,2]
  if(sum(scale<0)>0) { gpd.neglog.l <- 1e10}
  else {
    for(i in 1:n.site){
      above <- (data[,i]-thresh[i])[data[,i]>thresh[i]]
      if (abs(shape)>1e-5){
        above1 <- above[1 + shape*above/scale[i] > 0]
        neglog.ll <- sum(log(scale[i]) + (1+ 1/shape)*log(1+shape*above1/scale[i]))
      }
      if (abs(shape)<=1e-5){
        neglog.ll <- sum(log(scale[i]) + above/scale[i])
      }
      gpd.neglog.l <-gpd.neglog.l+neglog.ll
    }
  }
  return(gpd.neglog.l)
}


## Calcuate return level using GPD
gpd.est_rl<-function(scale,shape,thresh,zeta,time){ #time = t-year return level
  gpd.est.return <- thresh + ifelse(abs(shape)>1e-5, (scale/shape)*((time*zeta)^shape -1), scale*log(time*zeta)) 
  return(gpd.est.return)
}
