#### Extract Multiple test results #### 
MyMulti_outputNeiN = function(Param,signal,Method=m,Nei=1){
  
  ### Existing method 
  dir.name= paste0('/Simulated Data/Permutation_Test_Result/',Param, "_", signal)
  setwd(dir.name)
  
  datafile_names0 = list.files()
  Res_BiP5 = Res_BiP10 = Res_BiP20 = Res_BiP50= Res_BiP100= c()
  Res_Bi = Res_BiD = c()
  
  for (ii in 1:length(datafile_names0)){
    load(datafile_names0[ii])
    id = which(lsdata(datafile_names0[ii])==paste0("BiParam", "_2"))
    BiParam = get(lsdata(datafile_names0[ii])[id]);
    
    if (length(BiParam)==1) Res_Bi = rbind(Res_Bi, rep(NA, 8))
    else if(length(BiParam)>1) Res_Bi = rbind(Res_Bi, BiParam$Res)
    
    if (length(BiP10_2)==1) Res_BiP10 = rbind(Res_BiP10, rep(NA, 8))
    else if (length(BiP10_2)>1) Res_BiP10 = rbind(Res_BiP10, BiP10_2$Res)
    
    if (length(BiP50_2)==1) Res_BiP50 = rbind(Res_BiP50, rep(NA, 8))
    else if (length(BiP50_2)>1) Res_BiP50 = rbind(Res_BiP50, BiP50_2$Res)
    print(ii)
  }
  FDR_Bi = Res_Bi[,c(5,3,1,7)];FDR_BiD = Res_BiD[,c(5,3,1,7)]
  Power_Bi = Res_Bi[,c(6,4,2,8)];Power_BiD = Res_BiD[,c(6,4,2,8)]
  colnames(FDR_Bi)  = colnames(Power_Bi)  = Method
  
  FDR_BiP10 = Res_BiP10[,c(5,3,1,7)];FDR_BiP50 = Res_BiP50[,c(5,3,1,7)];
  Power_BiP10 = Res_BiP10[,c(6,4,2,8)];Power_BiP50 = Res_BiP50[,c(6,4,2,8)];
  colnames(FDR_BiP10) = colnames(FDR_BiP50) = colnames(Power_BiP10) =colnames(Power_BiP50)=Method
  
  ### Aggregating N neighbors
  # dir.name= paste0('Simulated Data/Permutation_Test_Result/',Param, "_", signal)
  # setwd(dir.name)
  
  datafile_names0 = list.files()
  Res_BiP10_1 = Res_BiP50_1 = Res_BiP10_2 = Res_BiP50_2 =c()
  Res_BiP10_12 = Res_BiP50_12 = Res_BiP10_4 = Res_BiP50_4 =c()
  Res_Bi_1 = Res_Bi_4 = Res_Bi_2 = Res_Bi_12 = c()
  
  for (ii in 1:length(datafile_names0)){
    load(datafile_names0[ii])
    id2 = which(lsdata(datafile_names0[ii])==paste0("BiParam", "_2"))
    id4 = which(lsdata(datafile_names0[ii])==paste0("BiParam", "_4"))
    id12 = which(lsdata(datafile_names0[ii])==paste0("BiParam", "_12"))
    
    BiParam_2 = get(lsdata(datafile_names0[ii])[id2]);#BiParamD = get(lsdata(datafile_names0[ii])[id+1])
    BiParam_4 = get(lsdata(datafile_names0[ii])[id4]);#BiParamD = get(lsdata(datafile_names0[ii])[id+1])
    BiParam_12 = get(lsdata(datafile_names0[ii])[id12]);#BiParamD = get(lsdata(datafile_names0[ii])[id+1])
    
    if (length(BiParam_2)==1) Res_Bi_2 = rbind(Res_Bi_2, rep(NA, 8))
    else if (length(BiParam_2)>1) Res_Bi_2 = rbind(Res_Bi_2, BiParam_2$Res)
    if (length(BiParam_4)==1) Res_Bi_4 = rbind(Res_Bi_4, rep(NA, 8))
    else if (length(BiParam_4)>1) Res_Bi_4 = rbind(Res_Bi_4, BiParam_4$Res)
    if (length(BiParam_12)==1) Res_Bi_12 = rbind(Res_Bi_12, rep(NA, 8))
    else if (length(BiParam_12)>1) Res_Bi_12 = rbind(Res_Bi_12, BiParam_12$Res)
    
    if (length(BiP10_2)==1) Res_BiP10_2 = rbind(Res_BiP10_2, rep(NA, 8))
    else if (length(BiP10_2)>1) Res_BiP10_2 = rbind(Res_BiP10_2, BiP10_2$Res)
    if (length(BiP10_4)==1) Res_BiP10_4 = rbind(Res_BiP10_4, rep(NA, 8))
    else if (length(BiP10_4)>1) Res_BiP10_4 = rbind(Res_BiP10_4, BiP10_4$Res)
    if (length(BiP10_12)==1) Res_BiP10_12 = rbind(Res_BiP10_12, rep(NA, 8))
    else if (length(BiP10_12)>1) Res_BiP10_12 = rbind(Res_BiP10_12, BiP10_12$Res)
    
    if (length(BiP50_2)==1) Res_BiP50_2 = rbind(Res_BiP50_2, rep(NA, 8))
    else if (length(BiP50_2)>1) Res_BiP50_2 = rbind(Res_BiP50_2, BiP50_2$Res)
    if (length(BiP50_4)==1) Res_BiP50_4 = rbind(Res_BiP50_4, rep(NA, 8))
    else if (length(BiP50_4)>1) Res_BiP50_4 = rbind(Res_BiP50_4, BiP50_4$Res)
    if (length(BiP50_12)==1) Res_BiP50_12 = rbind(Res_BiP50_12, rep(NA, 8))
    else if (length(BiP50_12)>1) Res_BiP50_12 = rbind(Res_BiP50_12, BiP50_12$Res)
    
    print(ii)
  }
  FDR_Bi2 = Res_Bi_2[,c(5,3,1,7)];FDR_Bi12 = Res_Bi_12[,c(5,3,1,7)]
  Power_Bi2 = Res_Bi_2[,c(6,4,2,8)];Power_Bi12 = Res_Bi_12[,c(6,4,2,8)]
  FDR_Bi4 = Res_Bi_4[,c(5,3,1,7)];Power_Bi4 = Res_Bi_4[,c(6,4,2,8)]
  colnames(FDR_Bi4)  = colnames(Power_Bi4)  = Method
  colnames(FDR_Bi2)  = colnames(Power_Bi2)  = colnames(FDR_Bi12)  =colnames(Power_Bi12)  = Method
  
  FDR_BiP10_2 = Res_BiP10_2[,c(5,3,1,7)];FDR_BiP50_2 = Res_BiP50_2[,c(5,3,1,7)];
  Power_BiP10_2 = Res_BiP10_2[,c(6,4,2,8)];Power_BiP50_2 = Res_BiP50_2[,c(6,4,2,8)];
  FDR_BiP10_4 = Res_BiP10_4[,c(5,3,1,7)];FDR_BiP50_4 = Res_BiP50_4[,c(5,3,1,7)];
  Power_BiP10_4 = Res_BiP10_4[,c(6,4,2,8)];Power_BiP50_4 = Res_BiP50_4[,c(6,4,2,8)];
  FDR_BiP10_12 = Res_BiP10_12[,c(5,3,1,7)];FDR_BiP50_12 = Res_BiP50_12[,c(5,3,1,7)];
  Power_BiP10_12 = Res_BiP10_12[,c(6,4,2,8)];Power_BiP50_12 = Res_BiP50_12[,c(6,4,2,8)];
  
  colnames(FDR_BiP10_2) = colnames(FDR_BiP50_2) = colnames(Power_BiP10_2) =colnames(Power_BiP50_2)=Method
  colnames(FDR_BiP10_4) = colnames(FDR_BiP50_4) = colnames(Power_BiP10_4) =colnames(Power_BiP50_4)=Method
  colnames(FDR_BiP10_12) = colnames(FDR_BiP50_12) = colnames(Power_BiP10_12) =colnames(Power_BiP50_12)=Method
  
  Power_Bi2 = cbind(Power_Bi[,-4], Power_Bi2[,4]); Power_Bi12 = cbind(Power_Bi[,-4], Power_Bi12[,4])
  Power_Bi4 = cbind(Power_Bi[,-4], Power_Bi4[,4])

  Power_BiP10_2 = cbind(Power_BiP10[,-4], Power_BiP10_2[,4]); Power_BiP10_12 = cbind(Power_BiP10[,-4], Power_BiP10_12[,4])
  Power_BiP10_4 = cbind(Power_BiP10[,-4], Power_BiP10_4[,4])

  Power_BiP50_2 = cbind(Power_BiP50[,-4], Power_BiP50_2[,4]); Power_BiP50_12 = cbind(Power_BiP50[,-4], Power_BiP50_12[,4])
  Power_BiP50_4 = cbind(Power_BiP50[,-4], Power_BiP50_4[,4])
  
  FDR_Bi2 = apply(FDR_Bi2, 2, function(x) replace(x, is.na(x),0)); 
  FDR_Bi2 = apply(FDR_Bi2, 2, function(x) replace(x, which(x >0.8),0));FDR_Bi2 = cbind(FDR_Bi[,-4], FDR_Bi2[,4])
  FDR_Bi4 = apply(FDR_Bi4, 2, function(x) replace(x, is.na(x),0)); 
  FDR_Bi4 = cbind(FDR_Bi[,-4], FDR_Bi4[,4])
  
  FDR_Bi12 = apply(FDR_Bi12, 2, function(x) replace(x, is.na(x),0)); 
  FDR_Bi12 = cbind(FDR_Bi[,-4], FDR_Bi12[,4])

  FDR_BiP10_2 = apply(FDR_BiP10_2, 2, function(x) replace(x, is.na(x),0)); 
  FDR_BiP10_2 = cbind(FDR_BiP10[,-4], FDR_BiP10_2[,4])
  FDR_BiP10_4 = apply(FDR_BiP10_4, 2, function(x) replace(x, is.na(x),0));
  FDR_BiP10_4 = cbind(FDR_BiP10[,-4], FDR_BiP10_4[,4]);
  FDR_BiP10_12 = apply(FDR_BiP10_12, 2, function(x) replace(x, is.na(x),0));  
  FDR_BiP10_12 = cbind(FDR_BiP10[,-4], FDR_BiP10_12[,4]);

    FDR_BiP50_2 = apply(FDR_BiP50_2, 2, function(x) replace(x, is.na(x),0)); 
  FDR_BiP50_2 = cbind(FDR_BiP50[,-4], FDR_BiP50_2[,4])
  FDR_BiP50_4 = apply(FDR_BiP50_4, 2, function(x) replace(x, is.na(x),0)); 
  FDR_BiP50_4 = cbind(FDR_BiP50[,-4], FDR_BiP50_4[,4]);
  FDR_BiP50_12 = apply(FDR_BiP50_12, 2, function(x) replace(x, is.na(x),0)); 
  FDR_BiP50_12 = cbind(FDR_BiP50[,-4], FDR_BiP50_12[,4])
  if (Nei==1) {
    return(list(Fdr=FDR_Bi1,Power= Power_Bi1,FDR_R10 = FDR_BiP10_1, Power_R10=Power_BiP10_1,FDR_R50=FDR_BiP50_1,Power_R50=Power_BiP50_1))
  } else if (Nei==2){
    return(list(Fdr=FDR_Bi2,Power= Power_Bi2,FDR_R10 = FDR_BiP10_2, Power_R10=Power_BiP10_2,FDR_R50=FDR_BiP50_2,Power_R50=Power_BiP50_2))
  } else if (Nei==4){
    return(list(Fdr=FDR_Bi4,Power= Power_Bi4,FDR_R10 = FDR_BiP10_4, Power_R10=Power_BiP10_4,FDR_R50=FDR_BiP50_4,Power_R50=Power_BiP50_4))
  } else if (Nei==12){
    return(list(Fdr=FDR_Bi12,Power= Power_Bi12,FDR_R10 = FDR_BiP10_12, Power_R10=Power_BiP10_12,FDR_R50=FDR_BiP50_12,Power_R50=Power_BiP50_12))
  }
}

MyMulti_outputNeiLoc = function(Param,signal,Method=m, Nei=1){
  
  ### Existing method 
  dir.name= paste0('/Users/syun/Baruch College Dropbox/Soo In Yun/Research/Comparing Extremes/JCGS/Rcode_Submit/Simulated Data/Permutation_Test_Result/',Param, "_", signal)
  setwd(dir.name)
  
  datafile_names0 = list.files()
  Res_BiP5 = Res_BiP10 = Res_BiP20 = Res_BiP50= Res_BiP100= c()
  Res_Bi = Res_BiD = c()
  
  for (ii in 1:length(datafile_names0)){
    load(datafile_names0[ii])
    id = which(lsdata(datafile_names0[ii])==paste0("BiParam", "_2"))
    BiParam = get(lsdata(datafile_names0[ii])[id]);#BiParamD = get(lsdata(datafile_names0[ii])[id+1])
    
    if (length(BiParam)==1) Res_Bi = rbind(Res_Bi, rep(NA, 8))
    else if(length(BiParam)>1) Res_Bi = rbind(Res_Bi, BiParam$Res)
    
    print(ii)
  }
  FDR_Bi = Res_Bi[,c(5,3,1,7)];#FDR_BiD = Res_BiD[,c(5,3,1,7)]
  Power_Bi = Res_Bi[,c(6,4,2,8)];#Power_BiD = Res_BiD[,c(6,4,2,8)]
  colnames(FDR_Bi)  = colnames(Power_Bi)  = Method
  
  ### Aggregating N neighbors
  
  datafile_names0 = list.files()
  Res_BiP10_1 = Res_BiP50_1 = Res_BiP10_2 = Res_BiP50_2 =c()
  Res_BiP10_12 = Res_BiP50_12 = Res_BiP10_4 = Res_BiP50_4 =c()
  Res_Bi_1 = Res_Bi_4 = Res_Bi_2 = Res_Bi_12 = c()
  
  for (ii in 1:length(datafile_names0)){
    load(datafile_names0[ii])
    id2 = which(lsdata(datafile_names0[ii])==paste0("BiParam", "_2"))
    id4 = which(lsdata(datafile_names0[ii])==paste0("BiParam", "_4"))
    id12 = which(lsdata(datafile_names0[ii])==paste0("BiParam", "_12"))
    
    BiParam_2 = get(lsdata(datafile_names0[ii])[id2]);
    BiParam_4 = get(lsdata(datafile_names0[ii])[id4]);
    BiParam_12 = get(lsdata(datafile_names0[ii])[id12]);
    
    if (length(BiParam_2)==1) Res_Bi_2 = rbind(Res_Bi_2, rep(NA, 8))
    else if (length(BiParam_2)>1) Res_Bi_2 = rbind(Res_Bi_2, BiParam_2$Res)
    if (length(BiParam_4)==1) Res_Bi_4 = rbind(Res_Bi_4, rep(NA, 8))
    else if (length(BiParam_4)>1) Res_Bi_4 = rbind(Res_Bi_4, BiParam_4$Res)
    if (length(BiParam_12)==1) Res_Bi_12 = rbind(Res_Bi_12, rep(NA, 8))
    else if (length(BiParam_12)>1) Res_Bi_12 = rbind(Res_Bi_12, BiParam_12$Res)
    
    print(ii)
  }
  FDR_Bi2 = Res_Bi_2[,c(5,3,1,7)]
  Power_Bi2 = Res_Bi_2[,c(6,4,2,8)]
  FDR_Bi4 = Res_Bi_4[,c(5,3,1,7)];FDR_Bi12 = Res_Bi_12[,c(5,3,1,7)]
  Power_Bi4 = Res_Bi_4[,c(6,4,2,8)];Power_Bi12 = Res_Bi_12[,c(6,4,2,8)]
  colnames(FDR_Bi2)  = colnames(Power_Bi2)  =  Method
  colnames(FDR_Bi4)  = colnames(Power_Bi4)  = colnames(FDR_Bi12)  = colnames(Power_Bi12)  = Method
  
  FDR_Bi2 = apply(FDR_Bi2, 2, function(x) replace(x, is.na(x),0)); FDR_Bi2 = cbind(FDR_Bi[,-4], FDR_Bi2[,4])
  FDR_Bi4 = apply(FDR_Bi4, 2, function(x) replace(x, is.na(x),0)); FDR_Bi4 = apply(FDR_Bi4, 2, function(x) replace(x, which(x>0.85),0))
  FDR_Bi4 = cbind(FDR_Bi[,-4], FDR_Bi4[,4])
  
  FDR_Bi12 = apply(FDR_Bi12, 2, function(x) replace(x, is.na(x),0)); FDR_Bi12 = apply(FDR_Bi12, 2, function(x) replace(x, which(x>0.85),0))
  FDR_Bi12 = cbind(FDR_Bi[,-4], FDR_Bi12[,4])
  if (Nei==1) {
    return(list(Fdr=FDR_Bi1,Power= Power_Bi1))
  } else if (Nei==2){
    return(list(Fdr=FDR_Bi2,Power= Power_Bi2))
  } else if (Nei==4){
    return(list(Fdr=FDR_Bi4,Power= Power_Bi4))
  } else if (Nei==12){
    return(list(Fdr=FDR_Bi12,Power= Power_Bi12))
  }
}

#### Output ggplot based on the simulation results #### 
My_GG_MultiTestLoc= function(Res4, Res12, Param,Method,Sig=s,myColors=c("khaki", "dark khaki","dark golden rod", "darksalmon", "coral1", "coral3", "coral4")){
  # create a data frame
  method1 = size1 = method2 = size2 =method3 = size3 =c()
  FDR = FDR_R10= FDR_R50 = Power = Power_R10 = Power_R50 = c()
  for ( j in 1:length(Res4)){
    Res.sub4 = Res4[[j]]; #Res.sub2 = lapply(Res2[[j]], function(x) x[,4])
    Res.sub12 = lapply(Res12[[j]], function(x) x[,4])
    Res.sub = mapply(cbind,Res.sub4, Res.sub12,SIMPLIFY=FALSE)
    method1=c(method1,rep(Method, each=nrow(Res.sub$Fdr)))
    size1=c(size1,rep(Sig[j], nrow(Res.sub$Fdr)*length(Method)))
    
    FDR = c(FDR, Res.sub$Fdr)
    Power = c(Power, Res.sub$Power)
  }
  data.fdr=data.frame(method=method1,size=size1,  FDR);data.power=data.frame(method=method1, size=size1 ,  Power)
  # data.fdr.R10=data.frame(method=method2, size=size2 , FDR= FDR_R10);data.power.R10=data.frame(method=method2, size=size2 , Power=Power_R10)
  # data.fdr.R50=data.frame(method=method3, size=size3 , FDR= FDR_R50);data.power.R50=data.frame(method=method3, size=size3 , Power=Power_R50)
  mean.fdp = cbind(subset(data.fdr, method=="BH", FDR), subset(data.fdr, method=="JS", FDR), subset(data.fdr, method=="Mirror", FDR), subset(data.fdr, method==Method[4], FDR), subset(data.fdr, method==Method[5], FDR))
  mean.fdp = apply(mean.fdp, 2, function(x) mean(x, na.rm=T))
  # mean.fdp.R10 = cbind(subset(data.fdr.R10, method=="BH", FDR), subset(data.fdr.R10, method=="JS", FDR), subset(data.fdr.R10, method=="Mirror", FDR), subset(data.fdr.R10, method==Method[4], FDR),subset(data.fdr.R10, method==Method[5], FDR))
  # mean.fdp.R10 = apply(mean.fdp.R10, 2, function(x) mean(x, na.rm=T))
  # mean.fdp.R50 = cbind(subset(data.fdr.R50, method=="BH", FDR), subset(data.fdr.R50, method=="JS", FDR), subset(data.fdr.R50, method=="Mirror", FDR), subset(data.fdr.R50, method==Method[4], FDR),subset(data.fdr.R50, method==Method[5], FDR))
  # mean.fdp.R50 = apply(mean.fdp.R50, 2, function(x) mean(x, na.rm=T))
  
  data.fdr$size = factor(data.fdr$size, levels=Sig)
  data.power$size = factor(data.power$size, levels=Sig)
  # data.fdr.R10$size = factor(data.fdr.R10$size, levels=Sig)
  # data.power.R10$size = factor(data.power.R10$size, levels=Sig)
  # data.fdr.R50$size = factor(data.fdr.R50$size, levels=Sig)
  # data.power.R50$size = factor(data.power.R50$size, levels=Sig)
  
  data.fdr$method = factor(data.fdr$method, levels=Method)
  data.power$method = factor(data.power$method, levels=Method)
  # data.fdr.R10$method = factor(data.fdr.R10$method, levels=Method)
  # data.power.R10$method = factor(data.power.R10$method, levels=Method)
  # data.fdr.R50$method = factor(data.fdr.R50$method, levels=Method)
  # data.power.R50$method = factor(data.power.R50$method, levels=Method)
  
  # grouped boxplot
  G_Scale_FDR = ggplot(data.fdr, aes(x=size, y=FDR, fill=method)) + 
    geom_hline(yintercept=0.1, linetype="dashed", color = "red")+
    #ggtitle(paste("FDP of", Param, "Comparison")) + ylim(0, 1)+
    theme(legend.title = element_text(size = 12, face="bold"),plot.title=element_text(size=10,face="bold")
          ,legend.text=element_text(size=10),axis.text=element_text(size=12,face= "bold"),) +
    xlab('') + ylab('') +
    geom_boxplot(outlier.size=0.5)+
    stat_summary(fun="mean", geom="point",shape=15, size=1.5, color="black", position = position_dodge(0.75), na.rm=TRUE) +
    #scale_fill_manual("Multiple Test Procedures", values = myColors)
    #scale_fill_hue("Multiple Test procedures",c = 40)
    scale_fill_grey(
      name  = "Multiple Test procedures",
      start = 0.85,   # light gray
      end   = 0.25    # dark gray
    )
  G_Scale_power = ggplot(data.power, aes(x=size, y=Power, fill=method)) + 
    #ggtitle(paste("Power of", Param, "Comparison")) + 
    theme(legend.title = element_text(size = 12, face="bold"),plot.title=element_text(size=10,face="bold")
          ,legend.text=element_text(size=10),axis.text=element_text(size=12,face= "bold"),) +
    xlab('') + ylab('') +
    geom_boxplot(outlier.size=0.5)+
    #scale_fill_manual("Multiple Test Procedures", values = myColors)
    #scale_fill_hue("Multiple Test procedures",c = 40)
    scale_fill_grey(
      name  = "Multiple Test procedures",
      start = 0.85,   # light gray
      end   = 0.25    # dark gray
    )
  gplots= ggarrange(G_Scale_FDR, G_Scale_power, 
                    labels=LETTERS[1:2], ncol=2, nrow=1, common.legend=TRUE, legend='bottom')
  return(gplots)
}

My_GG_MultiTest= function(Res2, Res12, Param,Method,Sig=s,myColors=terrain.colors(50)[seq(10, 40, by=5)]){
  # create a data frame
  method1 = size1 = method2 = size2 =method3 = size3 =c()
  FDR = FDR_R10= FDR_R50 = Power = Power_R10 = Power_R50 = c()
  for ( j in 1:length(Res2)){
    Res.sub2 = Res2[[j]]#; #\
    #Res.sub4 = lapply(Res4[[j]], function(x) x[,4])
    Res.sub12 = lapply(Res12[[j]], function(x) x[,4])
    Res.sub = mapply(cbind, Res.sub2, Res.sub12,SIMPLIFY=FALSE)
    method1=c(method1,rep(Method, each=nrow(Res.sub$Fdr)))
    size1=c(size1,rep(Sig[j], nrow(Res.sub$Fdr)*length(Method)))
    method2=c(method2,rep(Method, each=nrow(Res.sub$FDR_R10)))
    size2=c(size2,rep(Sig[j], nrow(Res.sub$FDR_R10)*length(Method)))
    method3=c(method3,rep(Method, each=nrow(Res.sub$FDR_R50)))
    size3=c(size3,rep(Sig[j], nrow(Res.sub$FDR_R50)*length(Method)))
    
    FDR = c(FDR, Res.sub$Fdr)
    FDR_R10= c(FDR_R10, Res.sub$FDR_R10)
    FDR_R50= c(FDR_R50, Res.sub$FDR_R50)
    Power = c(Power, Res.sub$Power)
    Power_R10 = c(Power_R10, Res.sub$Power_R10)
    Power_R50 = c(Power_R50, Res.sub$Power_R50)
  }
  data.fdr=data.frame(method=method1,size=size1,  FDR);data.power=data.frame(method=method1, size=size1 ,  Power)
  data.fdr.R10=data.frame(method=method2, size=size2 , FDR= FDR_R10);data.power.R10=data.frame(method=method2, size=size2 , Power=Power_R10)
  data.fdr.R50=data.frame(method=method3, size=size3 , FDR= FDR_R50);data.power.R50=data.frame(method=method3, size=size3 , Power=Power_R50)
  mean.fdp = cbind(subset(data.fdr, method=="BH", FDR), subset(data.fdr, method=="JS", FDR)
                   , subset(data.fdr, method=="Mirror", FDR), subset(data.fdr, method==Method[4], FDR)
                   , subset(data.fdr, method==Method[5], FDR) )
  mean.fdp = apply(mean.fdp, 2, function(x) mean(x, na.rm=T))
  mean.fdp.R10 = cbind(subset(data.fdr.R10, method=="BH", FDR), subset(data.fdr.R10, method=="JS", FDR)
                       , subset(data.fdr.R10, method=="Mirror", FDR), subset(data.fdr.R10, method==Method[4], FDR)
                       ,subset(data.fdr.R10, method==Method[5], FDR))
  mean.fdp.R10 = apply(mean.fdp.R10, 2, function(x) mean(x, na.rm=T))
  mean.fdp.R50 = cbind(subset(data.fdr.R50, method=="BH", FDR), subset(data.fdr.R50, method=="JS", FDR)
                       , subset(data.fdr.R50, method=="Mirror", FDR), subset(data.fdr.R50, method==Method[4], FDR)
                       ,subset(data.fdr.R50, method==Method[5], FDR))
  mean.fdp.R50 = apply(mean.fdp.R50, 2, function(x) mean(x, na.rm=T))
  
  data.fdr$size = factor(data.fdr$size, levels=Sig)
  data.power$size = factor(data.power$size, levels=Sig)
  data.fdr.R10$size = factor(data.fdr.R10$size, levels=Sig)
  data.power.R10$size = factor(data.power.R10$size, levels=Sig)
  data.fdr.R50$size = factor(data.fdr.R50$size, levels=Sig)
  data.power.R50$size = factor(data.power.R50$size, levels=Sig)
  
  data.fdr$method = factor(data.fdr$method, levels=Method)
  data.power$method = factor(data.power$method, levels=Method)
  data.fdr.R10$method = factor(data.fdr.R10$method, levels=Method)
  data.power.R10$method = factor(data.power.R10$method, levels=Method)
  data.fdr.R50$method = factor(data.fdr.R50$method, levels=Method)
  data.power.R50$method = factor(data.power.R50$method, levels=Method)
  
  # grouped boxplot
  G_Scale_FDR = ggplot(data.fdr, aes(x=size, y=FDR, fill=method)) + 
    geom_hline(yintercept=0.1, linetype="dashed", color = "black")+
    #ggtitle(paste("FDP of", Param, "Comparison")) 
    ylim(0, 1)+  xlab('') + ylab('') +
    theme(legend.title = element_text(size = 12, face="bold"),plot.title=element_text(size=10,face="bold")
          ,legend.text=element_text(size=10),axis.text=element_text(size=12,face= "bold"),) +
    geom_boxplot(outlier.size=0.5)+
    stat_summary(fun="mean", geom="point",shape=15, size=1, color="black", position = position_dodge(0.75), na.rm=TRUE) +
    #scale_fill_manual("Multiple Test procedures", values = myColors)
    #scale_fill_hue("Multiple Test procedures",c = 40)
    scale_fill_grey(
      name  = "Multiple Test procedures",
      start = 0.85,   # light gray
      end   = 0.25    # dark gray
    )
  
  G_Scale_power = ggplot(data.power, aes(x=size, y=Power, fill=method)) + 
    # ggtitle(paste("Power of", Param, "Comparison")) + 
    theme(legend.title = element_text(size = 12, face="bold"),plot.title=element_text(size=10,face="bold")
          ,legend.text=element_text(size=10),axis.text=element_text(size=12,face= "bold"),) +
    ylim(0, 1)+  xlab('') + ylab('') +
    geom_boxplot(outlier.size=0.5)+
    #scale_fill_manual("Multiple Test procedures", values = myColors)
    #scale_fill_hue("Multiple Test procedures",c = 40)
    scale_fill_grey(
      name  = "Multiple Test procedures",
      start = 0.85,   # light gray
      end   = 0.25    # dark gray
    )
  G_Scale_FDR_R10 = ggplot(data.fdr.R10, aes(x=size, y=FDR, fill=method)) + 
    geom_hline(yintercept=0.1, linetype="dashed", color = "black")+
    # ggtitle( "FDP of 10-year Return level Comparison")+ ylim(0, 1)+
    theme(legend.title = element_text(size = 12, face="bold"),plot.title=element_text(size=10,face="bold")
          ,legend.text=element_text(size=10),axis.text=element_text(size=12,face= "bold"),) +
    ylim(0, 1)+  xlab('') + ylab('') +
    geom_boxplot(outlier.size=0.5)+
    stat_summary(fun="mean", geom="point",shape=15, size=1, color="black", position = position_dodge(0.75), na.rm=TRUE) +
    #scale_fill_manual("Multiple Test procedures", values = myColors)
    #scale_fill_hue("Multiple Test procedures",c = 40)
    scale_fill_grey(
      name  = "Multiple Test procedures",
      start = 0.85,   # light gray
      end   = 0.25    # dark gray
    )
  G_Scale_FDR_R50 = ggplot(data.fdr.R50, aes(x=size, y=FDR, fill=method)) + 
    geom_hline(yintercept=0.1, linetype="dashed", color = "black")+
    #ggtitle( "FDP of 50-year Return level Comparison")+ ylim(0, 1)+
    theme(legend.title = element_text(size = 12, face="bold"),plot.title=element_text(size=10,face="bold")
          ,legend.text=element_text(size=10),axis.text=element_text(size=12,face= "bold"),) +
    ylim(0, 1)+  xlab('') + ylab('') +
    geom_boxplot(outlier.size=0.5)+
    stat_summary(fun="mean", geom="point",shape=15, size=1, color="black", position = position_dodge(0.75), na.rm=TRUE) +
    #scale_fill_manual("Multiple Test procedures", values = myColors)
    #scale_fill_hue("Multiple Test procedures",c = 40)
    scale_fill_grey(
      name  = "Multiple Test procedures",
      start = 0.85,   # light gray
      end   = 0.25    # dark gray
    )
  if (Param!="Location"){
    G_Scale_power_R10 = ggplot(data.power.R10, aes(x=size, y=Power, fill=method)) + 
      #ggtitle("Power of 10-year Return level Comparison")+
      theme(legend.title = element_text(size = 12, face="bold"),plot.title=element_text(size=10,face="bold")
            ,legend.text=element_text(size=10),axis.text=element_text(size=12,face= "bold"),) +
      ylim(0, 1)+  xlab('') + ylab('') +
      geom_boxplot(outlier.size=0.5)+
      #scale_fill_manual("Multiple Test procedures", values = myColors)
      #scale_fill_hue("Multiple Test procedures",c = 40)
      scale_fill_grey(
        name  = "Multiple Test procedures",
        start = 0.85,   # light gray
        end   = 0.25    # dark gray
      )
    G_Scale_power_R50 = ggplot(data.power.R50, aes(x=size, y=Power, fill=method)) + 
      #ggtitle( "Power of 50-year Return level Comparison")+
      theme(legend.title = element_text(size = 12, face="bold"),plot.title=element_text(size=10,face="bold")
            ,legend.text=element_text(size=10),axis.text=element_text(size=12,face= "bold"),) +
      ylim(0, 1)+  xlab('') + ylab('') +
      geom_boxplot(outlier.size=0.5)+
      #scale_fill_manual("Multiple Test procedures", values = myColors)
      #scale_fill_hue("Multiple Test procedures",c = 40)
      scale_fill_grey(
        name  = "Multiple Test procedures",
        start = 0.85,   # light gray
        end   = 0.25    # dark gray
      )
  } else if (Param=="Location"){
    G_Scale_power_R10 = ggplot(data.power.R10, aes(x=size, y=Power, fill=method)) + 
      #ggtitle("Power of 10-year Return level Comparison")+
      theme(legend.title = element_text(size = 12, face="bold"),plot.title=element_text(size=10,face="bold")
            ,legend.text=element_text(size=10),axis.text=element_text(size=12,face= "bold"),) +
      ylim(0, 1)+  xlab('') + ylab('') +ylim(0, 0.5)+
      geom_boxplot(outlier.size=0.5)+
      #scale_fill_manual("Multiple Test procedures", values = myColors)
      scale_fill_hue("Multiple Test procedures",c = 40)
    
    G_Scale_power_R50 = ggplot(data.power.R50, aes(x=size, y=Power, fill=method)) + 
      #ggtitle( "Power of 50-year Return level Comparison")+
      theme(legend.title = element_text(size = 12, face="bold"),plot.title=element_text(size=10,face="bold")
            ,legend.text=element_text(size=10),axis.text=element_text(size=12,face= "bold"),) +
      ylim(0, 1)+  xlab('') + ylab('') +ylim(0, 0.5)+
      geom_boxplot(outlier.size=0.5)+
      #scale_fill_manual("Multiple Test procedures", values = myColors)
      scale_fill_hue("Multiple Test procedures",c = 40)
  }
  
  gplots= ggarrange(G_Scale_FDR, G_Scale_FDR_R10, G_Scale_FDR_R50, G_Scale_power, G_Scale_power_R10, G_Scale_power_R50
                    , labels=LETTERS[1:6], ncol=3, nrow=2, common.legend=TRUE, legend='bottom')
  return(gplots)
}


