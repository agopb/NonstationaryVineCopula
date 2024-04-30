
DynamicVineCopula_main<-function(data,family,matrix,par1,par2){
  
  data1<-data
  family1<-family
  matrix1<-matrix
  
  monthi<-monthi-2
  
  hyper1<-list(matrix(),matrix(),matrix())
  AIC1<-rep(0,3)
  BIC1<-rep(0,3)
  logLik1<-rep(0,3)
  
  mar1<-matrix1[3,monthi]
  
  for(jj in 1:2)
  {
    par02<-cbind(par1[jj,monthi],par2[jj,monthi])
    family02<-family1[jj,monthi]
    
    mar1<-matrix1[jj,monthi]
    
    data02<-cbind(data1[,mar1%/%10],data1[,mar1%%10])
    
    result1<-DynamicVineCopula(data02,par=par02,family=family02)
    hyper1[[jj]]<-result1$TVPDEP
    AIC1[jj]<-result1$AIC
    BIC1[jj]<-result1$BIC
    logLik1[jj]<-result1$logLik
    
  }
  
  jj=3
  par02<-cbind(par1[jj,monthi],par2[jj,monthi])
  family02<-family1[jj,monthi]
  
  mar1<-matrix1[jj,monthi]
  
  data_h1<-cbind(data1[,mar1%/%100],data1[,mar1%%10])
  data_h2<-cbind(data1[,(mar1%/%10)%%10],data1[,mar1%%10])
  
  datah10<-rep(0,nrow(data1))
  datah20<-rep(0,nrow(data1))
  
  hyper10<-hyper1[[1]]
  hyper20<-hyper1[[2]]
  
  for(kk in 1:nrow(data1)){
    
    if(family1[1,monthi]==2){
      par10<-hyper10[kk,1]
      par20<-hyper10[kk,2]
    }else{      
      par10<-hyper10[kk]
      par20<-0
    }
    datah10[kk]<-.C("Hfunc2_vec",                      # h(u1|u2)
                    as.integer(family1[1,monthi]),
                    as.integer(1),
                    as.double(data_h1[kk,1]),
                    as.double(data_h1[kk,2]),
                    as.double(par10),
                    as.double(par20),
                    as.double(rep(0, 1)),
                    PACKAGE = "VineCopula")[[7]]
  }
  
  for(kk in 1:nrow(data1)){
    
    if(family1[2,monthi]==2){
      par10<-hyper20[kk,1]
      par20<-hyper20[kk,2]
    }else{      
      par10<-hyper20[kk]
      par20<-0
    }
    datah20[kk]<-.C("Hfunc2_vec",                      # h(u1|u2)
                    as.integer(family1[2,monthi]),
                    as.integer(1),
                    as.double(data_h2[kk,1]),
                    as.double(data_h2[kk,2]),
                    as.double(par10),
                    as.double(par20),
                    as.double(rep(0, 1)),
                    PACKAGE = "VineCopula")[[7]]
  }
  
  data02<-cbind(datah10,datah20)
  
  
  result1<-DynamicVineCopula(data02,par=par02,family=family02)
  
  hyper1[[jj]]<-result1$TVPDEP
  AIC1[jj]<-result1$AIC
  BIC1[jj]<-result1$BIC
  logLik1[jj]<-result1$logLik
  
  result2<-list(TVPDEP=hyper1,AIC=AIC1,BIC=BIC1,logLik=logLik1)
  return(result2)
}

DynamicVineCopula<-function(data,par,family){
  if(family==1){
    print("gaussian")
    rhobar1<-par[1]
    theta01<-c(log((1+rhobar1)/(1-rhobar1)),0,0)
    
    x1<-dynamic_normal(data,theta01,rhobar1,plot=FALSE)
    tvtpdep1<-x1$tvtpdep
    aic1<-x1$AIC
    bic1<-x1$BIC
    logLik1<-x1$Loglikelihood
  }
  if(family==2){
    print("t")
    rhobar1<-par
    theta01<-c(log((1+rhobar1[1])/(1-rhobar1[1])),0,0,sqrt(rhobar1[2]-2),0,0)
    
    x1<-dynamic_T(data,theta01,rhobar1,plot=FALSE)
    tvtpdep1<-x1$tvtpdep
    aic1<-x1$AIC
    bic1<-x1$BIC
    logLik1<-x1$Loglikelihood
    
  }
  if(family==3){
    print("Clayton")
    rhobar1<-par[1]
    theta01<-c(sqrt(rhobar1),0,0)
    
    x1<-dynamic_clayton(data,theta01,rhobar1,plot=FALSE)
    tvtpdep1<-x1$tvtpdep
    aic1<-x1$AIC
    bic1<-x1$BIC
    logLik1<-x1$Loglikelihood
    
  }
  if(family==4){
    print("gumbel")
    rhobar1<-par[1]
    theta01<-c(sqrt(rhobar1-1),0,0)
    
    x1<-dynamic_gumbel(data,theta01,rhobar1,plot=FALSE)
    tvtpdep1<-x1$tvtpdep
    aic1<-x1$AIC
    bic1<-x1$BIC
    logLik1<-x1$Loglikelihood
    
  }
  if(family==5){
    print("frank")
    rhobar1<-par[1]
    theta01<-c(rhobar1,0,0)
    
    x1<-dynamic_frank(data,theta01,rhobar1,plot=FALSE)
    tvtpdep1<-x1$tvtpdep
    aic1<-x1$AIC
    bic1<-x1$BIC
    logLik1<-x1$Loglikelihood
    
  }
  
  result<-list(TVPDEP=tvtpdep1,AIC=aic1,BIC=bic1,logLik=logLik1)
  return(result)
}


DynamicVineCopula_Hfunction<-function(data1,matrix1,family1,i,monthi){
  mar1<-matrix1[3,monthi]
  
  data_h1<-cbind(data1[,mar1%/%100],data1[,mar1%%10])
  data_h2<-cbind(data1[,(mar1%/%10)%%10],data1[,mar1%%10])
  
  jj=1
  family02<-family1[jj,monthi]
  path1<-paste("E:\\phd time\\paper5\\data\\DynamicVineCopula_hyper\\TVPDEP-",t1$lon[i],"-",t1$lat[i],"-",monthi,"-",jj,".txt",sep="")
  hyper10<-read.table(path1)
  datah10<-rep(0,nrow(hyper10))
  
  for(kk in 1:nrow(hyper10)){
    
    if(family1[1,monthi]==2){
      par10<-hyper10[kk,1]
      par20<-hyper10[kk,2]
    }else{      
      par10<-hyper10[kk]
      par20<-0
    }
    datah10[kk]<-.C("Hfunc2_vec",                      # h(u1|u2)
                    as.integer(family1[jj,monthi]),
                    as.integer(1),
                    as.double(data_h1[1]),
                    as.double(data_h1[2]),
                    as.double(par10),
                    as.double(par20),
                    as.double(rep(0, 1)),
                    PACKAGE = "VineCopula")[[7]]
  }
  
  jj=2
  family02<-family1[jj,monthi]
  path1<-paste("E:\\phd time\\paper5\\data\\DynamicVineCopula_hyper\\TVPDEP-",t1$lon[i],"-",t1$lat[i],"-",monthi,"-",jj,".txt",sep="")
  hyper20<-read.table(path1)
  datah20<-rep(0,nrow(hyper20))
  
  for(kk in 1:nrow(hyper20)){
    
    if(family1[1,monthi]==2){
      par10<-hyper20[kk,1]
      par20<-hyper20[kk,2]
    }else{      
      par10<-hyper20[kk]
      par20<-0
    }
    datah20[kk]<-.C("Hfunc2_vec",                      # h(u1|u2)
                    as.integer(family1[jj,monthi]),
                    as.integer(1),
                    as.double(data_h2[1]),
                    as.double(data_h2[2]),
                    as.double(par10),
                    as.double(par20),
                    as.double(rep(0, 1)),
                    PACKAGE = "VineCopula")[[7]]
  }
  
  jj=3
  family02<-family1[jj,monthi]
  path1<-paste("E:\\phd time\\paper5\\data\\DynamicVineCopula_hyper\\TVPDEP-",t1$lon[i],"-",t1$lat[i],"-",monthi,"-",jj,".txt",sep="")
  hyper30<-read.table(path1)
  datah30<-rep(0,nrow(hyper30))
  
  if((mar1%/%10)%%10==1){data_h3<-cbind(datah20,datah10)
  }else if((mar1%/%10)%%10==2|(mar1%/%10)%%10==3){data_h3<-cbind(datah10,datah20)
  }
  
  for(kk in 1:nrow(hyper30)){
    
    if(family1[1,monthi]==2){
      par10<-hyper30[kk,1]
      par20<-hyper30[kk,2]
    }else{      
      par10<-hyper30[kk]
      par20<-0
    }
    datah30[kk]<-.C("Hfunc2_vec",                      # h(u1|u2)
                    as.integer(family1[jj,monthi]),
                    as.integer(1),
                    as.double(data_h3[kk,1]),
                    as.double(data_h3[kk,2]),
                    as.double(par10),
                    as.double(par20),
                    as.double(rep(0, 1)),
                    PACKAGE = "VineCopula")[[7]]
  }
  
  return(data30)
  
}