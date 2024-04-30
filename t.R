
T_pdf=function(u,v,k1,k2){
  pdf <- .C("PDF_seperate_vec",
            as.integer(2),
            as.integer(1),
            as.double(u),
            as.double(v),
            as.double(k1),
            as.double(k2),
            as.double(rep(0, 1)),
            PACKAGE = "VineCopula")[[7]]
  return(pdf)
}

T_tvtp=function(theta,data,rhobar,Call){
  
  # print(Call)

  x = qstd(data[,1])
  y = qstd(data[,2])
  n=length(x)
  T = nrow(data)
  
  u = data[,1]
  v = data[,2]
  w1=theta[1]
  a1=theta[2]
  b1=theta[3]
  
  w2=theta[4]
  a2=theta[5]
  b2=theta[6]
  
  # print(theta)
  kappa =matrix(-999,n,2)

  kappa[1,1] = rhobar[1]
  kappa[1,2] = rhobar[2]
  

  for (jj in  2:n){
    if (jj<=10){
      psi1 = w1+ a1*kappa[jj-1,1] + b1*(mean((x[1:(jj-1)]*y[1:(jj-1)])))
      psi2= w2+ a2*kappa[jj-1,2] + b2*(mean((x[1:(jj-1)]*y[1:(jj-1)])))
      
      
    }else{
      psi1 = w1+ (a1*kappa[jj-1,1]) + b1*(mean((x[(jj-10):(jj-1)]*y[(jj-10):(jj-1)])))
      psi2 = w2+ (a2*kappa[jj-1,2]) + b2*(mean((x[(jj-10):(jj-1)]*y[(jj-10):(jj-1)])))
    }
    kappa[jj,1] = 1.998/(1+exp(-psi1))-0.999;		# a modified logistic transformation
    kappa[jj,2] = 2.001+psi2*psi2
    
  }
  
  rhohat = kappa;  # time-path of conditional copula parameter
  
  
  n=T
  
  if(sum(is.infinite(rhohat))){
    CL<- -n*100
  }else if(sum(is.nan(rhohat))){
    CL<- -n*100
  }else{
    CL=rep(0,T)
    for (i in 1:T){
      CL[i]=T_pdf(u[i],v[i],rhohat[i,1],rhohat[i,2])
    }
    CL = sum(log(CL))
  }
  
  if (is.infinite(CL))  # control for optimization
    CL<--n*100
  if (is.nan(CL))  # control for optimization
    CL<--n*100
  
  # print(CL)

  CL1=-CL
  
  if (Call=="optim")
    return(CL1)
  
  if (Call=="filtering")
  {
    specOut<-list(LL=CL1 ,beta=c(theta),TVTP=rhohat)
    
    return(specOut)
  }
}

dynamic_T=function(data,theta,rhobar1, plot){
  lower =rep(-5,6)
  upper =rep(5,6)

  model<-fmincon(theta,T_tvtp,data=data,rhobar=rhobar1,Call="optim",tol = 1e-03,maxiter=100000)
  
  model$value<- -model$value
  
  # table of results
  n=nrow(data)
  coef<- model$par
  BIC= -2*model$value+ (log(n)*length(coef))
  AIC = -2*model$value + 2*length(coef)
  
  rhot=T_tvtp(coef,data=data,rhobar=rhobar1,Call="filtering")$TVTP
  # if ( plot==TRUE){
  #   plot(ts( rhot[5:n]))}
  
  output=list(
    AIC=AIC,
    BIC=BIC,
    Loglikelihood=model$value,
    tvtpdep=rhot
  )
  output
}
