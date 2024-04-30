
frank_pdf=function(u,v,k1){
  pdf <- .C("PDF_seperate_vec",
            as.integer(5),
            as.integer(1),
            as.double(u),
            as.double(v),
            as.double(k1),
            as.double(0),
            as.double(rep(0, 1)),
            PACKAGE = "VineCopula")[[7]]
  return(pdf)
}

frank_tvtp=function(theta,data,rhobar,Call){
  
  # print(Call)
  # 
  # print(theta)
  u = data[,1]
  v = data[,2]
  x=u
  y=v
  T = nrow(data)
  
  w=theta[1]
  a=theta[2]
  b=theta[3]
  
  kappa =rep(-0.99,T)
  kappa[1] = rhobar;			# this is the MLE of kappa in the time-invariant version of this model
  
  for (jj in  2:T){
    if (jj<=10){
      psi1 = w+ a*kappa[jj-1] + b*(mean(abs(x[1:jj-1]-y[1:jj-1])))
    }else{
      psi1 = w+ (a*kappa[jj-1]) + b*(mean(abs(x[(jj-10):(jj-1)]-y[(jj-10):(jj-1)])))
    }
    kappa[jj] = psi1;		# a modified logistic transformation
  }
  
  rhohat = kappa;  # time-path of conditional copula parameter
  
  n=T
  
  if(sum(is.infinite(rhohat))){
    CL<- -n*100
    # print(CL)
  }else if(sum(is.nan(rhohat))){
    CL<- -n*100
    # print(CL)
  }else{
    CL=rep(0,T)
    for (i in 1:T){
      CL[i]=frank_pdf(u[i],v[i],rhohat[i])
    }
    CL = sum(log(CL))
  }
  
  n=T
  if (is.infinite(CL))  # control for optimization
    CL<--n*100
  if (is.nan(CL))  # control for optimization
    CL<--n*100
  
  # cat("Sum of log Likelihood for normdynamic ->",sprintf("%4.4f",c(CL,kappa[n])),"\n")
  
  if (Call=="optim")
    return(CL)
  
  if (Call=="filtering")
  {
    specOut<-list(LL=CL ,beta=c(theta),TVTP=rhohat)
    
    return(specOut)
  }
}

dynamic_frank=function(data,theta,rhobar1, plot){
  lower =rep(-5,3)
  upper =rep(5,3)
  
  model <- optim(theta,fn=frank_tvtp,data=data,rhobar=rhobar1, Call="optim",
                 control = list(maxit=100000,fnscale=-1),method="L-BFGS-B",
                 lower =lower,upper =upper, hessian=TRUE)
  
  # table of results
  n=nrow(data)
  coef<- model$par
  model$se <- sqrt(-diag(solve(model$hessian)))
  S.E.= as.vector(model$se)
  (paramsWithTs = cbind (model$par , model$par/model$se ) )
  stat=coef/S.E.
  pvalue <- 2*(1 - pnorm(abs(stat)))
  confident=matrix(0,length(coef),2)
  for (i in 1:length(coef)){
    confident[i,]=coef[i]+c(-1,1)*S.E.[i]*qt(0.975,60)}
  
  result <- rbind(coef,S.E.,stat,pvalue)
  BIC= -2*model$value+ (log(n)*length(coef))
  AIC = -2*model$value + 2*length(coef)
  
  rhot=frank_tvtp(coef,data=data,rhobar=rhobar1,Call="filtering")$TVTP
  if ( plot==TRUE){
    plot(ts( rhot[5:n]))}
  
  output=list(
    result=result,
    AIC=AIC,
    BIC=BIC,
    Loglikelihood=model$value,
    tvtpdep=rhot
  )
  output
}
