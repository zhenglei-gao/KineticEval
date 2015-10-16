##' Standard error estimator for transformed parameters
##' 
##' This function use Monte Carlo simulations to 
##' calculate the standard error of the transformed parameters.
##'@param mu mean vector.
##'@param sigma sd vector.
##'@param transformation a vecotr of transformation functions
##'@param N default 10,000, number of samples.
##'@return standard error of the transformed parameters.
##'@export
MCstd <- function(mu,sigma,transformations,N=10000,fixed){
 ## f <- ForwardCalcFF(FF) ## the starting values that will be used in the optimization program so that the optimization does not need to deal with the sum ff==1 contraint.###
 ## sigma <- summary(Fit)$par[,2]
 ## mu <- Fit$par
 ## transformations <- mkinmodini$ff
 ## fixed<- Fit$fixed[,1]
  ## names(fixed) <- rownames(Fit$fixed)
if(all(is.na(sigma))){
  ## Fix issue von Simon #36
  v_names <- names(mu)
  ff_names = names(transformations)
  ff <- NULL
  
  fftrans <- function(v,fixed){
    ff <- NULL
    if(is.null(names(v))) names(v)<-v_names
    for (ff_name in ff_names)
    {
      ff[[ff_name]] =
        eval(parse(text = transformations[ff_name]), as.list(c(v,fixed)))
    }
    return(ff)
  }
  y0 <- fftrans(mu,fixed)
  res <- cbind(y0,NA,NA,NA)
  colnames(res) <- c("Estimate", "Std. Error",'Lower CI','Upper CI')
}else{
  mat <- replicate(N,rnorm(length(mu),mu,sigma))
  ## mat has N columns.
  v_names <- names(mu)
  ff_names = names(transformations)
  ff <- NULL
 
  fftrans <- function(v,fixed){
    ff <- NULL
    if(is.null(names(v))) names(v)<-v_names
    for (ff_name in ff_names)
    {
      ff[[ff_name]] =
        eval(parse(text = transformations[ff_name]), as.list(c(v,fixed)))
    }
    return(ff)
  }
  y <- apply(mat,2,fftrans,fixed=fixed)
  y0 <- fftrans(mu,fixed)
  #pvals <- sapply(1:length(y0),function(k){
  #  if(length(unique(y[k,]))==1) return(0) else{
      ## t.test(y[k,],mu=0,alternative="greater")$p.value
      ## cannot use this since there are too many samples.
      ## sum(y[k,]>0)/N
  #  }    
  #})
  ##
  if(is.null(dim(y))){
    #browser()
    res <- cbind(y0,sd(y),t(quantile(y,c(0.025,0.975))))
  }else res <- cbind(y0,apply(y,1,sd),t(apply(y,1,function(z)quantile(z,c(0.025,0.975))))
               )
  colnames(res) <- c("Estimate", "Std. Error",'Lower CI','Upper CI')#,
                 #"Pr(>t)")
}
  return(res)
}




##' Standard error estimator for transformed parameters
##' 
##' This function use Monte Carlo simulations to 
##' calculate the standard error of the transformed parameters.
##'@param mu mean vector.
##'@param sigma sd vector.
##'@param transformation a vecotr of transformation functions
##'@param N default 10,000, number of samples.
##'@return standard error of the transformed parameters.
##'@export
MCMCstd <- function(mu_mcmc,transformations,fixed){
  ## f <- ForwardCalcFF(FF) ## the starting values that will be used in the optimization program so that the optimization does not need to deal with the sum ff==1 contraint.###
  ## sigma <- summary(Fit)$par[,2]
  ## mu <- Fit$par
  ## transformations <- mkinmodini$ff
  ## fixed<- Fit$fixed[,1]
  ## names(fixed) <- rownames(Fit$fixed)
  
  
  ## mat has N columns.
  v_names <- colnames(mu_mcmc)
  ff_names = names(transformations)
  ff <- NULL
  
  fftrans <- function(v,fixed){
    ff <- NULL
    if(is.null(names(v))) names(v)<-v_names
    for (ff_name in ff_names)
    {
      ff[[ff_name]] =
        eval(parse(text = transformations[ff_name]), as.list(c(v,fixed)))
    }
    return(ff)
  }
  y <- apply(mu_mcmc,1,fftrans,fixed=fixed)
  
  if(is.null(dim(y))){
    #browser()
    res <- cbind(mean(y),sd(y),t(quantile(y,c(0.025,0.975))))
    rownames(res)<- ff_names[1]
  }else res <- cbind(apply(y,1,mean),apply(y,1,sd),t(apply(y,1,function(z)quantile(z,c(0.025,0.975)))))
  colnames(res) <- c("Estimate", "Std. Error",'Lower CI','Upper CI')#,
  #"Pr(>t)")
  return(res)
}