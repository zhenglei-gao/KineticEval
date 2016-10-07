##' Fit a kinetic model using the NLS, IRLS and MCMC methods.
##'
##' Instead of implicitly assuming equal error variances or giving arbitrary weights decided by the researcher as in the NLS algorithm,  an iteratively reweighted least squares (IRLS) algorithm was implemented to obtain the maximum likelihood estimates of the kinetic model parameters.
##' @title Fit a kinetic model using the IRLS algorithm.
##' @param mkinmodini A list of class \code{\link{mkinmod.full}}, containing the kinetic model to be fitted to the data, and the initial parameter values, the observed data.
##' @param evalMethod Choice of evaluation method, NLS, IRLS, MCMC, directMLE(not implemented)
##' @param optimMethod Choise of optimization algorithm, the user should choose from LM and TRR.
##' @param eigen If TRUE,  the solution of the system of differential equation should be based on the spectral decomposition of the coefficient matrix in cases that this is possible.
##' @param plotfit If TRUE,the observed values and the numerical solutions should be plotted at each stage of the optimisation.
##' @param plotRes If TRUE,the observed values and the numerical solutions should be plotted only at the end of the optimisation.
##' @param plottitle The title of the plot for visualizing the optimization process.
##' @param quiet If TRUE, suppress printing out the current model cost after each(>1) improvement.
##' @param ctr A list of control values for the estimation algorithm to replace the default values including maximum iterations and absolute error tolerance.  Defaults to the output of \code{\link{kingui.control}}.
##' @param irls.control A list of control values for the estimation algorithm to replace the default values including the maximum number of iterations for the outer iteration and the error tolerance level for the error variance estimation updating.
##' @param MCMCoptions control parameters for MCMC evaluation method
##' @param update If not NULL, should be a list of starting values obtained from other optimization methods.
##' @param useHsolnp Whether to use the hessian matrix derived from the solnp optimization algorithm.
##' @param ... Further arguments that will be passed to \code{\link{kin_mod}} or \code{IRLSkinfit.full}.
##' @return A list with  "kingui" in the class attribute. A summary can be obtained by \code{\link{summary.kingui}}.
##' @author Zhenglei Gao
##' @examples
##' complex <- mkinmod.full(
##'  parent = list(type = "SFO", to = c("A1", "B1", "C1"), sink = FALSE),
##'  A1 = list(type = "SFO", to = "A2"),
##'  B1 = list(type = "SFO"),
##'  C1 = list(type = "SFO"),
##'  A2 = list(type = "SFO"),
##'  inpartri='default',
##'  outpartri='default',
##'  data=schaefer07_complex_case,
##'  weight=NULL)
##' Fit    <- KinEval(
##'            complex,
##'            evalMethod="IRLS",
##'            optimMethod="LM",
##'               plotfit      = TRUE,
##'               quiet     = TRUE,
##'               ctr       = kingui.control(),
##'            irls.control = list(
##'                              maxIter = 10,
##'                            tolerance = 0.001))
##' @keywords Kinetic-Evaluations
##' @import FME Rsolnp optimx BB minqa deSolve coda minpack.lm numDeriv
##' @note adapted from \code{mkinfit} and \code{}
##' @export
KinEval <- function(mkinmodini,
                    evalMethod=c('NLLS','IRLS','MCMC','directMLE'),
                    optimMethod=c("LM","TRR","Byrd","GENOUD","bobyqa","port","nls2","Nash","Marq", "Port", "Newton", "Nelder-Mead", "BFGS",
                                  "CG","L-BFGS-B", "SANN", "Pseudo",'trust',
                                  'spg','ucminf','nmk','Rcgmin','Rvmmin','deoptim','solnp'),
                    eigen = FALSE,
                    plotfit= TRUE, plotRes = FALSE, plottitle='',quiet = FALSE,
                    ctr=kingui.control(),irls.control=list(), MCMCoptions=MCMCcontrol(),
                    update=NULL,useHsolnp=FALSE,
                    ...)
{
  ## Must be simple and clear this time. :)
  
  #####################
  #### Control parameters ####
  if("logging" %in% loadedNamespaces()){
    logall <- TRUE
  }else{
    logall <- FALSE
  }
  method <- ctr$method
  odesolver <- ctr$odesolver
  atol <- ctr$atol
  rtol <- ctr$rtol
  control <- ctr$control
  marqctr <- ctr$marqctr
  goMarq <- ctr$goMarq
  submethod <- ctr$submethod
  Hmethod1 <- ctr$Hmethod1
  Hmethod2 <- ctr$Hmethod2
  runTRR <- ctr$runTRR
  if(ctr$runTRR4==FALSE) {
    runTRR4 <- (length(mkinmodini$map)<5)
    runTRR <- runTRR4
    if(logall & !runTRR) loginfo("Number of compartments greater or equal to 5, turn off the option to run Trust Regrion Algorithm!")
  }
  calcJacob <- ctr$calcJacob
  if(calcJacob){
    exactHessian<-FALSE
    exactJacob<-TRUE
  }else{
    exactHessian<-FALSE
    exactJacob<-FALSE
  }
  ## -----------------------------------
  ## Get the parametrization.
  inpartri <- mkinmodini$inpartri
  outpartri <- mkinmodini$outpartri
  ##
  
  ## mkinmodini is an object by mkinmod.full
  parms.ini <- mkinmodini$parms.ini
  state.ini <- mkinmodini$state.ini
  state.ini.orig <- mkinmodini$state.ini.orig
  lower <- mkinmodini$lower
  upper <- mkinmodini$upper
  fixed_parms <- mkinmodini$fixed_parms
  fixed_initials <- mkinmodini$fixed_initials
  mod_vars <- names(mkinmodini$diffs)
  observed <-  mkin_wide_to_long(mkinmodini$residue,time='time')
  observed$err <-c(as.matrix(mkinmodini$weightmat))
  ## Subset dataframe with mapped (modelled) variables
  observed <- subset(observed, name %in% names(mkinmodini$map))
  ## Get names of observed variables
  ## NOTE HERE: the order may not be the same as the input mkinmod.full differential equations list. 
  ## XXXXX TODO XXXX Reorder them maybe a good idea if the data is given from a data file while the mkinmod.full is defined not following the colnames order, 
  ## although it is already taken care of in the cost(P) function to reorder the odeini using mod_vars
  obs_vars = unique(as.character(observed$name))
  
  
  ## Name the parameters if they are not named yet ## usually they are already names
  if(is.null(names(parms.ini))) names(parms.ini) <- mkinmodini$parms
  
  ## Name the inital parameter values if they are not named yet
  if(is.null(names(state.ini))) names(state.ini) <- mod_vars
  
  ## Parameters to be optimised
  parms.fixed <- parms.ini[fixed_parms]
  optim_parms <- setdiff(names(parms.ini), fixed_parms)
  parms.optim <- parms.ini[optim_parms]
  
  ## # ### ### ### ### ###
  state.ini.fixed <- state.ini[fixed_initials]
  optim_initials <- setdiff(names(state.ini), fixed_initials)
  state.ini.optim <- state.ini[optim_initials]
  if(is.null(state.ini.orig)) state.ini.orig <- state.ini
  state.ini.orig.optim <- state.ini.orig[optim_initials]
  state.ini.optim.boxnames <- names(state.ini.optim)
  state.ini.fixed.boxnames <- names(state.ini.fixed)
  if(length(state.ini.optim) > 0) {
    names(state.ini.optim) <- paste('M0',names(state.ini.optim),  sep="_")
    names(state.ini.orig.optim) <- names(state.ini.optim)
    }
  if(length(state.ini.fixed) > 0) {
    names(state.ini.fixed) <- paste('M0',names(state.ini.fixed), sep="_")
  }
  
  eigen <- FALSE
  oldparms <- c(state.ini.optim,parms.optim)
  if (length(mkinmodini$map) == 1) {
    solution = "analytical"
  } else {
    if (is.matrix(mkinmodini$coefmat) & eigen) solution = "eigen"
    else solution = "deSolve"
  }
  ## always define mkindiff function since most of the time we will use it.
  mkindiff <- function(t, state, parms) {
    time <- t
    diffs <- vector()
    for (box in mod_vars)
    {
      diffname <- paste("d", box, sep="_")
      diffs[diffname] <- with(as.list(c(time,state, parms)),
                              eval(parse(text=mkinmodini$diffs[[box]])))
    }
    ##https://stat.ethz.ch/pipermail/r-sig-dynamic-models/2010q2/000031.html
    #bady <- (!is.finite(diffs))|(diffs<=0)
    #if(logall & sum(bady)>0) logwarn("NaN values simulated from solving the ODE, which were replaced with 0s.")
    #diffs[bady] <- 0 
    return(list(c(diffs)))
  }
  
  
  # -----------------------------------
  ## Get the evaluation and optimization method.
  evalMethod <- match.arg(evalMethod)
  optimMethod <- match.arg(optimMethod)
  runIRLS <- TRUE ## determin whther to run an IRLS step especially before the MCMC step
  ## if the parameter values are provided to the MCMC, then do not run
  if(evalMethod == 'MCMC' & (!is.null(update))) runIRLS <- FALSE
  if(evalMethod == "NLLS") runIRLS <- FALSE
  # ----------------------------------------
  # options(warn=-1) ## turn off the warning
  # ----------------------------------------
  
  if(length(parms.optim)+length(state.ini.optim)==0){
    ## no need to run optimization!
    fit <- mkinmodini
    fit$par <- c(state.ini.optim,parms.optim)
    environment(kin_mod_cost) <- environment()
    costout <- kin_mod_cost(P=c(state.ini.optim,parms.optim))
    
    fit$costout <- costout
    fit$out_predicted <- out_predicted
    class(fit) <- "mkinmod.full"
  }else{
    # -----------------------------------------
    ## Run the Optimization ##
    ## Always run a first step optimization.
    oldparms <- c(state.ini.optim,parms.optim)
    y <- observed$value
    pnames <- names(oldparms)
    environment(kin_mod) <- environment()
    environment(kin_mod_nls) <- environment()
    
    if(optimMethod %in% c("Marq", "Port", "Newton", "Nelder-Mead", "BFGS",
                          "CG","L-BFGS-B", "SANN", "Pseudo",'trust',
                          'spg','ucminf','nmk','Rcgmin','Rvmmin','deoptim','solnp')){
      ## Use the old IRLSkinfit, mkinfit, mcmckinfit
      if(evalMethod=="NLLS"){
        fit <- mkinfit.full(mkinmodini,
                            eigen = eigen,
                            plot = plotfit, plottitle=plottitle,quiet = quiet,
                            err = NULL, weight = "none", scaleVar = FALSE, ctr=ctr, ...)
      }else{
        if(evalMethod=="IRLS"){
          fit <- IRLSkinfit.full(mkinmodini,
                                 eigen = eigen,
                                 plot = plotfit, plottitle=plottitle,quiet = quiet,
                                 err = NULL, weight = "none", scaleVar = FALSE,
                                 ctr=ctr,irls.control=irls.control,
                                 update=update,useHsolnp=FALSE,...)  
        }else{
          attach(MCMCoptions)
          #fit <- mcmckinfit.full()
          stop("Please use either LM or TRR as the optimMethod for MCMC evalMethod")
        }
      }
      
      
    }else{
      cost.old <- 1e100
      calls <- 0
      out_predicted <- NA
      
      if(optimMethod=="port"){
        stop("Not fully implemented yet!!")
        res0 <- nls(y ~ kin_mod_nls(P,pnames=names(oldparms)),start=list(P=oldparms),lower=lower,upper=upper,algorithm="port")
        
      }
      
      if(plotfit==TRUE) x11()
      f <- function(P,plotfit=FALSE,plotRes=FALSE,...){
        temp <- kin_mod(P,inside=TRUE,...)
        
        if(plotfit==TRUE) res <- observed$value-temp else{
          if(plotRes==TRUE) {
            res <- observed$value-kin_mod(P,inside=FALSE,plot=TRUE,plottitle=plottitle,...)
          }else res <- observed$value-kin_mod(P,inside=FALSE,plot=FALSE,...)
        }
        id <- which(is.na(res))
        if(length(id)>0) return(res[-id]) else return(res)
      }
      f2<-function(P,plotfit=FALSE,plotRes=FALSE,...){
        sum(f(P,plotfit=plotfit,plotRes=plotRes,...)^2)
      }
      
      if(optimMethod=="LM") {
        res0 <- nls.lm(par=oldparms,lower=lower,upper=upper,fn=f,plotfit=plotfit,plotRes=plotRes,...)
        ## res0 <- nls.lm(par=oldparms,lower=lower,upper=upper,fn=f)
        ##browser()
        # renaming results for compatibility with other methods
        names(res0)[7] <- "iterations"  # called "niter" here
        names(res0)[9] <- "ssr"         # called "deviance" here
        ##names(res0)[3] <- "residuals"
        res0$residuals <- res0$fvec
        res0$hessian   <- 2*res0$hessian # returns 0.5*hessian!
        Diag          <- unlist(res0[6])
        res0$diag      <- NULL
        res0$diag      <- Diag ### ? What is diag here?
        ## if any on the boundary, we should use a trust region method
        if(atBoundary(res0$par,lower,upper)){
          if(runTRR){
            if(logall==TRUE){
              loginfo("Some paramter estimates are at boundary, running an STIR to step aside the boundary problem.")
            }
            print("It is not freezing, running a STIR to step aside the boundary problem. Please wait a while.")
            res00 <- res0
            res0 <- try(lsqnonlin(f,xstart=res0$par,l=lower,u=upper,plotfit=plotfit,plotRes=FALSE,...),silent=FALSE)
            ##browser()
            if(class(res0)=="try-error"){
              res0 <- res00
              rm(res00)
              if(logall==TRUE){
                logwarn("STIR step failed. Please report this case.")
                
                loginfo("Restore the results from LM optimization algorithm.")
              }
            }else{
              if(logall==TRUE){
                loginfo("STIR step succeeded.")
              }
              res0$par <- as.vector(res0$xcurr)
              names(res0$par) <- pnames
              res0$residuals <- res0$fvec
              res0$hessian <- 2*t(res0$JACOB)%*%(res0$JACOB)
            }
          }else{
            warning("Some Parameters on the boundary, please check the results 
                    and make sure they are the same as using other optimization algorithms 
                    like TRR!")
            if(logall==TRUE) logwarn("Some Parameters on the boundary, please check the results 
                    and make sure they are the same as using other optimization algorithms 
                    like TRR!")
            runByrd <- FALSE
            if(runByrd){
              print("not freezing!")
              res0 <- optim(res0$par,f2,gr=NULL,plotfit=plotfit,method="L-BFGS-B",lower=lower,upper=upper,hessian=TRUE)
              names(res0)[2] <- "ssr"
            }
          }
        }
      }
      if(optimMethod=="Byrd"){
        
        res0 <- optim(oldparms,f2,gr=NULL,plotfit=plotfit,method="L-BFGS-B",lower=lower,upper=upper,hessian=TRUE)
        names(res0)[2] <- "ssr"
      }
      if(optimMethod=="GENOUD"){
        stop("GENOUD is too slow, not fully tested and implemented yet!")
        
        replaceINF <- function(x){
          x[x==Inf] <- 10^5
        }
        res0 <- genoud(fn=f2,nvars=length(oldparms),
                       Domains=cbind(lower,replaceINF(upper)),hessian=TRUE,plotfit=plotfit)
      }
      if(optimMethod=="bobyqa"){
        ## 
        stop("not good result, too many function evaluations!")
        res1 <- bobyqa(par=res0$par,fn=f2,lower=lower,upper=upper,plotfit=plotfit)
      }
      if(optimMethod=="TRR") {
        res0 <- lsqnonlin(f,xstart=oldparms,l=lower,u=upper,plotfit=plotfit,plotRes=plotRes,...)
        ## Renaming and calculating more statistics to produce summary text and figures.
        res0$par <- as.vector(res0$xcurr)
        names(res0$par) <- pnames
        res0$residuals <- res0$fvec
        res0$hessian <- 2*t(res0$JACOB)%*%(res0$JACOB)
        ##res0$diag <- ?
        
      }
      if(optimMethod== "Nash") {
        ## stop("Not fully implemented yet")
        ## res0 <- nlxb(y ~ kin_mod(P,pnames=names(oldparms)),start=list(P=oldparms),lower=lower,upper=upper)
        res0 <- nlfb(start=oldparms, resfn=f,lower=lower,upper=upper)
      }
      if(optimMethod == "nls2"){
        stop("Not fully implemented yet")
        P <- oldparms
        res0 <- nls2(y ~ kin_mod(P,pnames=names(oldparms)),start=t(oldparms),lower=lower,upper=upper)
      }
      
      ## if(evalMethod='MCMC'), stop here.
      if(evalMethod=='IRLS' & runIRLS ){ ## when IRLS, need an iterative step:
        if(missing(irls.control) || length(irls.control)==0) irls.control <- list(tol=1e-3,maxIter=3)
        if(is.null(irls.control$tol)) irls.control$tol <- 1e-05
        if(is.null(irls.control$maxIter)) irls.control$maxIter <- 5
        
        y <- res0$fvec
        id <- which(!is.na(observed$value))
        observedComplete <- observed[id,]
        errstd <- tapply(y,observed$name[id],sd)
        errstd.old <- rep(1000,length(errstd))
        diffsigma <- KineticEval:::norm(errstd-errstd.old)
        niter <- 0
        while(diffsigma > irls.control$tol && niter <= irls.control$maxIter){
          cat("IRLS iteration at",niter, "; Diff in error variance ", diffsigma,"\n")
          errstd.old <- errstd
          ERR <- errstd[as.character(observed$name)]
          observed$err <- ERR
          f1 <- function(P,...){
            ##browser()
            res <- observed$value-kin_mod(P,plot=FALSE,...)
            res <- res/ERR
            res <- as.vector(res)
            id <- which(is.na(res))
            if(length(id)>0) return(res[-id]) else return(res)
          }
          
          f12 <- function(P,...){
            return(sum(f1(P,...)^2))
          }
          if(optimMethod=="LM") res0 <- nls.lm(par=res0$par,lower=lower,upper=upper,fn=f1)
          if(optimMethod=="TRR") {
            
            res0 <- lsqnonlin(f1,xstart=res0$par,l=lower,u=upper,...)
          }
          if(optimMethod=="Byrd"){
            res0 <- optim(res0$par,f12,method="L-BFGS-B",lower=lower,upper=upper,hessian=TRUE)
            res0$fvec <- f1(res0$par)
          }
          ERR1 <- as.vector(errstd[as.character(observedComplete$name)])
          y <- res0$fvec*ERR1
          errstd <- tapply(y,observed$name[id],sd)
          niter <- niter+1
          diffsigma <- KineticEval:::norm(errstd-errstd.old)
        }
        ## Renaming some list items to be compatabile.
        if(optimMethod=="LM"){
          names(res0)[7] <- "iterations"  # called "niter" here
          names(res0)[9] <- "ssr"         # called "deviance" here
          ##names(res0)[3] <- "residuals"
          res0$residuals <- res0$fvec
          res0$hessian   <- 2*res0$hessian # returns 0.5*hessian!
          Diag          <- unlist(res0[6])
          res0$diag      <- NULL
          res0$diag      <- Diag ### ? What is diag here?
        }
        if(optimMethod=="TRR"){
          ## Renaming and calculating more statistics to produce summary text and figures.
          res0$par <- as.vector(res0$xcurr)
          names(res0$par) <- pnames
          res0$residuals <- res0$fvec
          res0$hessian <- 2*t(res0$JACOB)%*%(res0$JACOB)
          ##res0$diag <- ?
        }
        if(optimMethod=="Byrd"){
          names(res0)[2]<- "ssr"
        }
        
      }
      #### Now deal with hessian!!!!!!
      np <- length(res0$par)
      ##browser()
      if(exactHessian) {
        ##browser()
        if(logall) loginfo("Calcluating Numerical Hessian using numDeriv.")
        if(runIRLS) res0$numH <- try(hessian(f12,res0$par)) else{
          res0$numH <- try(hessian(f2,res0$par))
        }
        res0$numCov <-  try(ginv(0.5*res0$numH), silent = TRUE)
        res0$covar <- res0$numCov
      }else{
        if(exactJacob){
          if(evalMethod=="NLLS") res0$covar <- calcCovar(res0$par,f=f,fval=res0$residuals,lb=lower,ub=upper,...)
          if(evalMethod=="IRLS") res0$covar <- calcCovar(res0$par,f=f1,fval=res0$residuals,lb=lower,ub=upper,...)
        }else{
          res0$covar <-  try(solve(0.5*res0$hessian), silent = TRUE)
          
          if(!is.numeric(res0$covar)){
            ## try once again!
            #browser()
            if(logall==TRUE) {
              loginfo("Convergence criteria met but final hessian is not positive definite. Trying generalized inverse.")
            }
            ## Make diagonal value 0 being NA instead of 
            junk <- qr(res0$hessian)
            
            res0$covar <-  try(ginv(0.5*res0$hessian), silent = TRUE)
            if(class(res0$covar)=="try-error"){
              if(logall==TRUE) loginfo("Hessian not generalized invertable, please report the case!")
              if(evalMethod=="NLLS") res0$covar <- calcCovar(res0$par,f=f,fval=res0$residuals,lb=lower,ub=upper,...)
              if(evalMethod=="IRLS") res0$covar <- calcCovar(res0$par,f=f1,fval=res0$residuals,lb=lower,ub=upper,...)
            }else{
              if(junk$rank < length(res0$par)){
                if(logall==TRUE) loginfo("Hessian is not full rank.")
                junkR <- diag(res0$covar)
                problem <- which(junkR < .Machine$double.eps)
                res0$covar[problem,] <- NA
                res0$covar[,problem] <- NA
              }
            }
          }
          
          
        }
      }
      if(!is.numeric(res0$covar)){
        print("Not able to estimate the covariance matrix due to non-positive definite hessian")
        res0$covar <- matrix(data = NA, nrow = np, ncol = np)
        ############################
      }
      npar <- length(res0$par)
      nobs <- length(res0$fvec)
      res0$df.residual <- (nobs-npar)
      
      ## Deal with fit$_var_ms, var_ms_unscaled, ms_unweighted
      
      ############################
      if(evalMethod=='MCMC'){ ## when MCMC, always use the different
        ## 'sigma' setting and start with the fitted parameters from IRLS/NLS results
        ## have some redundancy here. #########NEED TO MODIFY HERE!!!!!!
        ## calculate sigma
        y <- res0$fvec
        if(!exists("id")) id <- which(!is.na(observed$value))
        errstd <- tapply(y,observed$name[id],sd)
        ERR <- errstd[as.character(observed$name)]
        observed$err <- ERR
        environment(kin_mod_cost) <- environment()
        tmp <- kin_mod_cost(res0$par)
        ## Run MCMC
        niter <- MCMCoptions$niter
        attach(MCMCoptions)
        npar <- length(res0$par)
        nobs <- length(y)
        cov0 <- res0$covar*res0$ssr/(nobs-npar)*2.4^2/npar
        var0 <- tmp$var$SSR.unweighted/tmp$var$N 
        names(var0) <- tmp$var$name
		if(any(is.na(cov0))){
			res0 <- modMCMC(kin_mod_cost,res0$par,...,jump=NULL,lower=lower,upper=upper,prior=prior,
                        var0=var0,wvar0=wvar0,niter=niter,outputlength = outputlength,
                        burninlength = burninlength, updatecov = updatecov,
                        ntrydr=ntrydr,drscale=drscale,verbose=verbose)
		}else{
			res0 <- modMCMC(kin_mod_cost,res0$par,...,jump=cov0,lower=lower,upper=upper,prior=prior,
                        var0=var0,wvar0=wvar0,niter=niter,outputlength = outputlength,
                        burninlength = burninlength, updatecov = updatecov,
                        ntrydr=ntrydr,drscale=drscale,verbose=verbose)
		}
        detach(MCMCoptions)
        err1 <- sqrt(apply(res0$sig,2,mean))
        observed$err <- err1[as.character(paste('var_',observed$name,sep=''))]
        ## #
        
        res0$par <- apply(res0$pars,2,mean)
        ## names(res0$par) <- pnames
        ## ####################################
        ccend <- kin_mod_cost(res0$par)### using the estimated parameter to calculate the values.
        res0$residuals <-ccend$residual$res
        ##predicted_long <- mkin_wide_to_long(out_predicted, time = "time")
        ##res0$predicted <- out_predicted
        
      }
      
      # -----------------------------------------
      ## Result Summary Section for the 2 newly implemented method##
      # -----------------------------------------
      ## We need to return some more information for summary and plotting in the GUI wrap.
      if(!is.null(res0$message)) {
        res0$stopmess <- res0$message
        res0$message <- NULL
      }
      fit <- res0
      fit$solution <- solution
      if (solution == "eigen") {
        fit$coefmat <- mkinmodini$coefmat
      }
      if (solution == "deSolve") {
        fit$mkindiff <- mkindiff
      }
      fit$map <- mkinmodini$map
      fit$diffs <- mkinmodini$diffs
      fit$observed <- mkinmodini$residue
      predicted_long <- mkin_wide_to_long(out_predicted, time = "time") ## TODO??????????? check if we have our_predicted here
      fit$predicted <- out_predicted
      
      
      ## Collect initial parameter values in two dataframes
      ##
      if(outpartri=='default'){
        if(length(state.ini.optim)>0){
          fit$start0 <- data.frame(initial=state.ini.orig.optim,type=rep("state", length(state.ini.optim)),lower=lower[1:length(state.ini.optim)],upper=upper[1:length(state.ini.optim)])
          }else{
            fit$start0 <- data.frame(initial=state.ini.optim,type=rep("state", length(state.ini.optim)),lower=numeric(0),upper=numeric(0))
          }
        ##fit$start0 <- data.frame(initial=state.ini.optim,type=rep("state", length(state.ini.optim)),lower=lower[1:length(state.ini.optim)],upper=upper[1:length(state.ini.optim)])
        start0 <- mkinmodini$start[mkinmodini$start$fixed==0,]
        fit$start0 <- rbind(fit$start0,data.frame(initial=start0$initial,type=start0$type,lower=start0$lower,upper=start0$upper,row.names =rownames(start0)))
        fit$fixed0 <- data.frame(value = state.ini.fixed,type=rep("state", length(state.ini.fixed)),by=rep("user", length(state.ini.fixed)))
        fixed0 <- mkinmodini$start[mkinmodini$start$fixed==1,]
        if(nrow(fixed0)>0) fit$fixed0 <- rbind(fit$fixed0,data.frame(value=fixed0$initial,type=fixed0$type,by=rep('user',nrow(fixed0)),row.names=rownames(fixed0)))
      }else{
        ## for inpartri being defaul with formation fraction and 
        ## outpartri being water-sediment, the program cannot fix anything unless k*ff are both fixed.
        ## there is also no original start0 and fixed0
        fit$start0 <- NULL
        fit$fixed0 <- NULL
      }
      fit$start <- data.frame(initial = c(state.ini.orig.optim, parms.optim))
      fit$start$type = c(rep("state", length(state.ini.optim)), rep("deparm", length(parms.optim)))
      fit$start$lower <- lower
      fit$start$upper <- upper
      
      fit$fixed <- data.frame(
        value = c(state.ini.fixed, parms.fixed))
      fit$fixed$type = c(rep("state", length(state.ini.fixed)), rep("deparm", length(parms.fixed)))
      ## XXXXXX determine whether it is fixed by user or by KinGUII XXXXXXXXXXXXX
      fit$fixed$by <- c(rep("user", length(state.ini.fixed)), mkinmodini$fixed_flag)
      
      
      ## Calculate chi2 error levels according to FOCUS (2006)
      ## 0 values at sample time 0 should not be used.
      observed1 <- observed
      observed1 <- observed[!(observed$time==0 & observed$value==0),]
      means <- aggregate(value ~ time + name, data = observed1, mean, na.rm=TRUE)##using the mean of repeated measurements.
      ##browser()
      ##errdata <- merge(means, predicted_long, observed,by = c("time", "name"), suffixes = c("_mean", "_pred",'_obs'))
      errdata <- merge(means, predicted_long, by = c("time", "name"), suffixes = c("_mean", "_pred"))
      ## !!!here is the problem!!! observed has two values, thus not be able to really merge!!!!
      ## errdata <- merge(errdata, observed,by = c("time", "name"))
      ## names(errdata)[5] <- 'value_obs'
      errobserved <- merge(observed, predicted_long, by = c("time", "name"), suffixes = c("_obs", "_pred"))
      errdata <- errdata[order(errdata$time, errdata$name), ]
      
      errmin.overall <- chi2err(errdata, length(parms.optim) + length(state.ini.optim),errobserved)
      errmin <- data.frame(err.min = errmin.overall$err.min,
                           n.optim = errmin.overall$n.optim, df = errmin.overall$df,
                           err.sig = errmin.overall$err.sig,RMSE=errmin.overall$RMSE,
                           EF=errmin.overall$EF,R2=errmin.overall$R2)
      rownames(errmin) <- "All data"
      for (obs_var in obs_vars)
      {
        errdata.var <- subset(errdata, name == obs_var)
        errobserved.var <- subset(errobserved, name == obs_var)
        if(outpartri=='default'){
          ##n.k.optim <- (paste("k", obs_var, sep="_")) %in% (names(parms.optim))+length(grep(paste("f", obs_var,'to', sep="_"), names(parms.optim)))
          n.k.optim <- (paste("k", obs_var, sep="_")) %in% (names(parms.optim))+length(grep(paste("f",'.*','to',obs_var,sep="_"), names(parms.optim)))
        }
        if(outpartri=='water-sediment'){
          n.k.optim <- length(grep(paste("k_", obs_var, '_',sep=""), names(parms.optim)))
        }
        n.initials.optim <- as.numeric((paste('M0_',obs_var, sep="")) %in% (names(state.ini.optim)))#n.initials.optim <- length(grep(paste('M0_',obs_var, sep=""), names(state.ini.optim)))
        n.optim <- n.k.optim + n.initials.optim
        ## ## added
        k1name <- paste("k1", obs_var,  sep="_")
        k2name <- paste("k2", obs_var,  sep="_")
        gname <- paste("g", obs_var,  sep="_")
        tbname <- paste("tb", obs_var,  sep="_")
        alphaname <- paste("alpha", obs_var,  sep="_")
        betaname <- paste("beta", obs_var,  sep="_")
        ## #
        ## if ("alpha" %in% names(parms.optim)) n.optim <- n.optim + 1
        ## if ("beta" %in% names(parms.optim)) n.optim <- n.optim + 1
        ## if ("k1" %in% names(parms.optim)) n.optim <- n.optim + 1
        ## if ("k2" %in% names(parms.optim)) n.optim <- n.optim + 1
        ## if ("g" %in% names(parms.optim)) n.optim <- n.optim + 1
        ## if ("tb" %in% names(parms.optim)) n.optim <- n.optim + 1
        ## #
        if (alphaname %in% names(parms.optim)) n.optim <- n.optim + 1
        if (betaname %in% names(parms.optim)) n.optim <- n.optim + 1
        if (k1name %in% names(parms.optim)) n.optim <- n.optim + 1
        if (k2name %in% names(parms.optim)) n.optim <- n.optim + 1
        if (gname %in% names(parms.optim)) n.optim <- n.optim + 1
        if (tbname %in% names(parms.optim)) n.optim <- n.optim + 1
        ## Change the ouput order for k1, k2, g
        orderDFOP <- FALSE
        if(orderDFOP){
          if (k1name %in% names(parms.optim) | k2name %in% names(parms.optim)){
            k1 <- parms.optim[k1name]
            k2 <- parms.optim[k2name]
            if(k1 > k2) {
              parms.optim[k1name] <- k2
              parms.optim[k2name] <- k1
              parms.optim[gname] <- 1 - parms.optim[gname]
            }
          }
        }
        ##                             #errmin.tmp <- mkinerrmin(errdata.var, n.optim)
        errmin.tmp <- chi2err(errdata.var, n.optim,errobserved.var)
        errmin[obs_var, c("err.min", "n.optim", "df",'err.sig','RMSE','EF','R2')] <- errmin.tmp
      }
      fit$errmin <- errmin
      
      ## Calculate degradation(not dissipation) times DT50 and DT90 and formation fractions
    
      parms.all = c(fit$par, parms.fixed)
      fit$distimes <- data.frame(DT50 = rep(NA, length(obs_vars)), DT90 = rep(NA, length(obs_vars)),Kinetic=rep(NA,length(obs_vars)),row.names = obs_vars)
      #########################################################################
      if(mkinmodini$outpartri=='default'){
        ## Now deals with formation fractions in case the outpartri is 'default'##
        ## Keep the original IRLSkinfit.gui codes ##
        fit$ff <- vector()
        ff_names = names(mkinmodini$ff)
        for (ff_name in ff_names)
        {
          fit$ff[[ff_name]] =
            eval(parse(text = mkinmodini$ff[ff_name]), as.list(parms.all))
        }
        ## ###
        
        for (obs_var in obs_vars) {
          
          if(mkinmodini$sinkT[obs_var]){
            f_tot <- grep(paste(obs_var, "_",sep=''), names(fit$ff), value=TRUE)
            f_exp <- grep(paste(obs_var, "to",obs_var,sep='_'), names(fit$ff), value=TRUE)
            f_exp1 <- grep(paste(obs_var, "to",'sink',sep='_'), names(fit$ff), value=TRUE)
            fit$ff[[paste(obs_var,'to', "sink", sep="_")]] = 1 - sum(fit$ff[f_tot])+sum(fit$ff[f_exp])+sum(fit$ff[f_exp1])
            
          }else{
            ## when SINK is turned off by the user, then there is no need for calculation
            fit$ff[[paste(obs_var,'to', "sink", sep="_")]] = 0
          }
          type = names(mkinmodini$map[[obs_var]])[1]
          k1name <- paste("k1", obs_var,  sep="_")
          k2name <- paste("k2", obs_var,  sep="_")
          gname <- paste("g", obs_var,  sep="_")
          tbname <- paste("tb", obs_var,  sep="_")
          alphaname <- paste("alpha", obs_var,  sep="_")
          betaname <- paste("beta", obs_var,  sep="_")
          if (type == "SFO") {
              k_name <- paste("k", obs_var,sep="_")#grep(paste("k", obs_var,sep="_"), names(parms.all), value=TRUE)
              k_tot <- parms.all[k_name]
            
            DT50 = log(2)/k_tot
            DT90 = log(10)/k_tot
            ## for (k_name in k_names)
            ## {
            ##   fit$ff[[sub("k_", "", k_name)]] = parms.all[[k_name]] / k_tot
            ## }
          }
          if (type == "FOMC") {
            ## alpha = parms.all["alpha"]
            ## beta = parms.all["beta"]
            alpha = parms.all[alphaname]
            beta = parms.all[betaname]
            DT50 = beta * (2^(1/alpha) - 1)
            DT90 = beta * (10^(1/alpha) - 1)
            ## ff_names = names(mkinmodini$ff)
            ## for (ff_name in ff_names)
            ## {
            ##   fit$ff[[paste(obs_var, ff_name, sep="_")]] =
            ##     eval(parse(text = mkinmodini$ff[ff_name]), as.list(parms.all))
            ## }
            ## fit$ff[[paste(obs_var, "sink", sep="_")]] = 1 - sum(fit$ff)
          }
          if (type == "DFOP") {
            ## k1 = parms.all["k1"]
            ## k2 = parms.all["k2"]
            ## g = parms.all["g"]
           
            k1 = parms.all[k1name]
            k2 = parms.all[k2name]
            g = parms.all[gname]
            f <- function(t, x) {
              ((g * exp( - k1 * t) + (1 - g) * exp( - k2 * t)) - (1 - x/100))^2
            }
            fDT <- function(t, x) {
              ((g * exp( - k1 * t) + (1 - g) * exp( - k2 * t)) - (1 - x/100))
            }
            
            DTmax <- 1000
          
            DTmax1 <- log(2)/min(k1,k2)
            DTmin <- log(2)/max(k1,k2)
            if(DTmax1==Inf) {
              ## This is the case when one of the k is 0.
              ## The whole proportion won't degrade at all.
              nondeg <- (k1==0)*g+(k2==0)*(1-g)
              kdeg <- max(k1,k2)
              gdeg <- 1-nondeg
              if(nondeg > 0.5){
                DT50 <- Inf
                DT90 <- Inf
              }else if(nondeg>0.1){
                DT90 <- Inf
                DT50 <- log(100*gdeg/(100*gdeg-50))/kdeg
              }else{
                DT50 <- log(100*gdeg/(100*gdeg-50))/kdeg
                DT90 <- log(100*gdeg/(100*gdeg-90))/kdeg
              }
            }else{
              resUniroot <- uniroot(fDT,c(DTmin,DTmax1),x=50)
              DT50.o <- resUniroot$root
              if(resUniroot$f.root > 1e-4) print("Uniroot calculation for DFOP DT50 might be wrong!Do an optimization, check if the DT50 in the output make sense!")
              resOptim <- optimize(f, c(DTmin, DTmax1), x=50)
              DT50.o1 <- resOptim$minimum
              if(resOptim$objective > 1e-4) print("One dimentional optimization calculation for DFOP DT50 might be wrong!Doing a uniroot calculation, check if the DT50 in the output make sense!")
              DT50.o <- ifelse(resUniroot$f.root>resOptim$objective, DT50.o1,DT50.o)
              ##DT50 = ifelse(DTmax - DT50.o < 0.1, NA, DT50.o)
              DT50 <- DT50.o   
              DTmax1 <- log(10)/min(k1,k2)
              DTmin <- log(10)/max(k1,k2)
              DT90.o <- uniroot(fDT,c(DTmin,DTmax1),x=90)$root
              DT90.o1 <- optimize(f, c(DTmin, DTmax1), x=90)$minimum
              DT90.o <- ifelse(f(DT90.o,90)>f(DT90.o1,90), DT90.o1,DT90.o)
              DT90 <- DT90.o
            }
            
            }
          if (type == "HS") {
            k1 = parms.all[k1name]
            k2 = parms.all[k2name]
            tb = parms.all[tbname]
            HSendpoints <- "analytical"
            if(HSendpoints=="iterative"){
              f <- function(t, x) {
                fraction = ifelse(t <= tb, exp(-k1 * t), exp(-k1 * tb) * exp(-k2 * (t - tb)))
                (fraction - (1 - x/100))^2
              }
              DTmax <- 1000
              hso1 <- nlminb(0.0001,f, x=50)
              hso2 <- nlminb(tb,f, x=50)
              DT50.o <- ifelse(hso1$objective<=hso2$objective,hso1$par,hso2$par)
              DT50 = ifelse(DTmax - DT50.o < 0.1, NA, DT50.o)
              
              hso1 <- nlminb(0.0001,f, x=90)
              hso2 <- nlminb(tb,f, x=90)
              DT90.o <- ifelse(hso1$objective<=hso2$objective,hso1$par,hso2$par)
              DT90 = ifelse(DTmax - DT90.o < 0.1, NA, DT90.o)
            }else{## HSendpoints=="analytical"
              calcDTx <- function(k1,k2,tb,x){
                DTx1 <- log(100/(100-x))/k1
                DTx2 <- tb+(log(100/(100-x))-k1*tb)/k2
                
                if(DTx1 > tb){
                  return(DTx2)
                }else{
                  return(DTx1)
                  #                 if(exp(-k1*DTx1)>=(1-x/100)){
                  #                   ## fraction greater than x% then return DTx1
                  #                   return(DTx1)
                  #                 }else{
                  #                   return(DTx2)
                  #                 }
                }
              }
              DT50 <- calcDTx(k1,k2,tb,x=50)
              DT90 <- calcDTx(k1,k2,tb,x=90)
            }
            
            
          }
          ## if (type %in% c("DFOP", "HS")) {
          ##     DTmax <- 1000
          ##     DT50.o <- optimize(f, c(0.001, DTmax), x=50)$minimum
          ##     DT50 = ifelse(DTmax - DT50.o < 0.1, NA, DT50.o)
          ##     DT90.o <- optimize(f, c(0.001, DTmax), x=90)$minimum
          ##     DT90 = ifelse(DTmax - DT90.o < 0.1, NA, DT90.o)
          
          ## }
          if (type == "SFORB") {
            # FOCUS kinetics (2006), p. 60 f
            k_out_names = grep(paste("k", obs_var, "free", sep="_"), names(parms.all), value=TRUE)
            k_out_names = setdiff(k_out_names, paste("k", obs_var, "free", "bound", sep="_"))
            k_1output = sum(parms.all[k_out_names])
            k_12 = parms.all[paste("k", obs_var, "free", "bound", sep="_")]
            k_21 = parms.all[paste("k", obs_var, "bound", "free", sep="_")]
            
            sqrt_exp = sqrt(1/4 * (k_12 + k_21 + k_1output)^2 + k_12 * k_21 - (k_12 + k_1output) * k_21)
            b1 = 0.5 * (k_12 + k_21 + k_1output) + sqrt_exp
            b2 = 0.5 * (k_12 + k_21 + k_1output) - sqrt_exp
            
            SFORB_fraction = function(t) {
              ((k_12 + k_21 - b1)/(b2 - b1)) * exp(-b1 * t) +
                ((k_12 + k_21 - b2)/(b1 - b2)) * exp(-b2 * t)
            }
            f_50 <- function(t) (SFORB_fraction(t) - 0.5)^2
            max_DT <- 1000
            DT50.o <- optimize(f_50, c(0.01, max_DT))$minimum
            if (abs(DT50.o - max_DT) < 0.01) DT50 = NA else DT50 = DT50.o
            f_90 <- function(t) (SFORB_fraction(t) - 0.1)^2
            DT90.o <- optimize(f_90, c(0.01, 1000))$minimum
            if (abs(DT90.o - max_DT) < 0.01) DT90 = NA else DT90 = DT90.o
            for (k_out_name in k_out_names)
            {
              fit$ff[[sub("k_", "", k_out_name)]] = parms.all[[k_out_name]] / k_1output
            }
          }
          #fit$distimes[obs_var, ] = c(ifelse(is.na(DT50),NA,as.numeric(formatC(DT50,4,format='f'))), 
          #                           ifelse(is.na(DT90),NA,as.numeric(formatC(DT90,4,format='f'))),type)#c(DT50, DT90,type)
          fit$distimes[obs_var,1:2] <- c(DT50,DT90)
          fit$distimes[obs_var,3] <- type
        }
      }
      if(mkinmodini$outpartri=='water-sediment'){
        ## Now deals with multiple ks and additionally calculate formation fractions in case the outpartri is 'water-sediment'##
        fit$ff <- vector() ## for FOMC, HS, DFOP, they might still use formation fractions.
        ##ff_names = names(mkinmodini$ff)
        ##for (ff_name in ff_names)
        ##{
        ##   fit$ff[[ff_name]] =
        ##      eval(parse(text = mkinmodini$ff[ff_name]), as.list(parms.all))
        ##}
        ## ###
        ##browser()
        for (obs_var in obs_vars) {
          ##f_tot <- grep(paste(obs_var, "_",sep=''), names(fit$ff), value=TRUE)
          ##f_exp <- grep(paste(obs_var, "to",obs_var,sep='_'), names(fit$ff), value=TRUE)
          ##f_exp1 <- grep(paste(obs_var, "to",'sink',sep='_'), names(fit$ff), value=TRUE)
          ##fit$ff[[paste(obs_var,'to', "sink", sep="_")]] = 1 - sum(fit$ff[f_tot])+sum(fit$ff[f_exp])+sum(fit$ff[f_exp1])
          type = names(mkinmodini$map[[obs_var]])[1]
          k1name <- paste("k1", obs_var,  sep="_")
          k2name <- paste("k2", obs_var,  sep="_")
          gname <- paste("g", obs_var,  sep="_")
          tbname <- paste("tb", obs_var,  sep="_")
          alphaname <- paste("alpha", obs_var,  sep="_")
          betaname <- paste("beta", obs_var,  sep="_")
          if (type == "SFO") {
            k_names = grep(paste("k_", obs_var,'_', sep=""), names(parms.all), value=TRUE)
            cat("Seperating transportation rate from the degradation rates")
            ## check whether there are 2 arrows,
            ## 
            if(logall){loginfo("For water-sediment settings, degradation DT50 are reported.")}
            k_ins <- grep(paste0("k.*_to_",obs_var), names(parms.all), value=TRUE)
            rb <- sub(paste0("_to_",obs_var),"",k_ins)
            rb <- sub("k_","",rb)
            k_transport <- sapply(rb,function(x)grep(paste0(obs_var,"_to_",x),k_names,value=TRUE))
            k_deg <- setdiff(k_names,k_transport)
            k_tot = sum(parms.all[k_deg])
            DT50 = log(2)/k_tot
            DT90 = log(10)/k_tot
            for (k_name in k_names){
              fit$ff[[sub("k_", "", k_name)]] = parms.all[[k_name]] / sum(parms.all[k_name])
            }
          }
          if (type == "FOMC") {
            alpha = parms.all[alphaname]
            beta = parms.all[betaname]
            DT50 = beta * (2^(1/alpha) - 1)
            DT90 = beta * (10^(1/alpha) - 1)
            ff_names = names(mkinmodini$ff)
            for (ff_name in ff_names){
              fit$ff[[paste(obs_var, ff_name, sep="_")]] =
                eval(parse(text = mkinmodini$ff[ff_name]), as.list(parms.all))
            }
            fit$ff[[paste(obs_var, "sink", sep="_")]] = 1 - sum(fit$ff)
          }
          if (type == "DFOP") {
		 
            k1 = parms.all[k1name]
            k2 = parms.all[k2name]
            g = parms.all[gname]
            f <- function(t, x) {
              ((g * exp( - k1 * t) + (1 - g) * exp( - k2 * t)) - (1 - x/100))^2
            }
            
            fDT <- function(t, x) {
              ((g * exp( - k1 * t) + (1 - g) * exp( - k2 * t)) - (1 - x/100))
            }
            
            DTmax <- 1000
            
            DTmax1 <- log(2)/min(k1,k2)
            DTmin <- log(2)/max(k1,k2)
		  if(DTmax1==Inf) {
		    ## This is the case when one of the k is 0.
		    ## The whole proportion won't degrade at all.
		    nondeg <- (k1==0)*g+(k2==0)*(1-g)
		    kdeg <- max(k1,k2)
		    gdeg <- 1-nondeg
		    if(nondeg > 0.5){
		      DT50 <- Inf
		      DT90 <- Inf
		    }else if(nondeg>0.1){
		      DT90 <- Inf
		      DT50 <- log(100*gdeg/(100*gdeg-50))/kdeg
		    }else{
		      DT50 <- log(100*gdeg/(100*gdeg-50))/kdeg
		      DT90 <- log(100*gdeg/(100*gdeg-90))/kdeg
		    }
		  }else{
		    resUniroot <- uniroot(fDT,c(DTmin,DTmax1),x=50)
		    DT50.o <- resUniroot$root
		    if(resUniroot$f.root > 1e-4) print("Uniroot calculation for DFOP DT50 might be wrong!Do an optimization, check if the DT50 in the output make sense!")
		    resOptim <- optimize(f, c(DTmin, DTmax1), x=50)
		    DT50.o1 <- resOptim$minimum
		    if(resOptim$objective > 1e-4) print("One dimentional optimization calculation for DFOP DT50 might be wrong!Doing a uniroot calculation, check if the DT50 in the output make sense!")
		    DT50.o <- ifelse(resUniroot$f.root>resOptim$objective, DT50.o1,DT50.o)
		    ##DT50 = ifelse(DTmax - DT50.o < 0.1, NA, DT50.o)
		    DT50 <- DT50.o            
		    DT90.o <- uniroot(fDT,c(DTmin,DTmax1),x=90)$root
		    DT90.o1 <- optimize(f, c(DTmin, DTmax1), x=90)$minimum
		    DT90.o <- ifelse(f(DT90.o,90)>f(DT90.o1,90), DT90.o1,DT90.o)
		    DT90 <- DT90.o
		  }
            
          }
          if (type == "HS") {
            k1 = parms.all[k1name]
            k2 = parms.all[k2name]
            tb = parms.all[tbname]
            HSendpoints <- "analytical"
            if(HSendpoints=="iterative"){
              f <- function(t, x) {
                fraction = ifelse(t <= tb, exp(-k1 * t), exp(-k1 * tb) * exp(-k2 * (t - tb)))
                (fraction - (1 - x/100))^2
              }
              DTmax <- 1000
              hso1 <- nlminb(0.0001,f, x=50)
              hso2 <- nlminb(tb,f, x=50)
              DT50.o <- ifelse(hso1$objective<=hso2$objective,hso1$par,hso2$par)
              DT50 = ifelse(DTmax - DT50.o < 0.1, NA, DT50.o)
              
              hso1 <- nlminb(0.0001,f, x=90)
              hso2 <- nlminb(tb,f, x=90)
              DT90.o <- ifelse(hso1$objective<=hso2$objective,hso1$par,hso2$par)
              DT90 = ifelse(DTmax - DT90.o < 0.1, NA, DT90.o)
            }else{## HSendpoints=="analytical"
              calcDTx <- function(k1,k2,tb,x){
                DTx1 <- log(100/(100-x))/k1
                DTx2 <- tb+(log(100/(100-x))-k1*tb)/k2
                
                if(DTx1 > tb){
                  return(DTx2)
                }else{
                  return(DTx1)
                  #                 if(exp(-k1*DTx1)>=(1-x/100)){
                  #                   ## fraction greater than x% then return DTx1
                  #                   return(DTx1)
                  #                 }else{
                  #                   return(DTx2)
                  #                 }
                }
              }
              DT50 <- calcDTx(k1,k2,tb,x=50)
              DT90 <- calcDTx(k1,k2,tb,x=90)
            }
          }
          
          if (type == "SFORB") {
            # FOCUS kinetics (2006), p. 60 f
            k_out_names = grep(paste("k", obs_var, "free", sep="_"), names(parms.all), value=TRUE)
            k_out_names = setdiff(k_out_names, paste("k", obs_var, "free", "bound", sep="_"))
            k_1output = sum(parms.all[k_out_names])
            k_12 = parms.all[paste("k", obs_var, "free", "bound", sep="_")]
            k_21 = parms.all[paste("k", obs_var, "bound", "free", sep="_")]
            
            sqrt_exp = sqrt(1/4 * (k_12 + k_21 + k_1output)^2 + k_12 * k_21 - (k_12 + k_1output) * k_21)
            b1 = 0.5 * (k_12 + k_21 + k_1output) + sqrt_exp
            b2 = 0.5 * (k_12 + k_21 + k_1output) - sqrt_exp
            
            SFORB_fraction = function(t) {
              ((k_12 + k_21 - b1)/(b2 - b1)) * exp(-b1 * t) +
                ((k_12 + k_21 - b2)/(b1 - b2)) * exp(-b2 * t)
            }
            f_50 <- function(t) (SFORB_fraction(t) - 0.5)^2
            max_DT <- 1000
            DT50.o <- optimize(f_50, c(0.01, max_DT))$minimum
            if (abs(DT50.o - max_DT) < 0.01) DT50 = NA else DT50 = DT50.o
            f_90 <- function(t) (SFORB_fraction(t) - 0.1)^2
            DT90.o <- optimize(f_90, c(0.01, 1000))$minimum
            if (abs(DT90.o - max_DT) < 0.01) DT90 = NA else DT90 = DT90.o
            for (k_out_name in k_out_names)
            {
              fit$ff[[sub("k_", "", k_out_name)]] = parms.all[[k_out_name]] / k_1output
            }
          }
          #fit$distimes[obs_var, ] = c(ifelse(is.na(DT50),NA,formatC(DT50,4,format='f')), 
          #                            ifelse(is.na(DT90),NA,formatC(DT90,4,format='f')),
          #                            type)#c(DT50, DT90,type)
          fit$distimes[obs_var,1:2] <- c(DT50,DT90)
          fit$distimes[obs_var,3] <- type
        }
        
      }
      
      ## Collect observed, predicted and residuals
      observed0 <-  mkin_wide_to_long(mkinmodini$data0,time='time')
      observed0$err <- observed$err
      data <- merge(observed, predicted_long, by = c("time", "name"))
      data0 <- merge(observed0, predicted_long, by = c("time", "name"))
      ##names(data) <- c("time", "variable", "observed", "predicted")
      names(data) <- c("time", "variable", "observed","err-std", "predicted")
      ##browser()
      names(data0) <- c("time", "variable", "observed","err-std", "predicted")
      data$residual <- data$observed - data$predicted
      data0$residual <- data$residual
      data$variable <- ordered(data$variable, levels = obs_vars)
      data0$variable <- data$variable
      tmpid <- is.na(data0$residual) & !is.na(data0$observed)
      data0$'err-std'[tmpid] <- 0
      fit$data <- data[order(data$variable, data$time), ]
      fit$data0 <- data0[order(data0$variable, data0$time), ]
      fit$atol <- atol
      fit$inpartri <- inpartri
      fit$outpartri <- outpartri
      fit$optimmethod <- optimMethod
      fit$evalMethod <- evalMethod
      fit$modelmess <- mkinmodini$modelmess
      fit$trff <- mkinmodini$ff
      if(evalMethod=="MCMC")  class(fit) <- c('mcmckingui',"modMCMC") else{
        class(fit) <- c('kingui')
      }
      
    }## end if(optimMethod %in% ...) else
    
  }
  # ----------------------------------------
  return(fit)
}
##' Summary function for Kinetic models
##' 
##' @title summary method for mkinmod.full
##' @param mkinmodini Kinetic Model or fit results with all parameters fixed
##' @return a list containing DT50 and other information
##' @author Zhenglei Gao
##' @S3method summary mkinmod.full
##' @rdname summary.mkinmod.full
summary.mkinmod.full<-function(mkinmodini,ctr=kingui.control(),version=NULL){
 
  ## KineticEval::summary.mkinmod.full
  odesolver <- ctr$odesolver
  atol <- ctr$atol
  rtol <- ctr$rtol
  ## -----------------------------------
  ## Get the parametrization.
  inpartri <- mkinmodini$inpartri
  outpartri <- mkinmodini$outpartri
  ##
  
  ## mkinmodini is an object by mkinmod.full
  parms.ini <- mkinmodini$parms.ini
  state.ini <- mkinmodini$state.ini
  lower <- mkinmodini$lower
  upper <- mkinmodini$upper
  
  
  
  mod_vars <- names(mkinmodini$diffs)
  
  fixed_parms <-mkinmodini$parms##fixed_parms <- mkinmodini$fixed_parms
  fixed_initials <- mod_vars##mkinmodini$fixed_initials
  
  
  observed <-  mkin_wide_to_long(mkinmodini$residue,time='time')
  observed$err <-c(as.matrix(mkinmodini$weightmat))
  ## Subset dataframe with mapped (modelled) variables
  observed <- subset(observed, name %in% names(mkinmodini$map))
  ## Get names of observed variables
  ## NOTE HERE: the order may not be the same as the input mkinmod.full differential equations list. 
  ## XXXXX TODO XXXX Reorder them maybe a good idea if the data is given from a data file while the mkinmod.full is defined not following the colnames order, 
  ## although it is already taken care of in the cost(P) function to reorder the odeini using mod_vars
  obs_vars = unique(as.character(observed$name))
  
  
  ## Name the parameters if they are not named yet ## usually they are already names
  if(is.null(names(parms.ini))) names(parms.ini) <- mkinmodini$parms
  
  ## Name the inital parameter values if they are not named yet
  if(is.null(names(state.ini))) names(state.ini) <- mod_vars
  
  ## Parameters to be optimised no optimised parameters
  parms.fixed <- parms.ini## parms.fixed <- parms.ini[fixed_parms]
  optim_parms <- character(0)##setdiff(names(parms.ini), fixed_parms)
  parms.optim <- parms.ini[optim_parms]
  
  
  ## # ### ### ### ### ###
  state.ini.fixed <- state.ini[fixed_initials]
  optim_initials <- character(0)##setdiff(names(state.ini), fixed_initials)
  state.ini.optim <- state.ini[optim_initials]
  state.ini.optim.boxnames <- names(state.ini.optim)
  state.ini.fixed.boxnames <- names(state.ini.fixed)
  if(length(state.ini.optim) > 0) {
    names(state.ini.optim) <- paste('M0',names(state.ini.optim),  sep="_")
  }
  if(length(state.ini.fixed) > 0) {
    names(state.ini.fixed) <- paste('M0',names(state.ini.fixed), sep="_")
  }
  
  eigen <- FALSE
  oldparms <- c(state.ini.optim,parms.optim)
  if (length(mkinmodini$map) == 1) {
    solution = "analytical"
  } else {
    if (is.matrix(mkinmodini$coefmat) & eigen) solution = "eigen"
    else solution = "deSolve"
  }
  ## always define mkindiff function since most of the time we will use it.
  mkindiff <- function(t, state, parms) {
    time <- t
    diffs <- vector()
    for (box in mod_vars)
    {
      diffname <- paste("d", box, sep="_")
      diffs[diffname] <- with(as.list(c(time,state, parms)),
                              eval(parse(text=mkinmodini$diffs[[box]])))
    }
    ##https://stat.ethz.ch/pipermail/r-sig-dynamic-models/2010q2/000031.html
    #bady <- (!is.finite(diffs))|(diffs<=0)
    #diffs[bady] <- 0 
    return(list(c(diffs)))
  }
  
  if(is.null(mkinmodini$costout)){
    environment(kin_mod_cost) <- environment()
    costout <- kin_mod_cost(P=c(state.ini.optim,parms.optim),plot=TRUE)
  }else{
    costout <- mkinmodini$costout
    out_predicted <- mkinmodini$out_predicted
  }
  predicted_long <- mkin_wide_to_long(out_predicted, time = "time")
  y <- costout$residuals$mod
  #browser()
  id <- which(!is.na(observed$value))
  observedComplete <- observed[id,]
  errstd <- tapply(y[id],observed$name[id],sd)################## TO MODIFY ####################
  ERR <- errstd[as.character(observed$name)]
  observed$err <- ERR
  ## return fixed parameters
  out <- list()
  out$version <- version
  out$map <- mkinmodini$map
  out$diffs <- mkinmodini$diffs
  out$model <- costout$model
  out$var <- out$var
  
  out$fixed <- data.frame(value = c(state.ini.fixed, parms.fixed))
  out$fixed$type = c(rep("state", length(state.ini.fixed)), rep("deparm", length(parms.fixed)))
  
  
  ## Calculate chi2 error levels according to FOCUS (2006)
  ## 0 values at sample time 0 should not be used.
  observed1 <- observed
  observed1 <- observed[!(observed$time==0 & observed$value==0),]
  means <- aggregate(value ~ time + name, data = observed1, mean, na.rm=TRUE)
  ##using the mean of repeated measurements.!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ##errdata <- merge(means, predicted_long, observed,by = c("time", "name"), suffixes = c("_mean", "_pred",'_obs'))
  errdata <- merge(means, predicted_long, by = c("time", "name"), suffixes = c("_mean", "_pred"))
  ## !!!here is the problem!!! observed has two values, thus not be able to really merge!!!!
  ## errdata <- merge(errdata, observed,by = c("time", "name"))
  ## names(errdata)[5] <- 'value_obs'
  errobserved <- merge(observed, predicted_long, by = c("time", "name"), suffixes = c("_obs", "_pred"))
  errdata <- errdata[order(errdata$time, errdata$name), ]
  errmin.overall <- chi2err(errdata, length(parms.optim) + length(state.ini.optim),errobserved)
  errmin <- data.frame(err.min = errmin.overall$err.min,
                       n.optim = errmin.overall$n.optim, df = errmin.overall$df,
                       err.sig = errmin.overall$err.sig,RMSE=errmin.overall$RMSE,
                       EF=errmin.overall$EF,R2=errmin.overall$R2)
  rownames(errmin) <- "All data"
  
  for (obs_var in obs_vars)
  {
    errdata.var <- subset(errdata, name == obs_var)
    errobserved.var <- subset(errobserved, name == obs_var)
    if(outpartri=='default'){
      ##n.k.optim <- (paste("k", obs_var, sep="_")) %in% (names(parms.optim))+length(grep(paste("f", obs_var,'to', sep="_"), names(parms.optim)))
      n.k.optim <- (paste("k", obs_var, sep="_")) %in% (names(parms.optim))+length(grep(paste("f",'.*','to',obs_var,sep="_"), names(parms.optim)))
    }
    if(outpartri=='water-sediment'){
      n.k.optim <- length(grep(paste("k_", obs_var, '_',sep=""), names(parms.optim)))
    }
    n.initials.optim <- as.numeric((paste('M0_',obs_var, sep="")) %in% (names(state.ini.optim)))#n.initials.optim <- length(grep(paste('M0_',obs_var, sep=""), names(state.ini.optim)))
    n.optim <- n.k.optim + n.initials.optim
    ## ## added
    k1name <- paste("k1", obs_var,  sep="_")
    k2name <- paste("k2", obs_var,  sep="_")
    gname <- paste("g", obs_var,  sep="_")
    tbname <- paste("tb", obs_var,  sep="_")
    alphaname <- paste("alpha", obs_var,  sep="_")
    betaname <- paste("beta", obs_var,  sep="_")

    if (alphaname %in% names(parms.optim)) n.optim <- n.optim + 1
    if (betaname %in% names(parms.optim)) n.optim <- n.optim + 1
    if (k1name %in% names(parms.optim)) n.optim <- n.optim + 1
    if (k2name %in% names(parms.optim)) n.optim <- n.optim + 1
    if (gname %in% names(parms.optim)) n.optim <- n.optim + 1
    if (tbname %in% names(parms.optim)) n.optim <- n.optim + 1
    
    ##                             #errmin.tmp <- mkinerrmin(errdata.var, n.optim)
    errmin.tmp <- chi2err(errdata.var, n.optim,errobserved.var)
    errmin[obs_var, c("err.min", "n.optim", "df",'err.sig','RMSE','EF','R2')] <- errmin.tmp
  }
  out$errmin <- errmin
  ## return true formation fractions
  
  parms.all <- c(parms.ini,state.ini)
  out$ff <- vector()
  ff_names = names(mkinmodini$ff)
  for (ff_name in ff_names)
  {
    out$ff[[ff_name]] =
      eval(parse(text = mkinmodini$ff[ff_name]), as.list(parms.all))
  }
  ## ###
  ## return Chi2 error levels in percent
  
  ## Return DT50 and DT90s
  out$distimes <- data.frame(DT50 = rep(NA, length(obs_vars)), DT90 = rep(NA, length(obs_vars)),Kinetic=rep(NA,length(obs_vars)),row.names = obs_vars)
  parms.all <- parms.fixed
  if(mkinmodini$outpartri=='default'){
    ## Now deals with formation fractions in case the outpartri is 'default'##
    ## Keep the original IRLSkinfit.gui codes ##
   
    ## ###
    ##browser()
    for (obs_var in obs_vars) {
      type = names(mkinmodini$map[[obs_var]])[1]
      k1name <- paste("k1", obs_var,  sep="_")
      k2name <- paste("k2", obs_var,  sep="_")
      gname <- paste("g", obs_var,  sep="_")
      tbname <- paste("tb", obs_var,  sep="_")
      alphaname <- paste("alpha", obs_var,  sep="_")
      betaname <- paste("beta", obs_var,  sep="_")
      if (type == "SFO") {
        ##k_names = grep(paste("k", obs_var, sep="_"), names(parms.all), value=TRUE)
        ##k_tot = sum(parms.all[k_names])
        k_name <- paste("k", obs_var,sep="_")#grep(paste("k", obs_var,sep="_"), names(parms.all), value=TRUE)
        k_tot <- parms.all[k_name]
        DT50 = log(2)/k_tot
        DT90 = log(10)/k_tot
        ## for (k_name in k_names)
        ## {
        ##   fit$ff[[sub("k_", "", k_name)]] = parms.all[[k_name]] / k_tot
        ## }
      }
      if (type == "FOMC") {
        ## alpha = parms.all["alpha"]
        ## beta = parms.all["beta"]
        alpha = parms.all[alphaname]
        beta = parms.all[betaname]
        DT50 = beta * (2^(1/alpha) - 1)
        DT90 = beta * (10^(1/alpha) - 1)
        ## ff_names = names(mkinmodini$ff)
        ## for (ff_name in ff_names)
        ## {
        ##   fit$ff[[paste(obs_var, ff_name, sep="_")]] =
        ##     eval(parse(text = mkinmodini$ff[ff_name]), as.list(parms.all))
        ## }
        ## fit$ff[[paste(obs_var, "sink", sep="_")]] = 1 - sum(fit$ff)
      }
      if (type == "DFOP") {
       
        ## k1 = parms.all["k1"]
        ## k2 = parms.all["k2"]
        ## g = parms.all["g"]
        k1 = parms.all[k1name]
        k2 = parms.all[k2name]
        g = parms.all[gname]
        f <- function(t, x) {
          ((g * exp( - k1 * t) + (1 - g) * exp( - k2 * t)) - (1 - x/100))^2
        }
        ## browser() ## debug calc DT50 or DT90
        DTmax <- 1000
        DT50.o <- optimize(f, c(0.0001,DTmax), x=50)$minimum
        DTmax1 <- log(2)/min(k1,k2)
        if(DTmax1==Inf) DTmax1 <- .Machine$double.xmax
        DT50.o1 <- optimize(f, c(0, DTmax1), x=50)$minimum
        DT50.o <- ifelse(f(DT50.o,50)>f(DT50.o1,50), DT50.o1,DT50.o)
        ##DT50 = ifelse(DTmax - DT50.o < 0.1, NA, DT50.o)
        DT50 <- DT50.o
        DT90.o <- optimize(f, c(0.001, DTmax), x=90)$minimum
        DTmax1 <- log(10)/min(k1,k2)
        if(DTmax1==Inf) DTmax1 <- .Machine$double.xmax
        DT90.o1 <- optimize(f, c(0, DTmax1), x=90)$minimum
        ##DT90.o <- ifelse(f(DT90.o,90)>f(DT90.o1,90), DT90.o1,DT90.o)
        DT90 <- DT90.o
        ##DT90 = ifelse(DTmax - DT90.o < 0.1, NA, DT90.o)
      }
      if (type == "HS") {
        k1 = parms.all[k1name]
        k2 = parms.all[k2name]
        tb = parms.all[tbname]
        HSendpoints <- "analytical"
        if(HSendpoints=="iterative"){
          f <- function(t, x) {
            fraction = ifelse(t <= tb, exp(-k1 * t), exp(-k1 * tb) * exp(-k2 * (t - tb)))
            (fraction - (1 - x/100))^2
          }
          DTmax <- 1000
          hso1 <- nlminb(0.0001,f, x=50)
          hso2 <- nlminb(tb,f, x=50)
          DT50.o <- ifelse(hso1$objective<=hso2$objective,hso1$par,hso2$par)
          DT50 = ifelse(DTmax - DT50.o < 0.1, NA, DT50.o)
          
          hso1 <- nlminb(0.0001,f, x=90)
          hso2 <- nlminb(tb,f, x=90)
          DT90.o <- ifelse(hso1$objective<=hso2$objective,hso1$par,hso2$par)
          DT90 = ifelse(DTmax - DT90.o < 0.1, NA, DT90.o)
        }else{## HSendpoints=="analytical"
          calcDTx <- function(k1,k2,tb,x){
            DTx1 <- log(100/(100-x))/k1
            DTx2 <- tb+(log(100/(100-x))-k1*tb)/k2
            
            if(DTx1 > tb){
              return(DTx2)
            }else{
              return(DTx1)
              #                 if(exp(-k1*DTx1)>=(1-x/100)){
              #                   ## fraction greater than x% then return DTx1
              #                   return(DTx1)
              #                 }else{
              #                   return(DTx2)
              #                 }
            }
          }
          DT50 <- calcDTx(k1,k2,tb,x=50)
          DT90 <- calcDTx(k1,k2,tb,x=90)
        }
        
        
      }
      
      if (type == "SFORB") {
        # FOCUS kinetics (2006), p. 60 f
        k_out_names = grep(paste("k", obs_var, "free", sep="_"), names(parms.all), value=TRUE)
        k_out_names = setdiff(k_out_names, paste("k", obs_var, "free", "bound", sep="_"))
        k_1output = sum(parms.all[k_out_names])
        k_12 = parms.all[paste("k", obs_var, "free", "bound", sep="_")]
        k_21 = parms.all[paste("k", obs_var, "bound", "free", sep="_")]
        
        sqrt_exp = sqrt(1/4 * (k_12 + k_21 + k_1output)^2 + k_12 * k_21 - (k_12 + k_1output) * k_21)
        b1 = 0.5 * (k_12 + k_21 + k_1output) + sqrt_exp
        b2 = 0.5 * (k_12 + k_21 + k_1output) - sqrt_exp
        
        SFORB_fraction = function(t) {
          ((k_12 + k_21 - b1)/(b2 - b1)) * exp(-b1 * t) +
            ((k_12 + k_21 - b2)/(b1 - b2)) * exp(-b2 * t)
        }
        f_50 <- function(t) (SFORB_fraction(t) - 0.5)^2
        max_DT <- 1000
        DT50 <- optimize(f_50, c(0.01, max_DT))$minimum
        f_90 <- function(t) (SFORB_fraction(t) - 0.1)^2
        DT90 <- optimize(f_90, c(0.01, 1000))$minimum
      }
      #fit$distimes[obs_var, ] = c(ifelse(is.na(DT50),NA,as.numeric(formatC(DT50,4,format='f'))), 
      #                           ifelse(is.na(DT90),NA,as.numeric(formatC(DT90,4,format='f'))),type)#c(DT50, DT90,type)
      out$distimes[obs_var,1:2] <- c(DT50,DT90)
      out$distimes[obs_var,3] <- type
    }
  }
  if(mkinmodini$outpartri=='water-sediment'){
    ## Now deals with multiple ks and additionally calculate formation fractions in case the outpartri is 'water-sediment'##
    for (obs_var in obs_vars) {
      ##f_tot <- grep(paste(obs_var, "_",sep=''), names(fit$ff), value=TRUE)
      ##f_exp <- grep(paste(obs_var, "to",obs_var,sep='_'), names(fit$ff), value=TRUE)
      ##f_exp1 <- grep(paste(obs_var, "to",'sink',sep='_'), names(fit$ff), value=TRUE)
      ##fit$ff[[paste(obs_var,'to', "sink", sep="_")]] = 1 - sum(fit$ff[f_tot])+sum(fit$ff[f_exp])+sum(fit$ff[f_exp1])
      type = names(mkinmodini$map[[obs_var]])[1]
      k1name <- paste("k1", obs_var,  sep="_")
      k2name <- paste("k2", obs_var,  sep="_")
      gname <- paste("g", obs_var,  sep="_")
      tbname <- paste("tb", obs_var,  sep="_")
      alphaname <- paste("alpha", obs_var,  sep="_")
      betaname <- paste("beta", obs_var,  sep="_")
      if (type == "SFO") {
        k_names = grep(paste("k_", obs_var,'_', sep=""), names(parms.all), value=TRUE)
        k_tot = sum(parms.all[k_names])
        DT50 = log(2)/k_tot
        DT90 = log(10)/k_tot
      }
      if (type == "FOMC") {
        alpha = parms.all[alphaname]
        beta = parms.all[betaname]
        DT50 = beta * (2^(1/alpha) - 1)
        DT90 = beta * (10^(1/alpha) - 1)
        
      }
      if (type == "DFOP") {
       
        k1 = parms.all[k1name]
        k2 = parms.all[k2name]
        g = parms.all[gname]
        f <- function(t, x) {
          ((g * exp( - k1 * t) + (1 - g) * exp( - k2 * t)) - (1 - x/100))^2
        }
        DTmax1 <- log(2)/min(k1,k2)
        DTmax <- 1000
        DT50.o <- optimize(f, c(0, DTmax), x=50)$minimum
        DT50.o1 <- optimize(f, c(0, DTmax1), x=50)$minimum
        DT50.o <- ifelse(f(DT50.o,50)>f(DT50.o1,50), DT50.o1,DT50.o)
        ##DT50 = ifelse(DTmax - DT50.o < 0.1, NA, DT50.o)
        DT50 <- DT50.o
        DT90.o <- optimize(f, c(0, DTmax), x=90)$minimum
        DTmax1 <- log(10)/min(k1,k2)
        DT90.o1 <- optimize(f, c(0, DTmax1), x=90)$minimum
        DT90.o <- ifelse(f(DT90.o,90)>f(DT90.o1,90), DT90.o1,DT90.o)
        ##DT90 = ifelse(DTmax - DT90.o < 0.1, NA, DT90.o)
        DT90 <- DT90.o
        
      }
      if (type == "HS") {
        k1 = parms.all[k1name]
        k2 = parms.all[k2name]
        tb = parms.all[tbname]
        HSendpoints <- "analytical"
        if(HSendpoints=="iterative"){
          f <- function(t, x) {
            fraction = ifelse(t <= tb, exp(-k1 * t), exp(-k1 * tb) * exp(-k2 * (t - tb)))
            (fraction - (1 - x/100))^2
          }
          DTmax <- 1000
          hso1 <- nlminb(0.0001,f, x=50)
          hso2 <- nlminb(tb,f, x=50)
          DT50.o <- ifelse(hso1$objective<=hso2$objective,hso1$par,hso2$par)
          DT50 = ifelse(DTmax - DT50.o < 0.1, NA, DT50.o)
          
          hso1 <- nlminb(0.0001,f, x=90)
          hso2 <- nlminb(tb,f, x=90)
          DT90.o <- ifelse(hso1$objective<=hso2$objective,hso1$par,hso2$par)
          DT90 = ifelse(DTmax - DT90.o < 0.1, NA, DT90.o)
        }else{## HSendpoints=="analytical"
          calcDTx <- function(k1,k2,tb,x){
            DTx1 <- log(100/(100-x))/k1
            DTx2 <- tb+(log(100/(100-x))-k1*tb)/k2
            
            if(DTx1 > tb){
              return(DTx2)
            }else{
              return(DTx1)
              #                 if(exp(-k1*DTx1)>=(1-x/100)){
              #                   ## fraction greater than x% then return DTx1
              #                   return(DTx1)
              #                 }else{
              #                   return(DTx2)
              #                 }
            }
          }
          DT50 <- calcDTx(k1,k2,tb,x=50)
          DT90 <- calcDTx(k1,k2,tb,x=90)
        }
      }
      
      if (type == "SFORB") {
        # FOCUS kinetics (2006), p. 60 f
        k_out_names = grep(paste("k", obs_var, "free", sep="_"), names(parms.all), value=TRUE)
        k_out_names = setdiff(k_out_names, paste("k", obs_var, "free", "bound", sep="_"))
        k_1output = sum(parms.all[k_out_names])
        k_12 = parms.all[paste("k", obs_var, "free", "bound", sep="_")]
        k_21 = parms.all[paste("k", obs_var, "bound", "free", sep="_")]
        
        sqrt_exp = sqrt(1/4 * (k_12 + k_21 + k_1output)^2 + k_12 * k_21 - (k_12 + k_1output) * k_21)
        b1 = 0.5 * (k_12 + k_21 + k_1output) + sqrt_exp
        b2 = 0.5 * (k_12 + k_21 + k_1output) - sqrt_exp
        
        SFORB_fraction = function(t) {
          ((k_12 + k_21 - b1)/(b2 - b1)) * exp(-b1 * t) +
            ((k_12 + k_21 - b2)/(b1 - b2)) * exp(-b2 * t)
        }
        f_50 <- function(t) (SFORB_fraction(t) - 0.5)^2
        max_DT <- 1000
        DT50.o <- optimize(f_50, c(0.01, max_DT))$minimum
        if (abs(DT50.o - max_DT) < 0.01) DT50 = NA else DT50 = DT50.o
        f_90 <- function(t) (SFORB_fraction(t) - 0.1)^2
        DT90.o <- optimize(f_90, c(0.01, 1000))$minimum
        if (abs(DT90.o - max_DT) < 0.01) DT90 = NA else DT90 = DT90.o
      }
      #fit$distimes[obs_var, ] = c(ifelse(is.na(DT50),NA,formatC(DT50,4,format='f')), 
      #                            ifelse(is.na(DT90),NA,formatC(DT90,4,format='f')),
      #                            type)#c(DT50, DT90,type)
      out$distimes[obs_var,1:2] <- c(DT50,DT90)
      out$distimes[obs_var,3] <- type
    }
    
  }
  
  ## Collect observed, predicted and residuals
  observed0 <-  mkin_wide_to_long(mkinmodini$data0,time='time')
  observed0$err <- observed$err
  data <- merge(observed, predicted_long, by = c("time", "name"))
  data0 <- merge(observed0, predicted_long, by = c("time", "name"))
  ##names(data) <- c("time", "variable", "observed", "predicted")
  names(data) <- c("time", "variable", "observed","err-std", "predicted")
  ##browser()
  names(data0) <- c("time", "variable", "observed","err-std", "predicted")
  data$residual <- data$observed - data$predicted
  data0$residual <- data$residual
  data$variable <- ordered(data$variable, levels = obs_vars)
  data0$variable <- data$variable
  tmpid <- is.na(data0$residual) & !is.na(data0$observed)
  data0$'err-std'[tmpid] <- 0
  out$data <- data[order(data$variable, data$time), ]
  out$data0 <- data0[order(data0$variable, data0$time), ]
  ## Return Additional Statistics
  
  ## Return Data
  out$inpartri <- inpartri
  out$outpartri <- outpartri
  if(out$outpartri=="default") out$fixed0 <- out$fixed
  out$evalMethod <- NULL
  out$mess <- "The fitting of the specified kinetic model is not done."
  out$modelmess <- mkinmodini$modelmess
  class(out) <- "summary.mkinmod.full"
  return(out)
}
##################
##' Control parameters for MCMC method 
##'
##' Details see \code{modMCMC}
##' @title Generate MCMC control parameters
##' @param ... a list that could be in the list of MCMC control parameters
##' @return a list of control parameters for MCMC evaluation method
##' @author Zhenglei Gao
##' @export
MCMCcontrol <- function(...){
  niter <- 1000
  default <- list(jump = NULL,prior = NULL,
                  wvar0=0.1,niter = niter,
                  outputlength = niter, burninlength = 0, updatecov = niter,
                  ntrydr = 1, drscale = NULL, verbose = TRUE)
  spec <- list(...)
  for(nam in names(spec))
    default[[nam]] <- spec[[nam]]
  return(default)
}


##################
##' Calculating Covariance when LM returns non-positive definite Hessian
##'
##' Details
##' @title Estimate Covariance Matrix with Boundaries
##' @param par parameter vector
##' @param f the residual function
##' @param fval current value
##' @param lb lower bounds
##' @param ub upper bounds
##' @return a covariance matrix
##' @author Zhenglei Gao
##' @export
calcCovar <- function(par,f,fval,lb,ub,...){
  ## like in lsqnonlin function, they calculated the Jacobian with boundary adjustment.
  JACOB <- fdjacobian(func=f,par,fval=fval,lb=lb,ub=ub,...)
  hessian <- t(JACOB)%*%(JACOB)
  covar <-  try(solve(hessian), silent = TRUE)
  if(class(covar)=="try-error") covar <- try(ginv(hessian), silent = TRUE)
  return(covar)
}

##' check boundary reached 
##' 
##' @title check boundary reached
##' @param par fitted parameter vector
##' @param lower lower boundary
##' @param upper upper boundary
##' @return TRUE or FALSE
##' @author Zhenglei Gao
atBoundary <- function(par,lower,upper){
  return(any(par==lower | par==upper))
}

##' Lists model equations, the chi2 error levels
##' calculated according to FOCUS guidance (2006)
##' and optionally the data, consisting of observed, predicted and residual values, the
##' correlation matrix and so on.
##'
##' Expanded from \code{\link{summary.modFit}} and \code{\link{summary.mkinfit}}
##' @title S3 summary method for class 'kingui'
##' @method summary kingui
##' @S3method summary kingui
##' @param object A fitted object of class 'kingui' from the result of
##' \code{\link{mkinfit.full}} or  \code{\link{IRLSkinfit.full}}
##' @param data   If TRUE, include in the returned values a data frame containing the observed
##' and predicted values with residuals and estimated standard deviations or weights.
##' @param distimes If TRUE, DT50 and DT90 values should be included.
##' @param ff  If TRUE, the formation fraction should be calculated from the estimated
##' transformed parameters.
##' @param cov If TRUE, parameter covariances should be calculated.
##' @param version A version number indicating which version of the fit function has been used.
##' @param ... Optional arguments passed to methods like 'print'
##' @return The summary function returns a list of components from the results of the optimizat
##' ion routine used.
##' \item{par}{The fitted parameter values with corresponding standard error and t-test p
##' values.}
##' \item{version}{The version of function used.}
##' \item{pnames}{Parameter names.}
##' \item{diffs}{The differential equations used in the model.}
##' \item{data}{A data frame that will be printed.}
##' \item{start}{The starting values and bounds for the parameters.}
##' \item{fixed}{The fixed values for model parameters. }
##' \item{start0}{The starting values for the transformed parameters.}
##' \item{fixed0}{The fixed values for the transformed parameters.}
##' \item{errmin}{Chi2 error level and other related statistics for goodness-of-fit assessment.
##' }
##' \item{distimes}{The DT50 and DT90 values.}
##' \item{ff}{The formation fractions calculated from the fitted parameters.}
##' \item{outpartri}{Output parameterization.}
##' \item{optimmethod}{Optimization method used in the final fitting step and in calculating
##' Hessian. }
##' \item{residuals}{The residuals.}
##' \item{residualVariance}{The variance of the residuals.}
##' \item{sigma}{The standard deviation of the residuals.}
##' \item{df}{Degree of freedom for the residuals}
##' \item{cov.unscaled}{Unscaled Covariance}
##' \item{cov.scaled}{Scaled Covariance.}
##' \item{info}{Information inherated from the fitting procedure}
##' \item{niter}{Number of iterations.}
##' \item{stopmess}{Warning message.}
##' @author Zhenglei Gao
##' @rdname summary.kingui
summary.kingui <- function(object, data = TRUE, distimes = TRUE, ff = TRUE, cov = FALSE,version=NULL,...) {
  options(warn=-1)
  param  <- object$par
  pnames <- names(param)
  p      <- length(param)
  covar  <- object$covar
  #rownames(covar) <- colnames(covar) <-pnames
  rdf    <- object$df.residual
  resvar <- object$ssr / rdf
  se     <- sqrt(diag(covar) * resvar)
  #browser()
  lci <- param-qnorm(0.975)*se
  uci <- param+qnorm(0.975)*se
  names(se) <- pnames
  tval      <- param / se
  modVariance <- object$ssr / length(object$residuals)
  
  if (!all(object$start$lower >=0)) {
    message <- "Note that the one-sided t-test may not be appropriate if
      parameter values below zero are possible."
        warning(message)
  } else message <- "ok"
  
  param <- cbind(param, se, lci,uci, pt(tval, rdf, lower.tail = FALSE))
  ## adding one line when there is fixed by KinGUII object
  pnames1 <- pnames
  
#  fid <- which(object$fixed$by=='KinGUII')
#   if(length(fid)>0){
#     ## There is some reason that KinGUII fixed the parameters!!!!!!!
#     
#     nfid <- length(fid)
#     param <- rbind(param,(cbind(object$fixed$value[fid],matrix(NA,nfid,4))))
#     pnames1 <- c(pnames1, rownames(object$fixed)[fid])
#   }
  dimnames(param) <- list(pnames1, c("Estimate", "Std. Error",'Lower CI','Upper CI',
                                     "Pr(>t)"))
  ##browser()
  if(cov)
    ans <- list(residuals = object$residuals,
                residualVariance = resvar,
                sigma = sqrt(resvar),
                modVariance = modVariance,
                df = c(p, rdf), cov.unscaled = covar,
                cov.scaled = covar * resvar,
                info = object$info, niter = object$iterations,
                stopmess = message,
                par = param,version=version,pnames=pnames)
  else
    ans <- list(residuals = object$residuals,
                residualVariance = resvar,
                sigma = sqrt(resvar),
                modVariance = modVariance,
                df = c(p, rdf),
                info = object$info, niter = object$iterations,
                stopmess = message,
                par = param,version=version,pnames=pnames)
  
  ans$diffs <- object$diffs
  if(data) {
    if(!is.null(object$data0)) ans$data <- object$data0 else ans$data <- object$data
  }
  ans$start <- object$start
  ans$fixed <- object$fixed
  ans$start0 <- object$start0
  ans$fixed0 <- object$fixed0
  ans$errmin <- object$errmin
  if(distimes) ans$distimes <- object$distimes
  if(ff) {## of course we should always report formation fractions.
    ans$ff <- object$ff ## if report formation fractions.
    ##browser()
    fixed <- object$fixed[,1]
    names(fixed) <- rownames(object$fixed)
    if(length(object$trff)>0){
      tff <- MCstd(mu=object$par,sigma=se,transformations=object$trff,fixed=fixed,N=10000)
      ext <- setdiff(names(ans$ff),rownames(tff))
      ans$ff <- rbind(formatC(tff,digits=4,format='f'),cbind(formatC(ans$ff[ext],digits=4),matrix("",nrow=length(ext),ncol=3)))
    }else{
      ans$ff <- cbind(formatC(ans$ff,digits=4),matrix("",nrow=length(ans$ff),ncol=3))
      colnames(ans$ff) <- c("Estimate", "Std. Error",'Lower CI','Upper CI')
    }
    
  }
  ans$outpartri <- object$outpartri
  ans$optimmethod <- object$optimmethod
  ans$evalMethod <- object$evalMethod
  ans$mess <- object$mess
  ans$modelmess <- object$modelmess
  class(ans) <- c('summary.kingui', "summary.modFit")
  return(ans)
}
######################
##' Formatting a data frame by controlling the number of digits after the decimal point
##'
##' If the data frame contains numeric columns, formmating those columns.
##' @title Formating a data frame
##' @param x A data frame
##' @param digits How many digits should be printed after the decimal point.
##' @param ... Other parameters to be passed into  \code{\link{format}}
##' @return A formatted data frame
##' @author Zhenglei Gao
##' @keywords format
##' @export
myformat <- function(x,digits=4,...)
{
  ## x is a data frame
  for(i in 1:ncol(x))
  {
    if(is.numeric(x[,i]))
    {
      x[,i] <- formatC(x[,i],digits=digits,format='f')
    }else x[,i] <- format(x[,i],digits=digits,...)
  }
  return(x)
}

##' Print method for \code{\link{summary.kingui}}
##'
##' Expanded from \code{\link{print.summary.modFit}} and \code{\link{print.summary.mkinfit}}
##' @title Print method for \code{\link{summary.kingui}}
##' @method print summary.kingui
##' @param x An object of class \code{summary.kingui}
##' @param digits How many digits should be printed after the decimal point.
##' @param detailed Ignore for now.
##' @param ... Other parameters to be passed into \code{\link{summary.kingui}}
##' @return \code{NULL}
##' @author Zhenglei Gao
##' @S3method print summary.kingui
##' @rdname print.summary.kingui
print.summary.kingui <- function(x, digits = max(3, getOption("digits") - 3),detailed=FALSE, ...) {
  
  cat(paste("KineticEval Package Version:",packageVersion("KineticEval"),"\n"))
  if(!is.null(x$version)) cat(paste('Version:',x$version,'\n'))
  
  cat(paste('\n',sessionInfo()$R.version$version.string,'\n',sep=""))
  cat(paste("\nMethod:",x$evalMethod,"\n"))
  if(x$outpartri=='water-sediment')  cat("\nStudy:Water-Sediment\n")
  
  xx <- x[["diffs"]]
  
  cat("\nOptimization Algorithms Used:\n")
  print(x$optimmethod)
  cat("\n##Comment:\n")
  print(ifelse(is.null(x$mess),"",(x$mess)))
  cat("\nEquations:\n")
  for(i in 1:length(xx)) print(noquote(as.character(xx[i])))
  df  <- x$df
  rdf <- df[2]
  
  cat("\nStarting values for optimised parameters:\n")
  if(x$outpartri=='default') print(x$start0) else print(x$start)
  
  cat("\nFixed parameter values:\n")
  if(x$outpartri=='default') {
    if(length(x$fixed$value) == 0) cat("None\n")
    else print(x$fixed)
  }else{
    if(length(x$fixed$value) == 0) cat("None\n")
    else print(x$fixed)
  }
  
  fid <- which(x$fixed$by=='KinGUII')
  if(length(fid)>0){
    ## There is some reason that KinGUII fixed the parameters!!!!!!!
    #cat("\nWarning:Parameters fixed by KinGUII. Please check your log file and model set up. One or more of the SINKs are turned off but you did not fix the formation fractions to be 1.\n")
    cat("\nWarning: One or more compartments without sink. Parameters fixed by KinGUII.\n")
    if(!is.null(x$modelmess)) {
      cat("\nModel Message:\n")
      for(i in 1:length(x$modelmess)) cat(c(x$modelmess[i],"\n"))
    }
  }
  
  
  cat("\nOptimised parameters:\n")
  printCoefmat(x$par, digits = digits, ...)
  
  cat("\nResidual standard error:",
      format(signif(x$sigma, digits)), "on", rdf, "degrees of freedom\n")
  
  printff <- !is.null(x$ff)
  if(printff){
    cat("\nEstimated formation fractions:\n")
    print(as.data.frame(x$ff), digits=digits,...)
  }
  
  
  
  
  cat("\nChi2 error levels in percent :\n")
  print(x$errmin[,1:3], digits=digits,...)
  
  if(any(is.nan(x$errmin[,1]))){
    cat("\nNaNs in Chi2 errors: Replicates should be averaged before the calculation of the Chi2 err value. However, there are fewer time points compared to the number of parameters.Degree of freedom is smaller than 0, therefore no Chi2 error can be calculated.\n")
  }
  
  printdistimes <- !is.null(x$distimes)
  options("scipen"=10)
  if(printdistimes){
    cat("\nEstimated disappearance times:\n")
    print(x$distimes, digits=digits,...)
  }
  options("scipen"=0)
  
  
  
  
  printcor <- (!is.null(x$cov.unscaled))
  if (printcor){
    cat("\nAdditional Statistics:\n")
    print(x$errmin[,4:ncol(x$errmin)], digits=digits,...)
    Corr <- try(cov2cor(x$cov.unscaled),silent=TRUE)
    if(!is.numeric(Corr))
    {
      warning('Covariance matrix cannot be calculated.')
      np <- nrow(x$cov.unscaled)
      Corr <- matrix(data = NA, nrow = np, ncol = np)
    }
    rownames(Corr) <- colnames(Corr) <- x$pnames#rownames(x$par)
    cat("\nParameter correlation:\n")
    print(Corr, digits = digits, ...)
  }
  
  
  printdata <- !is.null(x$data)
  if (printdata){
    cat("\nData:\n")
    print(myformat(x$data, digits = digits, scientific=F,...), row.names = FALSE)
    #sprintf("%.4f",x$data)
  }
  
  invisible(x)
}
########################################
## utility functions
########################################

##' Calculating the Model Efficiency according to Focus Kinetics Guidelines.
##'
##' For EF>0, the values gives an indication of fraction of the total variance of the data
##' that can be explained by the model.
##' @title Calculating the Model Efficiency
##' @param o Observed values.
##' @param c Predicted values using model.
##' @return Model efficiency(EF) ranges from negative infinity to 1.
##' @author Zhenglei Gao
##' @export
EF <- function(o,c)
{
  naid <- !is.na(o)## Remove all the NA values
  c <- c[naid]
  o <- o[naid]
  if(length(c)>0) out <- 1-(sum((c-o)^2))/(sum((o-mean(o))^2)) else out <- NA
  out
}
##' Calculating the scaled root mean squared error according to Focus Kinetics Guidelines.
##'
##'
##' @title Calculating the root mean squared error.
##' @param obs Observed values.
##' @param pred Predicted values using model.
##' @return The root mean squared error.
##' @author Zhenglei Gao
##' @export
rmse <- function(obs, pred) sqrt(mean((obs-pred)^2,na.rm=TRUE))

##' This function calculates the smallest relative error resulting
##' in passing the chi-squared test as defined in the FOCUS kinetics
##' report from 2006.
##'
##'  This function is used internally by \code{\link{mkinfit.full}},
##' \code{\link{IRLSkinfit.full}}, and  \code{\link{mcmckinfit.full}}.
##' @title Calculate model error at which chi-square test is passed.
##' @usage chi2err(errdata, n.parms, errobserved, alpha = 0.05)
##' @param errdata  A data frame with mean observed values in column named 'value_mean'
##' and predicted values in column 'value_pred'.
##' @param n.parms The number of optimized parameters to be taken into account for
##' the data series.
##' @param errobserved A data frame with the true observed values instead of the
##' mean values for repeated measurements.
##' @param alpha The confidence level chosen for the chi-squared test.
##' @return A list of component to assess the model goodness-of-fit including:
##' \item{err.min}{Chi2 error level}
##' \item{n.optim}{Number of parameters}
##' \item{df}{Degree of freedom}
##' \item{RMSE}{The root mean squared error.}
##' \item{EF}{Model efficiency}
##' \item{R2}{Coefficient of determination}
##' @author Zhenglei Gao
##' @export
chi2err <- function(errdata, n.parms, errobserved,alpha = 0.05)
{
  if("logging" %in% loadedNamespaces()){
    logall <- TRUE
  }else{
    logall <- FALSE
  }
  means.mean <- mean(errdata$value_mean, na.rm = TRUE)
  df = length(errdata$value_mean) -sum(is.na(errdata$value_mean)) - n.parms
  C <- errdata$value_pred
  O <- errdata$value_mean
  b <-  sum((C-O)^2,na.rm=TRUE)
  a <- b/(means.mean)^2
  #a <- sum((C-O)^2/(means.mean)^2,na.rm=TRUE)
  if(df<=0){
    if(logall) logwarn("Replicates should be averaged before the calculation of the Chi2 err value. 
                       However, there are fewer time points compared to the number of parameters.
                       Degree of freedom is smaller than 0.")
    err.sig <- NaN
    err.min <- NaN
  }else{
    err.sig <-sqrt(b/qchisq(1-alpha,df))
    err.min <- sqrt(1/qchisq(1-alpha,df)*a)
  }
  
  O1 <- errobserved$value_obs##O1 <- errdata$value_obs
  C1 <- errobserved$value_pred#########
  ef <- EF(O1,C1)
  RMSE <- rmse(O1,C1)
  R2 <- cor(O1,C1,use='na.or.complete')^2
  return(list(err.min = 100*err.min, n.optim = n.parms, df = df,err.sig=err.sig,RMSE=RMSE,EF=ef,R2=R2))
}
