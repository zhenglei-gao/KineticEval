##' Function to update kinetic model and do a comparison of two parameter estimations.
##'
##' Inspired by Andrew's data example where he gets different results from KinGUI and KinGUII. Not a good sign!
##' Only compare two models! not more than two!! so either newparms is null or Fit is null.
##' @title update kinetic model with better parameter.
##' @param mkinmodini a list of class mkinmod
##' @param newparms a vector of parameter values
##' @param Fit an object of class kingui
##' @param compareChi2 logical, whether to compare chi2.
##' @param eigen logical,whether to use eigen method
##' @param odesolver which odesolver to use
##' @param atol the absolute tolerance level of the ode solver
##' @param rtol the relative tolerance lever of the ode solver.
##' @return a matrix of comparisons.
##' @export
##' @author Zhenglei Gao
update_kinmod <- function(mkinmodini,newparms=NULL,Fit=NULL,compareChi2=FALSE,
                            eigen=FALSE,odesolver='lsoda',atol=1e-9,rtol=1e-10)
{
  environment(kin_mod_cost) <- environment()
  ###############################################

  if(!is.null(Fit) & !is.null(newparms)){
    if(compareFit==FALSE) stop("Function only compares two models at once. Plase set compareFit being TRUE or
try to compare your models subsequently by updating the initial model at every step!")
  }
 ## --------------------------------------------------
  ## Get the parametrization.
  inpartri <- mkinmodini$inpartri
  outpartri <- mkinmodini$outpartri
  ##

  ## mkinmodini is an object by mkinmod.full
  parms.ini <- mkinmodini$parms.ini
  state.ini <- mkinmodini$state.ini
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
  ## NOTE HERE: the order may not be the same as the input mkinmod.full differential equations list. ## XXXXX TODO XXXX Reorder them maybe a good idea if the data is given from a data file while the mkinmod.full is defined not following the colnames order, although it is already taken care of in the cost(P) function to reorder the odeini using mod_vars
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
  state.ini.optim.boxnames <- names(state.ini.optim)
  state.ini.fixed.boxnames <- names(state.ini.fixed)
  if(length(state.ini.optim) > 0) {
    names(state.ini.optim) <- paste('M0',names(state.ini.optim),  sep="_")
  }
  if(length(state.ini.fixed) > 0) {
    names(state.ini.fixed) <- paste('M0',names(state.ini.fixed), sep="_")
  }
  
  
  oldparms <- c(state.ini.optim,parms.optim)
  if(!is.null(newparms)){
    if(length(newparms)!=length(oldparms)){
      stop('The provided parameter vector length is not right! Please check your Kinetic
model settings.')
    }
    
    if(is.null(names(newparms))){
      names(newparms) <- names(oldparms)
    }
  }
  # -----------------------------------------------------------------------
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
    return(list(c(diffs)))
  }

  # -----------------------------------------------------------------------

  oldpar <- par(no.readonly = TRUE)
  if(!is.null(newparms) || !is.null(Fit)) par(mfrow=c(1,2))
  mC0 <- kin_mod_cost(c(state.ini.optim, parms.optim),inside=FALSE,plot=TRUE,plottitle='Initial Model')
  if(!is.null(newparms)) {
    mC1 <- kin_mod_cost(newparms,inside=FALSE,plot=TRUE,plottitle='With New Paramters')
  }
  if(!is.null(Fit)) {
    mC1 <- kin_mod_cost(Fit$par,inside=FALSE,plot=TRUE,plottitle='Fitted Model')
    newparms <- Fit$par
  }

  comparison <- NULL


 ## --------------------------------------------------
  new_mkinmod <- mkinmodini
  if(is.null(newparms) & is.null(Fit)){
    ## No new parameters. Only do the part of the calculations for the initial model.
    cat("The initial model cost is:", mC0$model)

  }else{

    ## Should only report SSR or so.
    comparison$SSR <- cbind(c(mC0$model,mC1$model),rbind(t(mC0$var[c("SSR")]),t(mC1$var[c("SSR")])))
    rownames(comparison$SSR) <- c("Initial Model","Newly Fitted Model")
    colnames(comparison$SSR) <- c("ALL",as.character(mC0$var$name))

    #############
    comparison$SSR.unweighted <- cbind(c(NA,NA),rbind(t(mC0$var[c("SSR.unweighted")]),
                                                      t(mC1$var[c("SSR.unweighted")])))
    rownames(comparison$SSR.unweighted) <- c("Initial Model","Newly Fitted Model")
    colnames(comparison$SSR.unweighted) <- c("ALL",as.character(mC0$var$name))
    if(length(obs_vars==1)){
      comparison$SSR.unweighted[,1] <- comparison$SSR.unweighted[,2]
    }else{
      comparison$SSR.unweighted[,1]<- apply(comparison$SSR.unweighted[,2:ncol(comparison$SSR.unweighted)],1,sum)
    }
    if(compareChi2==TRUE){
      ###########
      ## Need to compare the chi-square error!

    }
    if(mC0$model > mC1$model){
      ##browser()
      new_mkinmod$state.ini[state.ini.optim.boxnames] <- as.numeric(newparms[1:length(state.ini.optim)])
      new_mkinmod$parms.ini[optim_parms] <- newparms[(length(state.ini.optim)+1):length(newparms)]
      #new_mkinmod$start$initial <-  new_mkinmod$parms.ini
    }
  }
  par(oldpar)
  return(list(new_mkinmod=new_mkinmod,comparison=comparison))
}
