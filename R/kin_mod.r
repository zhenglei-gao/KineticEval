## Define the model cost function, should be the same for all kinds of optimization ##
##' calculte the predicted/estimated residue values for kinetic models
##'
##' @title predicted residue values using defined kinetic model
##' @param P a vector of parameters
##' @param inside logical, whether the function is evaluated inside a fit function.
##' @param plot logical, whether to make a plot
##' @param plottitle the title of the plot
##' @param pnames0 parameter names nls function will forget them!
##' @param ... other parameters to be passed into plot
##' @return the calculated model values
##' @author Zhenglei Gao
##' @export
kin_mod <- function(P,inside=FALSE,plot=TRUE,plottitle,pnames0=NULL, ...)
{
  ##print(names(P))
  if(length(P)==0) warning("No parameters need to be optimized?")
  if(is.null(names(P))) {
   if(exists("pnames")) pnames0 <- pnames
   names(P) <- pnames0
  }
  if("logging" %in% loadedNamespaces()){
    logall <- TRUE
  }else{
    logall <- FALSE
  }
  ## names(P) <- pnames
  ## P is the parameter vector with names.
  if(exists('calls')) calls <<- calls+1
  if(length(state.ini.optim) > 0) {
    odeini <- c(P[1:length(state.ini.optim)], state.ini.fixed)
    names(odeini) <- c(state.ini.optim.boxnames,state.ini.fixed.boxnames)
  } else {
    odeini <- state.ini.fixed
    names(odeini) <- c( state.ini.fixed.boxnames)
  }
  ## has to change the odeini order since it is different from the mod_vars order.
  odeini <- odeini[mod_vars]
  if(length(parms.optim)>0) {
    odeparms <- c(P[(length(state.ini.optim) + 1):length(P)], parms.fixed)
  }else{
    odeparms <- parms.fixed
  }

  outtimes = unique(observed$time)
  evalparse <- function(string)
  {
    eval(parse(text=string), as.list(c(odeparms, odeini)))
  }

  ## Solve the system
  if (solution == "analytical") {
    parent.type = names(mkinmodini$map[[1]])[1]
    parent.name = names(mkinmodini$diffs)[[1]]
    o <- switch(parent.type,
                SFO = SFO.solution(outtimes,
                                   evalparse(parent.name),
                                   evalparse(paste("k", parent.name,  sep="_"))),
                FOMC = FOMC.solution(outtimes,
                                     evalparse(parent.name),
                                     evalparse(paste("alpha", parent.name,  sep="_")), evalparse(paste("beta", parent.name,  sep="_"))),
                DFOP = DFOP.solution(outtimes,
                                     evalparse(parent.name),
                                     evalparse(paste("k1", parent.name,  sep="_")), evalparse(paste("k2", parent.name,  sep="_")),evalparse(paste("g", parent.name,  sep="_"))),
                HS = HS.solution(outtimes,
                                 evalparse(parent.name),
                                 evalparse(paste("k1", parent.name,  sep="_")),evalparse(paste("k2", parent.name,  sep="_")),evalparse(paste("tb", parent.name,  sep="_"))),
                SFORB = SFORB.solution(outtimes,
                                       evalparse(parent.name),
                                       evalparse(paste("k", parent.name, "bound", sep="_")),
                                       evalparse(paste("k", sub("free", "bound", parent.name), "free", sep="_")),
                                       evalparse(paste("k", parent.name, "sink", sep="_")))
    )
    out <- cbind(outtimes, o)
    dimnames(out) <- list(outtimes, c("time", sub("_free", "", parent.name)))
  }
  if (solution == "eigen") {
    coefmat.num <- matrix(sapply(as.vector(mkinmodini$coefmat), evalparse),
                          nrow = length(mod_vars))
    e <- eigen(coefmat.num)
    zz.ev <- e$values
    if(min(zz.ev)[1]<0){
      ## in this case cannot use the eigen solution
      warning("\'coefmat is not positive definite!\n")
      solution <- 'deSolve' ## switch to deSolve methods,define the mkindiff function
      if(solution == "deSolve") {
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
      }
    }else{
      cc <- solve(e$vectors, odeini)
      f.out <- function(t) {
        e$vectors %*% diag(exp(e$values * t), nrow=length(mod_vars)) %*% cc
      }
      o <- matrix(mapply(f.out, outtimes),
                  nrow = length(mod_vars), ncol = length(outtimes))
      dimnames(o) <- list(mod_vars, outtimes)
      out <- cbind(time = outtimes, t(o))
    }
  }

  if (solution == "deSolve")
  {
    out <- ode(
      y = odeini,
      times = outtimes,
      func = mkindiff,
      parms = odeparms,
      atol = atol,
      rtol = rtol,
      method=odesolver
    )
    ## #sink("NUL")
    ## if(diagnostics(out)$istate[1]!=2)
    ##     {
    ##         #browser()
    ##          out <- ode(
    ##                     y = odeini,
    ##                     times = outtimes,
    ##                     func = mkindiff,
    ##                     parms = odeparms,
    ##                     atol = atol,
    ##                     rtol = rtol,
    ##                     method='ode45'
    ##                     )

    ##     }
    ## #sink()
  }

  ## Output transformation for models with unobserved compartments like SFORB
  out_transformed <- data.frame(time = out[,"time"])
  for (var in names(mkinmodini$map)) {
    if((length(mkinmodini$map[[var]]) == 1) || solution == "analytical") {
      out_transformed[var] <- out[, var]
    } else {
      out_transformed[var] <- rowSums(out[, mkinmodini$map[[var]]])
    }
  }
  assign("out_predicted", out_transformed, inherits=TRUE)
  if(sum(apply(out_transformed,2,function(x) sum(is.nan(x))>nrow(out_transformed)-2))>0)
  {
    ##browser()
    warning('Integration not completed, you may want to try some other solver!')
    out_transformed <- apply(out_transformed,2,function(x) {if(sum(is.nan(x))>nrow(out_transformed)-2) x <- rep(Inf,nrow(out_transformed)) else x <- x})
  }
  if(nrow(out_transformed)<length(outtimes))
  {
    tmp <- matrix(0,length(outtimes),ncol(out_transformed)-1)
    tmpnames <-names(out_transformed)
    out_transformed <- data.frame(time=outtimes,tmp)
    names(out_transformed) <- tmpnames
  }
  ##browser()
  ##mC <- modCost(out_transformed, observed, y = "value",
  ##              err = 'err', weight = "none", scaleVar = FALSE)
  yMod <- rep(NA,length(observed$value))
  Mod <- mkin_wide_to_long(as.data.frame(out_transformed),time="time")
  ## melt(as.data.frame(out_transformed),id.vars="time")
  names(Mod) <- c("name",  "time", "yMod")
  all <- merge(observed,Mod,by.x=c("name","time"),sort=FALSE)
  yMod <- all$yMod
  bady <- is.nan(yMod)
  if(sum(bady)>0){
    yMod[bady]<-0
    if(logall) logwarn("NaN values simulated from solving the ODE, which were replaced with 0s.")
  }
  #############
  
  if(missing(plottitle)) plottitle <- NULL
  ## Report and/or plot if the model is improved
  if(inside==TRUE){
    costcurr <- sum(((all$yMod-all$value)/all$err)^2,na.rm=TRUE)
    if (cost.old-costcurr > ctr$quiet.tol) {
      if(!quiet) cat("Model cost at call ", calls, ": ", costcurr, "\n")
      
      ## Plot the data and current model output if requested
      if(plot) {
        outtimes_plot = seq(min(observed$time), max(observed$time), length.out=100)
        if (solution == "analytical") {
          o_plot <- switch(parent.type,
                           SFO = SFO.solution(outtimes_plot,
                                              evalparse(parent.name),
                                              evalparse(paste("k", parent.name,  sep="_"))),
                           FOMC = FOMC.solution(outtimes_plot,
                                                evalparse(parent.name),
                                                evalparse(paste("alpha", parent.name,  sep="_")),evalparse(paste("beta", parent.name,  sep="_"))),
                           DFOP = DFOP.solution(outtimes_plot,
                                                evalparse(parent.name),
                                                evalparse(paste("k1", parent.name,  sep="_")),evalparse(paste("k2", parent.name,  sep="_")),evalparse(paste("g", parent.name,  sep="_"))),
                           HS = HS.solution(outtimes_plot,
                                            evalparse(parent.name),
                                            evalparse(paste("k1", parent.name,  sep="_")),evalparse(paste("k2", parent.name,  sep="_")),evalparse(paste("tb", parent.name,  sep="_"))),
                           SFORB = SFORB.solution(outtimes_plot,
                                                  evalparse(parent.name),
                                                  evalparse(paste("k", parent.name, "bound", sep="_")),
                                                  evalparse(paste("k", sub("free", "bound", parent.name), "free", sep="_")),
                                                  evalparse(paste("k", parent.name, "sink", sep="_")))
          )
          out_plot <- cbind(outtimes_plot, o_plot)
          dimnames(out_plot) <- list(outtimes_plot, c("time", sub("_free", "", parent.name)))
        }
        if(solution == "eigen") {
          o_plot <- matrix(mapply(f.out, outtimes_plot),
                           nrow = length(mod_vars), ncol = length(outtimes_plot))
          dimnames(o_plot) <- list(mod_vars, outtimes_plot)
          out_plot <- cbind(time = outtimes_plot, t(o_plot))
        }
        if (solution == "deSolve") {
          out_plot <- ode(
            y = odeini,
            times = outtimes_plot,
            func = mkindiff,
            parms = odeparms)
        }
        out_transformed_plot <- data.frame(time = out_plot[,"time"])
        for (var in names(mkinmodini$map)) {
          if((length(mkinmodini$map[[var]]) == 1) || solution == "analytical") {
            out_transformed_plot[var] <- out_plot[, var]
          } else {
            out_transformed_plot[var] <- rowSums(out_plot[, mkinmodini$map[[var]]])
          }
        }
        
        plot(0, type="n",
             xlim = range(observed$time), ylim = range(observed$value, na.rm=TRUE),
             xlab = "Time", ylab = "Observed",main=plottitle)
        col_obs <- pch_obs <- 1:length(obs_vars)
        names(col_obs) <- names(pch_obs) <- obs_vars
        for (obs_var in obs_vars) {
          points(subset(observed, name == obs_var, c(time, value)),
                 pch = pch_obs[obs_var], col = col_obs[obs_var])
        }
        matlines(out_transformed_plot$time, out_transformed_plot[-1])
        legend("topright", inset=c(0.05, 0.05), legend=obs_vars,
               col=col_obs, pch=pch_obs, lty=1:length(pch_obs))
      }
      
      assign("cost.old", costcurr, inherits=TRUE)
    }
  }else{
    ## Plot the data and current model output if requested
    if(plot) {
      outtimes_plot = seq(min(observed$time), max(observed$time), length.out=100)
      if (solution == "analytical") {
        o_plot <- switch(parent.type,
                         SFO = SFO.solution(outtimes_plot,
                                            evalparse(parent.name),
                                            evalparse(paste("k", parent.name,  sep="_"))),
                         FOMC = FOMC.solution(outtimes_plot,
                                              evalparse(parent.name),
                                              evalparse(paste("alpha", parent.name,  sep="_")),evalparse(paste("beta", parent.name,  sep="_"))),
                         DFOP = DFOP.solution(outtimes_plot,
                                              evalparse(parent.name),
                                              evalparse(paste("k1", parent.name,  sep="_")),evalparse(paste("k2", parent.name,  sep="_")),evalparse(paste("g", parent.name,  sep="_"))),
                         HS = HS.solution(outtimes_plot,
                                          evalparse(parent.name),
                                          evalparse(paste("k1", parent.name,  sep="_")),evalparse(paste("k2", parent.name,  sep="_")),evalparse(paste("tb", parent.name,  sep="_"))),
                         SFORB = SFORB.solution(outtimes_plot,
                                                evalparse(parent.name),
                                                evalparse(paste("k", parent.name, "bound", sep="_")),
                                                evalparse(paste("k", sub("free", "bound", parent.name), "free", sep="_")),
                                                evalparse(paste("k", parent.name, "sink", sep="_")))
        )
        out_plot <- cbind(outtimes_plot, o_plot)
        dimnames(out_plot) <- list(outtimes_plot, c("time", sub("_free", "", parent.name)))
      }
      if(solution == "eigen") {
        o_plot <- matrix(mapply(f.out, outtimes_plot),
                         nrow = length(mod_vars), ncol = length(outtimes_plot))
        dimnames(o_plot) <- list(mod_vars, outtimes_plot)
        out_plot <- cbind(time = outtimes_plot, t(o_plot))
      }
      if (solution == "deSolve") {
        out_plot <- ode(
          y = odeini,
          times = outtimes_plot,
          func = mkindiff,
          parms = odeparms)
      }
      out_transformed_plot <- data.frame(time = out_plot[,"time"])
      for (var in names(mkinmodini$map)) {
        if((length(mkinmodini$map[[var]]) == 1) || solution == "analytical") {
          out_transformed_plot[var] <- out_plot[, var]
        } else {
          out_transformed_plot[var] <- rowSums(out_plot[, mkinmodini$map[[var]]])
        }
      }
      
      plot(0, type="n",
           xlim = range(observed$time), ylim = range(observed$value, na.rm=TRUE),
           xlab = "Time", ylab = "Observed",main=plottitle)
      col_obs <- pch_obs <- 1:length(obs_vars)
      names(col_obs) <- names(pch_obs) <- obs_vars
      for (obs_var in obs_vars) {
        points(subset(observed, name == obs_var, c(time, value)),
               pch = pch_obs[obs_var], col = col_obs[obs_var])
      }
      matlines(out_transformed_plot$time, out_transformed_plot[-1])
      legend("topright", inset=c(0.05, 0.05), legend=obs_vars,
             col=col_obs, pch=pch_obs, lty=1:length(pch_obs))
    }
    
  }
  
  return(yMod)
  
}
