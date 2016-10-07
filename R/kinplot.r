##' Function to plot the fit results for the whole model or selected compartments.
##'
##' @title Plot the fit results
##' @param fit An object of class 'kingui'.
##' @param name Names of the compartmens that should be plotted.
##' @param xlab xlab name.
##' @param ylab ylab name.
##' @param xlim The plot range for x axis.
##' @param ylim The plot range for y axis.
##' @param legend whether to include legends.
##' @param pdf whether to put the plots into one pdf file.
##' @param ... optional arguments. Not used for now.
##' @return \code{NULL}
##' @note This function makes the plots.
##' @author Zhenglei Gao
##' @seealso \code{\link{kingraph}}
##' @export
##'
kinplot <- function(fit,name=NULL, xlab = "Time", ylab = "Observed", xlim = c(1,1.1)*range(fit$data$time), ylim = c(1,1.1)*range(fit$data$observed, na.rm = TRUE), legend = TRUE,pdf=FALSE, ...)
{
    ## example usage:
    ## kinplot(fit0,name=c('MESO','M851','M459','M460','M095','M944'))

  if(is.null(fit$solution)) solution <- "deSolve" else solution <- fit$solution
  fixed <- fit$fixed$value
  names(fixed) <- rownames(fit$fixed)
  parms.all <- c(fit$par, fixed)
  ininames <- c(
    rownames(subset(fit$start, type == "state")),
    rownames(subset(fit$fixed, type == "state")))
  odeini <- parms.all[ininames]
  odeini <- odeini[paste0("M0_",names(fit$diffs))]
  names(odeini) <- names(fit$diffs)

  outtimes <- seq(xlim[1], xlim[2], length.out=100)
  
  odenames <- c(
    rownames(subset(fit$start, type == "deparm")),
    rownames(subset(fit$fixed, type == "deparm")))
  odeparms <- parms.all[odenames]

  # Solve the system
  evalparse <- function(string)
  {
    eval(parse(text=string), as.list(c(odeparms, odeini)))
  }
  if (solution == "analytical") {
    parent.type = names(fit$map[[1]])[1]
    parent.name = names(fit$diffs)[[1]]
    o <- switch(parent.type,
                SFO = SFO.solution(outtimes,
                                   evalparse(parent.name),
                                   evalparse(paste("k", parent.name, sep="_"))),
                FOMC = FOMC.solution(outtimes,
                                     evalparse(parent.name),
                                     evalparse(paste("alpha", parent.name, sep="_")), evalparse(paste("beta", parent.name, sep="_"))),
                DFOP = DFOP.solution(outtimes,
                                     evalparse(parent.name),
                                     evalparse(paste("k1", parent.name, sep="_")), evalparse(paste("k2", parent.name, sep="_")), evalparse(paste("g", parent.name, sep="_"))),
                HS = HS.solution(outtimes,
                                 evalparse(parent.name),
                                 evalparse(paste("k1", parent.name, sep="_")), evalparse(paste("k2", parent.name, sep="_")), evalparse(paste("tb", parent.name, sep="_"))),
                SFORB = SFORB.solution(outtimes,
                                       evalparse(parent.name),
                                       evalparse(paste("k", parent.name, "free_bound", sep="_")),
                                       evalparse(paste("k", parent.name, "bound_free", sep="_")),
                                       evalparse(paste("k", parent.name, "free_sink", sep="_")))
                )
    out <- cbind(outtimes, o)
    dimnames(out) <- list(outtimes, c("time", parent.name))
  }
  if (solution == "eigen") {
    coefmat.num <- matrix(sapply(as.vector(fit$coefmat), evalparse),
                          nrow = length(odeini))
    e <- eigen(coefmat.num)
    c <- solve(e$vectors, odeini)
    f.out <- function(t) {
      e$vectors %*% diag(exp(e$values * t), nrow=length(odeini)) %*% c
    }
    o <- matrix(mapply(f.out, outtimes),
                nrow = length(odeini), ncol = length(outtimes))
    dimnames(o) <- list(names(odeini), NULL)
    out <- cbind(time = outtimes, t(o))
  }

  if (solution == "deSolve") {
    out <- ode(
      y = odeini,
      times = outtimes,
      func = fit$mkindiff,
      parms = odeparms,
      atol = fit$atol,
               method='lsode'
    )
  }

  # Output transformation for models with unobserved compartments like SFORB
  out_transformed <- data.frame(time = out[,"time"])
  for (var in names(fit$map)) {
    if(length(fit$map[[var]]) == 1) {
      out_transformed[var] <- out[, var]
    } else {
      out_transformed[var] <- rowSums(out[, fit$map[[var]]])
    }
  }
if(is.null(name)){
  # Plot the data and model output
  plot(0, type="n",
    xlim = xlim, ylim = ylim,
    xlab = xlab, ylab = ylab, ...)
  col_obs <- pch_obs <- 1:length(fit$map)
  names(col_obs) <- names(pch_obs) <- names(fit$map)
  for (obs_var in names(fit$map)) {
    points(subset(fit$data, variable == obs_var, c(time, observed)),
    pch = pch_obs[obs_var], col = col_obs[obs_var])
  }
  matlines(out_transformed$time, out_transformed[-1],lwd=2.5)
  if (legend == TRUE) {
    legend("topright", inset=c(0.05, 0.05), legend=names(fit$map),
      col=col_obs, pch=pch_obs, lty=1:length(pch_obs))
  }
}else{
    col_obs <- pch_obs <- 1:length(fit$map)
    names(col_obs) <- names(pch_obs) <- names(fit$map)
    if(pdf==F)
    {
        for (obs_var in name)
        {

            x11()
            plot(subset(fit$data, variable == obs_var, c(time, observed)),
                 pch = pch_obs[obs_var],xlim=xlim,ylim=c(1,1.1)*range(subset(fit$data, variable == obs_var, c( observed,predicted)),na.rm=TRUE), col = col_obs[obs_var],main=obs_var)
        #browser()
            lines(out_transformed$time,as.vector(out_transformed[obs_var][[1]]))

        }
    }else{
        for (obs_var in name)
        {
            plot(subset(fit$data, variable == obs_var, c(time, observed)),
                 pch = pch_obs[obs_var],xlim=xlim,ylim=c(1,1.1)*range(subset(fit$data, variable == obs_var, c( observed,predicted)),na.rm=TRUE), col = col_obs[obs_var],main=obs_var)
        #browser()
            lines(out_transformed$time,as.vector(out_transformed[obs_var][[1]]))
        }
    }
}

}
