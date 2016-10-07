##' Function to generate graph data for the GUI to make fit plot.
##'
##' A data table is generated and stored in file=filename for later usage.
##'
##' @param fit An object of class 'kingui'.
##' @param filename The file in which the
##' graph data will be stored.
##' @param xlab xlab
##' @param ylab ylab
##' @param legend logical, whether to have legend
##' @param length.out the number of points to be ploted
##' @param xlim The plot range for x axis.
##' @param ylim The plot range for y axis.
##' @param ... other parameters passed to plot.
##' @return \code{NULL}
##' @note This function only writes the data table for making the graph.
##' @author Zhenglei Gao
##' @seealso \code{\link{kinplot}}
##' @export
##'
kingraph <- function(fit, filename='graphdata.txt',xlab = "Time", ylab = "Observed", 
					 xlim = c(1,1.05)*range(fit$data$time), ylim = c(1,1.05)*range((summary(fit))$data$observed, na.rm = TRUE), 
					 legend = TRUE, length.out=100,...)
{
    ## function to generate the graph data.
    ## example usage:
    ## kinplot(fit0,name=c('MESO','M851','M459','M460','M095','M944'))

    ## change ylim if necessary. No need since we use summary in the arguments
    ## if(!is.null(fit$data0)) ylim = c(1,1.05)*range(fit$data0$observed, na.rm = TRUE)
  solution = fit$solution
  fixed <- fit$fixed$value
  names(fixed) <- rownames(fit$fixed)
  parms.all <- c(fit$par, fixed)
  ininames <- c(
    rownames(subset(fit$start, type == "state")),
    rownames(subset(fit$fixed, type == "state")))
  ## Note that 
  ## has to change the odeini order since it is different from the mod_vars order.
  ## although it is already taken care of in the cost(P) function
  # to reorder the odeini using mod_vars
  odeini <- parms.all[ininames]
  odeini <- odeini[paste0("M0_",names(fit$diffs))]
  names(odeini) <- names(fit$diffs)

  outtimes <- seq(xlim[1], xlim[2], length.out=length.out)

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
      atol = fit$atol
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

write.table(out_transformed,filename,sep='\t')

}
