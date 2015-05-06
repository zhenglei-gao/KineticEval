##' Fit a kinetic model using MCMC methods.
##'
##' This
##' function does kinetic evaluations using MCMC algorithm with function
##' \code{\link{modMCMC}} in the \code{\link{FME}} package.
##'
##'
##' @aliases mcmckinfit.gui
##' @param mkinmodini  A list of class
##' \code{\link{mkinmod.full}}, containing the kinetic model to be fitted to the
##' data, and the initial parameter values, the observed data.
##' @param eigen If TRUE, the solution of the
##' system of differential equations should be based on the spectral
##' decomposition of the coefficient matrix in cases that this is possible.
##' @param ctr Used when fitstart is TRUE. A
##' list of control values for the estimation algorithm to replac the default
##' values including maximum iterations and absolute error tolerance.  Defaults
##' to the output of \code{\link{kingui.control}}.
##' @param plottitle The title of the
##' plot for visualizing the optimization process.
##' @param plot When fitstart==TRUE, if
##' TRUE,the observed values and the numerical solutions should be plotted at
##' each stage of the first optimization step.
##' @param quiet If TRUE, suppress printing
##' out the current model cost after each(>1) improvement.
##' @param commonsigma If TRUE,the
##' error model has constant error variance and the NLS algorithm will be used
##' for the first optimizatio step.
##' @param jump Jump length, either a number,
##' a vector with length equal to the total number of parameters, a covariance
##' matrix, or a function that takes as input the current values of the
##' parameters and produces as output the perturbed parameters.  See details of
##' \code{\link{modMCMC}}.
##' @param prior Prior probability of the
##' parameters, either a function of the parameters or 'NULL'; in the latter
##' case a flat prior is used (i.e. all parameters are equally likely).
##' @param wvar0 "weight" for the initial
##' model variance. See details of \code{\link{modMCMC}}.
##' @param niter number of iterations for the
##' MCMC.
##' @param outputlength number of
##' iterations kept in the output.
##' @param burninlength number of
##' discarded initial iterations.
##' @param updatecov number of iterations
##' after which the parameter covariance matrix is (re)evaluated based on the
##' parameters kept thus far, and used to update the MCMC jumps.
##' @param ntrydr maximal number of tries
##' for the delayed rejection procedure.
##' @param drscale for each try during
##' delayed rejection, the cholesky decomposition of the proposal matrix is
##' scaled with this amount; if 'NULL', it is assumed to be 'c(0.2,0.25, 0.333,
##' 0.333, ...)'
##' @param verbose if 'TRUE': prints extra
##' output.
##' @param fitstart if 'TRUE': first
##' perform an optimization step and using the fitted parameters as the
##' starting values for MCMC.
##' @param update if not NULL, using the
##' values in the update as the starting values for MCMC.
##' @param err Either \code{NULL}, or the name
##' of the column with the \emph{error} estimates, used to weigh the residuals
##' (see details of \code{\link{modCost}}); if \code{NULL}, then the residuals
##' are not weighed.  In the GUI version, there is no need to consider this
##' argument since a default weight one matrix is setup in
##' \code{mkinmod.gui}. The err argument turned into 'err' automatically
##' in the codes.
##' @param weight only if
##' \code{err}=\code{NULL}: how to weigh the residuals, one of "none", "std",
##' "mean", see details of \code{\link{modCost}}.
##' @param scaleVar Will be passed to
##' \code{\link{modCost}}. Default is not to scale Variables according to the
##' number of observations.
##' @param \dots Further arguments that will
##' be passed to \code{\link{modFit}}.
##' @return  A list with "mcmckingui" and "modMCMC" in the class
##' attribute. A summary can be obtained by \code{\link{summary.mcmckingui}}.
##' @note \code{mcmckinfit.gui} is the deprecated version.
##' @author Zhenglei Gao
##' @seealso \code{\link{modMCMC}}
##' @keywords Kinetic-Evaluations
##' @export
##' @examples
##' \dontrun{
##' SFO_SFO_gui <- mkinmod.gui(Parent = list(type = "SFO", to = "Metab", sink = TRUE,
##'                            k = list(ini = 0.1,
##'                       fixed = 0,
##'                       lower = 0,
##'                       upper = Inf),
##'               M0 = list(ini = 195,
##'                       fixed = 0,
##'                       lower = 0,
##'                       upper = Inf),
##'                            FF = list(ini = c(.1),
##'                       fixed = c(0),
##'                       lower = c(0),
##'                       upper = c(1)),
##'                       time=c(0.0,2.8,   6.2,  12.0,  29.2,  66.8,  99.8,
##' 127.5, 154.4, 229.9, 272.3, 288.1, 322.9),
##'                     residue = c( 157.3, 206.3, 181.4, 223.0, 163.2,
##' 144.7,  85.0,  76.5,  76.4,  51.5,  45.5,  47.3, 42.7)),
##'                            Metab = list(type = "SFO",
##'                            k = list(ini = 0.1   ,
##'                       fixed = 0,
##'                       lower = 0,
##'                       upper = Inf),
##'               M0 = list(ini = 0,
##'                       fixed = 1,
##'                       lower = 0,
##'                       upper = Inf),
##'                     residue =c( 0.0,  0.0,  0.0,  1.6,  4.0, 12.3, 13.5,
##' 12.7, 11.4, 11.6, 10.9,  9.5,  7.6))                           )
##'
##' fit <- mcmckinfit.gui(SFO_SFO_gui)
##'
##' summary(fit)
##' plot(fit)
##' }
##'
mcmckinfit.full <- function (mkinmodini,eigen=FALSE,ctr=kingui.control(quiet.tol=1e-5,method="solnp"),plottitle='',
    plot = FALSE, quiet = FALSE, err = NULL, weight = "none",
    scaleVar = FALSE,commonsigma=FALSE,jump = NULL,prior = NULL,
    wvar0=0.1,niter = 1000,
    outputlength = niter, burninlength = 0, updatecov = niter,
    ntrydr = 1, drscale = NULL, verbose = TRUE,fitstart=TRUE,update=NULL, ...)
{
    ## Get the parametrization.
    inpartri <- mkinmodini$inpartri
    outpartri <- mkinmodini$outpartri

    ## options(warn=-1) ## turn off the warning
    #### Control parameters ####
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
    ##############

    ## mkinmodini is an object by mkinmod.gui
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
    obs_vars = unique(as.character(observed$name))

    ## Name the parameters if they are not named yet
    if(is.null(names(parms.ini))) names(parms.ini) <- mkinmodini$parms

    ## Name the inital parameter values if they are not named yet
    if(is.null(names(state.ini))) names(state.ini) <- mod_vars

    ## Parameters to be optimised
    parms.fixed <- parms.ini[fixed_parms]
    optim_parms <- setdiff(names(parms.ini), fixed_parms)
    parms.optim <- parms.ini[optim_parms]


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
    ## # If updating from previous fit.####
    if(!is.null(update))
    {
        parms.optim <- update$par[optim_parms]
        state.ini.optim <-update$par[names(state.ini.optim)]
         if(!is.null(update$ctr)){
             ctr <- update$ctr
             ## Control parameters ####
             method <- ctr$method
             odesolver <- ctr$odesolver
             atol <- ctr$atol
             rtol <- ctr$rtol
             control <- ctr$control
         }
    }
    ## ####################################################################
    if (length(mkinmodini$map) == 1) {
        solution = "analytical"
    } else {
        if (is.matrix(mkinmodini$coefmat) & eigen) solution = "eigen"
        else solution = "deSolve"
    }

    ## Create a function calculating the differentials specified by the model
    ## if necessary
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


    cost.old <- 1e100
    calls <- 0
    out_predicted <- NA

    ## Define the model cost function
    cost <- function(P)
    {
                                        #names(P) <- pnames
        assign("calls", calls+1, inherits=TRUE)
        if(length(state.ini.optim) > 0) {
            odeini <- c(P[1:length(state.ini.optim)], state.ini.fixed)
            names(odeini) <- c(state.ini.optim.boxnames,state.ini.fixed.boxnames)
        } else {
            odeini <- state.ini.fixed
            names(odeini) <- c( state.ini.fixed.boxnames)
        }
        odeini <- odeini[mod_vars]
        if(length(parms.optim)>0) {
          odeparms <- c(P[(length(state.ini.optim) + 1):length(P)], parms.fixed)
        }else{
          odeparms <- parms.fixed
        }##odeparms <- c(P[(length(state.ini.optim) + 1):length(P)], parms.fixed)

        outtimes = unique(observed$time)
        evalparse <- function(string)
        {
            eval(parse(text=string), as.list(c(odeparms, odeini)))
        }
                                        # Solve the system
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

                warning("\'coefmat is not positive definite!\n")
                solution <- 'deSolve' ## switch to deSolve methods
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
                                        #browser()
            warning('Integration not completed')
            out_transformed <- apply(out_transformed,2,function(x) {if(sum(is.nan(x))>nrow(out_transformed)-2) x <- rep(Inf,nrow(out_transformed)) else x <- x})
        }
        if(nrow(out_transformed)<length(outtimes))
        {
            tmp <- matrix(0,length(outtimes),ncol(out_transformed)-1)
            tmpnames <-names(out_transformed)
            out_transformed <- data.frame(time=outtimes,tmp)
            names(out_transformed) <- tmpnames
        }

        mC <- modCost(out_transformed, observed, y = "value",
                      err = 'err', weight = weight, scaleVar = scaleVar)

        ## Report and/or plot if the model is improved
        if (cost.old-mC$model >  ctr$quiet.tol) {
            if(!quiet) cat("Model cost at call ", calls, ": ", mC$model, "\n")

                                        # Plot the data and current model output if requested
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

            assign("cost.old", mC$model, inherits=TRUE)
        }
        return(mC)
    }
    if(plot) x11()

   #########     ##############################
        if(fitstart==TRUE)
        {
            ## ############# Get Initial Paramtervalues   #############
            ## Start with no weighting
            pnames=names(c(state.ini.optim, parms.optim))
            fn <- function(P){
                names(P) <- pnames
                FF<<-cost(P)
                return(FF$model)}
            a <- try(fit <- solnp(c(state.ini.optim, parms.optim),fun=fn,LB=lower,UB=upper,control=control),silent=TRUE)
            flag <- 1
            optimmethod <- 'solnp'
            if(class(a) == "try-error")
            {
                print('solnp fails, try PORT or other algorithm by users choice, might take longer time. Do something else!')
                warning('solnp fails, switch to  PORT or other algorithm by users choice')
                fit <- modFit1(cost, c(state.ini.optim, parms.optim), lower = lower, upper = upper, method=submethod,control=kingui.control(method=submethod,tolerance=ctr$control$tol)$control)
                optimmethod <- submethod
                flag <- 0
                ## now using submethod already

            }
            if(commonsigma==TRUE){
                ## no need for updating
                ## observed$err <- rep(1,nrow(observed))
                if(flag==1)## fit from solnp
                {
                    ## run a small loop with 'Marq' or some other method
                    fit <- modFit1(cost, fit$par, lower = lower, upper = upper, method=submethod,control=kingui.control(method=submethod,tolerance=ctr$control$tol)$control)
                    optimmethod <- submethod
                }
                cov0 <- summary(fit)$cov.scaled*2.4^2/length(fit$par)
                res <- modMCMC(cost,fit$par,...,jump=cov0,lower=lower,upper=upper,prior=prior,var0=var0,wvar0=wvar0,niter=niter,outputlength = outputlength, burninlength = burninlength, updatecov = updatecov,ntrydr=ntrydr,drscale=drscale,verbose=verbose)
            }else{
                ## ############## One reweighted estimation ############################
                ## Estimate the error variance(sd)
                olderr <- rep(1,length(mod_vars))
                ## # other list need to be attached to fit to give comparable results as in modFit.
                if(flag==1)## fit from solnp
                {
                    fit$ssr <- fit$values[length(fit$values)]
                    fit$residuals <-FF$residual$res
                    ## mean square per varaible
                    if (class(FF) == "modCost") {
                        names(fit$residuals)  <- FF$residuals$name
                        fit$var_ms            <- FF$var$SSR/FF$var$N
                        fit$var_ms_unscaled   <- FF$var$SSR.unscaled/FF$var$N
                        fit$var_ms_unweighted <- FF$var$SSR.unweighted/FF$var$N

                        names(fit$var_ms_unweighted) <- names(fit$var_ms_unscaled) <-
                            names(fit$var_ms) <- FF$var$name
                    } else fit$var_ms <- fit$var_ms_unweighted <- fit$var_ms_unscaled <- NA
                }
                err1 <- sqrt(fit$var_ms_unweighted)
                ERR <- err1[as.character(observed$name)]
                observed$err <-ERR
                diffsigma <- sum((err1-olderr)^2)
                ## At least do one iteration step to get a weighted LS using submethod
                fit <- modFit1(cost, fit$par, lower = lower, upper = upper, method=submethod,control=kingui.control(method=submethod,tolerance=ctr$control$tol)$control)
                optimmethod <- submethod
                ## ##########################################################
                ## Use this as the Input for MCMC algorithm
                ## ##########################
                cov0 <- summary(fit)$cov.scaled*2.4^2/length(fit$par)
                var0 <- fit$var_ms_unweighted
                if(any(is.na(cov0))){
                  res <- modMCMC(cost,fit$par,...,
                                 jump=NULL,lower=lower,upper=upper,prior=prior,var0=var0,wvar0=wvar0,niter=niter,outputlength = outputlength, burninlength = burninlength, updatecov = updatecov,ntrydr=ntrydr,drscale=drscale,verbose=verbose)
                }else{
                  res <- modMCMC(cost,fit$par,...,
                               jump=cov0,lower=lower,upper=upper,prior=prior,var0=var0,wvar0=wvar0,niter=niter,outputlength = outputlength, burninlength = burninlength, updatecov = updatecov,ntrydr=ntrydr,drscale=drscale,verbose=verbose)
                }
            }
        }else{
            ## Use either the starting value provided by the user or the default values.
            if(commonsigma==TRUE){
                res <- modMCMC(cost,c(state.ini.optim, parms.optim),...,jump=NULL,lower=lower,upper=upper,prior=prior,var0=var0,wvar0=wvar0,niter=niter,outputlength = outputlength, burninlength = burninlength, updatecov = updatecov,ntrydr=ntrydr,drscale=drscale,verbose=verbose)
                optimmethod <- NULL
            }else{
                fit <- modFit(cost, c(state.ini.optim, parms.optim), lower = lower,
                              upper = upper,method=method,control=control,...)
                optimmethod <- method
                tmpres <- fit$residuals
                oldERR <- observed$err
                err <- rep(NA,length(mod_vars))
                for(i in 1:length(mod_vars))
                {
                    box <- mod_vars[i]
                    ind <- which(names(tmpres)==box)
                    tmp <- tmpres[ind]
                    err[i] <- sd(tmp)
                }
                names(err) <- mod_vars
                ERR <- err[as.character(observed$name)]
                observed$err <-ERR
                olderr <- rep(1,length(mod_vars))
                ## At least do one iteration step to get a weighted LS
                fit <- modFit(cost, fit$par, lower = lower, upper = upper,method=submethod,control=control, ...)
                optimmethod <- submethod
                var0 <- fit$var_ms_unweighted
                                        #print(summary(fit))
                res <- modMCMC(cost,c(state.ini.optim, parms.optim),...,jump=NULL,lower=lower,upper=upper,prior=prior,var0=var0,wvar0=wvar0,niter=niter,outputlength = outputlength, burninlength = burninlength, updatecov = updatecov,ntrydr=ntrydr,drscale=drscale,verbose=verbose)

            }
        }
    fit <- res#### original mcmc output
    fit$optimmethod <- optimmethod
    ## calculate additional information for the output summary:
    ## Collect initial parameter values in two dataframes
    if(outpartri=='default'){
        if(length(state.ini.optim)>0){
            fit$start0 <- data.frame(initial=state.ini.optim,type=rep("state", length(state.ini.optim)),lower=lower[1:length(state.ini.optim)],upper=upper[1:length(state.ini.optim)])}else{
                fit$start0 <- data.frame(initial=state.ini.optim,type=rep("state", length(state.ini.optim)),lower=numeric(0),upper=numeric(0))
            }
        start0 <- mkinmodini$start[mkinmodini$start$fixed==0,]
        fit$start0 <- rbind(fit$start0,data.frame(initial=start0$initial,type=start0$type,lower=start0$lower,upper=start0$upper,row.names =rownames(start0)))
        fit$fixed0 <- data.frame(value = state.ini.fixed,type=rep("state", length(state.ini.fixed)),by=rep("user", length(state.ini.fixed)))
        fixed0 <- mkinmodini$start[mkinmodini$start$fixed==1,]
        if(nrow(fixed0)>0) fit$fixed0 <- rbind(fit$fixed0,data.frame(value=fixed0$initial,type=fixed0$type,by=rep('user',nrow(fixed0)),row.names=rownames(fixed0)))
    }else{
        fit$start0 <- NULL
        fit$fixed0 <- NULL
    }
        fit$start <- data.frame(initial = c(state.ini.optim, parms.optim))
        fit$start$type = c(rep("state", length(state.ini.optim)), rep("deparm", length(parms.optim)))
        fit$start$lower <- lower
        fit$start$upper <- upper

        fit$fixed <- data.frame(
                                value = c(state.ini.fixed, parms.fixed))
        fit$fixed$type = c(rep("state", length(state.ini.fixed)), rep("deparm", length(parms.fixed)))
        fit$fixed$by <- c(rep("user", length(state.ini.fixed)), mkinmodini$fixed_flag)

    ## ##############
    ## # Re-arrange observed and predicted data.
    err1 <- sqrt(apply(res$sig,2,mean))
    observed$err <- err1[as.character(paste('var_',observed$name,sep=''))]
    ## #
    fit$observed <- observed
    fit$diffs <- mkinmodini$diffs
    fit$par <- apply(res$pars,2,mean)
    ## ####################################
    ccend <- cost(fit$par)### using the estimated parameter to calculate the values.
    fit$residuals <-ccend$residual$res
    predicted_long <- mkin_wide_to_long(out_predicted, time = "time")
    fit$predicted <- out_predicted
    ## ###############################

    ## Calculate chi2 error levels according to FOCUS (2006)
    ## 0 values at sample time 0 should not be used.
    observed1 <- observed
    observed1 <- observed[!(observed$time==0 & observed$value==0),]
    means <- aggregate(value ~ time + name, data = observed1, mean, na.rm=TRUE)##using the mean of repeated measurements.
                                        #browser()
                                        #errdata <- merge(means, predicted_long, observed,by = c("time", "name"), suffixes = c("_mean", "_pred",'_obs'))
    errdata <- merge(means, predicted_long, by = c("time", "name"), suffixes = c("_mean", "_pred"))

    errdata <- errdata[order(errdata$time, errdata$name), ]
    errobserved <- merge(observed, predicted_long, by = c("time", "name"), suffixes = c("_obs", "_pred"))
    errmin.overall <- chi2err(errdata, length(parms.optim) + length(state.ini.optim),errobserved)#    errmin.overall <- chi2err(errdata, length(parms.optim) + length(state.ini.optim))
    errmin <- data.frame(err.min = errmin.overall$err.min,
                         n.optim = errmin.overall$n.optim, df = errmin.overall$df,err.sig = errmin.overall$err.sig,RMSE=errmin.overall$RMSE,EF=errmin.overall$EF,R2=errmin.overall$R2)
    rownames(errmin) <- "All data"
    for (obs_var in obs_vars)
    {
        errdata.var <- subset(errdata, name == obs_var)
        errobserved.var <- subset(errobserved, name == obs_var)
        if(outpartri=='default'){
            ##n.k.optim <- (paste("k", obs_var, sep="_")) %in% (names(parms.optim))+length(grep(paste("f", obs_var,'to', sep="_"), names(parms.optim)))#length(grep(paste("k", obs_var, sep="_"), names(parms.optim)))+length(grep(paste("f", obs_var, sep="_"), names(parms.optim)))
            n.k.optim <- (paste("k", obs_var, sep="_")) %in% (names(parms.optim))+length(grep(paste("f",'.*','to',obs_var,sep="_"), names(parms.optim)))
        }
        if(outpartri=='water-sediment'){
            n.k.optim <- length(grep(paste("k_", obs_var, '_',sep=""), names(parms.optim)))
        }
        n.initials.optim <- as.numeric((paste('M0_',obs_var, sep="")) %in% (names(state.ini.optim)))#length(grep(paste('M0_',obs_var, sep=""), names(state.ini.optim)))
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

                                        #errmin.tmp <- mkinerrmin(errdata.var, n.optim)
        ##errmin.tmp <- chi2err(errdata.var, n.optim)
        errmin.tmp <- chi2err(errdata.var, n.optim,errobserved.var)
        errmin[obs_var, c("err.min", "n.optim", "df",'err.sig','RMSE','EF','R2')] <- errmin.tmp
    }
    fit$errmin <- errmin

    ## Calculate dissipation times DT50 and DT90 and formation fractions
    parms.all = c(fit$par, parms.fixed)
    fit$distimes <- data.frame(DT50 = rep(NA, length(obs_vars)), DT90 = rep(NA, length(obs_vars)),Kinetic=rep(NA,length(obs_vars)),row.names = obs_vars)
    if(mkinmodini$outpartri=='default'){
        ## Now deals with formation fractions in case the outpartri is 'default'##
        ## Keep the original mcmckinfit.gui codes ##
        fit$ff <- vector()
        ff_names = names(mkinmodini$ff)
        for (ff_name in ff_names)
        {
            fit$ff[[ff_name]] =
                eval(parse(text = mkinmodini$ff[ff_name]), as.list(parms.all))
        }
        ## ###
        ## browser()
        for (obs_var in obs_vars) {
            f_tot <- grep(paste(obs_var, "_",sep=''), names(fit$ff), value=TRUE)
            f_exp <- grep(paste(obs_var, "to",obs_var,sep='_'), names(fit$ff), value=TRUE)
            f_exp1 <- grep(paste(obs_var, "to",'sink',sep='_'), names(fit$ff), value=TRUE)
            fit$ff[[paste(obs_var,'to', "sink", sep="_")]] = 1 - sum(fit$ff[f_tot])+sum(fit$ff[f_exp])+sum(fit$ff[f_exp1])
                                        #fit$ff[[paste(obs_var,'to', "sink", sep="_")]] = 1 - sum(fit$ff[f_tot])+sum(fit$ff[f_exp])
            type = names(mkinmodini$map[[obs_var]])[1]
            k1name <- paste("k1", obs_var,  sep="_")
            k2name <- paste("k2", obs_var,  sep="_")
            gname <- paste("g", obs_var,  sep="_")
            tbname <- paste("tb", obs_var,  sep="_")
            alphaname <- paste("alpha", obs_var,  sep="_")
            betaname <- paste("beta", obs_var,  sep="_")
            if (type == "SFO") {
                                        #k_names = grep(paste("k", obs_var, sep="_"), names(parms.all), value=TRUE)
                                        #k_tot = sum(parms.all[k_names])
                k_name <- paste("k", obs_var,sep="_")# k_name <- grep(paste("k", obs_var,sep="_"), names(parms.all), value=TRUE)
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
                DTmax <- 1000
                DT50.o <- optimize(f, c(0.0001,DTmax), x=50)$minimum
                DTmax1 <- log(2)/min(k1,k2)
                if(DTmax1==Inf) DTmax1 <- .Machine$double.xmax
                DT50.o1 <- optimize(f, c(0, DTmax1), x=50)$minimum
                DT50.o <- ifelse(f(DT50.o,50)>f(DT50.o1,50), DT50.o1,DT50.o)
                DT50 = ifelse(DTmax - DT50.o < 0.1, NA, DT50.o)
                DT90.o <- optimize(f, c(0.001, DTmax), x=90)$minimum
                DTmax1 <- log(10)/min(k1,k2)
                if(DTmax1==Inf) DTmax1 <- .Machine$double.xmax
                DT90.o1 <- optimize(f, c(0, DTmax1), x=90)$minimum
                DT90.o <- ifelse(f(DT90.o,90)>f(DT90.o1,90), DT90.o1,DT90.o)
                DT90 = ifelse(DTmax - DT90.o < 0.1, NA, DT90.o)
            }
            if (type == "HS") {
                ## k1 = parms.all["k1"]
                ## k2 = parms.all["k2"]
                ## tb = parms.all["tb"]
                k1 = parms.all[k1name]
                k2 = parms.all[k2name]
                tb = parms.all[tbname]
                f <- function(t, x) {
                    fraction = ifelse(t <= tb, exp(-k1 * t), exp(-k1 * tb) * exp(-k2 * (t - tb)))
                    (fraction - (1 - x/100))^2
                }
                ##DT50=1
                ##DT90=2
                DTmax <- 1000
                ##DT50 <- nlm(f, 0.0001, x=50)$estimate
                ##DT90 <- nlm(f, 0.0001, x=90)$estimate
                hso1 <- nlminb(0.0001,f, x=50)
                hso2 <- nlminb(tb,f, x=50)
                DT50.o <- ifelse(hso1$objective<=hso2$objective,hso1$par,hso2$par)
                DT50 = ifelse(DTmax - DT50.o < 0.1, NA, DT50.o)

                hso1 <- nlminb(0.0001,f, x=90)
                hso2 <- nlminb(tb,f, x=90)
                DT90.o <- ifelse(hso1$objective<=hso2$objective,hso1$par,hso2$par)
                DT90 = ifelse(DTmax - DT90.o < 0.1, NA, DT90.o)

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
            fit$distimes[obs_var, ] = c(ifelse(is.na(DT50),NA,formatC(DT50,4,format='f')), ifelse(is.na(DT90),NA,formatC(DT90,4,format='f')),type)#c(DT50, DT90,type)
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
                k_tot = sum(parms.all[k_names])
                DT50 = log(2)/k_tot
                DT90 = log(10)/k_tot
                for (k_name in k_names){
                    fit$ff[[sub("k_", "", k_name)]] = parms.all[[k_name]] / k_tot
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
                DTmax <- 1000
                DT50.o <- optimize(f, c(0.0001,DTmax), x=50)$minimum
                DTmax1 <- log(2)/min(k1,k2)
                DT50.o1 <- optimize(f, c(0, DTmax1), x=50)$minimum
                DT50.o <- ifelse(f(DT50.o,50)>f(DT50.o1,50), DT50.o1,DT50.o)
                DT50 = ifelse(DTmax - DT50.o < 0.1, NA, DT50.o)
                DT90.o <- optimize(f, c(0.001, DTmax), x=90)$minimum
                DTmax1 <- log(10)/min(k1,k2)
                DT90.o1 <- optimize(f, c(0, DTmax1), x=90)$minimum
                DT90.o <- ifelse(f(DT90.o,90)>f(DT90.o1,90), DT90.o1,DT90.o)
                DT90 = ifelse(DTmax - DT90.o < 0.1, NA, DT90.o)
            }
            if (type == "HS") {
                k1 = parms.all[k1name]
                k2 = parms.all[k2name]
                tb = parms.all[tbname]
                f <- function(t, x) {
                    fraction = ifelse(t <= tb, exp(-k1 * t), exp(-k1 * tb) * exp(-k2 * (t - tb)))
                    (fraction - (1 - x/100))^2
                }
                DTmax <- 1000
                ##DT50 <- nlm(f, 0.0001, x=50)$estimate
                ##DT90 <- nlm(f, 0.0001, x=90)$estimate
                hso1 <- nlminb(0.0001,f, x=50)
                hso2 <- nlminb(tb,f, x=50)
                DT50.o <- ifelse(hso1$objective<=hso2$objective,hso1$par,hso2$par)
                DT50 = ifelse(DTmax - DT50.o < 0.1, NA, DT50.o)

                hso1 <- nlminb(0.0001,f, x=90)
                hso2 <- nlminb(tb,f, x=90)
                DT90.o <- ifelse(hso1$objective<=hso2$objective,hso1$par,hso2$par)
                DT90 = ifelse(DTmax - DT90.o < 0.1, NA, DT90.o)


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
            fit$distimes[obs_var, ] = c(ifelse(is.na(DT50),NA,formatC(DT50,4,format='f')), ifelse(is.na(DT90),NA,formatC(DT90,4,format='f')),type)#c(DT50, DT90,type)
        }

    }
                                        #browser()
    ## Collect observed, predicted and residuals
    observed0 <-  mkin_wide_to_long(mkinmodini$data0,time='time')
    observed0$err <- observed$err
    data <- merge(observed, predicted_long, by = c("time", "name"))
    data0 <- merge(observed0, predicted_long, by = c("time", "name"))
    ##names(data) <- c("time", "variable", "observed", "predicted")
    names(data) <- c("time", "variable", "observed","err-std", "predicted")
    names(data0) <- c("time", "variable", "observed","err-std", "predicted")
    data$residual <- data$observed - data$predicted
    data$variable <- ordered(data$variable, levels = obs_vars)
    data0$residual <- data$residual
    data0$variable <- data$variable
    tmpid <- is.na(data0$residual) & !is.na(data0$observed)
    data0$'err-std'[tmpid] <- 0
    fit$data <- data[order(data$variable, data$time), ]
    fit$data0 <- data0[order(data0$variable, data0$time), ]
    fit$atol <- atol
    fit$solution <- solution
    if (solution == "eigen") {
        fit$coefmat <- mkinmodini$coefmat
    }
    if (solution == "deSolve") {
        fit$mkindiff <- mkindiff
    }
    fit$inpartri <- inpartri
    fit$outpartri <- outpartri
    fit$map <- mkinmodini$map
    ## #################################################


    class(fit) <- c('mcmckingui',"modMCMC")
    return(fit)
}

## -----------------------------------------------------------------------------
##' Lists model equations, the chi2 error levels
##' calculated according to FOCUS guidance (2006)
##' and optionally the data, consisting of observed, predicted and residual values, the
##' correlation matrix and also the same summary statistics for class \code{modMCMC}
##'
##' @title S3 summary method for class \code{mcmckingui}
##' @method summary mcmckingui
##' @param object A fitted object of class \code{mcmckingui} from the result of
##' \code{\link{mcmckinfit.full}}.
##' @param remove The iterations should be removed from summary statistics calculations.
##' @param data If TRUE, include in the returned values a data frame containing the observed
##' and predicted values with residuals and estimated standard deviations or weights.
##' @param distimes If TRUE, DT50 and DT90 values should be included.
##' @param ff If TRUE, the formation fraction should be calculated from the estimated
##' transformed parameters.
##' @param version  A version number indicating which version of the fit function
##' has been used.
##' @param ... Optional arguments passed to methods like \code{print}.
##' @return The summary function returns a list of components from the results of the optimizat
##' ion routine used.
##' @seealso \code{\link{summary.kingui}}.
##' @author Zhenglei Gao
##' @S3method summary mcmckingui
##' @rdname summary.mcmckingui
summary.mcmckingui <- function (object, remove = NULL, data=TRUE,distimes=TRUE,ff=TRUE,version="1.2011.922.1530",...) {
  options(warn=-1)
  mcmc <- object$pars
  ### drop more burnin period.
  if (! is.null (remove)) {
    if (max(remove) > nrow(mcmc))
      stop("too many runs should be removed from modMCMC object")
    if (min(remove) < 1)
      stop("cannot remove negative runs from modMCMC object")
    mcmc <- mcmc[-remove,]
  }
  ## Res <- data.frame(rbind(
  ##   mean = apply(mcmc, 2, mean),
  ##   sd   = apply(mcmc, 2, sd),
  ##   min  = apply(mcmc, 2, min),
  ##   max  = apply(mcmc, 2, max),
  ##   q025 = apply(mcmc, 2, quantile, probs = 0.25),
  ##   q050 = apply(mcmc, 2, quantile, probs = 0.5),
  ##   q075 = apply(mcmc, 2, quantile, probs = 0.75)
  ## ))

  ## if (!is.null(object$sig))
  ##   Res <- data.frame(Res,
  ##     sig = rbind(
  ##       mean = apply(object$sig, 2, mean),
  ##       sd   = apply(object$sig, 2, sd),
  ##       min  = apply(object$sig, 2, min),
  ##       max  = apply(object$sig, 2, max),
  ##       q025 = apply(object$sig, 2, quantile, probs = 0.25),
  ##       q050 = apply(object$sig, 2, quantile, probs = 0.5),
  ##       q075 = apply(object$sig, 2, quantile, probs = 0.75)
  ##   ))
  ###########################
    rdf <- (length(object$residuals)-length(object$par))

    resvar <- sum((object$residuals)^2) / (length(object$residuals)-length(object$par))
    param  <- object$par
    pnames <- names(param)
    p      <- length(param)
                                        #rownames(covar) <- colnames(covar) <-pnames

    se     <- apply(mcmc,2,sd)
    ##browser()
    lci <-apply(mcmc,2,quantile,probs=0.025)
    uci <-apply(mcmc,2,quantile,probs=0.975)
    names(se) <- pnames
    param <- cbind(param, se, lci,uci)
     ## adding one line when there is fixed by KinGUII object
    pnames1 <- pnames

    fid <- which(object$fixed$by=='KinGUII')
    if(length(fid)>0){
        nfid <- length(fid)
        param <- rbind(param,(cbind(object$fixed$value[fid],matrix(NA,nfid,3))))
        pnames1 <- c(pnames1, rownames(object$fixed)[fid])
    }
    dimnames(param) <- list(pnames1, c("Estimate", "Std. Error",'Lower CI','Upper CI'))
                                        #,"Pr(>t)"))

    #dimnames(param) <- list(pnames, c("Estimate", "Std. Error",'Lower CI','Upper CI'))
    ## tval      <- param / se
    ## modVariance <- object$ssr / length(object$residuals)

    ## if (!all(object$start$lower >=0)) {
    ##   message <- "Note that the one-sided t-test may not be appropriate if
    ##     parameter values below zero are possible."
    ##   warning(message)
    ## } else message <- "ok"

    ## param <- cbind(param, se, lci,uci, pt(tval, rdf, lower.tail = FALSE))
    ## dimnames(param) <- list(pnames, c("Estimate", "Std. Error",'Lower CI','Upper CI',
    ##                                    "Pr(>t)"))

    ans <- list(residuals = object$residuals,
                                        #residualVariance = resvar,
                sigma = sqrt(resvar),
                                        #modVariance = modVariance,
                df = c(p, rdf),
                Corr=cor(mcmc),
                info = object$info, niter = object$iterations,
                stopmess = message,
                par = param,version=version,pnames=pnames)

    ans$diffs <- object$diffs
    if(data)  {
        if(!is.null(object$data0)) ans$data <- object$data0 else ans$data <- object$data
    }
    ans$start <- object$start
    ans$fixed <- object$fixed
    ans$start0 <- object$start0
    ans$fixed0 <- object$fixed0
    ans$errmin <- object$errmin
    ans$optimmethod <- object$optimmethod
    if(distimes) ans$distimes <- object$distimes
    if(ff) {## actually always report ff.
      ans$ff <- object$ff
      fixed <- object$fixed[,1]
      names(fixed) <- rownames(object$fixed)
      if(length(object$trff)>0){
        tff <- MCMCstd(mu=object$pars,transformations=object$trff,fixed=fixed)
        ext <- setdiff(names(ans$ff),rownames(tff))
        ans$ff <- rbind(formatC(tff,digits=4,format='f'),cbind(formatC(ans$ff[ext],digits=4),matrix("",nrow=length(ext),ncol=3)))
      }else{
        ans$ff <- cbind(formatC(ans$ff,digits=4),matrix("",nrow=length(ans$ff),ncol=3))
        colnames(ans$ff) <- c("Estimate", "Std. Error",'Lower CI','Upper CI')
      }

    }
    ans$outpartri <- object$outpartri
    class(ans) <- c('summary.mcmckingui')
    options(warn=0)
    return(ans)
    ## #########################

}
##' S3 print method for class \code{summary.mcmckingui}
##'
##' @title S3 print method for class \code{summary.mcmckingui}
##' @method print summary.mcmckingui
##' @param x An object of class \code{summary.mcmckingui}.
##' @param digits Number of digits to be printed after the decimal point.
##' @param detailed Not used.
##' @param ... Optional arguments.
##' @return \code{NULL}
##' @author Zhenglei Gao
##' @S3method print summary.mcmckingui
##' @rdname print.summary.mcmckingui
print.summary.mcmckingui <- function(x,digits = max(3, getOption("digits") - 3),detailed=FALSE, ...) {
   cat(paste("KineticEval Package Version:",packageVersion("KineticEval"),"\n"))
    cat(paste('Version:',x$version,'\n'))
    ##cat('\nR version: 2.12.2 (2011-02-25)\n ')
    cat(paste('\n',sessionInfo()$R.version$version.string,'\n',sep=""))
    cat("\nMethods:MCMC\n")
    if(x$outpartri=='water-sediment')  cat("\nStudy:Water-Sediment\n")
    xx <- x[["diffs"]]
    cat("\nOptimization Algorithms Used:\n")
    print(x$optimmethod)
    cat("\nEquations:\n")
    for(i in 1:length(xx)) print(noquote(as.character(xx[i])))

    cat("\nStarting values for optimised parameters:\n")
    if(x$outpartri=='default') print(x$start0) else print(x$start)

    cat("\nFixed parameter values:\n")
    if(x$outpartri=='default') {
        if(length(x$fixed0$value) == 0) cat("None\n")
        else print(x$fixed0)
    }else{
        if(length(x$fixed$value) == 0) cat("None\n")
        else print(x$fixed)
    }


    cat("\nOptimised parameters:\n")
    printCoefmat(x$par, digits = digits, ...)

    cat("\nResidual standard error:",
        format(signif(x$sigma, digits)), "on", x$df[2], "degrees of freedom\n")

    printff <- !is.null(x$ff)
    if(printff){
        cat("\nEstimated formation fractions:\n")
        print(data.frame(ff = x$ff), digits=digits,...)
    }




    cat("\nChi2 error levels in percent :\n")
    print(x$errmin[,1:3], digits=digits,...)

    printdistimes <- !is.null(x$distimes)
    if(printdistimes){
        cat("\nEstimated disappearance times:\n")
        print(x$distimes, digits=digits,...)
    }





    printcor <- (!is.null(x$Corr))
    if (printcor){
        cat("\nAdditional Statistics:\n")
        print(x$errmin[,4:ncol(x$errmin)], digits=digits,...)
        Corr <- x$Corr
        rownames(Corr) <- colnames(Corr) <- x$pnames## rownames(x$par)
        cat("\nParameter correlation:\n")
        print(Corr, digits = digits, ...)
    }


    printdata <- !is.null(x$data)
    if (printdata){
        cat("\nData:\n")
        print(myformat(x$data, digits = digits, scientific = FALSE,...), row.names = FALSE)
    }

    invisible(x)
}
########################################
##' S3 plot method to plot for calss 'mcmckingui'
##'
##' Make the density correlation and trace plots of the generated Markov Chains.
##'

##' @param x An object of class
##' 'mcmckingui'
##' @param y NULL
##' @param fname1 The file name of the
##' density plot.
##' @param fname2 The file name of the
##' correlation plot.
##' @param fname3 The file name of the trace plot
##' @param pch What kind of points to use in
##' the plots.
##' @param device The plot device to be
##' used.
##' @param \dots Other arguments to be passed
##' to 'plot'.
##' @return Density and Correlation plots of the sampled parameters in 'wmf' or
##' other format.
##' @note \code{plot.mcmckingui0} is a deprecated version to be used with the object
##' returned by \code{mcmckinfit.gui}
##' @author Zhenglei Gao
##' @method plot mcmckingui
##' @S3method plot mcmckingui
##' @rdname plot.mcmckingui
plot.mcmckingui <- function(x,y=NULL,fname1='density',fname2='correlation',fname3='trace',pch=1,device='wmf',...)
{
    ## Make the density plot.
    ## Make the correlation plot
    ## Make the trace plot
	object <- x
	if(!is.null(y)) stop("do not try that!")
    if(device=='wmf') win.metafile(filename=paste(fname1,'wmf',sep='.'),width = 7,height = 7, pointsize = 12)
    if(device=='jpeg') jpeg(filename=paste(fname1,'jpeg',sep='.'),width = 700,height = 700, pointsize = 12)
    pnames <- colnames(object$pars)
    a <- autolayout(ncol(object$pars))
    par(mfrow=c(a$m,a$n))
    for(i in 1:ncol(object$pars))
        {plot(density(object$pars[,i]),main='')
    title(paste('Density of',pnames[i]))}
    dev.off()
    if(device=='wmf') win.metafile(filename=paste(fname2,'wmf',sep='.'),width = 7,height = 7, pointsize = 12)
    if(device=='jpeg') jpeg(filename=paste(fname2,'jpeg',sep='.'),width = 700,height = 700, pointsize = 12)
    pairs(object$pars,pch=pch)
    dev.off()

    if(device=='wmf') win.metafile(filename=paste(fname3,'wmf',sep='.'),width = 7,height = 7, pointsize = 12)
    if(device=='jpeg') jpeg(filename=paste(fname3,'jpeg',sep='.'),width = 700,height = 700, pointsize = 12)
    par(mfrow=c(a$m,a$n))
    traceplot(as.mcmc(object$pars))
    dev.off()
}

##' Utility function to set the layout when making multiple plots on the same page.
##'
##' @title Set the layout when making multiple plots on the same page.
##' @param n number of plots to make
##' @return number of rows and columns in the layout
##' @author Zhenglei Gao
##' @export
autolayout <- function(n)
{
    m <- ceiling(sqrt(n))
    n <- ceiling(n/m)
    return(list(m=m,n=n))
}
