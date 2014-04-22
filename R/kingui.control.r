#require(optimx)
##' Choose the algorithm to use and related control parameters in kinetic
##' evaluations using NLS and IRLS.
##'
##' @param method The method to be used, one
##' of "Rsolnp", "minqa","spg", "Marq", "Port", "Newton", "Nelder-Mead",
##' "BFGS", "CG", "L-BFGS-B", "SANN", "Pseudo" - see details.
##' @param maxIter Numerber of maximum
##' iterations
##' @param tolerance A positive numeric
##' value specifying the tolerance level for the relative offset convergence
##' criterion.
##' @param odesolver The integration
##' routines used. see \code{\link{ode}}
##' @param atol absolute error tolerance, either a scalar or an array as long
##' as 'y'. See \code{\link{ode}}
##' @param rtol Relative error tolerance,
##' either a scalar or an array as long as 'y'.See \code{\link{ode}}.
##' @param rhobeg Initial trust region
##' radius for method 'bobyqa'. See details of see \code{\link{bobyqa}}.
##' @param iprint Control the amount of
##' printing by setting iprint to 0, 1 2, or 3.  See details of see
##' \code{\link{bobyqa}}.
##' @param trace A logical variable
##' (TRUE/FALSE). If 'TRUE', information on the progress of optimization is
##' printed. See details of see \code{\link{spg}}.
##' @param goMarq If TRUE, using "Marq" for
##' the iterations after the first step in IRLS.
##' @param delta Control parameters for
##' 'solnp'. See details of see \code{\link{solnp}}.
##' @param rho Control parameters for 'solnp'.
##' See details of see \code{\link{solnp}}.
##' @param submethod If the method chosen
##' failed to produce results, run the optimization using a substitute method.
##' @param Hmethod1 The method used to calculate Hessian matrix
##' @param Hmethod2 The substitute method if \code{Hmethod1} fails.
##' @param quiet.tol The improvement level at which the objective function(cost
##' value) and the plot should be updated and reported.
##' @param runTRR4 if the number of compartments exceeds 4 in the system, do not
##' run TRR when reach boundary
##' @param ... Other characteristics of
##' different optimizer for the users to play with.
##' @return  A list of control parameters for the ODE solver and the
##' optimization routines.
##' @note Should have a control function for \code{\link{mcmckinfit.full}}
##' @author Zhenglei Gao
##' @keywords Kinetic-Evaluations
##' @examples
##'
##' kingui.control()
##'
##' @export
kingui.control <- function (method='LM',maxIter=100,tolerance=1e-8,odesolver='lsoda',atol=1e-9,
                            rtol=1e-10,rhobeg=0.05,iprint=1,trace=0,goMarq=0,delta=1e-6,rho=1,
                            submethod='Port',Hmethod1='LM',Hmethod2='L-BFGS-B',quiet.tol=1,
                            automate=TRUE,calcJacob=TRUE,runTRR4=FALSE,...)
{
    if(automate==TRUE){
      runTRR <- TRUE
      autofix <- TRUE
    }else{
      runTRR <- FALSE
      autofix <- FALSE
    }
    marqctr <- list(maxiter=maxIter,ftol=tolerance,ptol=tolerance,...)
    if (method %in% c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"))##using optim
    {
        if(method=='L-BFGS-B')
            {
                factr <- tolerance*1e+15
                return(list(method=method,control=list(maxit=maxIter,factr=factr,trace=trace,...),odesolver=odesolver,atol=atol,rtol=rtol,marqctr= marqctr,goMarq=goMarq,submethod=submethod,Hmethod1='L-BFGS-B',Hmethod2='Marq',quiet.tol=quiet.tol,runTRR=runTRR,autofix=autofix,calcJacob=calcJacob,runTRR4=runTRR4))
            }else{
                return(list(method=method,control=list(maxit=maxIter,reltol=tolerance,trace=trace,...),odesolver=odesolver,atol=atol,rtol=rtol,marqctr= marqctr,goMarq=goMarq,submethod=submethod,Hmethod1='L-BFGS-B',Hmethod2='Marq',quiet.tol=quiet.tol,runTRR=runTRR,autofix=autofix,calcJacob=calcJacob,runTRR4=runTRR4 ))
            }
    }else if (method == "Port") {
        return(list(method=method,control=list(iter.max=(maxIter+50),rel.tol=tolerance*1e-1,...),odesolver=odesolver,atol=atol,rtol=rtol,marqctr= marqctr,goMarq=goMarq,submethod=submethod,Hmethod1='L-BFGS-B',Hmethod2='Marq',quiet.tol=quiet.tol,runTRR=runTRR,autofix=autofix,calcJacob=calcJacob,runTRR4=runTRR4))
    }else if (method == "Marq") {
        return(list(method=method,control=list(maxiter=maxIter,ftol=tolerance,ptol=tolerance,...),odesolver=odesolver,atol=atol,rtol=rtol,marqctr= marqctr,goMarq=goMarq,submethod=submethod,Hmethod1='L-BFGS-B',Hmethod2='Marq',quiet.tol=quiet.tol,runTRR=runTRR,autofix=autofix,calcJacob=calcJacob,runTRR4=runTRR4 ))
    } else if (method == "LM") {
      return(list(method=method,control=list(maxiter=maxIter,ftol=tolerance,ptol=tolerance,...),odesolver=odesolver,atol=atol,rtol=rtol,marqctr= marqctr,goMarq=goMarq,submethod=submethod,Hmethod1='L-BFGS-B',Hmethod2='Marq',quiet.tol=quiet.tol,runTRR=runTRR,autofix=autofix,calcJacob=calcJacob,runTRR4=runTRR4 ))
    }else if (method == "TRR") {
      return(list(method=method,control=list(maxiter=maxIter,ftol=tolerance,ptol=tolerance,...),odesolver=odesolver,atol=atol,rtol=rtol,marqctr= marqctr,goMarq=goMarq,submethod=submethod,Hmethod1='L-BFGS-B',Hmethod2='Marq',quiet.tol=quiet.tol,runTRR=runTRR,autofix=autofix,calcJacob=calcJacob,runTRR4=runTRR4 ))
    }else if (method == "Pseudo") {
        return(list(method=method,control=list(numiter=maxIter,varleft=tolerance,...),odesolver=odesolver,atol=atol,rtol=rtol,marqctr= marqctr,goMarq=goMarq ,submethod=submethod,Hmethod1='L-BFGS-B',Hmethod2='Marq',quiet.tol=quiet.tol,runTRR=runTRR,autofix=autofix,calcJacob=calcJacob,runTRR4=runTRR4))

    }else if(method=='trust'){
        return(list(method=method,control = list(rhobeg=rhobeg,rhoend=tolerance/rhobeg,iprint=iprint,...),odesolver=odesolver,atol=atol,rtol=rtol,marqctr= marqctr,goMarq=goMarq,submethod=submethod ,Hmethod1='L-BFGS-B',Hmethod2='Marq',quiet.tol=quiet.tol,runTRR=runTRR,autofix=autofix,calcJacob=calcJacob,runTRR4=runTRR4))
    }else if(method=='spg'){
        return(list(method=method,control=list(maxit=maxIter,ftol=tolerance,trace=trace,...),odesolver=odesolver,atol=atol,rtol=rtol,marqctr= marqctr,goMarq=goMarq,submethod=submethod,Hmethod1='L-BFGS-B',Hmethod2='Marq',quiet.tol=quiet.tol,runTRR=runTRR,autofix=autofix,calcJacob=calcJacob,runTRR4=runTRR4 ))
    }else if(method=='nmk'){
        return(list(method=method,control=list(),odesolver=odesolver,atol=atol,rtol=rtol,marqctr= marqctr,goMarq=goMarq,submethod=submethod,Hmethod1='L-BFGS-B',Hmethod2='Marq',quiet.tol=quiet.tol,runTRR=runTRR,autofix=autofix,calcJacob=calcJacob,runTRR4=runTRR4 ))
    }else if(method=='dfoptim'){
        return(list(method=method,control=list(),odesolver=odesolver,atol=atol,rtol=rtol,marqctr= marqctr,goMarq=goMarq,submethod=submethod,Hmethod1='L-BFGS-B',Hmethod2='Marq',quiet.tol=quiet.tol,runTRR=runTRR,autofix=autofix,calcJacob=calcJacob,runTRR4=runTRR4))
    }else if(method=='DEoptim'){
        return(list(method=method,control=list(),odesolver=odesolver,atol=atol,rtol=rtol,marqctr= marqctr,goMarq=goMarq,submethod=submethod,Hmethod1='L-BFGS-B',Hmethod2='Marq',quiet.tol=quiet.tol,runTRR=runTRR,autofix=autofix,calcJacob=calcJacob,runTRR4=runTRR4 ))
    }else if(method=='Rvmmin'){
        return(list(method=method,control=list(),odesolver=odesolver,atol=atol,rtol=rtol,marqctr= marqctr,goMarq=goMarq,submethod=submethod,Hmethod1='L-BFGS-B',Hmethod2='Marq',quiet.tol=quiet.tol,runTRR=runTRR,autofix=autofix,calcJacob=calcJacob,runTRR4=runTRR4 ))
    }else if(method=='Rcgmin'){
        return(list(method=method,control=list(),odesolver=odesolver,atol=atol,rtol=rtol,marqctr= marqctr,goMarq=goMarq,submethod=submethod,Hmethod1='L-BFGS-B',Hmethod2='Marq',quiet.tol=quiet.tol,runTRR=runTRR,autofix=autofix,calcJacob=calcJacob,runTRR4=runTRR4 ))
        }else if(method=='ucminf'){
        return(list(method=method,control=list(),odesolver=odesolver,atol=atol,rtol=rtol,marqctr= marqctr,goMarq=goMarq,submethod=submethod,Hmethod1='L-BFGS-B',Hmethod2='Marq',quiet.tol=quiet.tol,runTRR=runTRR,autofix=autofix,calcJacob=calcJacob,runTRR4=runTRR4 ))
    }else if(method=='solnp'){
        return(list(method=method,control=list(delta=delta,tol=tolerance,rho=rho),odesolver=odesolver,atol=atol,rtol=rtol,marqctr= marqctr,goMarq=goMarq,submethod=submethod,Hmethod1='L-BFGS-B',Hmethod2='Marq',quiet.tol=quiet.tol,runTRR=runTRR,autofix=autofix,calcJacob=calcJacob,runTRR4=runTRR4))
    }
}



