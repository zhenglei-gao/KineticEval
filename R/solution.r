##' a modification of the FOMC model solution function.
##' 
##' Modified function based on \code{mkin::FOMC.solution}
##' @param t time points
##' @param parent.0 Concentration or state at time 0
##' @param alpha parameter for FOMC
##' @param beta parameter for FOMC
##' @author Zhenglei Gao
##' @export
##' @seealso \code{mkin::FOMC.solution}
FOMC.solution <- function (t, parent.0, alpha, beta) 
{
  if(beta<0){
    parent <- rep(-Inf,length(t))
  }else {
    if(beta==0){
      parent <- rep(0,length(t))
      parent[t==0] <- -Inf
    }else parent <-  parent.0/(t/beta + 1)^alpha
  }
  return(parent)
}