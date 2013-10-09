## a collection of functions for the trust region methods
## Inspired by Matlab lsqnonlin function

##STIR_LS  Functions
#{
  ## Steps 1: Consider the problem of finding the trial step size s_i between the ith point x_i and the next
  ## point x_i+1 that minimize \Psi_i()s


  ## Step 2: Compute the 2 dimentional subspace V with 2 spanning vectors: one vector in the direction of the
  ## gradient at x_i, the other in the approximate Gauss-Newto direction

  ## Step 3: Solve for the trial step s_i of the 2 dimentional subproblem

  ## Step 4: if for a predefined constant \tau in (0,1), u(x_i+s_i) < u(x_i), then x_i+1 becomes the current point,
  ## otherwise x_i+1=x_i

  ## Step 5: update \Delat_i

  ## Step 6, if \partial u(s_i) is below a chosen tolerance, end, o.w, repeat and increment i


##}

# exitflag: Integer identifying the reason the algorithm terminated. The following lists the values of exitflag
# and the corresponding reasons the algorithm terminated:
#
#   1 Function converged to a solution x.
#
#   2 Change in x was less than the specified tolerance.
#
#   3 Change in the residual was less than the specified tolerance.
#
#   4 Magnitude of search direction was smaller than the specified tolerance.
#
#   0 Number of iterations exceeded options.MaxIter or number of function evaluations exceeded options.MaxFunEvals.
#
#   -1 Output function terminated the algorithm.
#
#   -2 Problem is infeasible: the bounds lb and ub are inconsistent.
#
#   -4 Line search could not sufficiently decrease the residual along the current search direction.


##' the 2-norm is the sum of the squares of the elements
##'
##' see also the R core function \code{norm}
##' @title Compute the 2-norm of a vector
##' @param x a vector
##' @param na.rm logical, whether to remove NAs.
##' @return the 2-norm of the vector
##' @author Zhenglei Gao
norm <- function(x,na.rm=TRUE) sqrt(sum(x^2,na.rm=na.rm))

##' Define the control parameters for the STIR algorithm
##'
##' details check the Matlab function \code{lsqnonlin}
##' @title Define the control parameters for the STIR algorithm
##' @param nopar number of parameters
##' @param gradflag gradient function flag
##' @param pcflags pc flag
##' @param MaxPCGIter maximum number of iteration for PCG(conjugate gradient)
##' @param TolPCG Tolerance for PCG
##' @param MaxFunEvals maximum number of function evaluations
##' @param MaxIter maximum number of iterations
##' @param TolFun tolerance of function value difference
##' @param TolX tolerance of parameter value(norm) differece
##' @return a list of control parameters
##' @author Zhenglei Gao
trr.options <- function(nopar=1,gradflag=0,pcflags=Inf,MaxPCGIter=1,TolPCG=0.1,MaxFunEvals=100,MaxIter=500,TolFun=1e-6,TolX=1e-6)
{
  MaxPCGIter <- max(MaxPCGIter,nopar/2)
  MaxFunEvals <- max(MaxFunEvals,100*nopar)
  return(list(gradflag=gradflag,pcflags=pcflags,kmax=MaxPCGIter,
              pcgtol=TolPCG,maxfunevals=MaxFunEvals,itb=MaxIter,
              tol1=TolFun,tol2=TolX))
}

##' Secular equation
##'
##' Secular equation is 
##' @title Secular equation
##' @param lambda a vector
##' @param eigval eigenval 
##' @param alpha alpha
##' @param delta delta
##' @return the value of the secular equation at lambda
##' @author Zhenglei Gao
seceqn <- function(lambda,eigval,alpha,delta)
{
  m <- length(lambda)
  n <- length(eigval)
  unn <- rep(1,n)
  unm <- rep(1,m)
  M <- eigval%*%t(unm)+unn%*%t(lambda)
  MC <- M
  MM <- alpha%*%t(unm)
  M[M!=0] <- MM[M!=0]/M[M!=0]
  M[MC==0] <- Inf
  M <- M*M
  value <- sqrt(unm/(t(M)%*%unn))
  value[is.nan(value)] <- 0
  value <- (1/delta)*unm-value
  return(value)
}
##' Finding the solution for f=0 on the right side.
##'
##' Not working and not used.
##' @title Finding the solution for f=0 on the right side.
##' @param f function to solve for 0
##' @param ... other parameters passed into f
##' @param lower lower bound of the optimization
##' @param tol tolerance level
##' @param maxiter maximum number of iterations
##' @return a list 
##' \itemize{
##' \item{b}{solution}
##' \item{c}{?}
##' \item{itfun}{function evalution iterations}
##' }
##' @author Zhenglei Gao
rfzero <- function(f,...,lower=0,tol=1e-6,maxiter=100)
{
  stop("cannot be used still! Error inside!!")
  itfun <- 0
  x <- lower
  if(x!=0) dx <- abs(x)/2 else dx <- 0.5
  a <- x
  c <- a
  fa <- f(a,...)
  itfun <- itfun+1
  b <- x + dx
  ##b <- x+1
  fb <- f(b,...)
  itfun <- itfun +1
  while((fa>0) * (fb>0) ==1){
    dx <- 2*dx
    if((fa>0) != (fb>0)){
      break
    } else{
      b <- x + dx
      ##b <- x+1
      fb <- f(b,...)
      itfun <- itfun +1
      if(itfun>maxiter) break
    }


  }
  fc <- fb
  while(fb!=0){
    if((fb>0)* (fc>0) ==1){
      c <- a
      fc <- fa
      d <- b-a
      e <-d
    }
    if(abs(fc)<abs(fb)){
      a<-b
      b<-c
      c<-a
      fa <- fb
      fb <- fc
      fc <- fa
    }

    if(itfun>itbnd) break
    m <- 0.5*(c-b)
    toler <- 2.0*tol*max(abs(b),1)
    if(abs(m)<=toler || fb==0) break
    if(abs(e)<toler || (abs(fa)<=abs(fb))){
      d <- m
      e <- m
    }else{
      s <- fb/fa
      if(a==c){
        p <- 2.0*m*s
        q <- 1.0-s
      }else{
        q <- fa/fc
        r <- fb/fc
        p <- s*(2*m*q*(q-r)-(b-a)*(r-1))
        q <- (q - 1.0)*(r - 1.0)*(s - 1.0)
      }
      if(p>0) q <- -q else p <- -p
      if((2.0*p < 3.0*m*q - abs(toler*q)) && (p < abs(0.5*e*q))){
        e <- d
        d <- p/q
      }else{
        d <- m
        e <- m
      }
    }
    a <- b
    fa <- fb
    if(abs(d)>toler) b <- b+d else{
      if(b>c) b <- b-toler else b <- b+toler
    }
    fb <- f(b,...)
    itfun <- itfun+1

  }
  #browser()
  return(list(b=b,c=c,itfun=itfun))
}

##' Exact solution of very small dimentioned trust region problems
##'
##' It is meant only for small dimentional problems
##' @title Exact solution of trust region subproblem
##' @param g gradient
##' @param H Hessian
##' @param delta the trust region radias
##' @return a list of items describing the solution of the problem
##' \itemize{
##' \item{"s"}{solution}
##' \item{"val"}{The value of the quadratic at the solution}
##' \item{"def"}{logical, if H is postive definite, 1}
##' \item{"count"}{number of evaluations of the secular equation}
##' \item{"lambda"}{the corresponding Lagrange multiplier}
##' }
##' @author Zhenglei Gao
trust.exact <- function(g,H,delta)
{

  #   %TRUST  Exact soln of trust region problem
  #   %
  #   % [s,val,posdef,count,lambda] = TRUST(g,H,delta) Solves the trust region
  #   % problem: min{g^Ts + 1/2 s^THs: ||s|| <= delta}. The full
  #   % eigen-decomposition is used; based on the secular equation,
  #   % 1/delta - 1/(||s||) = 0. The solution is s, the value
  #   % of the quadratic at the solution is val; posdef = 1 if H
  #   % is pos. definite; otherwise posdef = 0. The number of
  #   % evaluations of the secular equation is count, lambda
  #   % is the value of the corresponding Lagrange multiplier.
  #   %
  #   % Copyright 1990-2010 The MathWorks, Inc.
  #   %   $Revision: 1.1.6.3 $  $Date: 2011/05/09 01:16:32 $
  #   % TRUST is meant to be applied to very small dimensional problems.
  norm <- function(x) sqrt(sum(x^2))
  g <- as.vector(g)
  H <- as.matrix(H)
  tol <- 1e-12
  tol2 <- 1e-8
  key <- 0
  itbnd <- 50
  lambda <- 0
  n <- length(g)
  coeff <- rep(0,n)
  tmp <- eigen(H,symmetric=TRUE)
  V <- tmp$vectors
  D <- tmp$values
  count <- 0
  eigval <- D
  jmin <- which.min(eigval)
  mineig <- eigval[jmin]
  alpha <- as.vector(-t(V)%*%g)
  sig <- sign(alpha[jmin])+(alpha[jmin]==0)

  ## positive definite case
  if(mineig > 0){
    coeff <- alpha / eigval;
    lambda <- 0;
    s <- V%*%coeff;
    posdef <- 1;
    nrms <- norm(s);
    if(nrms <= 1.2*delta) key <-1 else laminit <- 0
  }else{
    laminit <- -mineig
    posdef <- 0
  }

  if(key==0){
    tmp <- seceqn(laminit,eigval,alpha,delta)
    if(tmp>0){
      ##tmp1 <- rfzero(f=seceqn,eigval=eigval,alpha=alpha,delta=delta,lower=laminit,tol=tol,maxiter=itbnd)
      ##tmp1 <- optim(par=laminit,fn=seceqn^2,eigval=eigval,alpha=alpha,delta=delta,lower=laminit,control=list(tol=tol,maxiter=itbnd))
      tmp1 <- uniroot(f=seceqn,eigval=eigval,alpha=alpha,delta=delta,lower=laminit,upper=1e100,tol=tol,maxiter=itbnd)
      b <- tmp1$root
      count <- tmp1$iter
      vval <- abs(seceqn(b,eigval,alpha,delta))
      if(vval <= tol2){
        lambda <- b
        key <- 2
        lam <- lambda*rep(1,n)
        w <- eigval + lam
        arg1 <- (w==0) & (alpha == 0);
        arg2 <- (w==0) & (alpha != 0);
        coeff[w !=0] <- alpha[w !=0] / w[w!=0];
        coeff[arg1] = 0;
        coeff[arg2] = Inf;
        coeff <- as.vector(coeff)
        coeff[is.nan(coeff)]=0;
        s <- V%*%coeff;
        nrms <- norm(s);
        if( (nrms > 1.2*delta) || (nrms < .8*delta)){
          key <- 5
          lambda <- -mineig
        }
      }
    }else{
      lambda <- -mineig
      key <- 4
    }
    lam <- lambda*rep(1,n)
    if(key>2){
      arg <- abs(eigval+lam) < 10*eps*max(abs(eigval),1)
      alpha[arg]<-0
    }
    w <- eigval+lam
    arg1 <- w==0 & alpha==0
    arg2 <- w==0 & alpha!=0
    coeff[w!=0]<-alpha[w!=0]/w[w!=0]
    coeff[arg1] <- 0
    coeff[arg2] <- Inf
    coeff <- as.vector(coeff)
    coeff[is.nan(coeff)] <- 0
    s <- V%*%coeff
    nrms <- norm(s)
    if(key>2 && (nrms < 0.8*delta)){
      beta <- sqrt(delta^2-nrms^2)
      s <- s+ beta*sig*V[,jmin]
    }
    if((key>2) && (nrms>1.2*delta)){
      tmp1 <- uniroot(f=seceqn,eigval=eigval,alpha=alpha,delta=delta,lower=laminit,upper=1e100,tol=tol,maxiter=itbnd)
      b <- tmp1$root
      count <- tmp1$iter
      lambda <- b
      lam <- lambda*rep(1,n)
      w <- eigval+lam
      arg1 <- w==0 & alpha==0
      arg2 <- w==0 & alpha!=0
      coeff[w!=0] <- alpha[w!=0]/w[w!=0]
      coeff[arg1] <- 0
      coeff[arg2] <- Inf
      coeff <- as.vector(coeff)
      coeff[is.nan(coeff)] <- 0
      s <- V%*%coeff
      nrms <- norm(s)
    }
  }
  val <- t(g)%*%s+t(0.5*s)%*%(H%*%s)
  return(list(s=s,val=val,posdef=posdef,count=count,lambda=lambda))
}

##' Change the control parameters of the STIR algorithm
##'
##' Change the control parameters of the STIR algorithm
##' @title Change the control parameters of the STIR algorithm
##' @param npar number of parameters
##' @param user user defined control parameters list
##' @param ... control parameters that should be modified
##' @return a list of control parameters
##' @author Zhenglei Gao
trrcontrol <- function(npar=1,user=NULL,...)
{
  if(!missing(...)){
      opt <- list(...)
  }else opt <- vector()


  defaultopt <- list(findiffmethod="simple",MaxIter=100,
                     MaxFunEvals=400,
                     MaxPCGIter=1,
                     ScaleProblem="none",PrecondBandWidth=Inf,
                     TolFun=1e-6,TolPCG=0.1,TolX=1e-6,
                     TypicalX=rep(1,npar),
                     Jacobian="off")
  if(length(opt)>0){
    namesopt <- names(opt)

    for(nam in namesopt){
      defaultopt[[nam]] <- opt[[nam]]
    }
  }
  if(!is.null(user)){
    namesopt <- names(user)

    for(nam in namesopt){
      defaultopt[[nam]] <- user[[nam]]
    }
  }
  defaultopt$MaxPCGIter <- max(defaultopt$MaxPCGIter,floor(npar/2))
  defaultopt$MaxFunEvals <- max(defaultopt$MaxFunEvals,npar*100)
  return(defaultopt)
}

##' The R clone of subspace trust region reflective algorithm.
##'
##' This is a modification of the Matlab lsqnonlin function, with
##' only implementation of STIR, no LM
##' @title Optimization of lease square problems using STIR
##' @param fun function that returns a vector of residuals
##' @param xstart starting values for the parameters
##' @param l lower bound
##' @param u upper bound
##' @param options a list of user defined control parameters
##' @param ... other input arguments to fun
##' @return a list of result items including par, ssr, etc. 
##' @author Zhenglei
##' @export
lsqnonlin <- function(fun,xstart,l,u,options=NULL,...)
{
  npar <- length(xstart)
  if(missing(l)) l <- rep(-Inf,npar)
  if(missing(u)) u <- rep(Inf,npar)
  options <- trrcontrol(npar=npar,user=options)
  mtxmpy <- atamult

  res <- snls(fun,xstart,l,u,options=options,fval=fun(xstart,...),
              JACval=NULL,
              computeLambda=0,
              mtxmpy=atamult,...)
  res$ssr <- as.numeric(crossprod(res$fvec))
  res$par <- res$xcurr
  return(res)
}

##' solve box constrained nonlinear least squares problem.
##'
##' See sparse nonlinear ls solver
##' @title solve box constrained nonlinear least squares problem.
##' @param funfcn function returning residuals
##' @param xstart starting value
##' @param l lower bounds 
##' @param u upper bounds
##' @param verb logical, verbose
##' @param options control parameters for STIR optimization
##' @param fval current fvalue
##' @param JACval Jacobian value
##' @param computeLambda logical, whether to compute Lambda
##' @param mtxmpy Jacobian matrix multiply function
##' @param ... other input arguments for funfcn
##' @return a list
##' @author Zhenglei Gao
snls <- function(funfcn,xstart,l,u,verb=FALSE,options,
                 fval,JACval,
                 computeLambda,mtxmpy,...)
{
  n <- length(xstart)
  xcurr <- xstart
  msgData <- NULL
  dnewt <- vector()
  iter <- 0
  numFunEvals <- 1 ## done in calling function lsqnonlin
  numGradEvals <- 0
  tol2 <- options$TolX
  tol1 <- options$TolFun
  itb <- options$MaxIter
  maxfunevals <- options$MaxFunEvals
  pcgtol <- options$TolPCG
  kmax <- options$MaxPCGIter
  msgdata <- NULL
  ############
  ex <-0
  posdef <- 1
  npcg <- 0
  pcgit <-0

  delta <- 10
  nrmsx <- 1
  ratio <- 0

  degen <- Inf
  dv <- rep(1,n)
  x <- xstart
  oval <- Inf
  nbnds <- 1
  Z <- vector()
  if(all(u==Inf) && all(l==-Inf)){
    nbnds <- 0
    degen <- -1
  }

  xcurr <- xstart
  fvec <- fval### ????????????
  ## browser()
  gradflag <- options$Jacobian=="on"
  pcflags <- options$PrecondBandWidth
  ## %   Evaluate F and J
  if(!gradflag){
    A <- fdjacobian(func=funfcn,xcurr,fval=fvec,lb=l,ub=u,...)
    findiffevals <-n ### TODO modify to the right numbers
  }else{
    J <- JACval
    findiffevals <- 0
  }
  numFunEvals <- numFunEvals+findiffevals

  delbnd <- max(100*norm(xstart),1)
  if(is.vector(fvec)) fvec <- matrix(fvec)
  if(is.array(fvec)) fvec <- as.matrix(fvec)
  size_fvec <- dim(fvec)
  mm <- size_fvec[1]
  pp <- size_fvec[2]
  if(mm < n) stop("Note enough equations for the parameters,M<N")
  if(pp==2) dnewt <- fvec[,2]
  ## %   Determine gradient of the nonlinear least squares function
  g <- mtxmpy(A,fvec[,1],flag=-1)

  ## %   Evaluate F (initial point)
  val <- t(fvec[,1])%*%fvec[,1]
  ## Create the fault tolerance structure: no need in our case

  ## Main Loop: Generates feasible sequence x(iter), s.t. ||F(x(iter))|| is decreasing.
  while(!ex){
    ## Update:
    vres <- definv(g,x,l,u)
    v <- vres$v
    dv <- vres$dv

    gopt <- v*g
    optnrm <- base::norm(matrix(gopt),type="I") ## get the infinity norm
    r <- abs(pmin(u-x,x-l))
    degen <- min(r+abs(g))
    if(!nbnds) degen <--1

    ## Display:(ignore first)

    ## Test for Convergence
    diff <- abs(oval-val)
    oval <- val
    if((optnrm<tol1) && (posdef==1)){
      ex <- 1
      EXITFLAG <- 1
      if(iter==1) msgFlag <-100 else msgFlag <- EXITFLAG
    }else{
      if ((nrmsx < .9*delta) && (ratio > .25) && (diff < tol1*(1+abs(oval)))){
        ex <- 2
        EXITFLAG <- 3
        msgData <- paste('snls')
      }else{
        if((iter > 1) && (nrmsx < tol2)){
          ex <- 3
          EXITFLAG <- 2
          msgData <- paste('snls')
        }else{
          if(iter > itb){
            ex <- 4
            EXITFLAG <-0
            msgData <- paste('snls')
          }else{
            if(numFunEvals > maxfunevals){
              ex <- 4
              EXITFLAG <- 0
            }
          }
        }
      }
    }
    ## % Reset fault tolerance structure I think no need in R
    ## -------------
    ## %     Continue if ex = 0 (i.e., not done yet)
    if(ex == 0){
      ## Determine the trust region correction ******** using trdog()
      dd <- abs(v)
      ##D <- Matrix(diag(sqrt(dd)))
      if(length(dd)>1) D <- (diag(sqrt(dd))) else D <- matrix(sqrt(dd))
      sx <- rep(0,n)
      theta <- max(.95,1-optnrm)
      oposdef <- posdef
      trdogRes <- trdog(x,g,H=A,D,delta,dv,mtxmpy,pcmtx=aprecon,
                        pcoptions=pcflags,tol=pcgtol,kmax,theta,l,u,
                        Z,dnewt,'jacobprecon',...)
      posdef <- ifelse(is.null(trdogRes$posdef),oposdef,trdogRes$posdef)
      sx <- trdogRes$s
      snod <- trdogRes$snod
      qp <- trdogRes$qpval
      Z <- trdogRes$Z
      nrmsx <- norm(trdogRes$snod)
      npcg <- npcg + trdogRes$pcgit
      newx <- x+sx
      ## Newton Reflection:
      pertRes <- perturbTrustRegionReflective(newx,l,u)
      pert <- pertRes$pert
      newx <- pertRes$x
    ## -------------
      xcurr <- as.vector(newx)
    ## Evaluate F and J
      ## %   Evaluate F and J
      if(!gradflag){
        newfvec <- funfcn(xcurr,...)
        newA <- fdjacobian(func=funfcn,xcurr,fval=newfvec,lb=l,ub=u,...)
        findiffevals <-n ### TODO modify to the right numbers
      }else{
        res <- funfcn(xcurr,...)
        newfvec <- res$funval
        newA <- res$grad
        findiffevals <- 0
        numGradEvals <- numGradEvals+1
      }
      numFunEvals <- numFunEvals+1+findiffevals
      if(is.vector(newfvec)) newfvec <- matrix(newfvec)
      if(is.array(newfvec)) newfvec <- as.matrix(newfvec)
      size_fvec <- dim(newfvec)
      mm <- size_fvec[1]
      pp <- size_fvec[2]

      newval <- as.numeric(crossprod(newfvec[,1]))

    ## update the fault tolerance structure: In R, there is no need
    ## to check whether it is NaN or Inf or  complex.

    ## update the trust region radius:
    if(any(is.nan(newfvec)) || any(is.na(newfvec)) || any(newfvec==Inf) || any(newfvec==-Inf)){
     ## % Shrink the trust region if any element of the function vector
      ##% at the new step is not defined (i.e. it's inf/NaN/complex)
      delta <- min(nrmsx/20,delta/20)
    }else{
      newgrad <- mtxmpy(newA,newfvec[,1],flag=-1)
      aug <- .5*t(snod)%*%(dv*abs(g)*snod)
      ratio <- ifelse(qp==0,-Inf,as.numeric((0.5*(newval-val)+aug)/qp))
      if(ratio>=0.75 && nrmsx>=0.9*delta){
        delta <- min(delbnd,2*delta)
      }else{
        if(ratio<=0.25) delta <- min(nrmsx/4,delta/4)
      }
    }
    ## Accept or reject the trial point:
    if( 1 && (newval<val)){
      x <- newx
      val <- newval
      g <- newgrad
      A <- newA
      Z <- vector()
      fvec <- newfvec
      if(pp==2) dnewt <- newfvec[,2]## Extract the newton direction

    }
    iter <- iter+1
    } ## end  if(ex == 0)
  }## end  while(!ex)

  JACOB <- A #as(A,"Matrix")
  xcurr <-x
  if(computeLambda){
    LAMBDA.lower <- rep(0,length(l))
    LAMBDA.upper <- rep(0,length(u))
    argl <- abs(x-l) < active_tol
    argu <- abs(x-u) < active_tol
    LAMBDA.lower[argl]<-g[argl]
    LAMBDA.upper[argu]<--g[argu]
  }else{LAMBDA <- vector()}
  return(list(xcurr=xcurr, fvec=fvec,LAMBDA=LAMBDA,JACOB=JACOB,exitflag=EXITFLAG,msgdata=msgdata))
}


##' Scaling vector and derivative
##'
##' calculate distance to the bounds according to the sign of gradient.
##' @title Scaling vector and derivative
##' @param g gradient
##' @param x current point
##' @param l lower bounds
##' @param u upper bounds
##' @return a list
##' \itemize{
##' \item{"v"}{scaling vector}
##' \item{"dv"}{derivative}
##' }
##' @author Zhenglei Gao
definv <- function(g,x,l,u)
{
  ## Scaling vector and derivative
#   [v,dv]= DEFINEV(g,x,l,u) returns v, distances to the
#   %   bounds corresponding to the sign of the gradient g, where
#   %   l is the vector of lower bounds, u is the vector of upper
#   %   bounds. Vector dv is 0-1 sign vector

  n <- length(x);
  v <- rep(0,n);
  dv <- rep(0,n)
  arg1 <- (g < 0)  & (u <  Inf );
  arg2 <- (g >= 0) & (l > -Inf);
  arg3 <- (g < 0)  & (u == Inf);
  arg4 <- (g >= 0) & (l == -Inf);
  v[arg1]= (x[arg1] - u[arg1]);
  dv[arg1] = 1;
  v[arg2]  = (x[arg2] - l[arg2]);
  dv[arg2] = 1;
  v[arg3]  = -1;
  dv[arg3] = 0;
  v[arg4]  = 1;
  dv[arg4] = 0;
  return(list(v=v,dv=dv))
}


##' Jacobian matrix multiply function
##'
##' compute either A'*A*Y or A*Y or A'*Y
##' @title Jacobian matrix multiply function
##' @param A a matrix
##' @param Y a vector
##' @param flag indicate with multiply of A and Y
##' @return either A'*A*Y or A*Y or A'*Y
##' @author Zhenglei Gao
atamult <- function(A,Y,flag=0)
{
  if(flag==0) V <- t(A)%*%(A%*%Y) else{
    if(flag < 0) V <- t(A)%*%Y else V <- A%*%Y
  }
  return(V)
}


##' shake the current point slight away from bounds
##'
##' perturbs the current point X slightly to shake it loose from tight (less than DEL away)
##' bounds U and L to be strictly feasible. 
##' @title Perturb point from bounds
##' @param x current point
##' @param l lower bounds
##' @param u upper bounds
##' @param del a small value away from bounds
##' @param y is the reflected point y
##' @param sigma perterb y w.r.t sigma
##' @return perturbed x(and y)
##' @author Zhenglei Gao
perturbTrustRegionReflective <- function(x,l,u,del=100*.Machine$double.eps,y=NULL,sigma=NULL)
{
#   perturbs the current
#   %   point X slightly to shake it loose from tight (less than DEL away)
#   %   bounds U and L to be strictly feasible.  also perturbs
#  %   the reflected point Y with respect to SIGMA,
  if(min(abs(u-x))<del | min(abs(x-l))<del){
    upperi <- (u-x) < del;
    loweri <- (x-l) < del;
    x[upperi] = x[upperi] - del;
    x[loweri] = x[loweri] + del;
    if(!is.null(y)){
      y[upperi] = y[upperi] - del*sigma[upperi];
      y[loweri] = y[loweri] + del*sigma[loweri];
    }
    pert<-TRUE
  }else pert <- FALSE
  return(list(x=x,y=y,pert=pert))
}

##' Banded preconditioner function for least-squares problems.
##' 
##' produces the sparse nonsingular upper triangular matrix RPCMTX such that
##' RPCMTX'*RPCMTX is a preconditioner for the
##' matrix M = DM*(A'*A)*DM + DG, where DM is a positive
##' diagonal matrix, DG is a non-negative diagonal matrix,
##' and A is sparse rectangular matrix with more rows than columns.
##' PPVEC is the associated permutation (row) vector and UPPERBW
##' specifies the upperbandwidth of RPCMTX. 
##' @title Banded preconditioner for LS problems
##' @param A a rectangular matrix with more rows than columns
##' @param upperbandw upper bandwidth
##' @param DM a positive diagonal matrix
##' @param DG a non negative diagonal matrix
##' @param ... other parameters.
##' @return a list with RPCMTX an upper triangular matrix such that RPCMTX'RPCMTX 
##' is a preconditioner matrix for M=DM*A*A'*DM+DG, and ppvec being permutation vector.
##' @author Zhenglei Gao
aprecon <- function(A=NULL,upperbandw=NULL,DM,DG,...)
{
  tmp <- dim(A)
  n <- ncol(DM)## ??????
  ##RPCMTX <- Matrix(diag(n))
  RPCMTX <- (diag(n))
  ppvec <- 1:n
  ## % In case "A" isn't really A, but something else to use with JacobMult function.
  if((!(is.numeric(A)))||(is.null(A))||(tmp[1]<tmp[2])||(tmp[2]!=n)){

    d1 <- diag(DM)
    d2 <- diag(DG)
    dd <- sqrt(d1*d1+abs(d2))
    ## RPCMTX <- matrix(0,n,n) ????????
    ## diag(RPCMTX) <- dd
    ##RPCMTX <- Matrix(diag(dd))
    if(length(dd)==1) RPCMTX <-matrix(dd) else RPCMTX <- (diag(dd))

  }else{
    epsi <- sqrt(.Machine$double.eps)
    if(is.null(upperbandw)) upperbandw <-0
    ## Form Matrix M
    TM <- A%*%DM
    if(upperbandw==0){
      ## M <- t(TM)%*%TM+DG
      M <- crossprod(TM)+DG
      dnrms <- sqrt(apply(M*M,1,sum)) ##??????
      d <- pmax(sqrt(dnrms),epsi)
      ##RPCMTX <- Matrix(diag(d))
      if(length(d)==1) RPCMTX <-matrix(d) else RPCMTX <- (diag(d))
      # ppvec <- 1:n ## no need to do it
      ## But what if RPCMTX is singular?
    }else{
      if(upperbandw >= (n-1)) {## Attempt sparse QR
        dgg <- sqrt(diag(DG))
        if(length(dgg)==1) DDG <- matrix(dgg) else DDG <- (diag(dgg))
        TM <- rbind(TM,DDG)##???
        ##p <- colamd(TM)????
        ##RPCMTX <- qr(TM[,p])
        ##RPCMTX <- RPCMTX[1:n,1:n]
        ##ppvec <- p
        tmp <- qr(TM)
        if(is(tmp,"sparseQR")){
          RPCMTX <- qr.R(tmp)
          ppvec <- order(tmp@q)
        }else{
          if(is(tmp,"qr")){
            RPCMTX <- qr.R(tmp)
            ppvec <- tmp$pivot
          }

        }

        if(n>1) RPCMTX <- RPCMTX[1:n,1:n]

        ## %    Modify for singularity?
        mdiag <- min(abs(diag(RPCMTX)))
        lambda <- 1
        while(lambda<Inf && mdiag < epsi){
          ##TM <- rBind(A%*%DM,DDG+Matrix(lambda*diag(n)))##???
          TM <- rbind(A%*%DM,DDG+(lambda*diag(n)))##???
          lambda <- 4*lambda
#           p <- colamd(TM)
#           RPCMTX <- qr(TM[,p])
#           RPCMTX <- RPCMTX[1:n,1:n]
#           ppvec <- p
          tmp <- qr(TM)
          if(is(tmp,"sparseQR")){
            RPCMTX <- qr.R(tmp)
            ppvec <- tmp@q
          }else{
            if(is(tmp,"qr")){
              RPCMTX <- qr.R(tmp)
              ppvec <- tmp$pivot
            }

          }

          RPCMTX <- RPCMTX[1:n,1:n]
          mdiag <- min(abs(diag(RPCMTX)))
        }
      }else{
        if(upperbandw>0 && (upperbandw < (n-1))){
          ## % Band approximation.
          M <- crossprod(TM)+DG
          p <- 1:n
          ## M <- tril(triu(M(p,p),-upperbandw),upperbandw)??????
          M <- band(M,-upperbandw,+upperbandw)
          RPCMTX <- matrix(0,n,n)
          tmp <- try(chol(M),silent=TRUE)
          lambda <- 1
          while(attr(tmp,"class")=="try-error"){
            M <- M+lambda*(diag(n))
            tmp <- try(chol(M),silent=TRUE)
            lambda <- 4*lambda
          }
          RPCMTX <- tmp
          ppvec <- p
        }else{
          stop("Invalid upperbandw")
        }
      }
    }


  }


  return(list(R=RPCMTX,permR=ppvec))
}



##' compute new direction using a Preconditioned conjugate gradient procedure.
##'
##' for q(p)=.5p'Mp+g'p, M=DM.H.DM+DG
##' @title Preconditioned conjugate gradients
##' @param DM a matrix
##' @param DG a matrix
##' @param g a gradient vector
##' @param kmax maximum number of CG iterations.
##' @param tol stopping tolerance 
##' @param mtxmpy function to compute products with the Hessian matrix H
##' @param H Hessian
##' @param R the cholesky factor of the preconditioner of M
##' @param pR a permutation vector so R'R=M(pR,pR)
##' @param callerflag whether to use "jacobprecon"
##' @param pcoptions Preconditioner option
##' @param ... other parameters passed to mtxmpy and others
##' @return a list
##' \itemize{
##' \item{"p"}{computed direction}
##' \item{"posdef"}{whether positive definite}
##' \item{"k}{number of iterations}
##' }
##' @author Zhenglei Gao
pcgr <- function(DM,DG,g,kmax,tol,mtxmpy=atamult,
                 H,R,pR,callerflag="jacobprecon",
                 pcoptions=Inf,...)
{## Preconditioned conjugate gradients
  n <- nrow(DG)
  r <- -g
  p <- rep(0,n)

  ## precondition
  z <- preproj(r,R,pR)
  znrm <- norm(z)
  stoptol <- tol*znrm
  inner2 <- 0
  inner1 <- t(r)%*%z
  inner1 <- inner1[1,1]
  posdef <- 1

  kmax <- max(kmax,1)

  for(k in 1:kmax){
    if(k==1){
      d <- z
    }else{
      beta <- inner1/inner2
      d <- z+beta*d
    }
    ww <- DM %*% d
    if(callerflag=="jacobprecon") w <- mtxmpy(H,ww)
    ww <- DM%*%w+DG%*%d
    denom <- t(d)%*%ww
    denom <- denom[1,1]
    if(denom <=0){
      if(norm(d)==0){
        p <- d
      }else p <- d/norm(d)
      posdef <-0
      break
    }else{
      alpha <- inner1/denom
      p <- p+alpha*d
      r <- r-alpha*ww
    }
    z <- preproj(r,R,pR)
    if(norm(z)<=stoptol) break
    inner2 <- inner1
    inner1 <- t(r)%*%z
    inner1 <- inner1[1,1]
  }

  if(pcoptions==Inf) k <- 0
  return(list(v1=p,posdef=posdef,pcgit=k))


}


##' Apply preconditioner to vector r
##'
##' The conceptual preconditioner is H(ppvec,ppvec) = RPCMTX'*RPCMTX
##' @title Apply preconditioner to vector r
##' @param r a vector
##' @param RPCMTX the preconditioner
##' @param ppvec the permutation vector
##' @return a preconditioned w
##' @author Zhenglei Gao
preproj <- function(r,RPCMTX=NULL,ppvec=NULL)
{
  ## %PREPROJ Apply preconditioner
  n <- length(r)
  if(is.null(ppvec)) ppvec <- 1:n
  if(is.null(RPCMTX)) RPCMTX <- (diag(n))
  ## wbar <- t(RPCMTX)\r[ppvec]
  wbar <- solve(t(RPCMTX),r[ppvec])
  ## w[ppvec,1] <- RPCMTX\wbar
  w <- rep(NA,n)
  w[ppvec] <- as.vector(solve(RPCMTX,wbar))
  return(w)
}


##' Calculate the jacobian by finite differece with lower and upper bounds
##'
##' modified according to numDeriv::jacobian. Calculate the jacobian by finite differece with lower and upper bounds
##' @title Jacobian by finite difference
##' @param func function 
##' @param x parameter vector
##' @param method see \code{jacobian}
##' @param method.args control parameters for the choosen method
##' @param fval function value at x
##' @param lb lower bound
##' @param ub upper bound
##' @param ... other parameters to be passed into func
##' @return the Jacobian of the function at x
##' @author Zhenglei Gao
##' @export
fdjacobian <- function(func,x,method="Richardson",method.args=list(eps=sqrt(.Machine$double.eps)),fval=NULL,lb=-Inf,ub=Inf,...)
{
  ## modified according to numDeriv::jacobian
  if(is.null(fval))f <- func(x, ...) else f <- fval
  
  n <- length(x)   #number of variables.
  ##browser()
  if(length(lb)<n) lb <- rep(lb,n)
  if(length(ub)<n) ub <- rep(ub,n)
  if(method=="simple"){
    #  very simple numerical approximation
    args <- list(eps=sqrt(.Machine$double.eps)) # default
    args[names(method.args)] <- method.args
    eps <- args$eps
    df <-matrix(NA, length(f), n)
    for (i in 1:n) {
      dx <- x
      tmp <- dx[i]+eps
      if(tmp<=ub[i] && tmp>=lb[i]){
        dx[i] <- tmp
        df[,i] <- (func(dx, ...)-f)/eps

      }else{
        if(tmp>ub[i]){
          tmp <- dx[i]-eps
          dx[i] <- tmp
          df[,i] <- (func(dx, ...)-f)/eps

        }else{
          tmp <- dx[i]+eps
          dx[i] <- tmp
          df[,i] <- (func(dx, ...)-f)/eps

        }

      }

    }
    return(df)
  }else
    if(method=="complex"){ # Complex step gradient
      # Complex step Jacobian
      eps <- .Machine$double.eps
      h0  <-  rep(0, n)
      h0[1] <- eps * 1i
      v <- try(func(x+h0, ...))
      if(inherits(v, "try-error"))
        stop("function does not accept complex argument as required by method 'complex'.")
      if(!is.complex(v))
        stop("function does not return a complex value as required by method 'complex'.")

      h0[1]  <- 0
      jac <- matrix(NA, length(v), n)
      jac[, 1] <- Im(v)/eps
      if (n == 1) return(jac)
      for (i in 2:n) {
        h0[i] <- eps * 1i
        jac[, i] <- Im(func(x+h0, ...))/eps
        h0[i]  <- 0
      }
      return(jac)
    } else
      if(method=="Richardson"){
        args <- list(eps=1e-4, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7),
                     r=4, v=2, show.details=FALSE) # default
        args[names(method.args)] <- method.args
        eps <- args$eps
        d <- args$d
        r <- args$r
        v <- args$v
        a <- array(NA, c(length(f),r, n) )

        h <- abs(d*x)+eps*(abs(x) < args$zero.tol)
        for(k in 1:r)  { # successively reduce h
          for(i in 1:n)  {
            ##if(i==2) browser()
            a[,k,i] <- (func(x + h*(i==seq(n)), ...) -
                          func(x - h*(i==seq(n)), ...))/(2*h[i])
            #if((k != 1)) a[,(abs(a[,(k-1),i]) < 1e-20)] <- 0 #some func are unstable near zero
          }
          h <- h/v     # Reduced h by 1/v.
        }
        for(m in 1:(r - 1)) {
          a <- (a[,2:(r+1-m),,drop=FALSE]*(4^m)-a[,1:(r-m),,drop=FALSE])/(4^m-1)
        }
        # drop second dim of a, which is now 1 (but not other dim's even if they are 1
        return(array(a, dim(a)[c(1,3)]))
      } else stop("indicated method ", method, "not supported.")
}


##' Calculate the reflective 2d trust region trial step with box constraints
##'
##' the trial step s is determined by the best of the 3 steps: sclaed gradient, 2d trust region solution, 
##' and the reflected 2d trust region solution. The 2d subspace is defined by the scaled gradient direction and
##' a conjugate gradient process(an approximate Newton step or a direction of the negative curvature)
##' @title Reflected 2d trust region step with box constraints
##' @param x current parameter value
##' @param g gradient
##' @param H Hessian
##' @param D a matrix
##' @param delta trust region 
##' @param dv a vector
##' @param mtxmpy Jacobian-matrix multiply function
##' @param pcmtx pcg multiply function
##' @param pcoptions pcg otpions
##' @param tol tolerance
##' @param kmax max k
##' @param theta default 0.95
##' @param l lower bound
##' @param u upper bound
##' @param Z a matrix
##' @param dnewt default NULL, 
##' @param preconflag PCG flag
##' @param ... other parameters that passes to functions
##' @return a list 
##' @author Zhenglei Gao
trdog <- function(x,g,H,D,delta,dv,
                  mtxmpy=atamult,pcmtx=aprecon,pcoptions=Inf,
                  tol=0.1,kmax=2,theta=0.95,
                  l=-Inf,u=Inf,Z=NULL,
                  dnewt=NULL,preconflag="jacobprecon",...)
{
  ## a clone of the trdog.m function in Matlab
#   %TRDOG Reflected (2-D) trust region trial step (box constraints)
#   %
#   % [s,snod,qpval,posdef,pcgit,Z] = TRDOG(x,g,H,D,delta,dv,...
#      %                 mtxmpy,pcmtx,pcoptions,tol,theta,l,u,Z,dnewt,preconflag);
#   %
#   %   Determine the trial step `s', an approx. trust region solution.
# %   `s' is chosen as the best of 3 steps: the scaled gradient
# %   (truncated to  maintain strict feasibility),
#   %   a 2-D trust region solution (truncated to remain strictly feas.),
#   %   and the reflection of the 2-D trust region solution,
#   %   (truncated to remain strictly feasible).
#   %
#   %   The 2-D subspace (defining the trust region problem) is defined
#   %   by the scaled gradient direction and a CG process (returning
#    %   either an approximate Newton step of a direction of negative curvature.
#      Driver functions are: SNLS, SFMINBX
#        SNLS actually calls TRDOG with the Jacobian matrix (and a special
#      Jacobian-matrix multiply function in MTXMPY).
#
#     Copyright 1990-2011 The MathWorks, Inc.
#      $Revision: 1.1.6.3 $  $Date: 2011/05/09 01:16:31 $
#


#       Initialization
  ##require(Matrix)
  eps <- .Machine$double.eps
  n = length(g);  # g is the gradient
  pcgit = 0;
  grad = D%*%g;    ##
  DM = D;
  #%DG = sparse(1:n,1:n,full(abs(g).*dv)); ## Matlab implementation
  # DG <- diag(abs(g)*dv)  ## R implementation without sparse matrix definition
  ##DG <- Matrix(0, nrow = n, ncol = n, sparse = TRUE)
  DG <- matrix(0, nrow = n, ncol = n)
  diag(DG) <- abs(g)*dv
  posdef = 1;
  pcgit = 0;
  tol2 = sqrt(eps);
  v1 = dnewt;
  qpval1 = Inf;
  qpval2 = Inf;
  qpval3 = Inf;
  ##browser()
  ## % DETERMINE A 2-DIMENSIONAL SUBSPACE
  if(is.null(Z) || length(Z)==0){
    if(is.null(v1) || length(v1)==0){
      if(preconflag=="jacobprecon") {
        tmp <- pcmtx(H,pcoptions,DM,DG,...)
        R <- tmp$R
        permR <- tmp$permR
      }else{
        stop("Not implemented yet!!")
      }
      if(tol <= 0) tol <- 0.1
      tmp <- pcgr(DM,DG,grad,kmax,tol,mtxmpy,H,R,permR,preconflag,pcoptions,...)
      v1 <- tmp$v1
      posdef <- tmp$posdef
      pcgit <- tmp$pcgit
    }

    nrmv1 <- norm(v1)
    if(nrmv1>0) v1 <- v1/nrmv1
    if(n>1) Z <- matrix(NA,n,2) else Z <- matrix(NA,n,1)

    Z[,1] <- v1
    if(n>1){
      if(posdef<1){
        v2 <- D%*%sign(grad)
        nrmv2 <- norm(v2)
        if(nrmv2>0) v2 <- v2/nrmv2
        v2 <- v2 - v1*as.numeric(t(v1)%*%v2)
        nrmv2 <- norm(v2)
        if(nrmv2>tol2) {
          v2 <- as.vector(v2/nrmv2)
          Z[,2] <- v2
        }
      }else{
        nrmgrad <- norm(grad)
        if(nrmgrad>0) v2 <- grad/nrmgrad else v2 <- grad
        v2 <- v2 - as.numeric((t(v1)%*%v2))*v1
        nrmv2 <- norm(v2)
        if(nrmv2>tol2) {
          v2 <- as.vector(v2/nrmv2)
          Z[,2] <- v2
        }
      }
    }
  }## end  if(is.null(Z))
  W <- DM%*%Z
  if(preconflag=="jacobprecon") {
    WW <- mtxmpy(H,W,0)
  }else stop("not implemented yet!!")
  W <- DM%*%WW
  MM <- t(Z)%*%W+t(Z)%*%DG%*%Z
  rhs <- t(Z)%*%grad
  tmp <- trust.exact(rhs,MM,delta)
  st <- tmp$s
  qpval <- tmp$val
  po <- tmp$posdef
  fcnt <- tmp$count
  lambda <- tmp$lambda
  ss <- Z%*%st
  s <- abs(diag(D))*ss
  ssave <- s
  sssave <- ss
  stsave <- st

  ## check direction for NaNs
  if(any(is.nan(as.vector(s)))) stop("erros in s, NaNs")
  ## Truncate the TR solution
  arg <- as.vector(abs(s)>0)
  if(!any(arg)){
    alpha <- 1
    mmdis <- 1
    ipt <- vector()
  }else{
    dis <- pmax((u[arg]-x[arg])/s[arg], (l[arg]-x[arg])/s[arg])
    ipt <- which.min(dis)
    mmdis <- dis[ipt]
    mdis <- theta * mmdis
    alpha <- min(1,mdis)
  }
  s <- alpha*s
  st <- alpha*st
  ss <- alpha*ss
  qpval1 <- as.numeric(t(rhs)%*%st+t(.5*st)%*%MM%*%st)
  if(n>1){
    ##%   Evaluate along the reflected direction?
    qpval3 <- Inf
    ssssave <- mmdis*sssave
    ## here in this case, it could be that a Inf*0
    if(any(is.nan(ssssave))) normssssave <-Inf else normssssave <- norm(ssssave)
    if(normssssave<.9*delta){
      r <- mmdis*ssave
      nx <- x+r
      stsave <- mmdis*stsave
      qpval0 <- as.numeric(t(rhs)%*%stsave+t(0.5*stsave)%*%MM%*%stsave)
      if(preconflag=="jacobprecon") ng <- mtxmpy(H,r,0)
      ng <- ng+g
      ngrad <- D%*%ng
      ngrad <- ngrad +DG%*%ssssave
      ## %      nss is the reflected direction
      nss <- sssave
      nss[ipt]<- -nss[ipt]
      ZZ <-as.vector(nss/(norm(nss)))
      W<-DM%*%ZZ
      if(preconflag=="jacobprecon") WW <- mtxmpy(H,W,0)
      W <- DM%*%WW
      MM <- t(ZZ)%*%W+t(ZZ)%*%DG%*%ZZ
      nrhs <- t(ZZ)%*%ngrad
      tmp1 <- quad1d(nss,ssssave,delta)
      nss <- tmp1$nx
      tau <- tmp1$tau
      nst <- tau/norm(nss)
      ns <- as.vector(abs(diag(D))*nss)
      ## check direction for NaNs
      if(any(is.nan(ns))) stop("erros in ns, NaNs")
      ## % Truncate the reflected direction?
      arg <- as.vector((abs(ns) > 0))
      if(!(any(arg))) alpha <- 1 else{
        dis <- pmax((u[arg]-nx[arg])/ns[arg], (l[arg]-nx[arg])/ns[arg]);
        mdis <- min(dis);
        mdis <- theta*mdis;
        alpha <- min(1,mdis)
      }
      ns <- alpha*ns
      nst <- alpha*nst
      nss <- alpha*nss
      qpval3 <- qpval0+as.numeric(t(nrhs)%*%nst+t(0.5*nst)%*%MM%*%nst)
    }
    gnorm <- norm(grad)
    ZZ <- grad/(gnorm+(gnorm==0))## % Protect against norm of 0
    W <- DM%*%ZZ
    if(preconflag=="jacobprecon") WW <- mtxmpy(H,W,0)
    W <- DM%*%WW
    MM <- t(ZZ)%*%W+t(ZZ)%*%DG%*%ZZ
    rhs <- t(ZZ)%*%grad
    tmp <- trust.exact(rhs,MM,delta)
    st <- tmp$s
    qpval <- as.numeric(tmp$val)
    po <- tmp$posdef
    fcnt <- tmp$count
    lambda <- tmp$lambda

    ssg <- ZZ%*%st
    sg <- as.vector(abs(diag(D))*ssg)
    if(any(is.nan(sg))) stop("erros in sg, NaNs")
    arg <- as.vector((abs(sg) > 0))
    if(!(any(arg))) alpha <- 1 else{
      dis = pmax((u[arg]-x[arg])/sg[arg], (l[arg]-x[arg])/sg[arg])
      mdis <- min(dis)
      mdis <- theta*mdis;
      alpha <- min(1,mdis)
    }
    sg <- alpha*sg
    st <- alpha*st
    ssg <- alpha*ssg
    qpval2 <- as.numeric(t(rhs)%*%st+t(0.5*st)%*%MM%*%st)
  }

  ## Choose the best of s, sg, ns.
  if(qpval2 <= min(qpval1,qpval3)){
    qpval <- qpval2
    s <- sg
    snod <- ssg
  }else{
    if(qpval1 <= min(qpval2,qpval3)){
      qpval <- qpval1
      snod <- ss
    }else{
      qpval <- qpval3
      s <- ns+r
      snod <- nss+ssssave
    }

  }
  return(list(s=s,snod=snod,qpval=qpval,posdef=posdef,pcgit=pcgit,Z=Z))

}



##' 1-d quadratic zero finder for trust region step
##'
##' 1-d quadratic of the form ay^2+by+c 
##' @title quadratic zero finder for 1D
##' @param x a vector
##' @param ss a vector
##' @param delta a number
##' @return a list of new x value and the tau(min(1, step to zero))
##' @author Zhenglei Gao
quad1d <- function(x,ss,delta)
{
  #   %QUAD1D  1D quadratic zero finder for trust region step
  #   %
  #   % [nx,tau] = quad1d(x,ss,delta) tau is min(1,step-to-zero)
  #   % of a 1-D quadratic ay^2 + b*y + c.
  #   % a = x'*x; b = 2*(ss'*x); c = ss'*ss-delta^2). nx is the
  # % new x value, nx = tau*x;
  #
  # % Algorithm:
  # % numer = -(b + sign(b)*sqrt(b^2-4*a*c));
  # % root1 = numer/(2*a);
  # % root2 = c/(a*root1);   % because root2*root1 = (c/a);

  a <- as.numeric(crossprod(x))
  b <- as.numeric(2*(t(ss)%*%x))
  c <- as.numeric(crossprod(ss))-delta^2

  numer <- -(b + sign(b)*sqrt(b^2-4*a*c))
  r1 <- numer/(2*a)
  r2 <- c/(a*r1)

  tau <- max(r1,r2)
  tau <- min(1,tau)
  if(tau <= 0) stop("error in quad1d")
  nx <- tau*x
  return(list(nx=nx,tau=tau))
}

## utilities
##' calculate the gradient at the current parameter.
##'
##' using the chain rule of lease square problem.
##' @title Derive the gradient 
##' @param modelfun model function
##' @param par parameter vector
##' @param obs observed values
##' @param ... other parameters to be passed into modelfun
##' @return the gradient at par
##' @author Zhenglei Gao
gradfun_ls <- function(modelfun,par,obs,...)
{
  if(is.null(names(par))) names(par) <- pnames
  np <- length(par)
  resfun <- function(par,obs, ...){
    modelfun(par,...)-obs
  }
  res <- resfun(par,obs=obs)
  Jmat <- jacobian(modelfun,x=par)
  g <- t(Jmat)%*%res
  return(g)
}

##' Derive the objective function including the gradient and heissian 
##'
##' using the numDeriv functions jacobian and genD
##' @title derive the objective function
##' @param modelfun model function
##' @param par parameter values
##' @param obs observed values
##' @param Q whether to calculate the hessian matrix
##' @param pnames parameter names
##' @param ... other parameters to be passed into modelfun
##' @author Zhenglei Gao
##' @examples
##' \dontrun{
##' farb<- function(x=c(3,1),t=1:10,...){
##'   x1 <- x[1]
##'   x2 <- x[2]
##'   x[1]+x[2]*t
##'   }
##'   y <- farb(x=c(3,1),t=1:10)+rnorm(10,0,2)
##'   t <- 1:10
##'   x <- c(3,1)
##'   objfun <- objfun_ls(farb,c(3,1),obs=y)
##'   }
objfun_ls <- function(modelfun,par,obs,Q=FALSE,pnames=NULL,...)
{
  np <- length(par)
  if(is.null(names(par))) names(par) <- pnames
  resfun <- function(par,obs, ...){
    modelfun(par,...)-obs
  }
  objfun <- function(par,obs,...)
  {
    0.5*sum((modelfun(par,...)-obs)^2)
  }

  res <- resfun(par,obs=obs)
  ## calculate the jacobian of model function(returning a vector)
  id <- which(is.na(res))

  if(Q==TRUE)
  {
    D <- genD(modelfun,x=par)$D
    Jmat <- D[,1:np]
    g <- t(Jmat[-id,])%*%res[-id]
    Hmat1 <- D[,(np+1):(np+np*(np+1)/2)]
    LT <- res[-id]%*%Hmat1[-id,]
    M <- matrix(NA,np,np)
    k <- 0
    for(i in 1:np)
      for(j in 1:i)
      {
        k <- k+1
        M[i,j] <- LT[k]
      }
    ## M is symmetric
    ##browser()
    M[upper.tri(M)] <- t(M)[upper.tri(M)]
    B <- t(Jmat[-id,])%*%Jmat[-id,]+M
  }else {
    Jmat <- jacobian(modelfun,x=par)
    g <- t(Jmat[-id,])%*%res[-id]

    B <- t(Jmat[-id,])%*%Jmat[-id,]
  }
  ## using numerical approximation : grad(objfun,x=c(3,1),obs=y)
  ## using numerical approximation : Hmat <- hessian(objfun,x=par,obs=obs)

  list(value=0.5*sum(res^2,na.rm=TRUE),gradient=g,hessian=B)
}

