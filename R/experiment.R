## experimenting Scripts for cleaner coding!
DefDiff <- function(t, state, parms,mod_vars,diffeqns) {
  ## an updated version of mkindiff
  ## @example DefDiff(t,state,parms, mod_vars, diffeqns=mkinmodini$diffs)
  time <- t
  diffs <- vector()
  for (box in mod_vars)
  {
    diffname <- paste("d", box, sep="_")
    diffs[diffname] <- with(as.list(c(time,state, parms)),
                            eval(parse(text=diffeqns[[box]])))
  }
  ##https://stat.ethz.ch/pipermail/r-sig-dynamic-models/2010q2/000031.html
  #bady <- (!is.finite(diffs))|(diffs<=0)
  #diffs[bady] <- 0 
  return(list(c(diffs)))
}

compDiff <- function(t,state,parms){
  
}

## Calculate DTx (Formation fractions are calcualted in the evaluation function)
calcDT <- function(type="SFO",...){
  if(type=="SFO"){
    DT50 = log(2)/k_tot
    DT90 = log(10)/k_tot
  }
  if (type == "FOMC") {
    DT50 = beta * (2^(1/alpha) - 1)
    DT90 = beta * (10^(1/alpha) - 1)
  }
  if (type == "DFOP") {
    f <- function(t, x) {
      ((g * exp( - k1 * t) + (1 - g) * exp( - k2 * t)) - (1 - x/100))^2
    }
    fDT <- function(t, x) {
      ((g * exp( - k1 * t) + (1 - g) * exp( - k2 * t)) - (1 - x/100))
    }
    
    DTmax <- 1000
    DTmax1 <- log(2)/min(k1,k2)
    DTmin <- log(2)/max(k1,k2)
    if(DTmax1==Inf) DTmax1 <- .Machine$double.xmax
    resUniroot <- uniroot(fDT,c(DTmin,DTmax1),x=50)
    DT50.o <- resUniroot$root
    if(resUniroot$f.root > 1e-4) print("Uniroot calculation for DFOP DT50 might be wrong!Do an optimization, check if the DT50 in the output make sense!")
    resOptim <- optimize(f, c(DTmin, DTmax1), x=50)
    DT50.o1 <- resOptim$minimum
    if(resOptim$objective > 1e-4) print("One dimentional optimization calculation for DFOP DT50 might be wrong!Doing a uniroot calculation, check if the DT50 in the output make sense!")
    DT50.o <- ifelse(resUniroot$f.root>resOptim$objective, DT50.o1,DT50.o)
    ##DT50 = ifelse(DTmax - DT50.o < 0.1, NA, DT50.o)
    DT50 <- DT50.o            
    DTmin <- log(10)/max(k1,k2)
    DTmax1 <- log(10)/min(k1,k2)
    if(DTmax1==Inf) DTmax1 <- .Machine$double.xmax
    DT90.o <- uniroot(fDT,c(DTmin,DTmax1),x=90)$root
    DT90.o1 <- optimize(f, c(DTmin, DTmax1), x=90)$minimum
    DT90.o <- ifelse(f(DT90.o,90)>f(DT90.o1,90), DT90.o1,DT90.o)
    DT90 <- DT90.o
  }
  if (type == "HS") {
    
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
  return(c(DT50,DT90))
}

## set up differential equation -- replace some parameters with some values/strings.
## The question is where to implement the parameter transformation.
## Keep the differential equation the same but transform back to the ode parms
## from the optim parms. 
## Or transform the objective functions into the original scale?

## Second approach is better!! since we want to use the transformed parms in the
## optimization and original scale in final Hessian/Jacobian calculation.The hessian is calculated always relative to the
## optim parm scales.

## for transformations using log, if you fix your parameter at 0 or 1, then don't
## do the transformation!!

## flexible control of each parameter?? if(met$k$transform==NULL){ ## no need of this step
##  if(!is.null(KinOptions$transformk)){
##    k_opt <- KinOptions$transformk(k)
##    f_opt <- KinOptions$transformff(ff)
##    g_opt <- # KinOptions$transformff(ff)
## }
## put them into the two sets differential equations.
## }

## NONONO, maybe we don't need that since P is only a vector with named elements.
setDiff <- function(){
}

RHS <- function(){
  ## should be an updated version of kin_mod
}



KinOptions <- function(...){
  default <- kingui.control()
  select <- list(...)
  for(nam in names(select)){
    default[[nam]] <- select[[nam]]
  }
  return(default)
}