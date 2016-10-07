

## Full function to define the model and the observed data.##
##' Full function to set up a kinetic model with one or more compartments and the data.
##'
##' GUI version of \code{\link{mkinmod}}. The function takes a specification,
##' consisting of a list of the compartments in the data. Each compartment is
##' again represented by a list, specifying the kinetic model type, reaction or
##' transfer to other observed compartments, the initial parameter values,
##' lower and upper bounds, fixed or not, and observed data.
##'
##'
##' @param ...  Each list cell represents a
##' comparment which contains a list of comonents including 'type'(kinetic
##' reaction type, single first order kinetics "SFO" are implemented for all
##' compartments, while "FOMC", "DFOP" and "HS" can additionally be chosen for
##' the first variable which is assumed to be the source compartment), each
##' parameter name(a list of 'ini','fixed','lower','upper'),'residue'(measured
##' concentrations),'time'(sampling time),'weight'(weights to be used, defaul
##' 1),'sink'( Default TRUE, tranformation to unspecified compartments.),'to'(a
##' vector of compartment names that the source compartment will be transferred
##' to).
##' @param inpartri Input parameterization.
##' @param outpartri Output parameterization.
##' @param data Optional. Data can be read from a data file.
##' @param weight Optional. General weighting schemes using a weight matrix.
##' @return A list of class 'mkinmod.full' for use with
##' \code{mkinfit.full},\code{\link{IRLSkinfit.full}} and
##' \code{\link{mcmckinfit.full}} containing:
##'
##' \item{diffs}{ A vector of string representations of differential equations,
##' one for each modelling compartment.}
##' \item{parms}{ A vector of parameter
##' names occurring in the differential equations.}
##' \item{map}{A list
##' containing named character vectors for each compartments in the model.}
##' \item{parms.ini}{Initial values for all kinetic parameters in the model.}
##' \item{state.ini}{Initial state values for all compartments in the model.}
##' \item{lower}{Lower bounds for the parameters(including state variables) to
##' be optimized.}
##' \item{upper}{upper bounds for the parameters(including state
##' variables) to be optimized.}
##' \item{fixed_parms}{The names of the kinetic
##' parameters that are fixed during optimization.}
##' \item{fixed_initials}{ The
##' names of the initial states that are fixed during optimization.}
##' \item{residue}{The observed data matrix with a time column.}
##' \item{weightmat}{The weights matrix.}
##' \item{ff}{A vector of string
##' representations of the transformation between the formation fractions in
##' the model and the transfomed formation fractions in the optimization
##' process.}
##' @note \code{mkinmod} is a deprecated version doing the same for the use with
##' \code{mkinfit.gui},\code{IRLSkinfit.gui} and
##' \code{mcmckinfit.gui}.
##' @author Zhenglei Gao
##' @seealso \code{\link{mkinmod}}, \code{\link{completeCompound}}
##' @keywords Kinetic-Model
##' @return A list with  "mkinmod.full" in the class attribute. 
##' @examples
##'
##' SFO_SFO_full <- mkinmod.full(Parent = list(type = "SFO", to = "Metab", sink = TRUE,
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
##' 12.7, 11.4, 11.6, 10.9,  9.5,  7.6))
##' )
##' @keywords Kinetic-Models
##' @export
##' @exportClass mkinmod.full
mkinmod.full <- function(...,inpartri=c('default','water-sediment','advanced'),outpartri=c('default','water-sediment','advanced'),data=NULL,weight=NULL,autoInit=TRUE)
{## Example usage ##############################################################
 ## a <- mkinmod.full(
 ##    parent = list(type = "SFO", to = c("A1", "B1", "C1"), sink = FALSE),
 ##    A1 = list(type = "SFO", to = "A2"),
 ##    B1 = list(type = "SFO"),
 ##    C1 = list(type = "SFO"),
 ##    A2 = list(type = "SFO"),inpartri='water-sediment',outpartri='water-sediment',data=schaefer07_complex_case,weight=NULL)
 ##################################
 if("logging" %in% loadedNamespaces()){
   logall <- TRUE
 }else{
   logall <- FALSE
 }
 
 inpartri <-  match.arg(inpartri)
 outpartri <-  match.arg(outpartri)
 
 spec <- list(...)
 obs_vars <- names(spec)
 ncompart <- length(obs_vars)
 
 spec <- lapply(1:ncompart,function(x){completeCompound(spec[[x]],varname=obs_vars[x],first=(x==1),inpartri=inpartri,outpartri=outpartri,data=data,weight=weight,autoInit=autoInit)})
 names(spec) <- obs_vars
 parentname <- obs_vars[1]
 ## ########################
 
 ## ## initialize the output varibles,common for all parametrizations.
 parms <- vector()
 parms.ini <- vector()
 parms.lower <- vector()
 parms.upper <- vector()
 state.ini <- vector()
 state.ini.orig <- vector()
 state.lower <- vector()
 state.upper <- vector()
 diffs <- vector()
 ff <- vector()
 sinkT <- sapply(spec,function(x) x$sink)
 residue <- NULL
 data0 <- NULL
 weightmat <- NULL
 map <- list()
 fixed_parms <- NULL
 fixed_flag <- NULL
 fixed_initials <- NULL
 lower <- NULL
 upper <- NULL
 if(spec[[1]]$type %in% c("FOMC", "DFOP", "HS")) {
   mat = FALSE
 } else mat = TRUE
 if(mat){
   if(outpartri=='default') m <- matrix(nrow=ncompart, ncol=ncompart, dimnames=list(obs_vars, obs_vars))
 }
 if(is.null(data)){## in case data is not directly given as a matrix, then the residue component should be there
   time <- spec[[1]]$time
   residue <- time
   for (varname in obs_vars){
     residue <- cbind(residue,spec[[varname]]$residue)
   }
   data0 <- residue
 }else {
   residue <- data
   data0 <- data
 }
 colnames(residue) <- c('time',obs_vars)
 colnames(data0) <- c('time',obs_vars)
 if(is.null(weight)){ ## in case no weight is specified in the arguments.
   ## there should be a weight component or set at default 1.
   for (varname in obs_vars){
     if(!is.null(spec[[varname]]$weight)) {
       tmpweight <- spec[[varname]]$weight
       tmpid <- tmpweight<=0
       ##browser()
       if(sum(tmpid)>0) {
         residue[tmpid,varname] <- NA
         tmpweight[tmpid] <- 1
       }
       
       weightmat <- cbind(weightmat,tmpweight)
     }else {
       ## if weight being null, and data not null
       if(!is.null(data)) weightmat <- cbind(weightmat,rep(1,nrow(data)))
     }
   }
 }else{##XXXXXXXXXXXX TODO XXXXXXXXXXXX
   ## weight specified by the user should either be a function of the observation, or a fixed value computed using the observations.
 }
 ## check for weight being 0: since modCost does not accept 0 weights.
 
 ## ###########
 ## Establish list of differential equations. The compartments/compound has passed the error checks and are completed with paramters and initial values and new_parms information but not new_ff. These functions are the same for both parametrizations, the differences are taken care of in the completeCompound function. ##########
 for (varname in obs_vars)
 {
   
   if(varname==parentname)  {
     if(spec[[varname]]$M0$fixed==1 || is.null(spec[[varname]]$M0)) fixed_initials <- c(fixed_initials,varname)## add fixed initials
     state.ini <- c(state.ini,spec[[varname]]$M0$ini)
     if(!is.null(spec[[varname]]$M0.orig)){
       state.ini.orig <- c(state.ini.orig,spec[[varname]]$M0.orig$ini)
     }else state.ini.orig <- c(state.ini.orig,spec[[varname]]$M0$ini)
     names(state.ini)[length(state.ini)] <- varname
     state.lower <- c(state.lower,spec[[varname]]$M0$lower)
     state.upper <- c(state.upper,spec[[varname]]$M0$upper)
   }else{
     if(spec[[varname]]$type=="DFOP"){
       if(spec[[varname]]$M0_sub1$fixed==1 || is.null(spec[[varname]]$M0_sub1)) fixed_initials <- c(fixed_initials,paste(varname,"sub1",sep="_"))
       if(spec[[varname]]$M0_sub2$fixed==1 || is.null(spec[[varname]]$M0_sub2)) fixed_initials <- c(fixed_initials,paste(varname,"sub2",sep="_"))
    
       state.ini <- c(state.ini,spec[[varname]]$M0_sub1$ini,spec[[varname]]$M0_sub2$ini)
       state.ini.orig <- c(state.ini.orig,spec[[varname]]$M0_sub1$ini,spec[[varname]]$M0_sub2$ini)
       
       names(state.ini)[c(length(state.ini)-1,length(state.ini))] <- paste(varname,c("sub1","sub2"),sep="_")
       state.lower <- c(state.lower,spec[[varname]]$M0_sub1$lower,spec[[varname]]$M0_sub2$lower)
       state.upper <- c(state.upper,spec[[varname]]$M0_sub1$upper,spec[[varname]]$M0_sub2$upper)
     }
     if(spec[[varname]]$type!="DFOP"){
       if(spec[[varname]]$M0$fixed==1 || is.null(spec[[varname]]$M0)) fixed_initials <- c(fixed_initials,varname)## add fixed initials
       state.ini <- c(state.ini,spec[[varname]]$M0$ini)
       if(!is.null(spec[[varname]]$M0.orig)){
         state.ini.orig <- c(state.ini.orig,spec[[varname]]$M0.orig$ini)
       }else state.ini.orig <- c(state.ini.orig,spec[[varname]]$M0$ini)
       names(state.ini)[length(state.ini)] <- varname
       state.lower <- c(state.lower,spec[[varname]]$M0$lower)
       state.upper <- c(state.upper,spec[[varname]]$M0$upper)
     }
   }
   
   # New (sub)compartments (boxes) needed for the model type
   if(varname==parentname){
     new_boxes <- switch(spec[[varname]]$type,
                         SFO = varname,
                         FOMC = varname,
                         DFOP = varname,
                         HS = varname,
                         SFORB = paste(varname, c("free", "bound"), sep="_")
     )
   }else{
     new_boxes <- switch(spec[[varname]]$type,
                         SFO = varname,
                         FOMC = varname,
                         DFOP = paste(varname, c("sub1", "sub2"), sep="_"),
                         HS = varname,
                         SFORB = paste(varname, c("free", "bound"), sep="_")
     )
   }
   map[[varname]] <- new_boxes
   names(map[[varname]]) <- rep(spec[[varname]]$type, length(new_boxes)) ## does not make much sense for DFOP Metabolites.
   ## Start a new differential equation for each new box
   new_diffs <- paste("d_", new_boxes, " =", sep="") ### XXXXXXXXXXX TODOXXXXXXXXXXXX
   if(!is.null(spec[[varname]]$nonlinear_term))new_diffs[[1]] <- paste(new_diffs[[1]], "-", spec[[varname]]$nonlinear_term)
   if(!is.null(spec[[varname]]$sink_term))new_diffs[[1]] <- paste(new_diffs[[1]],spec[[varname]]$sink_term)
   ## In case of SFORB
   if(!is.null(spec[[varname]]$reversible_binding_terms)) new_diffs<- paste(new_diffs,spec[[varname]]$reversible_binding_terms)
   if(!is.null(spec[[varname]]$DFOP_terms)) new_diffs<- paste(new_diffs,"-",spec[[varname]]$DFOP_terms)
   ## Add variables to model
   parms <- c(parms, spec[[varname]]$new_parms)
   parms.ini <- c(parms.ini, spec[[varname]]$new_parms.ini)
   parms.lower <- c(parms.lower, spec[[varname]]$new_parms.lower)
   parms.upper <- c(parms.upper, spec[[varname]]$new_parms.upper)
   fixed_parms <- c(fixed_parms, spec[[varname]]$new_fixed)
   fixed_flag <- c(fixed_flag,spec[[varname]]$new_fixed_flag)
   names(new_diffs) <- new_boxes
   diffs <- c(diffs, new_diffs)
   if(spec[[varname]]$type %in% c("SFO") && mat) {
     if(outpartri=='default') m[varname,varname] <- paste("-k", varname, sep="_")
   }
 }## end for for (varname in obs_vars)
 ## ##################################################
 ## Add a start component for later report usage.
 start <- NULL
 ## first remove the fixed components.
 if(outpartri=='default'){
   start_parms <- setdiff(parms, fixed_parms)
   names(parms.ini) <- parms
   names(parms.lower) <- parms
   names(parms.upper) <- parms
   start_fixed <- rep(0,length(parms))
   names(start_fixed) <- parms
   start_fixed[fixed_parms] <- 1
   start <- data.frame(initial=parms.ini,lower=parms.lower,upper=parms.upper,fixed=start_fixed)
   ##start <-data.frame(initial=parms.ini[start_parms],lower=parms.lower[start_parms],upper=parms.upper[start_parms])
   ##rownames(start) <- start_parms
 }
 modelmess <- NULL
 ## ##################################################
 ## Transfer between compartments ## XXXXXXXXXXX TODO XXXXXXXXXXXX
 for (varname in obs_vars) {
   to <- spec[[varname]]$to
   if(inpartri=='default'){
     if(outpartri=='default'){
       FF <- spec[[varname]]$FF$ini  ## formation fraction for every compartment
       ## Initialize the mat(coefficient matrix/Jacobian matrix) fill in 0's
       if(is.null(to)){## in this case, sink cannot be false.
         for(k in obs_vars){
           if(k!=varname){
             if(mat){if(is.na(m[k,varname])) m[k,varname] <- 0}
           }
         }
         
       }
       if(!is.null(to)) {
         for(k in obs_vars){
           if(!(k %in% to)) {
             if(mat){if(is.na(m[k,varname])) m[k,varname] <- 0}
           }
         }
       }
       f <- ForwardCalcFF(FF) ## the starting values that will be used in the optimization program so that the optimization does not need to deal with the sum ff==1 contraint.###
       names(f) <-to
       origin_box <- switch(spec[[varname]]$type,
                            SFO = varname,
                            FOMC = varname,
                            DFOP = varname,
                            HS = varname,
                            SFORB = paste(varname, "free", sep="_"))
       fraction_left <- NULL
       nto <- length(to)
       
       for (target in to) {
         index <- match(target,to) ### find formation fraction to a
         ###### VERY IMPORTANT HERE, if it is in to and of type DFOP, then not possible to use the other 
         ###### formation of differential equations. so it must be 
         ##corresponding compartment
         if(target==parentname){## In this case the target is a parent compound, 
           ##there is no need to break DFOP into target_sub1 and target_sub2
           target_box <- switch(spec[[target]]$type,
                                SFO = target,
                                FOMC = target,
                                HS = target,
                                DFOP = target,
                                SFORB = paste(target, "free", sep="_"))
           
         }else{
           target_box <- switch(spec[[target]]$type,
                                SFO = target,
                                FOMC = target, ## dangerous
                                HS = target, ## dangerous
                                DFOP = paste(target,c("sub1","sub2"),sep="_"),
                                SFORB = paste(target, "free", sep="_"))
           
         }
         
         ## ##########################################################
         ## add starting values for the optimization with formation fractions
         start <-rbind(start,c(FF[index],spec[[varname]]$FF$lower[index],spec[[varname]]$FF$upper[index],spec[[varname]]$FF$fixed[index]))
         rownames(start)[nrow(start)] <- paste('ff',varname,'to',target,sep='_')
         ## if(spec[[varname]]$FF$fixed[index]==0) {
         ##     start <-rbind(start,c(FF[index],spec[[varname]]$FF$lower[index],spec[[varname]]$FF$upper[index]))
         ##     rownames(start)[nrow(start)] <- paste('ff',varname,target,sep='_')
         ## }else{
         ## }
         ## #########################################
         if(spec[[varname]]$type %in% c("SFORB")) {### should not be here since we don't use this kind of parametrization for SFORB
           stop('SFORB should use water-sediment settings')
         }
         if(spec[[target]]$type %in% c("SFO","FOMC","HS")){
           if(spec[[varname]]$type %in% c("SFO")) {
             fraction_to_target = paste('f',origin_box,'to', target, sep="_")
             fraction_not_to_target = paste("(1 - ", fraction_to_target, ")",
                                            sep="")
             if(is.null(fraction_left)) {
               fraction_really_to_target = fraction_to_target
               fraction_left = fraction_not_to_target
             } else {
               fraction_really_to_target = paste(fraction_left, " * ",
                                                 fraction_to_target, sep="")
               fraction_left = paste(fraction_left, " * ",
                                     fraction_not_to_target, sep="")
             }
             ff[paste(origin_box,'to', target, sep="_")] = fraction_really_to_target
             diffs[[target_box]] <- paste(diffs[[target_box]], "+",
                                          paste("k", origin_box, sep="_"),'*',
                                          ff[paste(origin_box,'to', target,
                                                   sep="_")], "*", origin_box)
             parms <- c(parms, fraction_to_target)
             parms.ini <- c(parms.ini,f[target])
             if(spec[[varname]]$FF$fixed[index]==1) {
               fixed_parms <- c(fixed_parms,fraction_to_target)
               fixed_flag <- c(fixed_flag,'user')
             }
             ## IF NO SINK, then fixed parms should be adding 1, and the last tranformed formation fraction should be fixed at 1!!!!!!!!!!!!
             #browser()
             if(spec[[varname]]$sink==FALSE){
               ## browser()
               ## special care needed, adding a flag for later usage!!!
               if(index==nto){
                 if(spec[[varname]]$FF$fixed[index]==0){
                   ## only when the last formation fraction is not fixed. we
                   ## need to add one fixed parameter, that is changing the
                   ## last transformed formation fraction being fixed at 1.
                   ## Actually it is let the sum of the formation fractions
                   ## in FF being 1, instead of fixing any original formation
                   ## fractions.
                   fixed_parms <- c(fixed_parms,fraction_to_target)
                   fixed_flag <- c(fixed_flag,'KinGUII')
                   modelmess <- c(modelmess,paste("Formation fraction from",varname,"to",target, "is fixed by KinGUII at", noquote("\"1-all formation fractions to other compartments\"")))
                   if(logall) logwarn(paste0("FF fractions fixed by KinGUII. Please check your model set up. SINK is turned off for ", varname, ", but you did not fix the formation fractions to be summed to 1."))
                 }else{
                   warning('You need to switch the order if the formation fraction for the last to compartment is fixed at a certain value but you turn off the sink compartment and you have multiple to compartments!')
                   # fixed_flag[length(fixed_flag)] <- 'KinGUII'
                   if(logall) logwarn("Fixed formation fractions have to be in the beginnig of the to compartments vector if you have multiple compartments!")
                   
                 }
                 parms.ini[length(parms.ini)] <-1
                 
               }
               
             }### end for if(spec[[varname]]$sink==FALSE)
             parms.lower <- c(parms.lower, spec[[varname]]$FF$lower[index])
             parms.upper <- c(parms.upper, spec[[varname]]$FF$upper[index])
             if(mat){
               if(is.na(m[target,varname]))
                 m[target,varname] <- paste(paste("k", origin_box, sep="_"),'*',
                                            fraction_really_to_target)
               else m[target,varname] <- paste(m[target,varname],'+',
                                               paste("k", origin_box, sep="_"),
                                               '*',fraction_really_to_target,
                                               sep='')
             }
             
           } ### end for  if(spec[[varname]]$type %in% c("SFO"))
           if(spec[[varname]]$type %in% c("DFOP","FOMC",  "HS")) {
             fraction_to_target = paste('f',origin_box,'to', target, sep="_")
             fraction_not_to_target = paste("(1 - ", fraction_to_target, ")",
                                            sep="")
             if(is.null(fraction_left)) {
               fraction_really_to_target = fraction_to_target
               fraction_left = fraction_not_to_target
             } else {
               fraction_really_to_target = paste(fraction_left, " * ",
                                                 fraction_to_target, sep="")
               fraction_left = paste(fraction_left, " * ",
                                     fraction_not_to_target, sep="")
             }
             ff[paste(origin_box,'to', target, sep="_")] = fraction_really_to_target
             if(spec[[varname]]$type!="DFOP"){
               diffs[[target_box]] <- paste(diffs[[target_box]], "+",
                                            fraction_really_to_target, "*",
                                            spec[[varname]]$nonlinear_term) 
             }
             if(spec[[varname]]$type=="DFOP"){
               if(varname!= parentname){
                 diffs[[target_box]] <- paste(diffs[[target_box]], "+",
                                              fraction_really_to_target, "*(",
                                              spec[[varname]]$DFOP_terms[1],"+",spec[[varname]]$DFOP_terms[2],")") 					  
               }else{## when the origin is the parent substance
                 diffs[[target_box]] <- paste(diffs[[target_box]], "+",
                                              fraction_really_to_target, "*",
                                              spec[[varname]]$nonlinear_term) 
               }					
               
             }
             parms <- c(parms, fraction_to_target)
             parms.ini <- c(parms.ini, f[index])
             if(spec[[varname]]$FF$fixed[index]==1) {
               fixed_parms <- c(fixed_parms,fraction_to_target)
               fixed_flag <- c(fixed_flag,'user')
             }
             ## IN NO SINK, then fixed parms should be adding 1!!!!!!!!!!!!
             #browser()
             if(spec[[varname]]$sink==FALSE){
               if(index==nto){
                 if(spec[[varname]]$FF$fixed[index]==0){
                   ## only when the last formation fraction is not fixed. we
                   ## need to add one fixed parameter, that is changing the
                   ## last transformed formation fraction being fixed at 1.
                   ## Actually it is let the sum of the formation fractions
                   ## in FF being 1, instead of fixing any original formation
                   ## fractions.
                   fixed_parms <- c(fixed_parms,fraction_to_target)
                   fixed_flag <- c(fixed_flag,'KinGUII')
                 }else{
                   warning('You need to switch the order if the formation fraction for the last to compartment is fixed at a certain value but you turn off the sink compartment and you have multiple to compartments!')
                   #fixed_flag[length(fixed_flag)] <- 'KinGUII'
                 }
                 parms.ini[length(parms.ini)] <-1
               }
               
             }### end for if(spec[[varname]]$sink==FALSE)
             
             parms.lower <- c(parms.lower, spec[[varname]]$FF$lower[index])
             parms.upper <- c(parms.upper, spec[[varname]]$FF$upper[index])
           }## end for  if(spec[[varname]]$type %in%  c("FOMC", "DFOP", "HS"))
           ###################################################################
         }## end for if(spec[[target]]$type %in% c("SFO","FOMC","HS"))
         if(spec[[target]]$type=="DFOP"){
           ## in this case, only one parent is possible, no children metabolites.
           fraction_to_target = paste('f',origin_box,'to', target, sep="_")
           fraction_not_to_target = paste("(1 - ", fraction_to_target, ")",
                                          sep="")
           if(is.null(fraction_left)) {
             fraction_really_to_target = fraction_to_target
             fraction_left = fraction_not_to_target
           } else {
             fraction_really_to_target = paste(fraction_left, " * ",
                                               fraction_to_target, sep="")
             fraction_left = paste(fraction_left, " * ",
                                   fraction_not_to_target, sep="")
           }
           ff[paste(origin_box,'to', target, sep="_")] = fraction_really_to_target
           
           if(spec[[varname]]$type %in% c("SFO")){
             ## in the case where the metabolites backtransformed into a DFOP parent!
             if(target==parentname){## In this case the target is a parent compound, 
               ##there is no need to break DFOP into target_sub1 and target_sub2
               
               diffs[[target_box]] <- paste(diffs[[target_box[1]]], "+",
                                            fraction_really_to_target, "*",
                                            paste("k", origin_box, sep="_"),'*',origin_box,sep="")
             }else{
               diffs[[target_box[1]]] <- paste(diffs[[target_box[1]]], "+",
                                               fraction_really_to_target, "*g_",target,"*",
                                               paste("k", origin_box, sep="_"),'*',origin_box,sep="")
               diffs[[target_box[2]]] <- paste(diffs[[target_box[2]]], "+",
                                               fraction_really_to_target, "*(1-g_",target,")*",
                                               paste("k", origin_box, sep="_"),'*',origin_box,sep="")
               
             }
            
           }
           if(spec[[varname]]$type %in% c("FOMC",  "HS")){
             #browser()
             diffs[[target_box[1]]] <- paste(diffs[[target_box[1]]], "+",
                                             fraction_really_to_target, "*g_",target,"*",
                                             spec[[varname]]$nonlinear_term,sep="")
             diffs[[target_box[2]]] <- paste(diffs[[target_box[2]]], "+",
                                             fraction_really_to_target, "*(1-g_",target,")*",
                                             spec[[varname]]$nonlinear_term,sep="")
           }
           if(spec[[varname]]$type %in% c("DFOP")){
             if(varname==parentname){
               diffs[[target_box[1]]] <- paste(diffs[[target_box[1]]], "+",
                                               fraction_really_to_target, "*g_",target,"*",
                                               spec[[varname]]$nonlinear_term,sep="")
               diffs[[target_box[2]]] <- paste(diffs[[target_box[2]]], "+",
                                               fraction_really_to_target, "*(1-g_",target,")*",
                                               spec[[varname]]$nonlinear_term,sep="")
               
             }else{
               diffs[[target_box[1]]] <- paste(diffs[[target_box[1]]], "+",
                                               fraction_really_to_target, "*g_",target,"*(",
                                               spec[[varname]]$DFOP_terms[1],"+",spec[[varname]]$DFOP_terms[2],")",sep="")
               diffs[[target_box[2]]] <- paste(diffs[[target_box[2]]], "+",
                                               fraction_really_to_target, "*(1-g_",target,")*(",
                                               spec[[varname]]$DFOP_terms[1],"+",spec[[varname]]$DFOP_terms[2],")",sep="")
             }
           }
           
           parms <- c(parms, fraction_to_target)
           parms.ini <- c(parms.ini, f[index])
           if(spec[[varname]]$FF$fixed[index]==1) {
             fixed_parms <- c(fixed_parms,fraction_to_target)
             fixed_flag <- c(fixed_flag,'user')
           }
           ## IN NO SINK, then fixed parms should be adding 1!!!!!!!!!!!!
           #browser()
           if(spec[[varname]]$sink==FALSE){
             if(index==nto){
               if(spec[[varname]]$FF$fixed[index]==0){
                 ## only when the last formation fraction is not fixed. we
                 ## need to add one fixed parameter, that is changing the
                 ## last transformed formation fraction being fixed at 1.
                 ## Actually it is let the sum of the formation fractions
                 ## in FF being 1, instead of fixing any original formation
                 ## fractions.
                 fixed_parms <- c(fixed_parms,fraction_to_target)
                 fixed_flag <- c(fixed_flag,'KinGUII')
               }else{
                 warning('You need to switch the order if the formation fraction for the last to compartment is fixed at a certain value but you turn off the sink compartment and you have multiple to compartments!')
                 fixed_flag[length(fixed_flag)] <- 'KinGUII'
               }
               parms.ini[length(parms.ini)] <-1
             }
             
           }### end for if(spec[[varname]]$sink==FALSE)
           
           parms.lower <- c(parms.lower, spec[[varname]]$FF$lower[index])
           parms.upper <- c(parms.upper, spec[[varname]]$FF$upper[index])
           
           
         }## end for if(spec[[target]]$type=="DFOP")
         
       }## end for for(target in to)
     }## end for  outpartri=='default' ||||||||||||inpartri=='default'
     if(outpartri=='water-sediment'){
       ## #### XXXXXXXXXXXXXXXX TODO XXXXXXXXXXXXXXXXXXXXXXX
       ##  if(outpartri=='water-sediment'){
       ## #### XXXXXXXXXXXXXXXX TODO XXXXXXXXXXXXXXXXXXXXXXX
       ##
       ##
       ## Right now only the SFO case is considered!!!!
       ## NEED TO CONSIDER coefficient matrix if appropriate. XXXX TODO XXXX
       if(!is.null(to)) {
         
         origin_box <- switch(spec[[varname]]$type,
                              SFO = varname,
                              FOMC = varname,
                              DFOP = varname,
                              HS = varname,
                              SFORB = paste(varname, "free", sep="_"))
         fraction_left <- NULL
         if(spec[[varname]]$type %in% c("SFO", "SFORB")) kFF <- spec[[varname]]$kFF$ini
         for (target in to) {
           index <- match(target,to)
           
           target_box <- switch(spec[[target]]$type,
                                SFO = target,
                                SFORB = paste(target, "free", sep="_"))
           ## if(spec[[varname]]$type %in% c("SFO", "SFORB")) {
           if(spec[[varname]]$type  %in% c("SFO", "SFORB")) {
             k_from_to <- paste("k", origin_box, 'to',target_box, sep="_")
             diffs[[origin_box]] <- paste(diffs[[origin_box]], "-",
                                          k_from_to, "*", origin_box)
             diffs[[target_box]] <- paste(diffs[[target_box]], "+",
                                          k_from_to, "*", origin_box)
             parms <- c(parms, k_from_to)
             parms.ini <- c(parms.ini, kFF[index])
             parms.lower <- c(parms.lower, spec[[varname]]$kFF$lower[index])
             parms.upper <- c(parms.upper, spec[[varname]]$kFF$upper[index])
             if(spec[[varname]]$kFF$fixed[index]==1) {
               fixed_parms <- c(fixed_parms,k_from_to)
               fixed_flag <- c(fixed_flag,'user')
             }
           }## end for: if(spec[[varname]]$type %in% c("SFO", "SFORB"))
           
         }
       }
       
     }##end for if(outpartri=='water-sediment')
     ##}
   }## end for inpartri=='default'
   if(inpartri=='water-sediment'){
     if(outpartri=='water-sediment'){
       ## #### XXXXXXXXXXXXXXXX TODO XXXXXXXXXXXXXXXXXXXXXXX
       ##
       ##
       ## Right now only the SFO case is considered!!!!
       ## NEED TO CONSIDER coefficient matrix if appropriate. XXXX TODO XXXX
       if(!is.null(to)) {
         
         origin_box <- switch(spec[[varname]]$type,
                              SFO = varname,
                              FOMC = varname,
                              DFOP = varname,
                              HS = varname,
                              SFORB = paste(varname, "free", sep="_"))
         fraction_left <- NULL
         if(spec[[varname]]$type %in% c("SFO", "SFORB")) kFF <- spec[[varname]]$kFF$ini
         for (target in to) {
           index <- match(target,to)
           target_box <- switch(spec[[target]]$type,
                                SFO = target,
                                SFORB = paste(target, "free", sep="_"))
           ## if(spec[[varname]]$type %in% c("SFO", "SFORB")) {
           if(spec[[varname]]$type  %in% c("SFO", "SFORB")) {
             k_from_to <- paste("k", origin_box, 'to',target_box, sep="_")
             diffs[[origin_box]] <- paste(diffs[[origin_box]], "-",
                                          k_from_to, "*", origin_box)
             diffs[[target_box]] <- paste(diffs[[target_box]], "+",
                                          k_from_to, "*", origin_box)
             parms <- c(parms, k_from_to)
             parms.ini <- c(parms.ini, kFF[index])
             parms.lower <- c(parms.lower, spec[[varname]]$kFF$lower[index])
             parms.upper <- c(parms.upper, spec[[varname]]$kFF$upper[index])
             if(spec[[varname]]$kFF$fixed[index]==1) {
               fixed_parms <- c(fixed_parms,k_from_to)
               fixed_flag <- c(fixed_flag,'user')
             }
           }## end for: if(spec[[varname]]$type %in% c("SFO", "SFORB"))
           
         }
       }
       
     }##end for if(outpartri=='water-sediment')
     if(outpartri=='default'){
       ## the user are not supposed to do it this way.
       warning('Not implemented yet!')
     }##end for if(outpartri=='default')
   }##end for if(inpartri=='water-sediment')
 } ## end ( for (varname in obs_vars) )#end of tranfer between compartments.
 ##### OUTPUT #####
 names(parms.ini) <- parms
 names(parms.lower) <- parms
 names(parms.upper) <- parms
 parms.fixed <- parms.ini[fixed_parms]
 optim_parms <- setdiff(names(parms.ini), fixed_parms)
 parms.optim <- parms.ini[optim_parms]
 parms.lower <- parms.lower[optim_parms]
 parms.upper <- parms.upper[optim_parms]
 names(state.ini.orig) <- names(state.ini)
 names(state.lower) <- names(state.ini)
 names(state.upper) <- names(state.ini)
 state.ini.fixed <- state.ini[fixed_initials]
 optim_initials <- setdiff(names(state.ini), fixed_initials)
 state.ini.optim <- state.ini[optim_initials]
 state.ini.orig.optim <- state.ini.orig[optim_initials]
 state.lower <- state.lower[optim_initials]
 state.upper <- state.upper[optim_initials]
 
 lower <- c(state.lower,parms.lower)
 upper <- c(state.upper,parms.upper)
 if(!is.null(weightmat)) colnames(weightmat) <- obs_vars
 if(!is.null(start)) start$type <-  rep("deparm", nrow(start))
 model <- list(diffs = diffs, parms = parms, map = map,parms.ini=parms.ini,state.ini=state.ini,state.ini.orig=state.ini.orig,lower=lower,upper=upper,fixed_parms=fixed_parms,fixed_flag=fixed_flag,fixed_initials=fixed_initials,residue=as.data.frame(residue),data0=as.data.frame(data0),weightmat=as.data.frame(weightmat),start=start,modelmess=modelmess)
 
 ## Create coefficient matrix if appropriate
 if (mat) {
   if(outpartri=='water-sediment'){
     boxes <- names(diffs)
     n <- length(boxes)
     m <- matrix(nrow=n, ncol=n, dimnames=list(boxes, boxes))
     for (from in boxes) {
       for (to in boxes) {
         if (from == to) {
           k.candidate = paste("k", from, 'to',c(boxes, "sink"), sep="_")
           k.candidate = sub("free.*bound", "free_bound", k.candidate)
           k.candidate = sub("bound.*free", "bound_free", k.candidate)
           k.effective = intersect(model$parms, k.candidate)
           m[from,to] = ifelse(length(k.effective) > 0,
                               paste("-", k.effective, collapse = " "), "0")
         } else {
           k.candidate = paste("k", from, 'to', to, sep="_")
           k.candidate = sub("free.*bound", "free_bound", k.candidate)
           k.candidate = sub("bound.*free", "bound_free", k.candidate)
           k.effective = intersect(model$parms, k.candidate)
           m[to, from] = ifelse(length(k.effective) > 0,
                                k.effective, "0")
         }
       }## end for  for (to in boxes)
     }## end for  for (from in boxes)
   }
   model$coefmat <- m
 }
 if (exists("ff")) model$ff = ff
 class(model) <- "mkinmod.full"
 model$inpartri <- inpartri
 model$outpartri <- outpartri
 model$sinkT <- sinkT
 invisible(model)
}
##' Auxiliary function including error checks for \code{mkinmod.full}
##'
##' @title Auxiliary function including error checks for \code{mkinmod.full}
##' @param compound  A list of properties for a single compound.
##' @param varname The compound name.
##' @param first  hether this compound is the parent compound
##' @param inpartri Input parameterization.
##' @param outpartri Output parameterization.
##' @param data If not NULL, The residue data frame.
##' @param weight If weight is NULL, check weight component.
##' @param update If not NULL, replace the components
##' in compound with the ones in the update.
##' @param ... Other optional arguments. Not used.
##' @return A list of components to be used in \code{\link{mkinmod.full}}. For the
##' differential functions related to other compounds, they are derived in the
##' \code{\link{mkinmod.full}} function.
##' @author Zhenglei Gao
##' @export
completeCompound <- function(compound=list(type='SFO',to='M1'),varname=NULL,first=FALSE,inpartri=c('default','water-sediment','advanced'),outpartri=c('default','water-sediment','advanced'),data=NULL,weight=NULL,autoInit=TRUE,update=NULL,...)
{
  ## ## auxiliary function(error checks) for mkinmod.full
  
  ## ## Arguments
  ## input compound is a list of properties.
  ## input first is whether this compound is the parent compound
  ## input advanced is whether the specification of all parameters initials are required, the time and residue values are needed in the list.
  
  ## if data null, check residue component.
  ## if weight null, check weight component.
  ## if !is.null update, replace the components in compound with the ones in the update.
  
  ## ## Return: A 'completed' compound with a list of properties. for the differential functions related to other compounds, they are derived in the mkinmod.gui function.
  if("logging" %in% loadedNamespaces()){
    logall <- TRUE
  }else{
    logall <- FALSE
  }
  
  if(is.null(varname)) varname <-  deparse(substitute(compound))
  compound$name <- varname
  inpartri <- match.arg(inpartri)
  outpartri <- match.arg(outpartri)
  if(is.null(compound$type)) {
    warning('no reaction type specified! assign to SFO')
    compound$type <-  "SFO"
  }
  if(!compound$type %in% c("SFO", "FOMC", "DFOP","HS", "SFORB")) stop("Available types are SFO, FOMC, DFOP, HS and SFORB only")
  if(is.null(compound$sink)) compound$sink <- TRUE # Turn on sink if not specified otherwise
  if(is.null(compound$to))
  {
    if(compound$sink==FALSE){
      warning('no to compartments, This has to be a metabolite and a SINK at the same time!')
      if(logall==TRUE) logwarn(paste("No sink defined for", varname, ". KinGUII fixed the degradation rate(s) to sink being 0." ))
      #compound$sink <- TRUE
      if(compound$type=="SFO") compound$k <- list(ini=0,fixed=1,lower=0,upper=Inf,flag="KinGUII")
      if(compound$type=="DFOP") {
        compound$k1 <- list(ini=0,fixed=1,lower=0,upper=Inf,flag="KinGUII")
        compound$k2 <- list(ini=0,fixed=1,lower=0,upper=Inf,flag="KinGUII")
        compound$g <- list(ini=1,fixed=1,lower=0,upper=Inf,flag="KinGUII")
        
        }
      if(compound$type=="HS" | compound$type=="FOMC"){
        if(logall==TRUE) logerror(paste("No sink defined for", varname, "And for Metabolite HS and FOMC are not defined." ))
        stop(paste("No sink defined for", varname, "And for Metabolite HS and FOMC are not defined." ))
      }
    }
  }else{### in case there are compartments that the compound transform to
    ##        ## shall we make the FF or kFF's into the new_parms component of the complete compound property list? XXXXXXXXXXXXXXXX TODO XXXXXXXXXXXXXXXXXX
    n <- length(compound$to)
    if(inpartri=='default'){
      if(outpartri=='default'){
        if(is.null(compound$FF)){
          compound$FF=list(ini=rep(0.1,n),fixed=rep(0,n),lower=rep(0,n),upper=rep(1,n))
          ## if sink==FALSE needs special care, which is done in the mkinmod.full function
        }else{
          ## if there are fixed formation fractions has to be in the first few
          ffind1 <- which(compound$FF$fixed==1)
          if(length(ffind1)>=1){
            if(n>1){
              ffind2 <- which(compound$FF$fixed==0)
              ffind <- c(ffind1,ffind2)
              compound$to <- compound$to[ffind]
              compound$FF=list(ini=compound$FF$ini[ffind],fixed=compound$FF$fixed[ffind],
                               lower=compound$FF$lower[ffind],upper=compound$FF$upper[ffind])
            }
          }
          
          
        }
      }
      if(outpartri=='water-sediment'){
        if(is.null(compound$FF)){
          compound$FF=list(ini=rep(0.1,n),fixed=rep(0,n),lower=rep(0,n),upper=rep(1,n))
          ## if sink==FALSE needs special care, which is done in the mkinmod.full function
        }
        ### we need 'compound$k'
        if(is.null(compound$k)){
          compound$k <- list(ini= 0.1,fixed = 0,lower = 0.0,upper = Inf)
        }
        if(sum(compound$FF$ini)>1) stop("Formation Fractions cannot add up to over 1")
        compound$kFF=list(ini=compound$FF$ini*compound$k$ini,fixed=compound$FF$fixed*compound$k$fixed,lower=rep(0,n),upper=rep(Inf,n))
        compound$k_compound_sink <- list(ini=(1-sum(compound$FF$ini))*compound$k$ini,fixed=prod(compound$FF$fixed),lower=0,upper=Inf)
      }
    }
    if(inpartri=='water-sediment'){
      ## there should be multiple k's for SFO called kFF since there are multiple to compartments
      if(compound$type=='SFO'){
        if(is.null(compound$kFF)){
          compound$kFF=list(ini=rep(0.1,n),fixed=rep(0,n),lower=rep(0,n),upper=rep(Inf,n))
        }
        
        
      }## end if(compound$type=='SFO')
    }
  }
  ## for different types and parametrization, check the parameter settings, if missing initials or missing lower uppers, fixed, assign default values.##
  ## For initial concentration, the same for all types and parametrizations.
  if(first==TRUE){
    if(is.null(compound$M0)) {## using default. common for all types parametrizations and kinetic models
      compound$M0 <- list(ini = 100,fixed = 0,lower = 0.0,upper = Inf)
      if(!is.null(compound$residue)){
        ## since it is first, there must be a time component!
        if(autoInit == TRUE) compound$M0$ini <- mean(compound$residue[compound$time==min(compound$time,na.rm=TRUE)],na.rm=TRUE)     ## here we can change to the fitted value instead if possible.     
        ## note that mean(numeric(0))=NaN
      }else{
        if(!is.null(data)){
          compound$M0$ini <- mean(data[data[,1]==min(data[,1],na.rm=T),2],na.rm=TRUE)     ## here we can change to the fitted value instead if possible.     
          if(logall) loginfo("Initial values for M0 changed to the average concentrations at time 0")
        }else{
          if(autoInit == TRUE) compound$M0 <- list(ini = 100,fixed = 0,lower = 0.0,upper = Inf) ## assign 0 for metabolites, 100 for parents.
          
        }
      }
      
      
    }else{
      ## check if reasonable!!
      if(!is.null(compound$residue)){
        ## tmp <- na.omit(compound$residue)
        if(compound$M0$fixed==0) {
          ##browser()
          compound$M0.orig <- compound$M0
          ##compound$M0$ini <- tmp[1]     ## here we can change to the fitted value instead if possible.     
          if(autoInit == TRUE) compound$M0$ini <- mean(compound$residue[compound$time==min(compound$time,na.rm=TRUE)],na.rm=TRUE)
          if(logall) loginfo("Initial values for M0 changed to the average of time 0")
        }
      }else{
        ## compound$M0 is set up by the user! and there is data available
        if(!is.null(data)){
          if(compound$M0$fixed==0) {
            compound$M0.orig <- compound$M0
            if(autoInit == TRUE) compound$M0$ini <- mean(data[data[,1]==0,2],na.rm=TRUE)     ## here we can change to the fitted value instead if possible.     
            
          }
        }
      }
      ### give a first dirty fit to the first compartment.
      ## should be done later.
    }
    if(is.na(compound$M0$ini)){
      compound$M0$ini <- 100
    }
  }else{
    if(compound$type=="SFO") {
      if(is.null(compound$M0)) compound$M0 <- list(ini = 0,fixed = 1,lower = 0.0,upper = Inf)
    }
    if(compound$type=="DFOP") {
      if(is.null(compound$M0_sub1))  compound$M0_sub1 <- list(ini = 0,fixed = 1,lower = 0.0,upper = Inf)
      if(is.null(compound$M0_sub2))  compound$M0_sub2 <- list(ini = 0,fixed = 1,lower = 0.0,upper = Inf)
    }
    
  }
  ## ###
  ## New (sub)compartments (boxes) needed for the model type. Same for both parametrization.
  if(first==TRUE){
    new_boxes <- switch(compound$type,
                        SFO = varname,
                        FOMC = varname,
                        DFOP = varname,
                        HS = varname,
                        SFORB = paste(varname, c("free", "bound"), sep="_")
    )}else{
      new_boxes <- switch(compound$type,
                          SFO = varname,
                          FOMC = varname,
                          DFOP = paste(varname,c("sub1","sub2"),sep="_"),
                          HS = varname,
                          SFORB = paste(varname, c("free", "bound"), sep="_")
      )
    }
  
  ## input parametrization ## ---->>>>>>>>> ## output parametrization ##
  ##
  if(inpartri=='default'){
    if(outpartri=='default'){
      if(compound$type=='FOMC'){
        ## Construct and add FOMC term and add FOMC parameters if needed
        ## From p. 53 of the FOCUS kinetics report
        alphaname<-paste("alpha", new_boxes[[1]],  sep="_")
        betaname<-paste("beta", new_boxes[[1]],  sep="_")
        nonlinear_term <- paste("(",alphaname,"/",betaname,") * ((time/",betaname,") + 1)^-1 *", new_boxes[[1]])
        compound$nonlinear_term <- nonlinear_term
        
        new_parms <- c(alphaname, betaname)
        if(is.null(compound$alpha)){
          compound$alpha <- list(ini= 0.1,fixed = 0,lower = 0.0,upper = Inf)
        }
        if(is.null(compound$beta)){
          compound$beta <- list(ini= 0.1,fixed = 0,lower = 0.0,upper = Inf)
        }
        new_parms.ini <- c(compound$alpha$ini,compound$beta$ini)
        new_parms.lower <- c(compound$alpha$lower,compound$beta$lower)
        new_parms.upper <- c(compound$alpha$upper,compound$beta$upper)
        new_fixed <- NULL
        new_fixed_flag <- NULL
        if (compound$alpha$fixed==1) {
          new_fixed <- c(new_fixed,alphaname)
          new_fixed_flag <- c(new_fixed_flag ,"user")
        }
        if (compound$beta$fixed==1) {
          new_fixed <- c(new_fixed,betaname)
          new_fixed_flag <-c(new_fixed_flag,"user")
        }
        ## # finish FOMC parameters without the formation fraction part. # ##
        ## # FOMC cannot using the parametrization without formation fraction # ##
      }
      if(compound$type=='DFOP'){
        if(first==FALSE) {
          ## print('Warning! using DFOP for metabolites are not fully tested!Use at your own risk!')
          ## warning(' using DFOP for metabolites are not fully tested!Use at your own risk!')
          k1name <-paste("k1", varname,  sep="_")
          k2name <- paste("k2", varname,  sep="_")
          gname <- paste("g", varname,  sep="_")
          ## nonlinear_term <- c(paste("-",k1name,"*",new_boxes[[1]],sep=""),paste("-",k2name,"*", new_boxes[[2]],sep=""))
          DFOP_terms <- c(paste(k1name,"*",new_boxes[[1]],sep=""),paste(k2name,"*", new_boxes[[2]],sep=""))
          ## Needs special care for the g parameter.
          compound$nonlinear_term <- NULL
          compound$DFOP_terms <- DFOP_terms
          
          new_parms <- c(k1name,k2name,gname)
          ## if nitial parameters are not specified##
          if(is.null(compound$k1)){
            compound$k1 <- list(ini= 1,fixed = 0,lower = 0.0,upper = Inf)
          }
          if(is.null(compound$k2)){
            compound$k2 <- list(ini= 0.1,fixed = 0,lower = 0.0,upper = Inf)
          }
          if(is.null(compound$g)){
            compound$g <- list(ini= 0.1,fixed = 0,lower = 0.0,upper = 1)
          }
          ####
          new_parms.ini <- c(compound$k1$ini,compound$k2$ini,compound$g$ini)
          new_parms.lower <- c(compound$k1$lower,compound$k2$lower,compound$g$lower)
          new_parms.upper <- c(compound$k1$upper,compound$k2$upper,compound$g$upper)
          new_fixed <- NULL
          new_fixed_flag <- NULL
          if (compound$k1$fixed==1) {
            new_fixed <- c(new_fixed,k1name)
            new_fixed_flag <- c(new_fixed_flag,ifelse((!is.null(compound$k1$flag)) && compound$k1$flag=="KinGUII","KinGUII","user"))
            
          }
          if (compound$k2$fixed==1) {
            new_fixed <- c(new_fixed,k2name)
            new_fixed_flag <- c(new_fixed_flag,ifelse((!is.null(compound$k2$flag)) && compound$k2$flag=="KinGUII","KinGUII","user"))
            
          }
          if (compound$g$fixed==1) {
            new_fixed <- c(new_fixed,gname)
            new_fixed_flag <- c(new_fixed_flag,ifelse((!is.null(compound$g$flag)) && compound$g$flag=="KinGUII","KinGUII","user"))
            
          }
          
        }else{
          k1name <-paste("k1", new_boxes[[1]],  sep="_")
          k2name <- paste("k2", new_boxes[[1]],  sep="_")
          gname <- paste("g", new_boxes[[1]],  sep="_")
          nonlinear_term <- paste("((",k1name,"*", gname ,"*", "exp(-",k1name, "* time) + ",k2name, "* (1 - ",gname,") * exp(-",k2name, "* time)) / (",gname, "* exp(-",k1name,"* time) + (1 - ",gname,") * exp(-",k2name, "* time))) *", new_boxes[[1]])
          compound$nonlinear_term <- nonlinear_term
          
          new_parms <- c(k1name,k2name,gname)
          ## if nitial parameters are not specified##
          if(is.null(compound$k1)){
            compound$k1 <- list(ini= 1,fixed = 0,lower = 0.0,upper = Inf)
          }
          if(is.null(compound$k2)){
            compound$k2 <- list(ini= 0.1,fixed = 0,lower = 0.0,upper = Inf)
          }
          if(is.null(compound$g)){
            compound$g <- list(ini= 0.1,fixed = 0,lower = 0.0,upper = 1)
          }
          ####
          new_parms.ini <- c(compound$k1$ini,compound$k2$ini,compound$g$ini)
          new_parms.lower <- c(compound$k1$lower,compound$k2$lower,compound$g$lower)
          new_parms.upper <- c(compound$k1$upper,compound$k2$upper,compound$g$upper)
          new_fixed <- NULL
          new_fixed_flag <- NULL
          if (compound$k1$fixed==1) {
            new_fixed <- c(new_fixed,k1name)
            new_fixed_flag <- c(new_fixed_flag,ifelse((!is.null(compound$k1$flag)) && compound$k1$flag=="KinGUII","KinGUII","user"))
            
          }
          if (compound$k2$fixed==1) {
            new_fixed <- c(new_fixed,k2name)
            new_fixed_flag <- c(new_fixed_flag,ifelse((!is.null(compound$k2$flag)) && compound$k2$flag=="KinGUII","KinGUII","user"))
            
          }
          if (compound$g$fixed==1) {
            new_fixed <- c(new_fixed,gname)
            new_fixed_flag <- c(new_fixed_flag,ifelse((!is.null(compound$g$flag)) && compound$g$flag=="KinGUII","KinGUII","user"))
            
          }
        }
      }## end of DFOP
      if(compound$type=='HS'){
        k1name<-paste("k1", new_boxes[[1]],  sep="_")
        k2name<-paste("k2", new_boxes[[1]],  sep="_")
        tbname<-paste("tb", new_boxes[[1]],  sep="_")
        nonlinear_term <- paste("ifelse(time <=",tbname,",",k1name,",", k2name,")", "*", new_boxes[[1]])
        compound$nonlinear_term <- nonlinear_term
        ## #
        ## if nitial parameters are not specified##
        if(is.null(compound$k1)){
          compound$k1 <- list(ini= 1,fixed = 0,lower = 0.0,upper = Inf)
        }
        if(is.null(compound$k2)){
          compound$k2 <- list(ini= 0.1,fixed = 0,lower = 0.0,upper = Inf)
        }
        if(is.null(compound$tb)){
          compound$tb <- list(ini= 0.1,fixed = 0,lower = 0.0,upper = Inf)
        }
        ## Dealing with the parameters for HS###########
        new_parms <- c(k1name, k2name, tbname)
        new_parms.ini <- c(compound$k1$ini,compound$k2$ini,compound$tb$ini)
        new_parms.lower <- c(compound$k1$lower,compound$k2$lower,compound$tb$lower)
        new_parms.upper <- c(compound$k1$upper,compound$k2$upper,compound$tb$upper)
        new_fixed <- NULL
        new_fixed_flag <- NULL
        if(compound$k1$fixed==1) {
          new_fixed <- c(new_fixed,k1name)
          new_fixed_flag <- c(new_fixed_flag,"user")
          
        }
        if(compound$k2$fixed==1) {
          new_fixed <- c(new_fixed,k2name)
          new_fixed_flag <- c(new_fixed_flag,"user")
        }
        if(compound$tb$fixed==1) {
          new_fixed <- c(new_fixed,tbname)
          new_fixed_flag <- c(new_fixed_flag,"user")
        }
      }
      if(compound$type=='SFORB'){
        stop('For SFORB model, it is not allowed to use a parameterization with formation fractions! SFORB is only implemented but not fully tested.')
      }
      if(compound$type=='SFO'){##
        k_compound_sink <- paste("k", new_boxes[[1]],  sep="_")
        sink_term <- paste("-", k_compound_sink, "*", new_boxes[[1]])
        if(is.null(compound$k)){
          compound$k <- list(ini= 0.1,fixed = 0,lower = 0.0,upper = Inf)
        }
        new_parms <- k_compound_sink
        new_parms.ini <-compound$k$ini
        new_parms.lower <- compound$k$lower
        new_parms.upper <- compound$k$upper
        new_fixed <- NULL
        new_fixed_flag <- NULL
        if(compound$k$fixed==1) {
          new_fixed<- c(new_fixed,k_compound_sink)
          new_fixed_flag <- c(new_fixed_flag,ifelse((!is.null(compound$k$flag)) && compound$k$flag=="KinGUII","KinGUII","user"))
          
        }
        compound$sink_term <- sink_term
      }
      ## common for all types of kinect models
      compound$new_parms <- new_parms
      compound$new_parms.ini <- new_parms.ini
      compound$new_parms.lower <- new_parms.lower
      compound$new_parms.upper <- new_parms.upper
      compound$new_fixed <- new_fixed
      compound$new_fixed_flag <- new_fixed_flag
      
    }## end for if(outpartri=='default') conditional on inputpartri=='default'
    if(outpartri=='water-sediment'){
      ## now only deal with SFO
      new_parms <-NULL
      new_parms.ini <- NULL
      new_parms.lower <- NULL
      new_parms.upper <- NULL
      new_fixed <- NULL
      new_fixed_flag <- NULL
      if(compound$sink==FALSE){
        if(is.null(compound$k_compoud_sink)){
          compound$k_compound_sink <-list(ini=0,fixed=1,lower=0,upper=Inf,flag="KinGUII")
        }else{
          if(compound$k_compound_sink$ini!=0 | compound$k_compound_sink$fixed!=1){
            compound$k_compound_sink <-list(ini=0,fixed=1,lower=0,upper=Inf,flag="KinGUII")
          }
        }
        
      }
      if(compound$sink==TRUE){
        if(is.null(compound$k_compoud_sink)){
          ## even if sink==FALSE no hurt for adding this
          compound$k_compound_sink <-list(ini=0.1,fixed=0,lower=0,upper=Inf)
          
        }
        k_compound_sink <- paste("k", new_boxes[[1]],'to', "sink", sep="_")
        sink_term <- paste("-", k_compound_sink, "*", new_boxes[[1]])
        compound$sink_term <- sink_term
        new_parms <-k_compound_sink
        new_parms.ini <-compound$k_compound_sink$ini
        new_parms.lower <-compound$k_compound_sink$lower
        new_parms.upper <-compound$k_compound_sink$upper
        if(compound$k_compound_sink$fixed==1) {
          new_fixed<- c(new_fixed,k_compound_sink)
          new_fixed_flag <- c(new_fixed_flag,ifelse((!is.null(compound$k_compound_sink$flag)) && compound$k_compound_sink$flag=="KinGUII","KinGUII","user"))
        }
      }###
      if(compound$type=='FOMC'){
        stop('For FOMC, it is not allowed to use a parameterization without formation fractions!')
      }
      if(compound$type=='DFOP'){## TODO
        stop('Not Implemented YET!! Please wait until the next update!')
      }
      if(compound$type=='HS'){## TODO
        stop('Not Implemented YET!! Please wait until the next update!')
      }
      if(compound$type=='SFO'){## xxxxxx  TODO xxxxxxxxxxxx
        ##############################################
        
        ##############################################
      }
      if(compound$type=='SFORB'){
        stop('For SFORB model, it is not allowed to use a parameterization with formation fractions! SFORB is only implemented but not fully tested.')
      }
      ## common for all types of kinect models
      compound$new_parms <- new_parms
      compound$new_parms.ini <- new_parms.ini
      compound$new_parms.lower <- new_parms.lower
      compound$new_parms.upper <- new_parms.upper
      compound$new_fixed <- new_fixed
      compound$new_fixed_flag <- new_fixed_flag
    }## end for if(outpartri=='water-sediment')
    
  }## end for if(inpartri=='default')  input parameter being the default fomat with FF inputs
  if(inpartri=='water-sediment'){
    if(outpartri=='water-sediment'){###XXX TODO
      ### **************need to add sink_term in any ways when thinking about transformation between compartments.*****************
      new_parms <-NULL
      new_parms.ini <- NULL
      new_parms.lower <- NULL
      new_parms.upper <- NULL
      new_fixed <- NULL
      new_fixed_flag <- NULL
      if(compound$sink==TRUE){
        if(is.null(compound$k_compoud_sink)){
          ## even if sink==FALSE no hurt for adding this
          compound$k_compound_sink <-list(ini=0.1,fixed=0,lower=0,upper=Inf)
          
        }
        k_compound_sink <- paste("k", new_boxes[[1]],'to', "sink", sep="_")
        sink_term <- paste("-", k_compound_sink, "*", new_boxes[[1]])
        compound$sink_term <- sink_term
        new_parms <-k_compound_sink
        new_parms.ini <-compound$k_compound_sink$ini
        new_parms.lower <-compound$k_compound_sink$lower
        new_parms.upper <-compound$k_compound_sink$upper
        if(compound$k_compound_sink$fixed==1) {
          new_fixed<- c(new_fixed,k_compound_sink)
          new_fixed_flag <- c(new_fixed_flag,ifelse((!is.null(compound$k_compound_sink$flag)) && compound$k_compound_sink$flag=="KinGUII","KinGUII","user"))         
        }
      }###
      if(compound$type=='SFORB'){
        k_free_bound <- paste("k", varname, "free", "bound", sep="_")
        k_bound_free <- paste("k", varname, "bound", "free", sep="_")
        if(is.null(compound$k_free_bound)){
          compound$k_free_bound <- list(ini=0.1,fixed=0,lower=0,upper=Inf)
        }
        if(is.null(compound$k_bound_free)){
          compound$k_bound_free <- list(ini=0.1,fixed=0,lower=0,upper=Inf)
        }
        reversible_binding_terms <- c(paste("-", k_free_bound, "*", new_boxes[[1]],
                                            "+", k_bound_free, "*", new_boxes[[2]]),
                                      paste("+", k_free_bound, "*", new_boxes[[1]],
                                            "-", k_bound_free, "*", new_boxes[[2]]))
        new_parms <- c(new_parms, k_free_bound, k_bound_free)
        new_parms.ini <-c(new_parms.ini, compound$k_free_bound$ini,compound$k_bound_free$ini )
        new_parms.lower <-c(new_parms.lower, compound$k_free_bound$lower,compound$k_bound_free$lower )
        new_parms.upper <-c(new_parms.upper, compound$k_free_bound$upper,compound$k_bound_free$upper )
        if(compound$k_free_bound$fixed==1) {
          new_fixed<- c(new_fixed,k_free_bound)
          new_fixed_flag <- c(new_fixed_flag,ifelse((!is.null(compound$k_free_bound$flag)) && compound$k_free_bound$flag=="KinGUII","KinGUII","user"))
          
        }
        if(compound$k_bound_free$fixed==1) {
          new_fixed<- c(new_fixed,k_bound_free)
          new_fixed_flag <- c(new_fixed_flag,ifelse((!is.null(compound$k_bound_free$flag)) && compound$k_bound_free$flag=="KinGUII","KinGUII","user"))
          
        }
        
      }##
      ## if(compound$type=='SFO'){
      ##     if(is.null(compound$kFF)){
      ##         compound$kFF <- list(ini= rep(0.1,nto),fixed=rep(0,nto), lower=rep(0,nto), upper=rep(Inf,nto))
      
      ##     }
      ## }
      if(compound$type=='FOMC'){
        stop('For FOMC, it is not allowed to use a parameterization without formation fractions!')
      }
      if(compound$type=='DFOP'){## TODO
        stop('Not Implemented YET!! Please wait until the next update!')
      }
      if(compound$type=='HS'){## TODO
        stop('Not Implemented YET!! Please wait until the next update!')
      }
      
      ## common for all types of kinect models
      compound$new_parms <- new_parms
      compound$new_parms.ini <- new_parms.ini
      compound$new_parms.lower <- new_parms.lower
      compound$new_parms.upper <- new_parms.upper
      compound$new_fixed <- new_fixed
      compound$new_fixed_fixed_flag <- new_fixed_flag
    }## end for if(outpartri=='water-sediment'){&& if(inpartri=='water-sediment')}
    if(outpartri=='default'){
      stop('Not Implemented YET!! Please wait until the next update!')
    }## end for if(outpartri=='default'){&& if(inpartri=='water-sediment')}
  }## end for if(inpartri=='water-sediment')
  
  if(!is.null(data)){## data is provided as a matrix/data.frame read from an external txt file.
    ## No need to consider parent compound has to have a time component, or the metabolite compound has to have a residue component.
    
  }else{## data is not provided, then the compound should have data assigned to it.
    
    if(first==TRUE){
      ## There should be a residue component and a time component.
      if(is.null(compound$time)) stop('There should be a time component for a parent compound')
    }
    if(is.null(compound$residue)){
      warning('There should be a residue component for a parent compound, setting them to be NA, a Ghost compartment')
    }
    #####
    
    
  }
  
  return(compound)
}

ForwardCalcFF <- function(trueff)
{
  ### From the formation faction provided by the user(true formation fractions), calculte the f used in the program
  l <- length(trueff)
  ff <- trueff
  for(i in 1:l)
  {
    if(i>1)
    {
      for(j in 1:(i-1))
        ff[i] <- ff[i]/(1-ff[j])
    }#else ff[i] <- trueff[i]
  }
  ff
}
fitParent <- function(parent){
  ## parent of a list with 0 and other stuff.
  ## A quick and dirty fit to get better initial values for parent compartment.
}