##'@export
KinReport <- function(mod,Fit,version = "2.2013.0923.1534",filename,exactHessian=FALSE){
  if(missing(filename)) filename <- deparse(substitute(mod))
  ## browser()
  filename1 <- paste0(filename,".kgo")
  filename2 <- paste0(filename,".kgg")
  
  if(class(mod)=="try-error"){
    ## when cannot set up the kinetic Model
    ## empty report and graph data are generated.
    logerror("Cannot set up the kinetic model, please check your input data!")
    list2ascii(emptyreport(version=version),filename1)
    emptygraph(filename2)
    loginfo("Empty report files are created.")
  }else{
    #
    if(class(Fit)[1]=="try-error"){
      ## Fit cannot be ran successfully.
      logerror("The kinetic model cannot be fit properly. Please check your data and your model!")
      ## write summary using the mod object instead.
      writesummary <- try(list2ascii(summary(mod, version = version),filename1),silent=TRUE)
      if(class(writesummary)=="try-error"){
        logerror("The kinetic model cannot be summarized. Please check your data and your model!")
        ## finally using an empty report
        list2ascii(emptyreport(version=version),filename1)
        loginfo("Empty .kgo report files are created.")
      }else{loginfo("Summary .kgo file from the model is created.")}
      ## write graph data using the mod object instead.
      writegraph <- try(datagraph(mod,filename2),silent=TRUE)
      if(class(writegraph)=="try-error"){
        logerror("The data cannot be extracted from the model. Please check your data and your model!")
        ## finally using an empty report
        emptygraph(filename2)
        loginfo("Empty .kgg report files are created.")
      }else{loginfo("Graph data (.kgg file) is generated using the kinetic model.")}
      
    }else{
      loginfo("Fit is completed successfully")
      if(exactHessian){
        if(class(Fit$numCov)!="try-error"){
          loginfo("Use the exact hessian calculated using numDeriv instead of the approximated one from the optimization!")
          Fit$covar <- Fit$numCov
        }
      }
      writesummary <- try(list2ascii(summary(Fit, cov = TRUE,version = version),filename1),silent=TRUE)
      if(class(writesummary)=="try-error"){
        logerror("The fitted information for the kinetic model cannot be summarized.")
        
        ## write summary using the mod object instead.
        writesummay <- try(list2ascii(summary(mod, version = version),filename1),silent=TRUE)
        if(class(writesummary)=="try-error"){
          ## finally using an empty report
          logerror("The kinetic model cannot be summarized. Please check your data and your model!")
          list2ascii(emptyreport(version=version),filename)
          loginfo("Empty .kgo report files are created.")
        }else{
          loginfo("Summary .kgo file is generated using the kinetic model.")
        }
      }else{
        loginfo("Summary of the fitted object (.kgo file) is generated.")
      }
      writegraph <- try(kingraph(Fit,filename2),silent=TRUE)
      
      if(class(writegraph)=="try-error"){
        ## write graph data using the mod object instead.
        logerror("The graph data of the fitted object cannot be generated.")
        writegraph <- try(datagraph(mod,filename2),silent=TRUE)
        if(class(writegraph)=="try-error"){
          logerror("The data cannot be extracted from the model. Please check your data and your model!")
          ## finally using an empty graph
          emptygraph(filename2)
          loginfo("Empty .kgg report files are created.")
        }else{
          loginfo("Graph data (.kgg file) is generated using the kinetic model.")
        }
      }else{
        loginfo("Graph data of the fitted object (.kgg file) is generated.")
        ## Generating the combined plot!
        win.metafile(paste0(filename,".fitplot.wmf"),pointsize=12)
        par(lwd=1)
        kinplot(Fit)
        dev.off()    
        
      }
      if("mcmckingui" %in% class(Fit)){
        ## MCMC results need extra Fit plot!
        plot(Fit,y=NULL,
             paste0(filename,'.Density'),
             paste0(filename,'.Correlation'),
             paste0(filename,'.Trace'),
             device='wmf')
      }
    }
  }
}
##'@export
emptyreport <- function(version="2.2013.0923.1534"){
  ## Produce empty report
  ## used when the mkinmod.full does not produce anything.
  ## example: empres <- emptyreport()
  res <- NULL
  res$version <- version
  res$mess <- "Wrong model setup"
  class(res) <- "emptyReport"
  return(res)
  
}
##'@export
emptygraph <- function(filename="mod.kgg"){
  out <- data.frame("time"=NA,"Parent"=NA)
  write.table(out,filename,sep="\t")
}
##'@export
datagraph <- function(mod,filename="mod.kgg"){
  write.table(mod$data,filename,sep="\t")
}
print.emptyReport <- function(x){
  cat(paste("KineticeEval Package Version:",packageVersion("KineticEval"),"\n"))
  cat(paste('Version:',x$version,'\n'))
  ##cat('\nR version: 2.12.2 (2011-02-25)\n ')
  cat(paste('\n',sessionInfo()$R.version$version.string,'\n',sep=""))
  if(!is.null(x$outpartri)){
    if(x$outpartri=='water-sediment')  cat("\nStudy:Water-Sediment\n")
  }
  cat("\nOptimization Algorithms Used:\n")
  print(ifelse(is.null(x$optimmethod),"",(x$optimmethod)))
  cat("\n##Comment:\n")
  print(ifelse(is.null(x$mess),"",(x$mess)))
  cat("\nEquations:\n")
  
  cat("\nStarting values for optimised parameters:\n")
  
  cat("\nFixed parameter values:\n")
  
  
  cat("\nOptimised parameters:\n")
  cat("\nEstimated formation fractions:\n")
  printff <- !is.null(x$ff)
  if(printff){
    print(data.frame(ff = x$ff), digits=digits,...)
  }
  
  
  
  
  cat("\nChi2 error levels in percent :\n")
  cat("\nEstimated disappearance times:\n")
  printdistimes <- !is.null(x$distimes)
  if(printdistimes){
    
    print(x$distimes, digits=digits,...)
  }
  
  
  
  cat("\nAdditional Statistics:\n")
  
  printcor <- (!is.null(x$cov.unscaled))
  if (printcor){
    
    np <- nrow(x$fixed)
    Corr <- matrix(data = NA, nrow = np, ncol = np)
    
    rownames(Corr) <- colnames(Corr) <- rownames(x$fixed)#rownames(x$par)
    cat("\nParameter correlation:\n")
    print(Corr, digits = digits, ...)
  }
  
  
  printdata <- !is.null(x$data)
  if (printdata){
    cat("\nData:\n")
    print(myformat(x$data, digits = digits, scientific=F,...), row.names = FALSE)
    #sprintf("%.4f",x$data)
  }
  
  invisible(x)
}

##' Print method for \code{\link{summary.mkinmod.full}}
##'
##' Expanded from \code{\link{print.summary.modFit}} and \code{\link{print.summary.mkinfit}}
##' @title Print method for \code{\link{summary.kingui}}
##' @method print summary.kingui
##' @param x An object of class \code{summary.kingui}
##' @param digits How many digits should be printed after the decimal point.
##' @param detailed Ignore for now.
##' @param ... Other parameters to be passed into \code{\link{summary.kingui}}
##' @return \code{NULL}
##' @author Zhenglei Gao
##' @S3method print summary.mkinmod.full
##' @rdname print.summary.mkinmod.full
print.summary.mkinmod.full <- function(x, digits = max(3, getOption("digits") - 3),detailed=FALSE, ...) {
  cat(paste("KineticeEval Package Version:",packageVersion("KineticEval"),"\n"))
  cat(paste('Version:',x$version,'\n'))
  ##cat('\nR version: 2.12.2 (2011-02-25)\n ')
  cat(paste('\n',sessionInfo()$R.version$version.string,'\n',sep=""))
  if(!is.null(x$outpartri)){
    if(x$outpartri=='water-sediment')  cat("\nStudy:Water-Sediment\n")
  }
  xx <- x[["diffs"]]
  cat("\nOptimization Algorithms Used:\n")
  print(ifelse(is.null(x$optimmethod),"",(x$optimmethod)))
  cat("\n##Comment:\n")
  print(ifelse(is.null(x$mess),"",(x$mess)))
  cat("\nEquations:\n")
  for(i in 1:length(xx)) print(noquote(as.character(xx[i])))
  df  <- x$df
  rdf <- df[2]
  
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
  ifelse(is.null(x$par),"",printCoefmat(x$par, digits = digits, ...))
  if(is.null(x$sigma)){
    
  }else{
    cat("\nResidual standard error:",
        format(signif(x$sigma, digits)), "on", rdf, "degrees of freedom\n")
  }
  
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
  
  
  
  cat("\nAdditional Statistics:\n")
  print(x$errmin[,4:ncol(x$errmin)], digits=digits,...)
  
  printcor <- (!is.null(x$cov.unscaled))
  if (printcor){
   
    np <- nrow(x$fixed)
    Corr <- matrix(data = NA, nrow = np, ncol = np)
    
    rownames(Corr) <- colnames(Corr) <- rownames(x$fixed)#rownames(x$par)
    cat("\nParameter correlation:\n")
    print(Corr, digits = digits, ...)
  }
  
  
  printdata <- !is.null(x$data)
  if (printdata){
    cat("\nData:\n")
    print(myformat(x$data, digits = digits, scientific=F,...), row.names = FALSE)
    #sprintf("%.4f",x$data)
  }
  
  invisible(x)
}