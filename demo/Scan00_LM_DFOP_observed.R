# Trial           : Scan00
## To compare the results using calcJacob=TRUE and not.

# Check inst. versions of 'KineticEval' and 'R'
if(!require(KineticEval)) stop("Package KineticEval is not available")
if( as.numeric(version$major)<2||(as.numeric(version$major)==2 & as.numeric(version$minor)<15.1) ) stop("Please install R version > 2.15.0")


# Init 'R'
data(versioninfo)
require(methods)
require(logging)


# Create log file
basicConfig(level='FINEST')
addHandler(writeToFile, file="Scan00 LM DFOP observed.log", level='DEBUG')


#
# Residue data, pathway and kinetics
#

 Scan00 <- try(mkinmod.full(
   observed = list(
       time = c(    0,     0,     3,     3,     7,     7,    14,    14,    30,    30,    50,    50,    70,    70,   101,   101),
    residue = c(100.5,  99.6,  93.8,    93,  90.9,  90.9,  88.3,  88.6,  63.8,    73,  52.5,  28.7,  12.1,  18.1,   9.3,   6.2),
     weight = c(    1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1),
                   sink  = TRUE,
       type = "DFOP",
         k1 = list(ini   = 0.1,
                   fixed = 0,
                   lower = 0.0,
                   upper = Inf),
         k2 = list(ini   = 0.01,
                   fixed = 0,
                   lower = 0.0,
                   upper = Inf),
          g = list(ini   = 0.5,
                   fixed = 0,
                   lower = 0.0,
                   upper = 1),
         M0 = list(ini   = 100.05,
                   fixed = 0,
                   lower = 0.0,
                   upper = Inf))),silent=TRUE)

#
# Fit and optimizer
#

  Fit0    <- try(KinEval(
              Scan00,
              evalMethod = 'IRLS',
             optimMethod = 'LM',
               plotfit   = TRUE,
               quiet     = TRUE,
               ctr       = kingui.control(
                              maxIter = 100,
                            tolerance = 1E-06,
                            odesolver = 'lsoda',
                            calcJacob=FALSE),
            irls.control = list(
                              maxIter = 10,
                            tolerance = 0.001)),silent=TRUE)

Fit    <- try(KinEval(
  Scan00,
  evalMethod = 'IRLS',
  optimMethod = 'LM',
  plotfit   = TRUE,
  quiet     = TRUE,
  ctr       = kingui.control(
    maxIter = 100,
    tolerance = 1E-06,
    odesolver = 'lsoda',
    calcJacob=TRUE),
  irls.control = list(
    maxIter = 10,
    tolerance = 0.001)),silent=TRUE)

#
# Output
#


 KinReport(Scan00,
            Fit,
            version = versioninfo[1,1],
            filename = "Scan00 LM DFOP observed")


#
# End of R-script
#

