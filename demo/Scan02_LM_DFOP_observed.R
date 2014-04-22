
# Check inst. versions of 'KineticEval' and 'R'
if(!require(KineticEval)) stop("Package KineticEval is not available")
if( as.numeric(version$major)<2||(as.numeric(version$major)==2 & as.numeric(version$minor)<15.1) ) stop("Please install R version > 2.15.0")


# Init 'R'
data(versioninfo)
require(methods)
require(logging)


# Create log file
basicConfig(level='FINEST')
addHandler(writeToFile, file="Scan02 LM DFOP observed.log", level='DEBUG')


#
# Residue data, pathway and kinetics
#

 Scan02 <- try(mkinmod.full(
   observed = list(
       time = c(    0,     0,     3,     3,     7,     7,    14,    14,    30,    30,    50,    50,    70,    70,   101,   101),
    residue = c(100.1,  97.1,  66.3,  71.2,  62.9,  65.4,  57.2,  61.1,  36.1,    39,  26.8,  27.6,  20.7,  19.6,  12.5,    15),
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
         M0 = list(ini   = 98.6,
                   fixed = 0,
                   lower = 0.0,
                   upper = Inf))),silent=TRUE)

#
# Fit and optimizer
#

  Fit    <- try(KinEval(
              Scan02,
              evalMethod = 'IRLS',
             optimMethod = 'LM',
               plotfit   = TRUE,
               quiet     = TRUE,
               ctr       = kingui.control(
                              maxIter = 100,
                            tolerance = 1E-06,
                            odesolver = 'lsoda'),
            irls.control = list(
                              maxIter = 10,
                            tolerance = 0.001)),silent=TRUE)

#
# Output
#


 KinReport(Scan02,
            Fit,
            version = versioninfo[1,1],
            filename = "Scan02 LM DFOP observed")


#
# End of R-script
#

