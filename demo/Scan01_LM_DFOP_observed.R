
# Check inst. versions of 'KineticEval' and 'R'
if(!require(KineticEval)) stop("Package KineticEval is not available")
if( as.numeric(version$major)<2||(as.numeric(version$major)==2 & as.numeric(version$minor)<15.1) ) stop("Please install R version > 2.15.0")


# Init 'R'
data(versioninfo)
require(methods)
require(logging)


# Create log file
basicConfig(level='FINEST')
addHandler(writeToFile, file="Scan01 LM DFOP observed.log", level='DEBUG')


#
# Residue data, pathway and kinetics
#

 Scan01 <- try(mkinmod.full(
   observed = list(
       time = c(    0,     0,     3,     3,     7,     7,    14,    14,    30,    30,    50,    50,    70,    70,   101,   101),
    residue = c(100.5,  99.6,    80,  80.3,  71.2,  69.5,  63.7,  64.2,  40.2,  47.8,  32.9,  14.8,   6.3,  10.5,   4.9,   2.9),
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

  Fit    <- try(KinEval(
              Scan01,
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


 KinReport(Scan01,
            Fit,
            version = versioninfo[1,1],
            filename = "Scan01 LM DFOP observed")


#
# End of R-script
#

