# Trial           : newWaterSediment03280258
# File name       : newWaterSediment03280258 LM SFO  ParentW.r
# Target path     : C:\Emod\KinGUIIv2.1\WorkingDirectory\TestcasesJW\new_WaterSediment
# Created         : on 28 Mar 2014
#                   at 03:45
#                   by ernom on ADEMONC7868(4CPUs)
# KinGUII version : 2.2014.224.1704
# Comments        : 
# Data File  : C:\Emod\KinGUIIv2.1\WorkingDirectory\TestcasesJW\new_WaterSediment\WaterSedSlow.txt
# #  kdegPW=0.05,
# #  kdegPS=0.05,
# #  kdegMW = 0.07,
# #  kdegMS = 0.02,
# #  kmetW = 0.2,
# #  kmetS = 0.3,
# #  kPWS = 0.1,
# #  kPSW = 0.2,
# #  kMWS = 0.4,
# #  kMSW = 0.3

# Check inst. versions of 'KineticEval' and 'R'
if(!require(KineticEval)) stop("Package KineticEval is not available")
if( as.numeric(version$major)<2||(as.numeric(version$major)==2 & as.numeric(version$minor)<15.1) ) stop("Please install R version > 2.15.0")


# Init 'R'
data(versioninfo)
require(methods)
require(logging)


# Create log file
basicConfig(level='FINEST')
addHandler(writeToFile, file="newWaterSediment03280258 LM SFO  ParentW.log", level='DEBUG')


#
# Residue data, pathway and kinetics
#

 newWaterSediment03280258 <- try(mkinmod.full(
    ParentW = list(
       time = c(         0,          1,          2,          3,          4,          5,          6,          7,          8,          9,         10,         11,         12,         13,         14,         15,         16,         17,         18,         19,         20,         21,         22,         23,         24,         25,         26,         27,         28,         29,         30),
    residue = c(  99.97374,   71.83075,   52.47262,   35.97579,   28.50663,   20.36662,   15.82219,   11.52977,   9.313581,   5.923083,   3.613877,   2.919582,   2.439954,   1.709784,   1.488776,   1.357079,   1.956343,  0.6525393,  0.1870977,          0,   0.883141,          0,  0.2994365, 0.03312699,   1.946525,          0,  0.8791396,          0,          0,          0,  0.2914969),
     weight = c(         1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1),
                      to = c("ParentS","MetW"),
         FF = list(ini   = c(0.1,0.1),
                   fixed = c(0,0),
                   lower = c(0.0,0.0),
                   upper = c(1.0,1.0)),
                   sink  = TRUE,
       type = "SFO",
          k = list(ini   = 0.1,
                   fixed = 0,
                   lower = 0.0,
                   upper = Inf),
         M0 = list(ini   = 0.1,
                   fixed = 0,
                   lower = 0.0,
                   upper = Inf)),
    ParentS = list(
       time = c(         0,          1,          2,          3,          4,          5,          6,          7,          8,          9,         10,         11,         12,         13,         14,         15,         16,         17,         18,         19,         20,         21,         22,         23,         24,         25,         26,         27,         28,         29,         30),
    residue = c(         0,    6.89447,   7.443442,   8.920093,   7.839042,   7.087455,   4.944009,   3.862302,   4.448714,   3.212410,   2.595926,   2.738703,   1.717933,          0,   1.636281,  0.8892295,   1.442022,  0.3111406,  0.9786097, 0.02926622,          0,          0,  0.3656540,  0.9979095,          0,          0, 0.08981933,          0,   1.320007,   1.363480,  0.3409676),
     weight = c(         1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1),
                      to = c("ParentW","MetS"),
         FF = list(ini   = c(0.1,0.1),
                   fixed = c(0,0),
                   lower = c(0.0,0.0),
                   upper = c(1.0,1.0)),
                   sink  = TRUE,
       type = "SFO",
          k = list(ini   = 0.1,
                   fixed = 0,
                   lower = 0.0,
                   upper = Inf),
         M0 = list(ini   = 0,
                   fixed = 1,
                   lower = 0.0,
                   upper = Inf)),
       MetW = list(
       time = c(       0,        1,        2,        3,        4,        5,        6,        7,        8,        9,       10,       11,       12,       13,       14,       15,       16,       17,       18,       19,       20,       21,       22,       23,       24,       25,       26,       27,       28,       29,       30),
    residue = c(       0, 13.97654, 19.89485, 21.35401, 25.47581, 25.98307, 25.06308, 25.41788, 24.36118, 23.39334,  23.9051, 24.16468, 22.28418, 23.03641, 19.06472, 20.35343, 20.84758,  19.7099, 17.55303, 16.85331, 16.50262, 15.69098, 16.10898, 14.17591, 12.64925, 14.20439, 13.10763, 12.62307, 13.46294, 11.84691, 10.85763),
     weight = c(       1,        1,        1,        1,        1,        1,        1,        1,        1,        1,        1,        1,        1,        1,        1,        1,        1,        1,        1,        1,        1,        1,        1,        1,        1,        1,        1,        1,        1,        1,        1),
                      to = c("MetS"),
         FF = list(ini   = c(0.1),
                   fixed = c(0),
                   lower = c(0.0),
                   upper = c(1.0)),
                   sink  = TRUE,
       type = "SFO",
          k = list(ini   = 0.1,
                   fixed = 0,
                   lower = 0.0,
                   upper = Inf),
         M0 = list(ini   = 0,
                   fixed = 1,
                   lower = 0.0,
                   upper = Inf)),
       MetS = list(
       time = c(       0,        1,        2,        3,        4,        5,        6,        7,        8,        9,       10,       11,       12,       13,       14,       15,       16,       17,       18,       19,       20,       21,       22,       23,       24,       25,       26,       27,       28,       29,       30),
    residue = c(       0,  2.74402,  10.6276, 17.68945, 21.46951, 27.29593, 30.70874, 31.42721, 32.09972, 31.83171, 32.35434, 31.09248,  32.5204, 30.24904, 29.90265, 29.83513, 28.15532, 26.91444, 26.46741, 25.63678, 22.59982, 23.96371, 21.79595, 21.39598, 19.93387, 17.82440, 19.12311, 18.45624, 17.53936, 16.54211, 16.90705),
     weight = c(       1,        1,        1,        1,        1,        1,        1,        1,        1,        1,        1,        1,        1,        1,        1,        1,        1,        1,        1,        1,        1,        1,        1,        1,        1,        1,        1,        1,        1,        1,        1),
                      to = c("MetW"),
         FF = list(ini   = c(0.1),
                   fixed = c(0),
                   lower = c(0.0),
                   upper = c(1.0)),
                   sink  = TRUE,
       type = "SFO",
          k = list(ini   = 0.1,
                   fixed = 0,
                   lower = 0.0,
                   upper = Inf),
         M0 = list(ini   = 0,
                   fixed = 1,
                   lower = 0.0,
                   upper = Inf)),
   inpartri = 'default',
  outpartri = 'water-sediment'),silent=TRUE)

#
# Fit and optimizer
#

  Fit    <- try(KinEval(
              newWaterSediment03280258,
              evalMethod = 'NLLS',
             optimMethod = 'TRR',
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


 KinReport(newWaterSediment03280258,
            Fit,
            version = versioninfo[1,1],
            filename = "newWaterSediment03280258 LM SFO  ParentW")


#
# End of R-script
#

