## Demo Of KineticEval
## library(roxygen2)
## library(devtools)
## roxygenize('KineticEval')
## build("KineticEval")
## install.packages("KineticEval_1.0-1.tar.gz",type="source")
library(KineticEval)


### Complicated cases:

data(BCS2)
ex2 <- mkinmod.full(Parent= list(type = "SFO",to = c( "Met1", "Met2","Met4", "Met5"),
                                 k = list(ini = 0.1,fixed = 0,lower = 0,upper = Inf),
                                 M0 = list(ini = 100,fixed = 0,lower = 0,upper = Inf),
                                 FF = list(ini = c(.1,.1,.1,.1),fixed = c(0,0,0,0),lower = c(0,0,0,0),upper = c(1,1,1,1))),
                    Met1 = list(type = "SFO",to = c("Met3", "Met4")),
                    Met2 = list(type = "SFO",to = c("Met3")),
                    Met3 = list(type = "SFO" ),
                    Met4 = list(type = "SFO", to = c("Met5")),
                    Met5 = list(type = "SFO"),
                    data=BCS2)
system.time(res_ex2_LM <- KinEval(ex2,evalMethod='NLLS',optimMethod="LM"))
## system.time(res_ex2_TRR <- KinEval(ex2,evalMethod='NLLS',optimMethod="TRR"))
## compare_multi_kinmod(ex2,rbind(t(res_ex2_LM$par),t(res_ex2_TRR$par)))


