## Demo Of KineticEval
## library(roxygen2)
## library(devtools)
## roxygenize('KineticEval')
## build("KineticEval")
## install.packages("KineticEval_1.0-1.tar.gz",type="source")
library(KineticEval)

## Datasets

data(BCS1)
BCS1
data(andrew)
andrew
data(ex1)
ex1
data(ex1_a)
ex1_a
## Dealing with Vesion
data(versioninfo)

## Compare the algorithms
mkinmodini <- mkinmod.full(Parent=list(type="SFO",to="Metab"),
                           Metab=list(type="SFO",M0=list(ini=0,fixed=0,lower=0,upper=Inf)),
                           data=andrew)
res0 <- KinEval(mkinmodini,evalMethod='NLLS',optimMethod="LM")
res0.1 <- KinEval(mkinmodini,evalMethod='IRLS',optimMethod="LM")


res1 <- KinEval(mkinmodini,evalMethod='NLLS',optimMethod="TRR")
res1.1 <- KinEval(mkinmodini,evalMethod='IRLS',optimMethod="TRR")


compare_multi_kinmod(mkinmodini,rbind(t(res0$par),t(res1$par)))

## Problematic cases
### Andrew
summary(res1)

### BCS1, with different starting values
res0 <- mkinfit.full(ex1,plot = TRUE, quiet= TRUE,ctr = kingui.control(method = "Marq",submethod = 'Port',maxIter = 100,tolerance = 1E-06, odesolver = 'lsoda'))
res0_a <- mkinfit.full(ex1_a,plot = TRUE, quiet= TRUE,ctr = kingui.control(method = "Marq",submethod = 'Port',maxIter = 100,tolerance = 1E-06, odesolver = 'lsoda'))
compare_multi_kinmod(ex1,rbind(t(res0$par),t(res0_a$par)))

res0 <- KinEval(ex1,evalMethod='NLLS',optimMethod="LM")
res0_a <- KinEval(ex1_a,evalMethod='NLLS',optimMethod="LM")
compare_multi_kinmod(ex1_a,rbind(t(res0$par),t(res0_a$par)))

res0 <- KinEval(ex1,evalMethod='NLLS',optimMethod="TRR")
res0_a <- KinEval(ex1_a,evalMethod='NLLS',optimMethod="TRR")
compare_multi_kinmod(ex1_a,rbind(t(res0$par),t(res0_a$par)))

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


