## Demo Of KineticEval

library(KineticEval)

## Datasets
pause <- function() invisible(readline())
data(BCS1)
BCS1
pause()
data(andrew)
andrew
pause()
## ex1 <- mkinmod.full(Parent=list(type="SFO",to="Metab"),Metab=list(type="SFO"),data=BCS1)
## ex1_a <- mkinmod.full(Parent=list(type="SFO",to="Metab",k=list(ini=0.0058,fixed=0,lower=0,upper=Inf)),Metab=list(type="SFO"),data=BCS1)
## save(ex1,file="KineticEval/data/ex1.rda")
## save(ex1_a,file="KineticEval/data/ex1_a.rda")
data(ex1)
ex1
pause()
data(ex1_a)
ex1_a
pause()
## Dealing with Vesion
data(versioninfo)
par(ask=FALSE)
## Compare the algorithms
mkinmodini <- mkinmod.full(Parent=list(type="SFO",to="Metab"),
                           Metab=list(type="SFO",M0=list(ini=0,fixed=0,lower=0,upper=Inf)),
                           data=andrew)
res0 <- KinEval(mkinmodini,evalMethod='NLLS',optimMethod="LM")

res0.1 <- KinEval(mkinmodini,evalMethod='IRLS',optimMethod="LM")


res1 <- KinEval(mkinmodini,evalMethod='NLLS',optimMethod="TRR")
pause()
res1.1 <- KinEval(mkinmodini,evalMethod='IRLS',optimMethod="TRR")
pause()

compare_multi_kinmod(mkinmodini,rbind(t(res0$par),t(res1$par)))
pause()
## Problematic cases
### Andrew
summary(res1)

### BCS1, with different starting values
res0 <- mkinfit.full(ex1,plot = TRUE, quiet= TRUE,ctr = kingui.control(method = "Marq",submethod = 'Port',maxIter = 100,tolerance = 1E-06, odesolver = 'lsoda'))
res0_a <- mkinfit.full(ex1_a,plot = TRUE, quiet= TRUE,ctr = kingui.control(method = "Marq",submethod = 'Port',maxIter = 100,tolerance = 1E-06, odesolver = 'lsoda'))
compare_multi_kinmod(ex1,rbind(t(res0$par),t(res0_a$par)))
pause()
res0 <- KinEval(ex1,evalMethod='NLLS',optimMethod="LM")
res0_a <- KinEval(ex1_a,evalMethod='NLLS',optimMethod="LM")
compare_multi_kinmod(ex1_a,rbind(t(res0$par),t(res0_a$par)))
pause()
res0 <- KinEval(ex1,evalMethod='NLLS',optimMethod="TRR")
res0_a <- KinEval(ex1_a,evalMethod='NLLS',optimMethod="TRR")
compare_multi_kinmod(ex1_a,rbind(t(res0$par),t(res0_a$par)))
pause()
### Complicated cases:

## system.time(res_ex2_TRR <- KinEval(ex2,evalMethod='NLLS',optimMethod="TRR"))
## compare_multi_kinmod(ex2,rbind(t(res_ex2_LM$par),t(res_ex2_TRR$par)))


