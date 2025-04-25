setwd('/home/369')
library(xtable)
library(forestplot)
library(forestploter)

logit <- function(x) log((1+x)/(1-x))
source('head.R')
source('generate.R')
est1=est2=est3=est4=est5=est6=estrho0=estzeta0=estrho0_p=estrho=estzeta=NULL
se1=se2=se3=se4=se5=se6=serho0=sezeta0=serho0_p=sezeta0_p=serho=sezeta=NULL
tm=NULL
bs_rho0 = bs_rho = bs_zeta0 = bs_zeta = NULL
for (i in 1:1000){
  fit = try(load(paste0('res',i,'.Rdata')), silent=TRUE)
  if ('try-error' %in% class(fit)) next
  est1 = rbind(est1, res1_con$est)
  est2 = rbind(est2, res2_con$est)
  est3 = rbind(est3, res1_ord$est)
  est4 = rbind(est4, res2_ord$est)
  est5 = rbind(est5, res1_sur$est)
  est6 = rbind(est6, res2_sur$est)
  se1 = rbind(se1, res1_con$se)
  se2 = rbind(se2, res2_con$se)
  se3 = rbind(se3, res1_ord$se)
  se4 = rbind(se4, res2_ord$se)
  se5 = rbind(se5, res1_sur$se)
  se6 = rbind(se6, res2_sur$se)
  tm = rbind(tm, times)
  estrho0 = rbind(estrho0, est_rho0)
  estzeta0 = rbind(estzeta0, est_zeta0)
  estrho0_p = rbind(estrho0_p, est_rho0_p)
  estrho = rbind(estrho, est_rho)
  estzeta = rbind(estzeta, est_zeta)
  serho0 = rbind(serho0, se_rho0)
  sezeta0 = rbind(sezeta0, se_zeta0)
  serho0_p = rbind(serho0_p, se_rho0)
  sezeta0_p = rbind(sezeta0_p, se_zeta0)
  serho = rbind(serho, se_rho)
  sezeta = rbind(sezeta, se_zeta)
}
# logit transformation for rho
# serho0 = 2/(1-estrho0^2)*serho0
# serho = 2/(1-estrho^2)*serho
# sezeta = 2/(1-estzeta^2)*sezeta
# estrho0 = logit(estrho0)
# estrho = logit(estrho)
# estzeta = logit(estzeta)

true1 = true1_con
true2 = true2_con
true3 = true1_ord
true4 = true2_ord
true5 = true1_sur
true6 = true2_sur
S1 = true1[1]^2+true1[2]+true1[6]
S2 = true2[1]^2+true2[2]+true2[6]
S3 = true3[1]^2+true3[2]+1
S4 = true4[1]^2+true4[2]+1
S5 = true5[1]^2+true5[2]+pi^2/6
S6 = true6[1]^2+true6[2]+pi^2/6
S1p = true1[1]^2+true1[2]
S2p = true2[1]^2+true2[2]
S3p = true3[1]^2+true3[2]
S4p = true4[1]^2+true4[2]
S5p = true5[1]^2+true5[2]
S6p = true6[1]^2+true6[2]
gamma12 = sqrt(true1[2]*true2[2])*0.5
gamma13 = sqrt(true1[2]*true3[2])*0.5
gamma14 = sqrt(true1[2]*true4[2])*0.5
gamma15 = sqrt(true1[2]*true5[2])*0.5
gamma16 = sqrt(true1[2]*true6[2])*0.5
gamma23 = sqrt(true2[2]*true3[2])*0.5
gamma24 = sqrt(true2[2]*true4[2])*0.5
gamma25 = sqrt(true2[2]*true5[2])*0.5
gamma26 = sqrt(true2[2]*true6[2])*0.5
gamma34 = sqrt(true3[2]*true4[2])*0.5
gamma35 = sqrt(true3[2]*true5[2])*0.5
gamma36 = sqrt(true3[2]*true6[2])*0.5
gamma45 = sqrt(true4[2]*true5[2])*0.5
gamma46 = sqrt(true4[2]*true6[2])*0.5
gamma56 = sqrt(true5[2]*true6[2])*0.5
rho1 = true1[2]/S1
rho2 = true2[2]/S2
rho3 = true3[2]/S3
rho4 = true4[2]/S4
rho5 = true5[2]/S5
rho6 = true6[2]/S6
zeta1 = true1[1]^2/S1
zeta2 = true2[1]^2/S2
zeta3 = true3[1]^2/S3
zeta4 = true4[1]^2/S4
zeta5 = true5[1]^2/S5
zeta6 = true6[1]^2/S6
rho1p = true1[2]/S1p
rho2p = true2[2]/S2p
rho3p = true3[2]/S3p
rho4p = true4[2]/S4p
rho5p = true5[2]/S5p
rho6p = true6[2]/S6p
rho12 = gamma12/sqrt(S1*S2)
rho13 = gamma13/sqrt(S1*S3)
rho14 = gamma14/sqrt(S1*S4)
rho15 = gamma15/sqrt(S1*S5)
rho16 = gamma16/sqrt(S1*S6)
rho23 = gamma23/sqrt(S2*S3)
rho24 = gamma24/sqrt(S2*S4)
rho25 = gamma25/sqrt(S2*S5)
rho26 = gamma26/sqrt(S2*S6)
rho34 = gamma34/sqrt(S3*S4)
rho35 = gamma35/sqrt(S3*S5)
rho36 = gamma36/sqrt(S3*S6)
rho45 = gamma45/sqrt(S4*S5)
rho46 = gamma46/sqrt(S4*S6)
rho56 = gamma56/sqrt(S5*S6)
zeta12 = true1[1]*true2[1]/sqrt(S1*S2)
zeta13 = true1[1]*true3[1]/sqrt(S1*S3)
zeta14 = true1[1]*true4[1]/sqrt(S1*S4)
zeta15 = true1[1]*true5[1]/sqrt(S1*S5)
zeta16 = true1[1]*true6[1]/sqrt(S1*S6)
zeta23 = true2[1]*true3[1]/sqrt(S2*S3)
zeta24 = true2[1]*true4[1]/sqrt(S2*S4)
zeta25 = true2[1]*true5[1]/sqrt(S2*S5)
zeta26 = true2[1]*true6[1]/sqrt(S2*S6)
zeta34 = true3[1]*true4[1]/sqrt(S3*S4)
zeta35 = true3[1]*true5[1]/sqrt(S3*S5)
zeta36 = true3[1]*true6[1]/sqrt(S3*S6)
zeta45 = true4[1]*true5[1]/sqrt(S4*S5)
zeta46 = true4[1]*true6[1]/sqrt(S4*S6)
zeta56 = true5[1]*true6[1]/sqrt(S5*S6)
true_rho0 = c(rho1,rho2,rho3,rho4,rho5,rho6)
true_zeta0 = c(zeta1,zeta2,zeta3,zeta4,zeta5,zeta6)
true_rho0p = c(rho1p,rho2p,rho3p,rho4p,rho5p,rho6p)
true_rho = c(rho12,rho13,rho14,rho15,rho16,rho23,rho24,rho25,rho26,
             rho34,rho35,rho36,rho45,rho46,rho56)
true_zeta = c(zeta12,zeta13,zeta14,zeta15,zeta16,zeta23,zeta24,zeta25,zeta26,
             zeta34,zeta35,zeta36,zeta45,zeta46,zeta56)
# true_rho0 = logit(true_rho0)
# true_rho = logit(true_rho)
# true_zeta = logit(true_zeta)

mean.est1 = colMeans(est1)
mean.est2 = colMeans(est2)
mean.est3 = colMeans(est3)
mean.est4 = colMeans(est4)
mean.est5 = colMeans(est5)
mean.est6 = colMeans(est6)
mean.se1 = colMeans(se1, na.rm=TRUE)
mean.se2 = colMeans(se2, na.rm=TRUE)
mean.se3 = colMeans(se3, na.rm=TRUE)
mean.se4 = colMeans(se4, na.rm=TRUE)
mean.se5 = colMeans(se5, na.rm=TRUE)
mean.se6 = colMeans(se6, na.rm=TRUE)
sd1 = apply(est1, 2, sd)
sd2 = apply(est2, 2, sd)
sd3 = apply(est3, 2, sd)
sd4 = apply(est4, 2, sd)
sd5 = apply(est5, 2, sd)
sd6 = apply(est6, 2, sd)
n = nrow(est1)
cv1 = colMeans((est1-1.96*se1-rep(1,n)%*%t(true1)<=0)*
               (est1+1.96*se1-rep(1,n)%*%t(true1)>=0), na.rm=TRUE)
cv2 = colMeans((est2-1.96*se2-rep(1,n)%*%t(true2)<=0)*
               (est2+1.96*se2-rep(1,n)%*%t(true2)>=0), na.rm=TRUE)
cv3 = colMeans((est3-1.96*se3-rep(1,n)%*%t(true3)<=0)*
               (est3+1.96*se3-rep(1,n)%*%t(true3)>=0), na.rm=TRUE)
cv4 = colMeans((est4-1.96*se4-rep(1,n)%*%t(true4)<=0)*
               (est4+1.96*se4-rep(1,n)%*%t(true4)>=0), na.rm=TRUE)
cv5 = colMeans((est5-1.96*se5-rep(1,n)%*%t(true5)<=0)*
               (est5+1.96*se5-rep(1,n)%*%t(true5)>=0), na.rm=TRUE)
cv6 = colMeans((est6-1.96*se6-rep(1,n)%*%t(true6)<=0)*
               (est6+1.96*se6-rep(1,n)%*%t(true6)>=0), na.rm=TRUE)

result1 = cbind(true1, mean.est1, sd1, mean.se1, cv1)
result2 = cbind(true2, mean.est2, sd2, mean.se2, cv2)
result3 = cbind(true3, mean.est3, sd3, mean.se3, cv3)
result4 = cbind(true4, mean.est4, sd4, mean.se4, cv4)
result5 = cbind(true5, mean.est5, sd5, mean.se5, cv5)
result6 = cbind(true6, mean.est6, sd6, mean.se6, cv6)
colnames(result1) = c('True','Est','SD','SE','Coverage')
colnames(result2) = c('True','Est','SD','SE','Coverage')
colnames(result3) = c('True','Est','SD','SE','Coverage')
colnames(result4) = c('True','Est','SD','SE','Coverage')
colnames(result5) = c('True','Est','SD','SE','Coverage')
colnames(result6) = c('True','Est','SD','SE','Coverage')

mean.estrho0 = colMeans(estrho0, na.rm=TRUE)
mean.estzeta0 = colMeans(estzeta0, na.rm=TRUE)
mean.estrho0p = colMeans(estrho0_p, na.rm=TRUE)
mean.estrho = colMeans(estrho, na.rm=TRUE)
mean.estzeta = colMeans(estzeta, na.rm=TRUE)
mean.serho0 = colMeans(serho0, na.rm=TRUE)
mean.sezeta0 = colMeans(sezeta0, na.rm=TRUE)
mean.serho0p = colMeans(serho0_p, na.rm=TRUE)
mean.serho = colMeans(serho, na.rm=TRUE)
mean.sezeta = colMeans(sezeta, na.rm=TRUE)
sd.rho0 = apply(estrho0, 2, sd)
sd.zeta0 = apply(estzeta0, 2, sd)
sd.rho0p = apply(estrho0_p, 2, sd)
sd.rho = apply(estrho, 2, sd)
sd.zeta = apply(estzeta, 2, sd)
low.rho0 = apply(estrho0, 2, function(x) quantile(x,.025))
up.rho0 = apply(estrho0, 2, function(x) quantile(x,.975))
low.zeta0 = apply(estzeta0, 2, function(x) quantile(x,.025))
up.zeta0 = apply(estzeta0, 2, function(x) quantile(x,.975))
low.rho = apply(estrho, 2, function(x) quantile(x,.025))
up.rho = apply(estrho, 2, function(x) quantile(x,.975))
low.zeta = apply(estzeta, 2, function(x) quantile(x,.025))
up.zeta = apply(estzeta, 2, function(x) quantile(x,.975))
bs_rho0 = estrho0 - rep(1,n)%*%t(true_rho0)
bs_zeta0 = estzeta0 - rep(1,n)%*%t(true_zeta0)
bs_rho0p = estrho0_p - rep(1,n)%*%t(true_rho0p)
bs_rho = estrho - rep(1,n)%*%t(true_rho)
bs_zeta = estzeta - rep(1,n)%*%t(true_zeta)
cv.rho0 = colMeans((estrho0-1.96*serho0-rep(1,n)%*%t(true_rho0)<=0)*
               (estrho0+1.96*serho0-rep(1,n)%*%t(true_rho0)>=0), na.rm=TRUE)
cv.zeta0 = colMeans((estzeta0-1.96*sezeta0-rep(1,n)%*%t(true_zeta0)<=0)*
                     (estzeta0+1.96*sezeta0-rep(1,n)%*%t(true_zeta0)>=0), na.rm=TRUE)
cv.rho0p = colMeans((estrho0_p-1.96*serho0_p-rep(1,n)%*%t(true_rho0p)<=0)*
               (estrho0_p+1.96*serho0_p-rep(1,n)%*%t(true_rho0p)>=0), na.rm=TRUE)
cv.rho = colMeans((estrho-1.96*serho-rep(1,n)%*%t(true_rho)<=0)*
               (estrho+1.96*serho-rep(1,n)%*%t(true_rho)>=0), na.rm=TRUE)
cv.zeta = colMeans((estzeta-1.96*sezeta-rep(1,n)%*%t(true_zeta)<=0)*
               (estzeta+1.96*sezeta-rep(1,n)%*%t(true_zeta)>=0), na.rm=TRUE)
result.rho0 = rbind(true_rho0, mean.estrho0, sd.rho0, mean.serho0, cv.rho0)
result.zeta0 = rbind(true_zeta0, mean.estzeta0, sd.zeta0, mean.sezeta0, cv.zeta0)
result.rho0p = rbind(true_rho0p, mean.estrho0p, sd.rho0p, mean.serho0p, cv.rho0p)
result.rho = rbind(true_rho, mean.estrho, sd.rho, mean.serho, cv.rho)
result.zeta = rbind(true_zeta, mean.estzeta, sd.zeta, mean.sezeta, cv.zeta)
rownames(result.rho0) = c('True','Est','SD','SE','Coverage')
rownames(result.zeta0) = c('True','Est','SD','SE','Coverage')
rownames(result.rho0p) = c('True','Est','SD','SE','Coverage')
rownames(result.rho) = c('True','Est','SD','SE','Coverage')
rownames(result.zeta) = c('True','Est','SD','SE','Coverage')

Type = c('(1,2)',
         '(1,3)',
         '(1,4)',
         '(1,5)',
         '(1,6)',
         '(2,3)',
         '(2,4)',
         '(2,5)',
         '(2,6)',
         '(3,4)',
         '(3,5)',
         '(3,6)',
         '(4,5)',
         '(4,6)',
         '(5,6)')
result_rho = data.frame(cbind(colMeans(na.omit(bs_rho)), colMeans(na.omit(serho))))
result_rho = cbind(Type, result_rho)
colnames(result_rho) = c('Type','Bias','SE')
result_rho = transform(result_rho, Lower=Bias-1.96*SE, Upper=Bias+1.96*SE)
result_zeta = data.frame(cbind(colMeans(na.omit(bs_zeta)), colMeans(na.omit(sezeta))))
result_zeta = cbind(Type, result_zeta)
colnames(result_zeta) = c('Type','Bias','SE')
result_zeta = transform(result_zeta, Lower=Bias-1.96*SE, Upper=Bias+1.96*SE)

xtable(rbind(result1,result.rho0[,1],result.zeta0[,1],result.rho0p[,1]), digits=3)
xtable(rbind(result2,result.rho0[,2],result.zeta0[,2],result.rho0p[,2]), digits=3)
xtable(rbind(result3,result.rho0[,3],result.zeta0[,3],result.rho0p[,3]), digits=3)
xtable(rbind(result4,result.rho0[,4],result.zeta0[,4],result.rho0p[,4]), digits=3)
xtable(rbind(result5,result.rho0[,5],result.zeta0[,5],result.rho0p[,5]), digits=3)
xtable(rbind(result6,result.rho0[,6],result.zeta0[,6],result.rho0p[,6]), digits=3)
colMeans(tm)/60
colnames(result.rho) = Type
colnames(result.zeta) = Type
xtable(t(result.rho), digits=3)
xtable(t(result.zeta), digits=3)
result.co = t(rbind(result.rho,result.zeta))
rownames(result.co) = Type
xtable(result.co, digits=3)

par(mfrow=c(1,2))
fig1 <- forestplot(
  title = 'Heritability',
  1:6,
  mean = mean.estrho0-true_rho0,
  lower = mean.estrho0-true_rho0-1.96*mean.serho0,
  upper = mean.estrho0-true_rho0+1.96*mean.serho0,
  zero = 0,
  boxsize = 0.1,
  graph.pos = 2,
  xticks = seq(-.3,.3,.1)
)
fig1
fig2 <- forestplot(
  title = 'Environmental effect',
  1:6,
  mean = mean.estzeta0-true_zeta0,
  lower = mean.estzeta0-true_zeta0-1.96*mean.sezeta0,
  upper = mean.estzeta0-true_zeta0+1.96*mean.sezeta0,
  zero = 0,
  boxsize = 0.1,
  graph.pos = 2,
  xticks = seq(-.3,.3,.1)
)
fig2
fig1 <- forestplot(
  title = 'Coheritability',
  result_rho[,1],
  mean = result_rho[,2],
  lower = result_rho[,4],
  upper = result_rho[,5],
  zero = 0,
  boxsize = 0.1,
  graph.pos = 2,
  xticks = seq(-.3,.3,.1)
)
fig1
fig2 <- forestplot(
  title = 'Environmental correlation',
  result_zeta[,1],
  mean = result_zeta[,2],
  lower = result_zeta[,4],
  upper = result_zeta[,5],
  zero = 0,
  boxsize = 0.1,
  graph.pos = 2,
  xticks = seq(-.3,.3,.1)
)
fig2

dt = data.frame(matrix(0,6,6))
colnames(dt) = c('Phe','rho','low1','up1','low2','up2')
dt$Phe = 1:6
dt$rho = mean.estrho0-true_rho0
low.estrho0 = apply(estrho0,2,function(x) quantile(x,.025))
up.estrho0 = apply(estrho0,2,function(x) quantile(x,.975))
dt$low1 = mean.estrho0-true_rho0-1.96*sd.rho0
dt$up1 = mean.estrho0-true_rho0+1.96*sd.rho0
dt$low2 = mean.estrho0-true_rho0-1.96*mean.serho0
dt$up2 = mean.estrho0-true_rho0+1.96*mean.serho0
dt$' ' = paste(rep(" ",15), collapse=" ")
tm <- forest_theme(
  legend_name = '',
  legend_value = c('SD','SE'))
rho0_p <- forest(dt[,c(1,7)],
                 est = list(dt$rho, dt$rho),
                 lower = list(dt$low1, dt$low2),
                 upper = list(dt$up1, dt$up2),
                 ci_column = 2,
                 xlim = c(-.4,.4),
                 vert_line = 0,
                 nudge_y = 0.2,
                 ticks_at = c(-.4,-.2,0,.2,.4),
                 title = 'Heritability',
                 theme = tm)
rho0_p
dt = data.frame(matrix(0,6,6))
colnames(dt) = c('Phe','zeta','low1','up1','low2','up2')
dt$Phe = 1:6
dt$zeta = mean.estzeta0-true_zeta0
low.estzeta0 = apply(estzeta0,2,function(x) quantile(x,.025))
up.estzeta0 = apply(estzeta0,2,function(x) quantile(x,.975))
dt$low1 = mean.estzeta0-true_zeta0-1.96*sd.zeta0
dt$up1 = mean.estzeta0-true_zeta0+1.96*sd.zeta0
dt$low2 = mean.estzeta0-true_zeta0-1.96*mean.sezeta0
dt$up2 = mean.estzeta0-true_zeta0+1.96*mean.sezeta0
dt$' ' = paste(rep(" ",15), collapse=" ")
tm <- forest_theme(
  legend_name = '',
  legend_value = c('SD','SE'))
zeta0_p <- forest(dt[,c(1,7)],
                 est = list(dt$zeta, dt$zeta),
                 lower = list(dt$low1, dt$low2),
                 upper = list(dt$up1, dt$up2),
                 ci_column = 2,
                 xlim = c(-.2,.2),
                 vert_line = 0,
                 nudge_y = 0.2,
                 ticks_at = c(-.2,-.1,0,.1,.2),
                 title = 'Environmental effect',
                 theme = tm)
zeta0_p

dt = data.frame(matrix(0,15,6))
colnames(dt) = c('Pair','rho','low1','up1','low2','up2')
dt$Pair = Type
dt$rho = mean.estrho-true_rho
low.estrho = apply(estrho,2,function(x) quantile(x,.025))
up.estrho = apply(estrho,2,function(x) quantile(x,.975))
dt$low1 =  mean.estrho-true_rho-1.96*sd.rho
dt$up1 =  mean.estrho-true_rho+1.96*sd.rho
dt$low2 = mean.estrho-true_rho-1.96*mean.serho
dt$up2 = mean.estrho-true_rho+1.96*mean.serho
dt$' ' = paste(rep(" ",20), collapse=" ")
tm <- forest_theme(
  legend_name = '',
  legend_value = c('SD','SE'))
rho_p <- forest(dt[,c(1,7)],
                 est = list(dt$rho, dt$rho),
                 lower = list(dt$low1, dt$low2),
                 upper = list(dt$up1, dt$up2),
                 xlim = c(-.2,.3),
                 ci_column = 2,
                 vert_line = 0,
                 nudge_y = 0.2,
                 ticks_at = c(-.2,-.1,0,.1,.2),
                 title = 'Coheritability',
                 theme = tm)
rho_p
dt = data.frame(matrix(0,15,6))
colnames(dt) = c('Pair','zeta','low1','up1','low2','up2')
dt$Pair = Type
dt$zeta = mean.estzeta-true_zeta
low.estzeta = apply(estzeta,2,function(x) quantile(x,.025))
up.estzeta = apply(estzeta,2,function(x) quantile(x,.975))
dt$low1 = mean.estzeta-true_zeta-1.96*sd.zeta
dt$up1 = mean.estzeta-true_zeta+1.96*sd.zeta
dt$low2 = mean.estzeta-true_zeta-1.96*mean.sezeta
dt$up2 = mean.estzeta-true_zeta+1.96*mean.sezeta
dt$' ' = paste(rep(" ",20), collapse=" ")
tm <- forest_theme(
  legend_name = '',
  legend_value = c('SD','SE'))
zeta_p <- forest(dt[,c(1,7)],
                  est = list(dt$zeta, dt$zeta),
                  lower = list(dt$low1, dt$low2),
                  upper = list(dt$up1, dt$up2),
                  xlim = c(-.2,.2),
                  ci_column = 2,
                  vert_line = 0,
                  nudge_y = 0.2,
                  ticks_at = c(-.2,-.1,0,.1,.2),
                  title = 'Environmental correlation',
                  theme = tm)
zeta_p

setwd('/home/default/370')
estrho0=estzeta0=estrho=estzeta=NULL
serho0=sezeta0=serho=sezeta=tm=NULL
for (i in 1:1000){
  fit = try(load(paste0('res',i,'.Rdata')), silent=TRUE)
  if ('try-error' %in% class(fit)) next
  estrho0 = rbind(estrho0, est_rho0)
  estzeta0 = rbind(estzeta0, est_zeta0)
  estrho = rbind(estrho, est_rho)
  estzeta = rbind(estzeta, est_zeta)
  serho0 = rbind(serho0, se_rho0)
  sezeta0 = rbind(sezeta0, se_zeta0)
  serho = rbind(serho, se_rho)
  sezeta = rbind(sezeta, se_zeta)
  tm = rbind(tm, times)
}
mean.estrho01 = colMeans(estrho0, na.rm=TRUE)
mean.estzeta01 = colMeans(estzeta0, na.rm=TRUE)
mean.estrho1 = colMeans(estrho, na.rm=TRUE)
mean.estzeta1 = colMeans(estzeta, na.rm=TRUE)
low.rho01 = apply(estrho0, 2, function(x) quantile(x,.025))
up.rho01 = apply(estrho0, 2, function(x) quantile(x,.975))
low.zeta01 = apply(estzeta0, 2, function(x) quantile(x,.025))
up.zeta01 = apply(estzeta0, 2, function(x) quantile(x,.975))
low.rho1 = apply(estrho, 2, function(x) quantile(x,.025))
up.rho1 = apply(estrho, 2, function(x) quantile(x,.975))
low.zeta1 = apply(estzeta, 2, function(x) quantile(x,.025))
up.zeta1 = apply(estzeta, 2, function(x) quantile(x,.975))
sd.rho01 = apply(estrho0, 2, sd)
sd.zeta01 = apply(estzeta0, 2, sd)
sd.rho1 = apply(estrho, 2, sd)
sd.zeta1 = apply(estzeta, 2, sd)
tm1 = mean(colSums(tm))

setwd('/home/default/369')
estrho0=estzeta0=estrho=estzeta=NULL
serho0=sezeta0=serho=sezeta=tm=NULL
for (i in 1:1000){
  fit = try(load(paste0('res',i,'.Rdata')), silent=TRUE)
  if ('try-error' %in% class(fit)) next
  estrho0 = rbind(estrho0, est_rho0)
  estzeta0 = rbind(estzeta0, est_zeta0)
  estrho = rbind(estrho, est_rho)
  estzeta = rbind(estzeta, est_zeta)
  serho0 = rbind(serho0, se_rho0)
  sezeta0 = rbind(sezeta0, se_zeta0)
  serho = rbind(serho, se_rho)
  sezeta = rbind(sezeta, se_zeta)
  tm = rbind(tm, times)
}
mean.estrho02 = colMeans(estrho0, na.rm=TRUE)
mean.estzeta02 = colMeans(estzeta0, na.rm=TRUE)
mean.estrho2 = colMeans(estrho, na.rm=TRUE)
mean.estzeta2 = colMeans(estzeta, na.rm=TRUE)
low.rho02 = apply(estrho0, 2, function(x) quantile(x,.025))
up.rho02 = apply(estrho0, 2, function(x) quantile(x,.975))
low.zeta02 = apply(estzeta0, 2, function(x) quantile(x,.025))
up.zeta02 = apply(estzeta0, 2, function(x) quantile(x,.975))
low.rho2 = apply(estrho, 2, function(x) quantile(x,.025))
up.rho2 = apply(estrho, 2, function(x) quantile(x,.975))
low.zeta2 = apply(estzeta, 2, function(x) quantile(x,.025))
up.zeta2 = apply(estzeta, 2, function(x) quantile(x,.975))
sd.rho02 = apply(estrho0, 2, sd)
sd.zeta02 = apply(estzeta0, 2, sd)
sd.rho2 = apply(estrho, 2, sd)
sd.zeta2 = apply(estzeta, 2, sd)
tm2 = mean(colSums(tm))

setwd('/default/371')
estrho0=estzeta0=estrho=estzeta=NULL
serho0=sezeta0=serho=sezeta=tm=NULL
for (i in 1:1000){
  fit = try(load(paste0('res',i,'.Rdata')), silent=TRUE)
  if ('try-error' %in% class(fit)) next
  estrho0 = rbind(estrho0, est_rho0)
  estzeta0 = rbind(estzeta0, est_zeta0)
  estrho = rbind(estrho, est_rho)
  estzeta = rbind(estzeta, est_zeta)
  serho0 = rbind(serho0, se_rho0)
  sezeta0 = rbind(sezeta0, se_zeta0)
  serho = rbind(serho, se_rho)
  sezeta = rbind(sezeta, se_zeta)
  tm = rbind(tm, times)
}
mean.estrho03 = colMeans(estrho0, na.rm=TRUE)
mean.estzeta03 = colMeans(estzeta0, na.rm=TRUE) 
mean.estrho3 = colMeans(estrho, na.rm=TRUE)
mean.estzeta3 = colMeans(estzeta, na.rm=TRUE)
low.rho03 = apply(estrho0, 2, function(x) quantile(x,.025))
up.rho03 = apply(estrho0, 2, function(x) quantile(x,.975))
low.zeta03 = apply(estzeta0, 2, function(x) quantile(x,.025))
up.zeta03 = apply(estzeta0, 2, function(x) quantile(x,.975))
low.rho3 = apply(estrho, 2, function(x) quantile(x,.025))
up.rho3 = apply(estrho, 2, function(x) quantile(x,.975))
low.zeta3 = apply(estzeta, 2, function(x) quantile(x,.025))
up.zeta3 = apply(estzeta, 2, function(x) quantile(x,.975))
sd.rho03 = apply(estrho0, 2, sd)
sd.zeta03 = apply(estzeta0, 2, sd)
sd.rho3 = apply(estrho, 2, sd)
sd.zeta3 = apply(estzeta, 2, sd)
tm3 = mean(colSums(tm))

dt = data.frame(matrix(0,6,10))
colnames(dt) = c('Phe','rho1','rho2','rho3','low1','up1','low2','up2','low3','up3')
dt$Phe = 1:6
dt$rho1 = mean.estrho01 - true_rho0
dt$rho2 = mean.estrho02 - true_rho0
dt$rho3 = mean.estrho03 - true_rho0
dt$low1 = mean.estrho01 - true_rho0 - 1.96*sd.rho01
dt$up1 = mean.estrho01 - true_rho0 + 1.96*sd.rho01
dt$low2 = mean.estrho02 - true_rho0 - 1.96*sd.rho02
dt$up2 = mean.estrho02 - true_rho0 + 1.96*sd.rho02
dt$low3 = mean.estrho03 - true_rho0 - 1.96*sd.rho03
dt$up3 = mean.estrho03 - true_rho0 + 1.96*sd.rho03
dt$' ' = paste(rep(" ",20), collapse=" ")
tm <- forest_theme(
  legend_name = 'Setting',
  legend_value = c('Small families','Moderate families','Large families'))
rho0_p <- forest(dt[,c(1,11)],
                 est = list(dt$rho1, dt$rho2, dt$rho3),
                 lower = list(dt$low1, dt$low2, dt$low3),
                 upper = list(dt$up1, dt$up2, dt$up3),
                 xlim = c(-.4,.4),
                 ci_column = 2,
                 vert_line = 0,
                 nudge_y = 0.2,
                 ticks_at = c(-.4,-.2,0,.2,.4),
                 title = 'Heritability',
                 theme = tm)
rho0_p
dt = data.frame(matrix(0,6,10))
colnames(dt) = c('Phe','zeta1','zeta2','zeta3','low1','up1','low2','up2','low3','up3')
dt$Phe = 1:6
dt$zeta1 = mean.estzeta01 - true_zeta0
dt$zeta2 = mean.estzeta02 - true_zeta0
dt$zeta3 = mean.estzeta03 - true_zeta0
dt$low1 = mean.estzeta01 - true_zeta0 - 1.96*sd.zeta01
dt$up1 = mean.estzeta01 - true_zeta0 + 1.96*sd.zeta01
dt$low2 = mean.estzeta02 - true_zeta0 - 1.96*sd.zeta02
dt$up2 = mean.estzeta02 - true_zeta0 + 1.96*sd.zeta02
dt$low3 = mean.estzeta03 - true_zeta0 - 1.96*sd.zeta03
dt$up3 = mean.estzeta03 - true_zeta0 + 1.96*sd.zeta03
dt$' ' = paste(rep(" ",20), collapse=" ")
tm <- forest_theme(
  legend_name = 'Setting',
  legend_value = c('Small families','Moderate families','Large families'))
zeta0_p <- forest(dt[,c(1,11)],
                  est = list(dt$zeta1, dt$zeta2, dt$zeta3),
                  lower = list(dt$low1, dt$low2, dt$low3),
                  upper = list(dt$up1, dt$up2, dt$up3),
                  xlim = c(-.2,.2),
                  ci_column = 2,
                  vert_line = 0,
                  nudge_y = 0.2,
                  ticks_at = c(-.2,-.1,0,.1,.2),
                  title = 'Environmental effect',
                  theme = tm)
zeta0_p

dt = data.frame(matrix(0,15,10))
colnames(dt) = c('Pair','rho1','rho2','rho3','low1','up1','low2','up2','low3','up3')
dt$Pair = Type
dt$rho1 = mean.estrho1 - true_rho
dt$rho2 = mean.estrho2 - true_rho
dt$rho3 = mean.estrho3 - true_rho
dt$low1 = mean.estrho1 - true_rho - 1.96*sd.rho1
dt$up1 = mean.estrho1 - true_rho + 1.96*sd.rho1
dt$low2 = mean.estrho2 - true_rho - 1.96*sd.rho2
dt$up2 = mean.estrho2 - true_rho + 1.96*sd.rho2
dt$low3 = mean.estrho3 - true_rho - 1.96*sd.rho3
dt$up3 = mean.estrho3 - true_rho + 1.96*sd.rho3
dt$' ' = paste(rep(" ",20), collapse=" ")
tm <- forest_theme(
  legend_name = 'Setting',
  legend_value = c('Small families','Moderate families','Large families'))
rho_p <- forest(dt[,c(1,11)],
                est = list(dt$rho1, dt$rho2, dt$rho3),
                lower = list(dt$low1, dt$low2, dt$low3),
                upper = list(dt$up1, dt$up2, dt$up3),
                xlim = c(-.2,.2),
                ci_column = 2,
                vert_line = 0,
                nudge_y = 0.2,
                ticks_at = c(-.2,-.1,0,.1,.2),
                title = 'Coheritability',
                theme = tm)
rho_p
dt = data.frame(matrix(0,15,10))
colnames(dt) = c('Pair','zeta1','zeta2','zeta3','low1','up1','low2','up2','low3','up3')
dt$Pair = Type
dt$zeta1 = mean.estzeta1 - true_zeta
dt$zeta2 = mean.estzeta2 - true_zeta
dt$zeta3 = mean.estzeta3 - true_zeta
dt$low1 = mean.estzeta1 - true_zeta - 1.96*sd.zeta1
dt$up1 = mean.estzeta1 - true_zeta + 1.96*sd.zeta1
dt$low2 = mean.estzeta2 - true_zeta - 1.96*sd.zeta2
dt$up2 = mean.estzeta2 - true_zeta + 1.96*sd.zeta2
dt$low3 = mean.estzeta3 - true_zeta - 1.96*sd.zeta3
dt$up3 = mean.estzeta3 - true_zeta + 1.96*sd.zeta3
dt$' ' = paste(rep(" ",20), collapse=" ")
tm <- forest_theme(
  legend_name = 'Setting',
  legend_value = c('Small families','Moderate families','Large families'))
zeta_p <- forest(dt[,c(1,11)],
                 est = list(dt$zeta1, dt$zeta2, dt$zeta3),
                 lower = list(dt$low1, dt$low2, dt$low3),
                 upper = list(dt$up1, dt$up2, dt$up3),
                 xlim = c(-.2,.2),
                 ci_column = 2,
                 vert_line = 0,
                 nudge_y = 0.2,
                 ticks_at = c(-.2,-.1,0,.1,.2),
                 title = 'Environmental correlation',
                 theme = tm)
zeta_p
c(tm1,tm2,tm3)