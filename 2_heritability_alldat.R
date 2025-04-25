## Estimating single-trait parameters

lid = Sys.getenv("SLURM_ARRAY_TASK_ID")
lid = as.numeric(lid)

setwd('/home/jobs')
source('head.R')
source('continuous.R')
source('probit.R')
source('ordinal.R')
source('survival.R')

K1 = 152
K2 = 97
K3 = 27
K4 = 14
K = K1+K2+K3+K4
#i = (lid-1)%%K+1
i = lid

setwd('/home/ukb')
if (i<=K1) load('condat_all.Rdata')
if (i>K1&i<=K1+K2) load('bindat_all.Rdata')
if (i>K1+K2&i<=K1+K2+K3) load('orddat_all.Rdata')
if (i>K1+K2+K3) load('surdat_all.Rdata')

# analysis
setwd('/home/ukb/results_alldat')

#N = n2*2+n3*3+n4*4
N = nrow(X)
X = as.matrix(X[1:N,])
familyID = familyID[1:N]
KM = KM[1:N,]

if (i <= K1){
  setwd('/home/ukb/results_fam_all')
  load(paste0('fit_con',i,'.Rdata'))
  par0 = fit$est
  setwd('/home/ukb/results_alldat')
  X0 = cbind(1,X)
  Y_i = con.dat[1:N,i]
  rm(con.dat)
  famlist = sort(unique(familyID))
  Ynew = NULL; familyIDnew = NULL; KMnew = NULL; Xnew = NULL
  for (j in 1:length(famlist)){
    famid = famlist[j]
    famsize = sum(familyID==famid)
    if (famsize==0) next
    if (famsize>0){
      Yi = Y_i[familyID==famid]
      KMi = as.matrix(KM[familyID==famid,1:famsize])
      Xi = as.matrix(X0[familyID==famid,])
      if (ncol(Xi)==1) Xi = t(Xi)
      obs = !is.nan(Yi)
      ni = sum(obs)
      if (ni==0) next
      KMj = matrix(0,ni,4)
      KMj[,1:ni] = KMi[obs,obs]
      Ynew = append(Ynew, Yi[obs])
      familyIDnew = append(familyIDnew, rep(famid,ni))
      KMnew = rbind(KMnew, KMj)
      Xnew = rbind(Xnew, Xi[obs,])
    }
  }
  Y_i = Ynew; familyID = familyIDnew; KM = KMnew; X0 = Xnew
  #fit0 = lm(Y_i~0+X0)
  t0 = Sys.time()
  fit = cal_con(Y_i,X0,familyID,KM,par0)
  t1 = Sys.time()
  fit$tm = as.numeric(difftime(t1,t0,units = "secs"))
  save(fit, file=paste0('fit_con',i,'.Rdata'))
}

if (i>K1&i<=K2+K1){
  i = i-K1
  setwd('/home/ukb/results_fam_all')
  load(paste0('fit_bin',i,'.Rdata'))
  par0 = fit$est
  setwd('/home/ukb/results_alldat')
  X_bin = X
  Y_i = bin.dat[1:N,i]
  rm(bin.dat)
  famlist = sort(unique(familyID))
  Ynew = NULL; familyIDnew = NULL; KMnew = NULL; Xnew = NULL
  for (j in 1:length(famlist)){
    famid = famlist[j]
    famsize = sum(familyID==famid)
    if (famsize==0) next
    if (famsize>0){
      Yi = Y_i[familyID==famid]
      KMi = as.matrix(KM[familyID==famid,1:famsize])
      Xi = as.matrix(X_bin[familyID==famid,])
      if (ncol(Xi)==1) Xi = t(Xi)
      obs = !is.nan(Yi)
      ni = sum(obs)
      if (ni==0) next
      KMj = matrix(0,ni,4)
      KMj[,1:ni] = KMi[obs,obs]
      Ynew = append(Ynew, Yi[obs])
      familyIDnew = append(familyIDnew, rep(famid,ni))
      KMnew = rbind(KMnew, KMj)
      Xnew = rbind(Xnew, Xi[obs,])
    }
  }
  Y_i = Ynew; familyID = familyIDnew; KM = KMnew; X_bin = Xnew
  #fit0 = glm(Y_i~1+X_bin, family='binomial'(link=probit))
  #par0 = rep(0.5,3+ncol(X_bin))
  #par0[2+1:ncol(X_bin)] = fit0$coefficients[-1]
  #par0[3+ncol(X_bin)] = -fit0$coefficients[1]
  t0 = Sys.time()
  fit = cal_bin(Y_i,X_bin,familyID,KM,par0)
  t1 = Sys.time()
  fit$tm = as.numeric(difftime(t1,t0,units = "secs"))
  save(fit, file=paste0('fit_bin',i,'.Rdata'))
}

if (i>K1+K2&i<=K1+K2+K3){
  i = i-K1-K2
  setwd('/home/ukb/results_fam_all')
  load(paste0('fit_ord',i,'.Rdata'))
  par0 = fit$est
  setwd('/home/ukb/results_alldat')
  X_ord = X
  Y_i = ord.dat[1:N,i]
  rm(ord.dat)
  L = max(Y_i, na.rm=TRUE)
  famlist = sort(unique(familyID))
  Ynew = NULL; familyIDnew = NULL; KMnew = NULL; Xnew = NULL
  for (j in 1:length(famlist)){
    famid = famlist[j]
    famsize = sum(familyID==famid)
    if (famsize==0) next
    if (famsize>0){
      Yi = Y_i[familyID==famid]
      KMi = as.matrix(KM[familyID==famid,1:famsize])
      Xi = as.matrix(X_ord[familyID==famid,])
      if (ncol(Xi)==1) Xi = t(Xi)
      obs = !is.nan(Yi)
      ni = sum(obs)
      if (ni==0) next
      KMj = matrix(0,ni,4)
      KMj[,1:ni] = KMi[obs,obs]
      Ynew = append(Ynew, Yi[obs])
      familyIDnew = append(familyIDnew, rep(famid,ni))
      KMnew = rbind(KMnew, KMj)
      Xnew = rbind(Xnew, Xi[obs,])
    }
  }
  Y_i = Ynew; familyID = familyIDnew; KM = KMnew; X_ord = Xnew
  #fit0 = lm(Y_i~1+X_ord)
  #coeff = fit0$coefficients[-1]
  #resids = sd(fit0$residuals)
  #par0 = rep(0.5,1+ncol(X_ord)+L)
  #par0[2+1:ncol(X_ord)] = coeff/resids
  #par0[1+ncol(X_ord)+2:L] = (2:L+fit0$coefficients[1]-1)/resids
  t0 = Sys.time()
  fit = cal_ord(Y_i,X_ord,familyID,KM,par0)
  t1 = Sys.time()
  fit$tm = as.numeric(difftime(t1,t0,units = "secs"))
  save(fit, file=paste0('fit_ord',i,'.Rdata'))
}

if (i>K1+K2+K3){
  #library(survival)
  i = i-K1-K2-K3
  setwd('/home/ukb/results_fam_all')
  load(paste0('fit_sur',i,'.Rdata'))
  par0 = fit$est
  setwd('/home/ukb/results_alldat')
  X_sur = X[,-2]
  T_sur = sur.T.dat[1:N,i]
  rm(sur.T.dat)
  D_sur = as.numeric(sur.D.dat[1:N,i])
  rm(sur.D.dat)
  #fit0 = coxph(Surv(T_sur,D_sur)~X_sur)
  #par0[2+1:ncol(X_sur)] = fit0$coefficients
  t0 = Sys.time()
  fit = cal_sur(T_sur,D_sur,X_sur,familyID,KM, par0)
  t1 = Sys.time()
  fit$tm = as.numeric(difftime(t1,t0,units = "secs"))
  save(fit, file=paste0('fit_sur',i,'.Rdata'))
}
