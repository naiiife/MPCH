# estimating coheritability

lid = Sys.getenv("SLURM_ARRAY_TASK_ID")
lid = as.numeric(lid)
i = lid
#i = (lid-1)%/%290+1
#j = (lid-1)%%290+1
#lid = (i-1)*290+j

setwd('/home/jobs')
source('head.R')
source('con_con.R')
source('con_ord.R')
source('con_sur.R')
source('ord_ord.R')
source('ord_sur.R')
source('sur_sur.R')

setwd('/home/ukb')
load('ukb_famdat.Rdata')
familyID = familyID[1:N]
eid = eid[1:N]
K1 = 152
K2 = 97
K3 = 27
K4 = 14
K = K1+K2+K3+K4

kinship = rep(1,length(familyID))
fID1 = sapply(1:nrow(pairs), function(l) familyID[which(eid==pairs[l,1])])
fID2 = sapply(1:nrow(pairs), function(l) familyID[which(eid==pairs[l,2])])
eli = fID1==fID2
pairs = pairs[eli,]
P1 = c(pairs[,1], pairs[,2])
P2 = c(pairs[,2], pairs[,1])
fIDco = fID1[eli]
fIDco = c(familyID, fIDco, fIDco)
kinco = c(kinship, pairs[,3], pairs[,3])
res = matrix(NA,nrow=K,ncol=12)

for (j in 1:K){
  if (i>=j) next
  setwd('/home/ukb/results_fam_all')
  fiti = fitj = NULL
  er = 0
  if (i<=K1) {
    tr = try(load(paste0('fit_con',i,'.Rdata')))
  }
  if (i>K1&i<=K1+K2) {
    tr = try(load(paste0('fit_bin',i-K1,'.Rdata')))
  }
  if (i>K1+K2&i<=K1+K2+K3) {
    tr = try(load(paste0('fit_ord',i-K1-K2,'.Rdata')))
  }
  if (i>K1+K2+K3) {
    tr = try(load(paste0('fit_sur',i-K1-K2-K3,'.Rdata')))
  }
  if (!'try-error' %in% class(tr)) {
    fiti = fit
  } else {
    er = 1
  }
  
  if (j<=K1) {
    tr = try(load(paste0('fit_con',j,'.Rdata')))
  }
  if (j>K1&j<=K1+K2) {
    tr = try(load(paste0('fit_bin',j-K1,'.Rdata')))
  }
  if (j>K1+K2&j<=K1+K2+K3) {
    tr = try(load(paste0('fit_ord',j-K1-K2,'.Rdata')))
  }
  if (j>K1+K2+K3) {
    tr = try(load(paste0('fit_sur',j-K1-K2-K3,'.Rdata')))
  }
  if (!'try-error' %in% class(tr)) {
    fitj = fit
  } else {
    er = 1
  }
  
  if (er==1) next
  setwd('/home/ukb/')
  if (i<=K1){
    load('condat.Rdata')
    Y_i = con.dat[1:N,i]
    rm(con.dat)
  }
  if (i>K1&i<=K1+K2){
    load('bindat.Rdata')
    Y_i = bin.dat[1:N,i-K1]+1
    rm(bin.dat)
  }
  if (i>K1+K2&i<=K1+K2+K3){
    load('orddat.Rdata')
    Y_i = ord.dat[1:N,i-K1-K2]
    rm(ord.dat)
  }
  if (i>K1+K2+K3){
    load('surdat.Rdata')
    T_i = sur.T.dat[1:N,i-K1-K2-K3]
    D_i = sur.D.dat[1:N,i-K1-K2-K3]
    rm(sur.T.dat)
    rm(sur.D.dat)
  }
  if (j<=K1){
    load('condat.Rdata')
    Y_j = con.dat[1:N,j]
    rm(con.dat)
  }
  if (j>K1&j<=K1+K2){
    load('bindat.Rdata')
    Y_j = bin.dat[1:N,j-K1]+1
    rm(bin.dat)
  }
  if (j>K1+K2&j<=K1+K2+K3){
    load('orddat.Rdata')
    Y_j = ord.dat[1:N,j-K1-K2]
    rm(ord.dat)
  }
  if (j>K1+K2+K3){
    load('surdat.Rdata')
    T_j = sur.T.dat[1:N,j-K1-K2-K3]
    D_j = sur.D.dat[1:N,j-K1-K2-K3]
    rm(sur.T.dat)
    rm(sur.D.dat)
  }
  
  setwd('/home/ukb/results_co/')
  X = as.matrix(X[1:N,])
  X0 = cbind(1,X)
  X_sur = X[,-2]
  p = ncol(X)
  
  if (i<=K1){
    if (j<=K1){
      Yi = c(Y_i, 
             sapply(1:length(P1), function(l) Y_i[which(eid==P1[l])]))
      Yj = c(Y_j, 
             sapply(1:length(P2), function(l) Y_j[which(eid==P2[l])]))
      Xi = rbind(X0, 
                 t(sapply(1:length(P1), function(l) X0[which(eid==P1[l]),])))
      Xj = rbind(X0, 
                 t(sapply(1:length(P2), function(l) X0[which(eid==P2[l]),])))
      mis = complete.cases(cbind(Yi,Yj))
      t0 = Sys.time()
      par0 = c(0,-sqrt(fiti$gamma*fitj$gamma),sqrt(fiti$gamma*fitj$gamma))
      resij = con_con(Yi[mis],Yj[mis],Xi[mis,],Xj[mis,],fIDco[mis],
                      fiti,fitj,kinco[mis],par0)
      t1 = Sys.time()
      resij$tm = as.numeric(difftime(t1,t0,units='secs'))
    }
    if (j>K1&j<=K1+K2){
      Yi = c(Y_i, 
             sapply(1:length(P1), function(l) Y_i[which(eid==P1[l])]))
      Yj = c(Y_j, 
             sapply(1:length(P2), function(l) Y_j[which(eid==P2[l])]))
      Xi = rbind(X0, 
                 t(sapply(1:length(P1), function(l) X0[which(eid==P1[l]),])))
      Xj = rbind(X, 
                 t(sapply(1:length(P2), function(l) X[which(eid==P2[l]),])))
      mis = complete.cases(cbind(Yi,Yj))
      t0 = Sys.time()
      par0 = c(0,-sqrt(fiti$gamma*fitj$gamma),sqrt(fiti$gamma*fitj$gamma))
      resij = con_ord(Yi[mis],Yj[mis],Xi[mis,],Xj[mis,],fIDco[mis],
                      fiti,fitj,kinco[mis],par0)
      t1 = Sys.time()
      resij$tm = as.numeric(difftime(t1,t0,units='secs'))
    }
    if (j>K1+K2&j<=K1+K2+K3){
      Yi = c(Y_i, 
             sapply(1:length(P1), function(l) Y_i[which(eid==P1[l])]))
      Yj = c(Y_j, 
             sapply(1:length(P2), function(l) Y_j[which(eid==P2[l])]))
      Xi = rbind(X0, 
                 t(sapply(1:length(P1), function(l) X0[which(eid==P1[l]),])))
      Xj = rbind(X, 
                 t(sapply(1:length(P2), function(l) X[which(eid==P2[l]),])))
      mis = complete.cases(cbind(Yi,Yj)) 
      t0 = Sys.time()
      par0 = c(0,-sqrt(fiti$gamma*fitj$gamma),sqrt(fiti$gamma*fitj$gamma))
      resij = con_ord(Yi[mis],Yj[mis],Xi[mis,],Xj[mis,],fIDco[mis],
                      fiti,fitj,kinco[mis],par0) 
      t1 = Sys.time()
      resij$tm = as.numeric(difftime(t1,t0,units='secs'))
    }
    if (j>K1+K2+K3){
      Yi = c(Y_i, 
             sapply(1:length(P1), function(l) Y_i[which(eid==P1[l])]))
      Tj = c(T_j, 
             sapply(1:length(P2), function(l) T_j[which(eid==P2[l])]))
      Dj = c(D_j, 
             sapply(1:length(P2), function(l) D_j[which(eid==P2[l])]))
      Xi = rbind(X0, 
                 t(sapply(1:length(P1), function(l) X0[which(eid==P1[l]),])))
      Xj = rbind(X_sur, 
                 t(sapply(1:length(P2), function(l) X_sur[which(eid==P2[l]),])))
      mis = complete.cases(cbind(Yi,Tj,Dj))
      fitj$h = fitj$h[1:(p+1)]
      t0 = Sys.time()
      par0 = c(0,-sqrt(fiti$gamma*fitj$gamma),sqrt(fiti$gamma*fitj$gamma))
      resij = con_sur(Yi[mis],Tj[mis],Dj[mis],Xi[mis,],Xj[mis,],fIDco[mis],
                      fiti,fitj,kinco[mis],par0)
      t1 = Sys.time()
      resij$tm = as.numeric(difftime(t1,t0,units='secs'))
    }
  }
  
  if (i>K1&i<=K1+K2){
    if (j>K1&j<=K1+K2){
      Yi = c(Y_i, 
             sapply(1:length(P1), function(l) Y_i[which(eid==P1[l])]))
      Yj = c(Y_j, 
             sapply(1:length(P2), function(l) Y_j[which(eid==P2[l])]))
      Xi = rbind(X, 
                 t(sapply(1:length(P1), function(l) X[which(eid==P1[l]),])))
      Xj = rbind(X, 
                 t(sapply(1:length(P2), function(l) X[which(eid==P2[l]),])))
      mis = complete.cases(cbind(Yi,Yj))
      t0 = Sys.time()
      par0 = c(0,-sqrt(fiti$gamma*fitj$gamma),sqrt(fiti$gamma*fitj$gamma))
      resij = ord_ord(Yi[mis],Yj[mis],Xi[mis,],Xj[mis,],fIDco[mis],
                      fiti,fitj,kinco[mis],par0) 
      t1 = Sys.time()
      resij$tm = as.numeric(difftime(t1,t0,units='secs'))
    }
    if (j>K1+K2&j<=K1+K2+K3){
      Yi = c(Y_i, 
             sapply(1:length(P1), function(l) Y_i[which(eid==P1[l])]))
      Yj = c(Y_j, 
             sapply(1:length(P2), function(l) Y_j[which(eid==P2[l])]))
      Xi = rbind(X, 
                 t(sapply(1:length(P1), function(l) X[which(eid==P1[l]),])))
      Xj = rbind(X, 
                 t(sapply(1:length(P2), function(l) X[which(eid==P2[l]),])))
      mis = complete.cases(cbind(Yi,Yj))
      t0 = Sys.time()
      par0 = c(0,-sqrt(fiti$gamma*fitj$gamma),sqrt(fiti$gamma*fitj$gamma))
      resij = ord_ord(Yi[mis],Yj[mis],Xi[mis,],Xj[mis,],fIDco[mis],
                      fiti,fitj,kinco[mis],par0)
      t1 = Sys.time()
      resij$tm = as.numeric(difftime(t1,t0,units='secs'))
    }
    if (j>K1+K2+K3){
      Yi = c(Y_i, 
             sapply(1:length(P1), function(l) Y_i[which(eid==P1[l])]))
      Tj = c(T_j, 
             sapply(1:length(P2), function(l) T_j[which(eid==P2[l])]))
      Dj = c(D_j, 
             sapply(1:length(P2), function(l) D_j[which(eid==P2[l])]))
      Xi = rbind(X, 
                 t(sapply(1:length(P1), function(l) X[which(eid==P1[l]),])))
      Xj = rbind(X_sur, 
                 t(sapply(1:length(P2), function(l) X_sur[which(eid==P2[l]),])))
      mis = complete.cases(cbind(Yi,Tj,Dj))
      fitj$h = fitj$h[1:(p+1)]
      t0 = Sys.time()
      par0 = c(0,-sqrt(fiti$gamma*fitj$gamma),sqrt(fiti$gamma*fitj$gamma))
      resij = ord_sur(Yi[mis],Tj[mis],Dj[mis],Xi[mis,],Xj[mis,],fIDco[mis],
                      fiti,fitj,kinco[mis],par0)
      t1 = Sys.time()
      resij$tm = as.numeric(difftime(t1,t0,units='secs'))
    }
  }
  
  if (i>K1+K2&i<=K1+K2+K3){
    if (j>K1+K2&j<=K1+K2+K3){
      Yi = c(Y_i, 
             sapply(1:length(P1), function(l) Y_i[which(eid==P1[l])]))
      Yj = c(Y_j, 
             sapply(1:length(P2), function(l) Y_j[which(eid==P2[l])]))
      Xi = rbind(X, 
                 t(sapply(1:length(P1), function(l) X[which(eid==P1[l]),])))
      Xj = rbind(X, 
                 t(sapply(1:length(P2), function(l) X[which(eid==P2[l]),])))
      mis = complete.cases(cbind(Yi,Yj))
      t0 = Sys.time()
      par0 = c(0,-sqrt(fiti$gamma*fitj$gamma),sqrt(fiti$gamma*fitj$gamma))
      resij = ord_ord(Yi[mis],Yj[mis],Xi[mis,],Xj[mis,],fIDco[mis],
                      fiti,fitj,kinco[mis],par0)
      t1 = Sys.time()
      resij$tm = as.numeric(difftime(t1,t0,units='secs'))
    }
    if (j>K1+K2+K3){
      Yi = c(Y_i, 
             sapply(1:length(P1), function(l) Y_i[which(eid==P1[l])]))
      Tj = c(T_j, 
             sapply(1:length(P2), function(l) T_j[which(eid==P2[l])]))
      Dj = c(D_j, 
             sapply(1:length(P2), function(l) D_j[which(eid==P2[l])]))
      Xi = rbind(X, 
                 t(sapply(1:length(P1), function(l) X[which(eid==P1[l]),])))
      Xj = rbind(X_sur, 
                 t(sapply(1:length(P2), function(l) X_sur[which(eid==P2[l]),])))
      mis = complete.cases(cbind(Yi,Tj,Dj))
      fitj$h = fitj$h[1:(p+1)]
      t0 = Sys.time()
      par0 = c(0,-sqrt(fiti$gamma*fitj$gamma),sqrt(fiti$gamma*fitj$gamma))
      resij = ord_sur(Yi[mis],Tj[mis],Dj[mis],Xi[mis,],Xj[mis,],fIDco[mis],
                      fiti,fitj,kinco[mis],par0)
      t1 = Sys.time()
      resij$tm = as.numeric(difftime(t1,t0,units='secs'))
    }
  }
  
  if (i>K1+K2+K3){
    if (j>K1+K2+K3){
      Ti = c(T_i, 
             sapply(1:length(P1), function(l) T_i[which(eid==P1[l])]))
      Di = c(D_i, 
             sapply(1:length(P1), function(l) D_i[which(eid==P1[l])]))
      Tj = c(T_j, 
             sapply(1:length(P2), function(l) T_j[which(eid==P2[l])]))
      Dj = c(D_j, 
             sapply(1:length(P2), function(l) D_j[which(eid==P2[l])]))
      Xi = rbind(X_sur, 
                 t(sapply(1:length(P1), function(l) X_sur[which(eid==P1[l]),])))
      Xj = rbind(X_sur, 
                 t(sapply(1:length(P2), function(l) X_sur[which(eid==P2[l]),])))
      mis = complete.cases(cbind(Ti,Di,Tj,Dj))
      fiti$h = fiti$h[1:(p+1)]
      fitj$h = fitj$h[1:(p+1)]
      t0 = Sys.time()
      par0 = c(0,-sqrt(fiti$gamma*fitj$gamma),sqrt(fiti$gamma*fitj$gamma))
      resij = sur_sur(Ti[mis],Di[mis],Tj[mis],Dj[mis],Xi[mis,],Xj[mis,],fIDco[mis],
                      fiti,fitj,kinco[mis],par0)
      t1 = Sys.time()
      resij$tm = as.numeric(difftime(t1,t0,units='secs'))
    }
  }
  
  rho.f = resij$rho12_p
  se.rho.f = resij$se.rho12_p
  rho = resij$rho12
  se.rho = resij$se.rho12
  zeta.f = resij$zeta12_p
  se.zeta.f = resij$se.zeta12_p
  zeta = resij$zeta12
  se.zeta = resij$se.zeta12
  gamma = resij$gamma12
  se.gamma = resij$se
  tm = resij$tm
  res[j,] = c(rho.f, se.rho.f, zeta.f, se.zeta.f, 
              rho, se.rho, zeta, se.zeta, gamma, se.gamma,
              sum(mis),tm)
  save(res, file=paste0('res',lid,'.Rdata'))
}

setwd('/home/ukb/results_co/')
save(res, file=paste0('res',lid,'.Rdata'))