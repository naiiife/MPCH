# estimating coheritability

lid = Sys.getenv("SLURM_ARRAY_TASK_ID")
lid = as.numeric(lid)
i = lid
i = (lid-1)%/%29+1 #1:290
b = (lid-1)%%29+1  #1:29
#lid = (i-1)*290+j

setwd('/home/jobs')
source('head.R')
source('con_con.R')
source('con_ord.R')
source('con_sur.R')
source('ord_ord.R')
source('ord_sur.R')
source('sur_sur.R')

K1 = 152
K2 = 97
K3 = 27
K4 = 14
K = K1+K2+K3+K4
size = 502155
p = 14
setwd('/home/ukb')
load('ukb_famdat.Rdata')
gamma = read.csv('gamma.csv')[,-1]
se = read.csv('segamma.csv')[,-1]
familyID = familyID[1:N]
familyID = c(familyID, max(familyID)+1:(size-N))
N = size

kinship = rep(1,N)
fID1 = sapply(1:nrow(pairs), function(l) familyID[which(eid==pairs[l,1])])
fID2 = sapply(1:nrow(pairs), function(l) familyID[which(eid==pairs[l,2])])
eli = fID1==fID2
pairs = pairs[eli,]
P1 = c(pairs[,1], pairs[,2])
P2 = c(pairs[,2], pairs[,1])
fIDco = fID1[eli]
fIDco = c(familyID, fIDco, fIDco)
kinco = c(kinship, pairs[,3], pairs[,3])

rg = 1:10 + (b-1)*10
rg = rg[rg<=290]
res = matrix(NA,nrow=10,ncol=10)

fiti = fitj = NULL
er = 0
if (i<=K1){
  load('condat_all.Rdata')
  X0 = cbind(1,as.matrix(X))
  Y_i = con.dat[,i]
  rm(con.dat)
  Yi = c(Y_i, 
         sapply(1:length(P1), function(l) Y_i[which(eid==P1[l])]))
  Xi = rbind(X0, 
             t(sapply(1:length(P1), function(l) X0[which(eid==P1[l]),])))
}
if (i>K1&i<=K1+K2){
  load('bindat_all.Rdata')
  Y_i = bin.dat[,i-K1]+1
  X = as.matrix(X)
  rm(bin.dat)
  Yi = c(Y_i, 
         sapply(1:length(P1), function(l) Y_i[which(eid==P1[l])]))
  Xi = rbind(X, 
             t(sapply(1:length(P1), function(l) X[which(eid==P1[l]),])))
}
if (i>K1+K2&i<=K1+K2+K3){
  load('orddat_all.Rdata')
  Y_i = ord.dat[,i-K1-K2]
  X = as.matrix(X)
  rm(ord.dat)
  Yi = c(Y_i, 
         sapply(1:length(P1), function(l) Y_i[which(eid==P1[l])]))
  Xi = rbind(X, 
             t(sapply(1:length(P1), function(l) X[which(eid==P1[l]),])))
}
if (i>K1+K2+K3){
  load('surdat_all.Rdata')
  T_i = sur.T.dat[,i-K1-K2-K3]
  D_i = sur.D.dat[,i-K1-K2-K3]
  X_sur = as.matrix(X)[,-2]
  rm(sur.T.dat)
  rm(sur.D.dat)
  Ti = c(T_i, 
         sapply(1:length(P1), function(l) T_i[which(eid==P1[l])]))
  Di = c(D_i, 
         sapply(1:length(P1), function(l) D_i[which(eid==P1[l])]))
  Xi = rbind(X_sur, 
             t(sapply(1:length(P1), function(l) X_sur[which(eid==P1[l]),])))
}
p = ncol(X)
setwd('/home/ukb/results_alldat')
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

for (jb in 1:10){
  j = rg[jb]
  if (i>=j) next
  setwd('/home/ukb/results_alldat')
  t0 = Sys.time()
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
  setwd('/home/ukb')
  if (j<=K1){
    load('condat_all.Rdata')
    Y_j = con.dat[,j]
    X0 = cbind(1,as.matrix(X))
    rm(con.dat)
    Yj = c(Y_j, 
           sapply(1:length(P2), function(l) Y_j[which(eid==P2[l])]))
    Xj = rbind(X0, 
               t(sapply(1:length(P2), function(l) X0[which(eid==P2[l]),])))
  }
  if (j>K1&j<=K1+K2){
    load('bindat_all.Rdata')
    Y_j = bin.dat[,j-K1]+1
    X = as.matrix(X)
    rm(bin.dat)
    Yj = c(Y_j, 
           sapply(1:length(P2), function(l) Y_j[which(eid==P2[l])]))
    Xj = rbind(X, 
               t(sapply(1:length(P2), function(l) X[which(eid==P2[l]),])))
  }
  if (j>K1+K2&j<=K1+K2+K3){
    load('orddat_all.Rdata')
    Y_j = ord.dat[,j-K1-K2]
    X = as.matrix(X)
    rm(ord.dat)
    Yj = c(Y_j, 
           sapply(1:length(P2), function(l) Y_j[which(eid==P2[l])]))
    Xj = rbind(X, 
               t(sapply(1:length(P2), function(l) X[which(eid==P2[l]),])))
  }
  if (j>K1+K2+K3){
    load('surdat_all.Rdata')
    T_j = sur.T.dat[,j-K1-K2-K3]
    D_j = sur.D.dat[,j-K1-K2-K3]
    X_sur = as.matrix(X[,-2])
    rm(sur.T.dat)
    rm(sur.D.dat)
    Tj = c(T_j, 
           sapply(1:length(P2), function(l) T_j[which(eid==P2[l])]))
    Dj = c(D_j, 
           sapply(1:length(P2), function(l) D_j[which(eid==P2[l])]))
    Xj = rbind(X_sur, 
               t(sapply(1:length(P2), function(l) X_sur[which(eid==P2[l]),])))
  }
  
  setwd('/home/ukb/results_alldat_co/')
  print(i)
  print(j)
  if (i<=K1){
    if (j<=K1){
      mis = complete.cases(cbind(Yi,Yj))
      g12 = sqrt(fiti$gamma*fitj$gamma)
      par0 = c(gamma[j,i],max(gamma[j,i]-3*se[j,i],-g12),min(gamma[j,i]+3*se[j,i],g12))
      resij = try(con_con(Yi[mis],Yj[mis],Xi[mis,],Xj[mis,],fIDco[mis],
                          fiti,fitj,kinco[mis],par0))
      if ('try-error'%in%class(resij)) next
    }
    if (j>K1&j<=K1+K2){
      mis = complete.cases(cbind(Yi,Yj))
      g12 = sqrt(fiti$gamma*fitj$gamma)
      par0 = c(gamma[j,i],max(gamma[j,i]-3*se[j,i],-g12),min(gamma[j,i]+3*se[j,i],g12))
      resij = try(con_ord(Yi[mis],Yj[mis],Xi[mis,],Xj[mis,],fIDco[mis],
                          fiti,fitj,kinco[mis],par0))
      if ('try-error'%in%class(resij)) next
    }
    if (j>K1+K2&j<=K1+K2+K3){
      mis = complete.cases(cbind(Yi,Yj))
      g12 = sqrt(fiti$gamma*fitj$gamma)
      par0 = c(gamma[j,i],max(gamma[j,i]-3*se[j,i],-g12),min(gamma[j,i]+3*se[j,i],g12))
      resij = try(con_ord(Yi[mis],Yj[mis],Xi[mis,],Xj[mis,],fIDco[mis],
                          fiti,fitj,kinco[mis],par0))
      if ('try-error'%in%class(resij)) next
    }
    if (j>K1+K2+K3){
      mis = complete.cases(cbind(Yi,Tj,Dj))
      fitj$h = fitj$h[1:(p+1)]
      g12 = sqrt(fiti$gamma*fitj$gamma)
      par0 = c(gamma[j,i],max(gamma[j,i]-3*se[j,i],-g12),min(gamma[j,i]+3*se[j,i],g12))
      resij = try(con_sur(Yi[mis],Tj[mis],Dj[mis],Xi[mis,],Xj[mis,],fIDco[mis],
                          fiti,fitj,kinco[mis],par0))
      if ('try-error'%in%class(resij)) next
    }
  }
  
  if (i>K1&i<=K1+K2){
    if (j>K1&j<=K1+K2){
      mis = complete.cases(cbind(Yi,Yj))
      g12 = sqrt(fiti$gamma*fitj$gamma)
      par0 = c(gamma[j,i],max(gamma[j,i]-3*se[j,i],-g12),min(gamma[j,i]+3*se[j,i],g12))
      resij = try(ord_ord(Yi[mis],Yj[mis],Xi[mis,],Xj[mis,],fIDco[mis],
                          fiti,fitj,kinco[mis],par0))
      if ('try-error'%in%class(resij)) next
    }
    if (j>K1+K2&j<=K1+K2+K3){
      mis = complete.cases(cbind(Yi,Yj))
      g12 = sqrt(fiti$gamma*fitj$gamma)
      par0 = c(gamma[j,i],max(gamma[j,i]-3*se[j,i],-g12),min(gamma[j,i]+3*se[j,i],g12))
      resij = try(ord_ord(Yi[mis],Yj[mis],Xi[mis,],Xj[mis,],fIDco[mis],
                          fiti,fitj,kinco[mis],par0))
      if ('try-error'%in%class(resij)) next
    }
    if (j>K1+K2+K3){
      mis = complete.cases(cbind(Yi,Tj,Dj))
      fitj$h = fitj$h[1:(p+1)]
      g12 = sqrt(fiti$gamma*fitj$gamma)
      par0 = c(gamma[j,i],max(gamma[j,i]-3*se[j,i],-g12),min(gamma[j,i]+3*se[j,i],g12))
      resij = try(ord_sur(Yi[mis],Tj[mis],Dj[mis],Xi[mis,],Xj[mis,],fIDco[mis],
                          fiti,fitj,kinco[mis],par0))
      if ('try-error'%in%class(resij)) next
    }
  }
  
  if (i>K1+K2&i<=K1+K2+K3){
    if (j>K1+K2&j<=K1+K2+K3){
      mis = complete.cases(cbind(Yi,Yj))
      g12 = sqrt(fiti$gamma*fitj$gamma)
      par0 = c(gamma[j,i],max(gamma[j,i]-3*se[j,i],-g12),min(gamma[j,i]+3*se[j,i],g12))
      resij = try(ord_ord(Yi[mis],Yj[mis],Xi[mis,],Xj[mis,],fIDco[mis],
                          fiti,fitj,kinco[mis],par0))
      if ('try-error'%in%class(resij)) next
    }
    if (j>K1+K2+K3){
      mis = complete.cases(cbind(Yi,Tj,Dj))
      fitj$h = fitj$h[1:(p+1)]
      g12 = sqrt(fiti$gamma*fitj$gamma)
      par0 = c(gamma[j,i],max(gamma[j,i]-3*se[j,i],-g12),min(gamma[j,i]+3*se[j,i],g12))
      resij = try(ord_sur(Yi[mis],Tj[mis],Dj[mis],Xi[mis,],Xj[mis,],fIDco[mis],
                          fiti,fitj,kinco[mis],par0))
      if ('try-error'%in%class(resij)) next
    }
  }
  if (i>K1+K2+K3){
    if (j>K1+K2+K3){
      mis = complete.cases(cbind(Ti,Di,Tj,Dj))
      fiti$h = fiti$h[1:(p+1)]
      fitj$h = fitj$h[1:(p+1)]
      g12 = sqrt(fiti$gamma*fitj$gamma)
      par0 = c(gamma[j,i],max(gamma[j,i]-3*se[j,i],-g12),min(gamma[j,i]+3*se[j,i],g12))
      resij = try(sur_sur(Ti[mis],Di[mis],Tj[mis],Dj[mis],Xi[mis,],Xj[mis,],fIDco[mis],
                          fiti,fitj,kinco[mis],par0))
      if ('try-error'%in%class(resij)) next
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
  t1 = Sys.time()
  tm = difftime(t1,t0, units='secs')
  res[jb,] = c(rho.f, se.rho.f, zeta.f, se.zeta.f, 
               rho, se.rho, zeta, se.zeta,
               sum(mis), tm)
  save(res, file=paste0('res',i,'_',b,'.Rdata'))
}

save(res, file=paste0('res',i,'_',b,'.Rdata'))
