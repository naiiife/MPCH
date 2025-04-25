## Continuous, additional simulation

lid = Sys.getenv("SLURM_ARRAY_TASK_ID")
if (lid=="") lid=0
set.seed(lid)
source('head.R')
source('generatedata.R')
source('con_con.R')

# original method for continuous outcomes
cal_con <- function(Y,X,familyID,KM,par0){
  familyID = familyID + 1
  N = length(Y)
  p = ncol(X)
  Fam_list = sort(unique(familyID))
  Fam_Num = length(Fam_list)
  Fam_size = sapply(Fam_list, function(x) sum(familyID==x))
  ngrid = length(w_list)
  if (p==1) X = t(X)
  HX = solve(t(X) %*% X) %*% t(X)
  N = length(Y)
  Fam_list = sort(unique(familyID))
  Fam_Num = length(Fam_list)
  Fam_size = sapply(Fam_list, function(x) sum(familyID==x))
  
  # Ee, Eeps, Ee2, Eeps2, Eeeps, EepsGeps are vectors of N*1
  theta_k = par0[1]
  gamma_kk = par0[2]
  alpha_k = par0[3:(2+p)]
  sigmau2 = par0[3+p]
  iterrun = 0
  while (TRUE){
    # E step
    Ee = Eeps = Eeeps = Ee2 = Eeps2 = EepsGeps = NULL
    theta_k0 = theta_k; gamma_kk0 = gamma_kk; alpha_k0 = alpha_k; sigmau20 = sigmau2
    for (i in 1:Fam_Num){
      ni = Fam_size[i]
      fi = Fam_list[i]
      Yi = Y[familyID==fi]
      Xi = X[familyID==fi,]
      if (ni==1) Xi=t(Xi)
      Gi = as.matrix(KM[familyID==fi,1:ni])
      Ji = matrix(1, ni, ni)
      Oi = rep(1, ni)
      Ii = diag(Oi)
      Giinv = ginv(Gi)
      Si = ginv(theta_k^2*Ji + gamma_kk*Gi + sigmau2*Ii)
      Y_Xb = Yi - Xi %*% alpha_k
      Eei = as.numeric(theta_k^2 * t(Oi) %*% Si %*% Y_Xb) # 1x1
      Eepsi = gamma_kk * Gi %*% Si %*% Y_Xb          # nix1
      Vei = as.numeric(theta_k^2 - theta_k^4 * t(Oi) %*% Si %*% Oi)   # 1x1
      Vepsi = gamma_kk * Gi - gamma_kk^2 * Gi %*% Si %*% Gi    # nixni
      Eei2 = Eei^2 + Vei                                       # 1x1
      Eepsi2 = as.numeric(sum(Eepsi^2) + sum(diag(Vepsi)))     # 1x1
      Eeepsi = Eei*Eepsi - theta_k^2 * gamma_kk * Gi %*% Si %*% Oi    # nix1
      EepsGepsi = as.numeric(t(Eepsi) %*% Giinv %*% Eepsi +
                               sum(diag(Giinv %*% Vepsi)))     # 1x1
      
      Ee = append(Ee, rep(Eei, ni))
      Eeps = append(Eeps, Eepsi)
      Ee2 = append(Ee2, Eei2)
      Eeps2 = append(Eeps2, Eepsi2)
      Eeeps = append(Eeeps, Eeepsi)
      EepsGeps = append(EepsGeps, EepsGepsi)
    }
    # M step
    Y_Xb = as.numeric(Y - X %*% alpha_k)
    sigmau2 = (sum(Y_Xb^2) - 2*sum(Y_Xb * (Ee+Eeps))) / N + 
      (sum(Ee2*Fam_size) + 2*sum(Eeeps) + sum(Eeps2)) / N
    Y_ee = Y - Ee - Eeps
    alpha_k = HX %*% Y_ee
    theta_k = sqrt(sum(Ee2) / Fam_Num)
    gamma_kk = sum(EepsGeps) / N
    tol = max(max(abs(alpha_k-alpha_k0)), sum(abs(theta_k-theta_k0)),
              sum(abs(gamma_kk-gamma_kk0)))
    iterrun = iterrun + 1
    if(iterrun>=5000) tol = 0
    if (tol < 0.0001) break
    print(c(theta_k,gamma_kk,alpha_k,sigmau2,tol))
  }
  est = c(theta_k,gamma_kk,alpha_k,sigmau2)
  phi = matrix(0,Fam_Num,p+4)
  phi[,1] = Fam_list
  res = list(theta=theta_k, gamma=gamma_kk, alpha=as.numeric(alpha_k), 
             sigmau2=sigmau2, iterrun=iterrun, est=est, phi=phi)
  return(res)
}

## no environmental effect
cal_con_ne <- function(Y,X,familyID,KM,par0){
  familyID = familyID + 1
  N = length(Y)
  p = ncol(X)
  Fam_list = sort(unique(familyID))
  Fam_Num = length(Fam_list)
  Fam_size = sapply(Fam_list, function(x) sum(familyID==x))
  ngrid = length(w_list)
  if (p==1) X = t(X)
  HX = solve(t(X) %*% X) %*% t(X)
  N = length(Y)
  Fam_list = sort(unique(familyID))
  Fam_Num = length(Fam_list)
  Fam_size = sapply(Fam_list, function(x) sum(familyID==x))
  
  # Ee, Eeps, Ee2, Eeps2, Eeeps, EepsGeps are vectors of N*1
  theta_k = 0
  gamma_kk = par0[2]
  alpha_k = par0[3:(2+p)]
  sigmau2 = par0[3+p]
  iterrun = 0
  while (TRUE){
    # E step
    Ee = Eeps = Eeeps = Ee2 = Eeps2 = EepsGeps = NULL
    theta_k0 = theta_k; gamma_kk0 = gamma_kk; alpha_k0 = alpha_k; sigmau20 = sigmau2
    for (i in 1:Fam_Num){
      ni = Fam_size[i]
      fi = Fam_list[i]
      Xi = X[familyID==fi,]
      if (ni==1) Xi=t(Xi)
      Yi = Y[familyID==fi]
      Gi = as.matrix(KM[familyID==fi,1:ni])
      Ji = matrix(1, ni, ni)
      Oi = rep(1, ni)
      Ii = diag(Oi)
      Giinv = ginv(Gi)
      Si = ginv(theta_k^2*Ji + gamma_kk*Gi + sigmau2*Ii)
      Y_Xb = Yi - Xi %*% alpha_k
      Eei = as.numeric(theta_k^2 * t(Oi) %*% Si %*% Y_Xb) # 1x1
      Eepsi = gamma_kk * Gi %*% Si %*% Y_Xb          # nix1
      Vei = as.numeric(theta_k^2 - theta_k^4 * t(Oi) %*% Si %*% Oi)   # 1x1
      Vepsi = gamma_kk * Gi - gamma_kk^2 * Gi %*% Si %*% Gi    # nixni
      Eei2 = Eei^2 + Vei                                       # 1x1
      Eepsi2 = as.numeric(sum(Eepsi^2) + sum(diag(Vepsi)))     # 1x1
      Eeepsi = Eei*Eepsi - theta_k^2 * gamma_kk * Gi %*% Si %*% Oi    # nix1
      EepsGepsi = as.numeric(t(Eepsi) %*% Giinv %*% Eepsi +
                               sum(diag(Giinv %*% Vepsi)))     # 1x1
      
      Ee = append(Ee, rep(Eei, ni))
      Eeps = append(Eeps, Eepsi)
      Ee2 = append(Ee2, Eei2)
      Eeps2 = append(Eeps2, Eepsi2)
      Eeeps = append(Eeeps, Eeepsi)
      EepsGeps = append(EepsGeps, EepsGepsi)
    }
    # M step
    Y_Xb = as.numeric(Y - X %*% alpha_k)
    sigmau2 = (sum(Y_Xb^2) - 2*sum(Y_Xb * (Ee+Eeps))) / N + 
      (sum(Ee2*Fam_size) + 2*sum(Eeeps) + sum(Eeps2)) / N
    Y_ee = Y - Ee - Eeps
    alpha_k = HX %*% Y_ee
    theta_k = 0
    gamma_kk = sum(EepsGeps) / N
    tol = max(max(abs(alpha_k-alpha_k0)), sum(abs(theta_k-theta_k0)),
              sum(abs(gamma_kk-gamma_kk0)))
    iterrun = iterrun + 1
    if(iterrun>=5000) tol = 0
    if (tol < 0.0001) break
    print(c(theta_k,gamma_kk,alpha_k,sigmau2,tol))
  }
  est = c(theta_k,gamma_kk,alpha_k,sigmau2)
  phi = matrix(0,Fam_Num,p+4)
  phi[,1] = Fam_list
  res = list(theta=theta_k, gamma=gamma_kk, alpha=as.numeric(alpha_k), 
             sigmau2=sigmau2, iterrun=iterrun, est=est, phi=phi)
  return(res)
}

## joint MLE
con_joint <- function(Y1,Y2,X,familyID,KM){
  familyID = familyID + 1
  N = length(Y1)
  X = as.matrix(X)
  p = ncol(X)
  Fam_list = sort(unique(familyID))
  Fam_Num = length(Fam_list)
  Fam_size = sapply(Fam_list, function(x) sum(familyID==x))
  ngrid = length(w_list)
  N = length(Y1)
  iterrun = 0
  alpha_1 = alpha_2 = rep(0,p)
  theta_1 = 0.8; theta_2 = 0.6
  gamma_11 = 0.8; gamma_22 = 0.49; rho = 0.5
  sigmau1 = 1; sigmau2 = 0.64
  gamma_12 = rho*sqrt(gamma_11*gamma_22)
  tol = 1
  est_all = rep(0,2*p+7)
  while (tol>0.0001){
    Xs = NULL; Ys = NULL
    est0_all = est_all
    for (i in 1:Fam_Num){
      ni = Fam_size[i]
      fi = Fam_list[i]
      Y1i = Y1[familyID==fi]
      Y2i = Y2[familyID==fi]
      Xi = as.matrix(X[familyID==fi,])
      Gi = as.matrix(KM[familyID==fi,1:ni])
      Ji = matrix(1, ni, ni)
      Ii = diag(rep(1,ni))
      S1 = theta_1^2*Ji+gamma_11*Gi+sigmau1*Ii
      S2 = theta_2^2*Ji+gamma_22*Gi+sigmau2*Ii
      S3 = theta_1*theta_2*Ji+gamma_12*Gi
      S = rbind(cbind(S1,S3),cbind(S3,S2))
      svds = svd(S)
      U = svds$u; B = t(U); D = svds$d
      B11 = B[1:ni,1:ni]; B22 = B[ni+1:ni,ni+1:ni]
      B12 = B[1:ni,ni+1:ni]; B21 = B[ni+1:ni,1:ni]
      Yis = B%*%c(Y1i,Y2i)/sqrt(D)
      X1s = rbind(B11%*%Xi,B21%*%Xi)/sqrt(D)
      X2s = rbind(B12%*%Xi,B22%*%Xi)/sqrt(D)
      Xis = cbind(X1s,X2s)
      Ys = append(Ys,Yis)
      Xs = rbind(Xs,Xis)
    }
    alpha = ginv(t(Xs)%*%Xs)%*%t(Xs)%*%Ys
    alpha_1 = alpha[1:p]
    alpha_2 = alpha[p+1:p]
    
    cal_loglik = function(par){
      theta_1 = par[1]; theta_2 = par[2]
      gamma_11 = par[3]; gamma_22 = par[4]; rho = par[5]
      sigmau1=par[6]; sigmau2=par[7]
      gamma_12 = rho*sqrt(gamma_11*gamma_22)
      loglik = 0
      for (i in 1:Fam_Num){
        ni = Fam_size[i]
        fi = Fam_list[i]
        Y1i = Y1[familyID==fi]
        Y2i = Y2[familyID==fi]
        Xi = as.matrix(X[familyID==fi,])
        Gi = as.matrix(KM[familyID==fi,1:ni])
        Ji = matrix(1, ni, ni)
        Ii = diag(rep(1,ni))
        S1 = theta_1^2*Ji+gamma_11*Gi+sigmau1*Ii
        S2 = theta_2^2*Ji+gamma_22*Gi+sigmau2*Ii
        S3 = theta_1*theta_2*Ji+gamma_12*Gi
        S = rbind(cbind(S1,S3),cbind(S3,S2))
        ui = c(Y1i-Xi%*%alpha_1,Y2i-Xi%*%alpha_2)
        detS = det(S)
        if (detS==0) detS = 0.000001
        li = -1/2*log(detS)-1/2*as.numeric(t(ui)%*%ginv(S)%*%ui)
        loglik = loglik + li
      }
      return(-loglik)
    }
    est = c(theta_1,theta_2,gamma_11,gamma_22,rho,sigmau1,sigmau2)
    opt = optim(est,cal_loglik,method='L-BFGS-B',
                lower=c(0,0,0,0,-1,0,0),upper=c(Inf,Inf,Inf,Inf,1,Inf,Inf),
                control=list(maxit=10))
    par = opt$par
    theta_1 = par[1]; theta_2 = par[2]
    gamma_11 = par[3]; gamma_22 = par[4]; rho = par[5]
    sigmau1=par[6]; sigmau2=par[7]
    gamma_12 = rho*sqrt(gamma_11*gamma_22)
    est_all = c(par, alpha)
    tol = max(abs(est_all-est0_all))
    iterrun = iterrun + 1
    print(c(par,tol))
    if(iterrun>=5000) break
  }
  
  S1 = theta_1^2+gamma_11+sigmau1
  S2 = theta_2^2+gamma_22+sigmau2
  zeta1 = theta_1^2/S1; zeta2 = theta_2^2/S2
  rho1 = gamma_11/S1; rho2 = gamma_22/S2
  zeta12 = theta_1*theta_2/sqrt(S1*S2)
  rho12 = gamma_12/sqrt(S1*S2)
  res = list(theta_1=theta_1, theta_2=theta_2, gamma_11=gamma_11, 
             gamma_22=gamma_22, gamma12=gamma_12, alpha_1=alpha_1,
             alpha_2=alpha_2,sigmau1=sigmau1,sigmau2=sigmau2, 
             rho1=rho1,rho2=rho2,zeta1=zeta1,zeta2=zeta2,
             rho12=rho12,zeta12=zeta12,
             iterrun=iterrun)
  return(res)
}

## HEc
con_hec <- function(Y1,Y2,X,familyID,KM){
  Y1.p = predict(lm(Y1~X), data.frame(Y1,X))
  f1 = Y1 - Y1.p
  var1 = var(f1, na.rm=TRUE)
  Y2.p = predict(lm(Y2~X), data.frame(Y2,X))
  f2 = Y2 - Y2.p
  var2 = var(f2, na.rm=TRUE)
  v12 = v1 = v2 = trk = 0
  Fam_list = unique(familyID)
  for (m in Fam_list){
    ind = (familyID==m)
    if (sum(ind)==0) next
    Ki = KM[ind,]
    Ki = Ki[,1:nrow(Ki)]
    f1i = f1[ind]
    f2i = f2[ind]
    diag(Ki) = 0
    v12 = v12 + t(f1i)%*% Ki %*%f2i
    v1 = v1 + t(f1i)%*% Ki %*%f1i
    v2 = v2 + t(f2i)%*% Ki %*%f2i
    trk = trk + sum(diag(Ki%*%Ki))
  }
  rho1 = as.numeric(v1/trk/var1)
  rho2 = as.numeric(v2/trk/var2)
  rho12 = as.numeric(v12/trk/sqrt(var1*var2))
  return(list(rho1=rho1,rho2=rho2,rho12=rho12))
}


t0 = Sys.time()
res1_con = cal_con(Y1_con,X,familyID,KM,true1_con)
res2_con = cal_con(Y2_con,X,familyID,KM,true2_con)
kinship = rep(1,length(familyID))
res12 = con_con(Y1_con,Y2_con,X,X,familyID,res1_con,res2_con,kinship)
t1 = Sys.time()
t.1 = as.numeric(difftime(t1, t0, units = "secs"))
res.mpch = c(res12$rho1,res12$rho2,res12$zeta1,res12$zeta2,res12$rho12,res12$zeta12,t.1)
res = con_joint(Y1_con,Y2_con,X,familyID,KM)
t2 = Sys.time()
t.2 = as.numeric(difftime(t2, t1, units = "secs"))
res.joint = c(res$rho1,res$rho2,res$zeta1,res$zeta2,res$rho12,res$zeta12,t.2)
res = con_hec(Y1_con,Y2_con,X,familyID,KM)
res.hec = c(res$rho1,res$rho2,res$rho12)
res1_con = cal_con_ne(Y1_con,X,familyID,KM,true1_con)
res2_con = cal_con_ne(Y2_con,X,familyID,KM,true2_con)
kinship = rep(1,length(familyID))
res12 = con_con(Y1_con,Y2_con,X,X,familyID,res1_con,res2_con,kinship)
t3 = Sys.time()
t.3 = as.numeric(difftime(t3, t2, units = "secs"))
res.ne = c(res12$rho1,res12$rho2,res12$rho12,t.3)

save(res.mpch, res.joint, res.hec, res.ne,
     file=paste0('res',lid,'.Rdata'))


## analysis
# library(forestploter)
# setwd('...')
# est1 = est2 = est3 = est4 = NULL
# for (i in 1:1000){
#   fit = try(load(paste0('res',i,'.Rdata')), silent=TRUE)
#   if ('try-error' %in% class(fit)) next
#   est1 = rbind(est1, res.mpch)
#   est2 = rbind(est2, res.joint)
#   est3 = rbind(est3, res.hec)
#   est4 = rbind(est4, res.ne)
# }
# mean.rho1.1 = mean(est1[,1], na.rm=TRUE)
# mean.rho2.1 = mean(est1[,2], na.rm=TRUE)
# mean.zeta1.1 = mean(est1[,3], na.rm=TRUE)
# mean.zeta2.1 = mean(est1[,4], na.rm=TRUE)
# mean.rho12.1 = mean(est1[,5], na.rm=TRUE)
# mean.zeta12.1 = mean(est1[,6], na.rm=TRUE)
# mean.rho1.2 = mean(est2[,1], na.rm=TRUE)
# mean.rho2.2 = mean(est2[,2], na.rm=TRUE)
# mean.zeta1.2 = mean(est2[,3], na.rm=TRUE)
# mean.zeta2.2 = mean(est2[,4], na.rm=TRUE)
# mean.rho12.2 = mean(est2[,5], na.rm=TRUE)
# mean.zeta12.2 = mean(est2[,6], na.rm=TRUE)
# mean.rho1.3 = mean(est3[,1], na.rm=TRUE)
# mean.rho2.3 = mean(est3[,2], na.rm=TRUE)
# mean.rho12.3 = mean(est3[,3], na.rm=TRUE)
# mean.rho1.4 = mean(est4[,1], na.rm=TRUE)
# mean.rho2.4 = mean(est4[,2], na.rm=TRUE)
# mean.rho12.4 = mean(est4[,3], na.rm=TRUE)
# sd.rho1.1 = sd(est1[,1], na.rm=TRUE)
# sd.rho2.1 = sd(est1[,2], na.rm=TRUE)
# sd.zeta1.1 = sd(est1[,3], na.rm=TRUE)
# sd.zeta2.1 = sd(est1[,4], na.rm=TRUE)
# sd.rho12.1 = sd(est1[,5], na.rm=TRUE)
# sd.zeta12.1 = sd(est1[,6], na.rm=TRUE)
# sd.rho1.2 = sd(est2[,1], na.rm=TRUE)
# sd.rho2.2 = sd(est2[,2], na.rm=TRUE)
# sd.zeta1.2 = sd(est2[,3], na.rm=TRUE)
# sd.zeta2.2 = sd(est2[,4], na.rm=TRUE)
# sd.rho12.2 = sd(est2[,5], na.rm=TRUE)
# sd.zeta12.2 = sd(est2[,6], na.rm=TRUE)
# sd.rho1.3 = sd(est3[,1], na.rm=TRUE)
# sd.rho2.3 = sd(est3[,2], na.rm=TRUE)
# sd.rho12.3 = sd(est3[,3], na.rm=TRUE)
# sd.rho1.4 = sd(est3[,1], na.rm=TRUE)
# sd.rho2.4 = sd(est3[,2], na.rm=TRUE)
# sd.rho12.4 = sd(est3[,3], na.rm=TRUE)
# sd.h1.1 = sd(est1[,1]+est1[,3], na.rm=TRUE)
# sd.h2.1 = sd(est1[,2]+est1[,4], na.rm=TRUE)
# sd.h12.1 = sd(est1[,5]+est1[,6], na.rm=TRUE)
# sd.h1.2 = sd(est2[,1]+est2[,3], na.rm=TRUE)
# sd.h2.2 = sd(est2[,2]+est2[,4], na.rm=TRUE)
# sd.h12.2 = sd(est2[,5]+est2[,6], na.rm=TRUE)
# t1 = mean(est1[,7],rm.na=TRUE)
# t2 = mean(est2[,7],rm.na=TRUE)
# true = c(0.44444444,0.52127660,0.11111111,0.09574468,0.240665,0.10314212)
# S1 = true1_con[1]^2+true1_con[2]+true1_con[length(true1_con)]]
# S2 = true2_con[1]^2+true2_con[2]+true2_con[length(true2_con)]]
# true[1] = true1_con[2]/S1
# true[2] = true2_con[2]/S2
# true[3] = true1_con[1]^2/S1
# true[4] = true2_con[1]^2/S2
# true[5] = true1_con[2]*true2_con[2]*0.5/sqrt(S1*S2)
# true[6] = true1_con[1]*true2_con[1]/sqrt(S1*S2)
# truej = c(true[1]+true[3],true[2]+true[4],true[5]+true[6])
# est1 = c(mean.rho1.1,mean.rho2.1,mean.zeta1.1,mean.zeta2.1,
#          mean.rho12.1,mean.zeta12.1) - true
# est2 = c(mean.rho1.2,mean.rho2.2,mean.zeta1.2,mean.zeta2.2,
#          mean.rho12.2,mean.zeta12.2) - true
# est3 = c(mean.rho1.3,mean.rho2.3,mean.rho12.3) - true[c(1,2,5)]
# est4 = c(mean.rho1.4,mean.rho2.4,mean.rho12.4) - true[c(1,2,5)]
# est5 = c(mean.rho1.1+mean.zeta1.1,mean.rho2.1+mean.zeta2.1,mean.rho12.1+mean.zeta12.1)-truej
# est6 = c(mean.rho1.2+mean.zeta1.2,mean.rho2.2+mean.zeta2.2,mean.rho12.2+mean.zeta12.2)-truej
# sd1 = c(sd.rho1.1,sd.rho2.1,sd.zeta1.1,sd.zeta2.1,
#          sd.rho12.1,sd.zeta12.1)
# sd2 = c(sd.rho1.2,sd.rho2.2,sd.zeta1.2,sd.zeta2.2,
#         sd.rho12.2,sd.zeta12.2)
# sd3 = c(sd.rho1.3,sd.rho2.3,sd.rho12.3)
# sd4 = c(sd.rho1.4,sd.rho2.4,sd.rho12.4)
# sd5 = c(sd.h1.1,sd.h2.1,sd.h12.1)
# sd6 = c(sd.h1.2,sd.h2.2,sd.h12.2)
# 
# dt = data.frame(matrix(0,6,7))
# colnames(dt) = c('Parameter','est1','est2','low1','up1','low2','up2')
# dt$Parameter = c('Heritability 1', 'Heritability 2',
#            'Environ eff 1', 'Environ eff 2',
#            'Coheritability', 'Environ corr')
# dt$est1 = est1
# dt$est2 = est2
# dt$low1 = dt$est1 - 1.96*sd1
# dt$up1 = dt$est1 + 1.96*sd1
# dt$low2 = dt$est2 - 1.96*sd2
# dt$up2 = dt$est2 + 1.96*sd2
# dt$' ' = paste(rep(" ",20), collapse=" ")
# tm <- forest_theme(
#   legend_name = 'Method',
#   legend_value = c('MPCH','Joint'))
# rhozeta_j <- forest(dt[,c(1,8)],
#                  est = list(dt$est1, dt$est2),
#                  lower = list(dt$low1, dt$low2),
#                  upper = list(dt$up1, dt$up2),
#                  xlim = c(-.2,.2),
#                  ci_column = 2,
#                  vert_line = 0,
#                  nudge_y = 0.2,
#                  ticks_at = c(-.2,-.1,0,.1,.2),
#                  title = 'Bias of estimates',
#                  theme = tm)
# rhozeta_j
# 
# dt = data.frame(matrix(0,3,10))
# colnames(dt) = c('Parameter','est1','est2','est3','low1','up1','low2','up2','low3','up3')
# dt$Parameter = c('Heritability 1', 'Heritability 2', 'Coheritability')
# est1 = est1[c(1,2,5)]; sd1 = sd1[c(1,2,5)]
# dt$est1 = est1
# dt$est2 = est4
# dt$est3 = est3
# dt$low1 = est1 - 1.96*sd1
# dt$up1 = est1 + 1.96*sd1
# dt$low2 = est4 - 1.96*sd4
# dt$up2 = est4 + 1.96*sd4
# dt$low3 = est3 - 1.96*sd3
# dt$up3 = est3 + 1.96*sd3
# dt$' ' = paste(rep(" ",20), collapse=" ")
# tm <- forest_theme(
#   legend_name = 'Method',
#   legend_value = c('MPCH','MPCH_NE','HEc'))
# h_e <- forest(dt[,c(1,11)],
#                  est = list(dt$est1, dt$est2, dt$est3),
#                  lower = list(dt$low1, dt$low2, dt$low3),
#                  upper = list(dt$up1, dt$up2, dt$up3),
#                  xlim = c(-.1,.7),
#                  ci_column = 2,
#                  vert_line = 0,
#                  nudge_y = 0.2,
#                  ticks_at = c(-.1,.1,.3,.5,.7),
#                  title = 'Bias of estimates',
#                  theme = tm)
# h_e