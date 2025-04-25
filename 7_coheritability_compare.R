setwd('/home/ukb')
load('condat.Rdata')
load('ukb_famdat.Rdata')
X0 = cbind(1,as.matrix(X[1:N,]))
n = n2+n3+n4
K1 = ncol(con.dat)
con.dat = con.dat[1:N,]
KM = KM[1:N,]
familyID = familyID[1:N]

#HEc
setwd('/home/ukb/family_data_all')
rho_hec = genc_hec = se.rho_hec = matrix(NA,K1,K1)
for (i in 1:K1){
  Y_i = con.dat[,i]
  Y1.p = predict(lm(Y_i~0+X0), data.frame(X0))
  f1 = Y_i - Y1.p
  var1 = var(f1, na.rm=TRUE)
  for (j in i:K1){
    if (i>=j) next
    Y_j = con.dat[,j]
    Y2.p = predict(lm(Y_j~0+X0), data.frame(X0))
    f2 = Y_j - Y2.p
    var2 = var(f2, na.rm=TRUE)
    v12 = v1 = v2 = tr = trk = 0
    for (m in 1:n-1){
      ind = (familyID==m)
      if (sum(ind)==0) next
      Ki = KM[ind,]
      Ki = Ki[,1:nrow(Ki)]
      f1i = f1[ind]
      f2i = f2[ind]
      ind = (!is.na(f1i))&(!is.na(f2i))
      diag(Ki) = 0
      v12 = v12 + t(f1i[ind])%*% Ki[ind,ind] %*%f2i[ind]
      v1 = v1 + t(f1i[ind])%*% Ki[ind,ind] %*%f1i[ind]
      v2 = v2 + t(f2i[ind])%*% Ki[ind,ind] %*%f2i[ind]
      trk = trk + sum(diag(Ki[ind,ind]%*%Ki[ind,ind]))
    }
    r = v12/trk/sqrt(var1*var2)
    rho_hec[j,i] = r
    genc_hec[j,i] = v12/sqrt(abs(v1*v2))
    mu = 0.5*log((1+r)/(1-r))
    se = 1/sqrt(abs(sqrt(abs(trk))-3))
    se.rho_hec[j,i] = 4*exp(2*mu)/(exp(2*mu)+1)^2*se
  }
  cat(i, '')
}
#save(rho_hec,se.rho_hec, file='HEc_coheritability.Rdata')
load('HEc_coheritability.Rdata')

phe1 = phe2 = NULL
rho12 = se.rho12 = zeta12 = se.zeta12 = obs = NULL
rho_hec12 = se.rho_hec12 = NULL
for (i in 1:K1){
  for (j in 1:K1){
    if (i>=j) next
    phe1 = append(phe1, phes[i])
    phe2 = append(phe2, phes[j])
    rho12 = append(rho12, rho[j,i])
    se.rho12 = append(se.rho12, se.rho[j,i])
    zeta12 = append(zeta12, zeta[j,i])
    se.zeta12 = append(se.zeta12, se.zeta[j,i])
    rho_hec12 = append(rho_hec12, rho_hec[j,i])
    se.rho_hec12 = append(se.rho_hec12, se.rho_hec[j,i])
    obs = append(obs, mis[j,i])
  }
  cat(i,'')
}

rhohec = round(cbind(rho12,se.rho12,zeta12,se.zeta12,rho_hec12,se.rho_hec12),4)
coh = data.frame(cbind(phe1,phe2,rhohec))
colnames(coh) = c('Phenotype 1', 'Phenotype 2',
                  'Coheritability', 'Std Err',
                  'Environmental correlation', 'Std Err',
                  'Coheritability (HEc)', 'Std Err')
write.csv(coh, file='Table5.csv')

# heat map
rho_hec[is.na(rho_hec)] = 0
rho_hec = rho_hec + t(rho_hec)
r_hec = table3[kept,ncol(table3)-1]
names(r_hec) = phes
diag(rho_hec) = r_hec[1:K1]
colnames(rho) = rownames(rho) = phes
colnames(rho_hec) = rownames(rho_hec) = phes[1:K1]
heatmap.2(rho[1:K1,1:K1], scale = "none", col = bluered(100), breaks=seq(-1,1,0.02),
          trace = "none", density.info = "none", #Rowv=FALSE, Colv=FALSE,
          main='Coheritability by MPCH, continuous phenotypes',margins=c(20,20))
heatmap.2(rho_hec[1:K1,1:K1], scale = "none", col = bluered(100), breaks=seq(-1,1,0.02),
          trace = "none", density.info = "none", #Rowv=FALSE, Colv=FALSE,
          main='Coheritability by HEc, continuous phenotypes',margins=c(20,20))

uid = table3[,'Field.UID']
sub = c(21001,50,4079,4080,3063,30690,30760,30780,30750,30870,
        30160,30150,30020,30120,30130,30140,30080,30090,30010,30000,30250)
nms1 = phes[sub]
uid1 = as.numeric(uid[1:K1])
sub = sapply(sub, function(k) which(uid1==k))
nms1 = c('BMI','Height','DBP','SBP','FEV1','Chol','HDL','LDL','HbA1c',
         'Trigly','Baso','Eosino','Haemo','Lymph','Mono','Neutro',
         'Plate','Crit','RBC','WBC','Reticu')
rho_hec_sub = rho_hec[sub,sub]
rho_sub = rho[sub,sub]
zeta_sub = zeta[sub,sub]
colnames(rho_sub) = rownames(rho_sub) = nms1
colnames(zeta_sub) = rownames(zeta_sub) = nms1
colnames(rho_hec_sub) = rownames(rho_hec_sub) = nms1
heatmap.2(rho_sub, scale = "none", col = bluered(100), breaks=seq(-1,1,0.02),
          trace = "none", density.info = "none", #Rowv=FALSE, Colv=FALSE,
          cellnote=round(rho_sub,2), notecol='black', notecex=0.4,
          main='(A) Coheritability by MPCH')
heatmap.2(rho_hec_sub, scale = "none", col = bluered(100), breaks=seq(-1,1,0.02),
          trace = "none", density.info = "none", #Rowv=FALSE, Colv=FALSE,
          cellnote=round(rho_hec_sub,2), notecol='black', notecex=0.4,
          main='(B) Coheritability by HEc')
heatmap.2(zeta_sub, scale = "none", col = bluered(100), breaks=seq(-1,1,0.02),
          trace = "none", density.info = "none", #Rowv=FALSE, Colv=FALSE,
          cellnote=round(zeta_sub,2), notecol='black', notecex=0.4,
          main='Env correlation')

# sort
nms_sort = c('Mono','Eosino','Lymph','WBC','Neutro','HbA1c','Plate','Crit',
             'Baso','Trigly','BMI','RBC','Reticu','DBP','SBP','Haemo',
             'Height','LDL','Chol','FEV1','HDL')
rho_sub = rho_sub[nms_sort,rev(nms_sort)]
rho_hec_sub = rho_hec_sub[nms_sort,rev(nms_sort)]
heatmap.2(rho_sub, scale = "none", col = bluered(100), breaks=seq(-1,1,0.02),
          trace = "none", density.info = "none", #Rowv=FALSE, Colv=FALSE,
          cellnote=round(rho_sub,2), notecol='black', notecex=0.4,
          main='(A) Coheritability by MPCH')
heatmap.2(rho_hec_sub, scale = "none", col = bluered(100), breaks=seq(-1,1,0.02),
          trace = "none", density.info = "none", Rowv=FALSE, Colv=FALSE,
          cellnote=round(rho_hec_sub,2), notecol='black', notecex=0.4,
          main='(B) Coheritability by HEc')
