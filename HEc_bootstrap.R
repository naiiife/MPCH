## bootstrap HEc

B = 200
rho.hec.b = matrix(NA,B,K1)
con.dat = con.dat[1:N,]
for (b in 1:B){
  ss = sample(unique(familyID), replace=TRUE)
  Y_b = X0_b = KM_b = familyID_b = NULL
  im = 0
  for (i in ss){
    Y_b = rbind(Y_b, con.dat[familyID==i,])
    X0_b = rbind(X0_b, X0[familyID==i,])
    KM_b = rbind(KM_b, KM[familyID==i,])
    familyID_b = append(familyID_b, rep(im,sum(familyID==i)))
    im = im + 1
  }
  for (k in 1:K1){
    Y1.p = predict(lm(Y_b[,k]~0+X0_b), data.frame(X0_b))
    f1 = Y_b[,k] - Y1.p
    var1 = var(f1, na.rm=TRUE)
    v1 = tr = 0
    for (m in unique(familyID_b)){
      ind = (familyID_b==m)
      if (sum(ind)==0) next
      Ki = KM_b[ind,]
      Ki = Ki[,1:nrow(Ki)]
      diag(Ki) = 0
      f1i = f1[ind]
      ind = !is.na(f1i)
      if (sum(ind)==0) next
      v1 = v1 + t(f1i[ind])%*% Ki[ind,ind] %*%f1i[ind]
      tr = tr + sum(diag(Ki[ind,ind]%*%Ki[ind,ind]))
    }
    rho.hec.b[b,k] = v1/tr/var1
    cat(paste0(b, '-', k), '')
  }
}
se.rho.hec = apply(rho.hec.b,2,sd)
ci.rho.hec.b = sapply(1:K1, function(k) quantile(rho.hec.b[,k],c(.025,.975)))
rho.hec.low = ci.rho.hec.b[1,]
rho.hec.up = ci.rho.hec.b[2,]