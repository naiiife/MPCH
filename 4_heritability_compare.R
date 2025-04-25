setwd('/home/ukb')
load('condat.Rdata')
n = n2+n3+n4
N = n2*2+n3*3+n4*4
con.dat = con.dat[1:N,]
X0 = as.matrix(cbind(1,X[1:N,]))
familyID = familyID[1:N]
KM = KM[1:N,]
K1 = ncol(con.dat)

uid = table3[,'Field.UID']
phes = table3[,'Phenotype']
types = table3[,'Data.type']
sub = c(21001,50,4079,4080,3063,30690,30760,30780,30750,30870,
        30160,30150,30020,30120,30130,30140,30080,30090,30010,30000,30250)
nms1 = phes[1:K1]
uid1 = as.numeric(uid[1:K1])
sub = sapply(sub, function(k) which(uid1==k))

## Proposed
rho = rho.low = rho.up = rep(NA,K)
names(rho) = names(rho.low) = names(rho.up) = phes
for (k in 1:K){
  if (k<=K1){
    fit = fit_con[[k]]
    if (is.null(fit)) next
  }
  if (k>K1&k<=K1+K2){
    fit = fit_bin[[k-K1]]
    if (is.null(fit)) next
    fit$sigmau21 = 1
  }
  if (k>K1+K2&k<=K1+K2+K3){
    fit = fit_ord[[k-K1-K2]]
    if (is.null(fit)) next
    fit$sigmau21 = 1
  }
  if (k>K1+K2+K3&k<=K1+K2+K3+K4){
    fit = fit_sur[[k-K1-K2-K3]]
    if (is.null(fit)) next
    fit$sigmau21 = pi^2/6
  }
  D1 = rep(0,ncol(fit$phi))
  theta_1 = fit$theta
  gamma_1 = fit$gamma
  sigmau21 = fit$sigmau2
  sigma_1 = sqrt(theta_1^2+gamma_1+sigmau21)
  rho_1 = gamma_1/sigma_1^2
  D1[2] = -2*theta_1/(theta_1^2+sigmau21)
  D1[3] = 1/gamma_1
  if (k<=K1) D1[length(D1)] = -1/(theta_1^2+sigmau21)
  IF_rho_1 = as.numeric(fit$phi%*%D1)
  rho[k] = rho_1
  rho.low[k] = expit(logit(rho_1)-1.96*sd(IF_rho_1)/sqrt(nrow(fit$phi)))
  rho.up[k] = expit(logit(rho_1)+1.96*sd(IF_rho_1)/sqrt(nrow(fit$phi)))
}


## HEc
setwd('/home/ukb/family_data_all/')
rho.hec = rho.hec.low = rho.hec.up = rep(0,K1)
for (i in 1:K1){
  Y_i = con.dat[,i]
  Y1.p = predict(lm(Y_i~0+X0), data.frame(Y_i,X0))
  f1 = Y_i - Y1.p
  var1 = var(f1, na.rm=TRUE)
  v1 = tr = 0
  for (m in unique(familyID)){
    ind = (familyID==m)
    if (sum(ind)==0) next
    Ki = KM[ind,]
    Ki = Ki[,1:nrow(Ki)]
    diag(Ki) = 0
    f1i = f1[ind]
    ind = !is.na(f1i)
    if (sum(ind)==0) next
    v1 = v1 + t(f1i[ind])%*% Ki[ind,ind] %*%f1i[ind]
    tr = tr + sum(diag(Ki[ind,ind]%*%Ki[ind,ind]))
  }
  rho.hec[i] = v1/tr/var1
  cat(i, '')
}
set.seed(2024)
source('HEc_bootstrap.R')
names(rho.hec) = nms1
save(rho.hec,rho.hec.low,rho.hec.up,se.rho.hec, file='HEc_heritability.Rdata')
#load('HEc_heritability.Rdata')

setwd('/home/ukb/')
her = cbind(round(rho.hec,4), round(se.rho.hec,4))
her = rbind(her, matrix(NA,nrow(table3)-K1,2))
colnames(her) = c('Heritability (HEc)', 'Std Err')
table4 = cbind(table3, her)
rownames(table4) = 1:nrow(table4)
write.csv(table4, file='Table3.csv')

df = data.frame(value=c(rho[1:K1],rho.hec),
                Method=c(rep('MPCH',K1),rep('HEc',K1)))
df$method = as.factor(df$Method)
ggplot(df, aes(x=value,color=Method,fill=Method)) + xlim(0,1) +
  geom_histogram(aes(y=..density..), position="identity", alpha=0.4)+
  geom_density(alpha=0, lwd=1)+
  theme(legend.position="bottom") +
  labs(title='Genetic heritability + environmental effect\n(K1=152 phenotypes)', 
       x='Proportion')

# phenotypes with high heritability
dt = data.frame(matrix(0,20,6))
colnames(dt) = c('Phenotype','Type','rho','low1','up1','ci')
highest = c('50','30790','5096','5099','20261',
            '20160','PHE401','4631','PHE240','PHE490',
            '1747','1717','20487','2060','5164',
            '3786','42016','2976','42000','42002')
highest = sapply(highest, function(k) which(uid==k))
dt$Phenotype = phes[highest]
dt$Type = c(rep('(Continuous)',5),rep('(Binary)',5),
            rep('(Ordinal)',5),rep('(Time to event)',5))
dt$rho = rho[highest]
dt$low1 = rho.low[highest]
dt$up1 = rho.up[highest]
for (i in 1:20){
  dt$ci[i] = paste0(round(dt$rho[i],3),' (',
                    round(dt$low1[i],3),', ',round(dt$up1[i],3),')')
}
colnames(dt)[6] = 'Estimate (95% CI)'
dt$' ' = paste(rep(" ",30), collapse=" ")
herithigh_p <- forest(dt[,c(1,2,7,6)],
                      est = list(dt$rho),
                      lower = list(dt$low1),
                      upper = list(dt$up1),
                      ci_column = 3,
                      vert_line = 0,
                      nudge_y = 0.2,
                      ticks_at = c(0,.2,.4,.6,.8,1))
herithigh_p

# compare with HEc
dt = data.frame(matrix(0,K1,9))
colnames(dt) = c('Phenotype','rho','rho.hec','low1','up1','low2','up2','ci','ci.hec')
dt$Phenotype = paste0(nms1,'  ')
dt$rho = rho[1:K1]
dt$rho.hec = rho.hec
dt$low1 = rho.low[1:K1]
dt$up1 = rho.up[1:K1]
dt$low2 = rho.hec.low
dt$up2 = rho.hec.up
for (i in 1:K1){
  dt$ci[i] = paste0(round(dt$rho[i],3),' (',
                    round(dt$low1[i],3),', ',round(dt$up1[i],3),')')
  dt$ci.hec[i] = paste0(round(dt$rho.hec[i],3),' (',
                        round(dt$low2[i],3),', ',round(dt$up2[i],3),')')
}
colnames(dt)[8] = 'Est (95% CI) MPCH'
colnames(dt)[9] = 'Est (95% CI) HEc'
dt$' ' = paste(rep(" ",30), collapse=" ")
dt = dt[sub,]
tm <- forest_theme(
  legend_name = 'Method',
  legend_value = c('MPCH','HEc'))
herit_p <- forest(dt[,c(1,10,8,9)],
                  est = list(dt$rho, dt$rho.hec),
                  lower = list(dt$low1, dt$low2),
                  upper = list(dt$up1, dt$up2),
                  ci_column = 2,
                  vert_line = 0,
                  nudge_y = 0.2,
                  ticks_at = c(0,.2,.4,.6,.8,1),
                  theme = tm)
herit_p
