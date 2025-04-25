setwd('/home/ukb/results_alldat_co/')
K = 290
rho = se.rho = zeta = se.zeta = gamma = mis = tm = matrix(NA,K,K)
fracgen = se.fracgen = fracenv = se.fracenv = segamma = matrix(NA,K,K)
for (i in 1:K){
  cat(i,'')
  for (b in 1:29){
    res = matrix(NA,10,10)
    tr = try(load(paste0('res',i,'_',b,'.Rdata')))
    if ('try-error' %in% tr) next
    rg = (b-1)*10+1:10
    fracgen[rg,i] = res[,1]
    se.fracgen[rg,i] = res[,2]
    fracenv[rg,i] = res[,3]
    se.fracenv[rg,i] = res[,4]
    rho[rg,i] = res[,5]
    se.rho[rg,i] = res[,6]
    zeta[rg,i] = res[,7]
    se.zeta[rg,i] = res[,8]
    mis[rg,i] = res[,9]
    tm[rg,i] = res[,10]
  }
}

setwd('/home/ukb/')

kept = (1:K)[!(1:K)%in%c(193,195)]
fracgen = fracgen[kept,kept]
fracenv = fracenv[kept,kept]
rho = rho[kept,kept]
zeta = zeta[kept,kept]
se.fracgen = se.fracgen[kept,kept]
se.fracenv = se.fracenv[kept,kept]
se.rho = se.rho[kept,kept]
se.zeta = se.zeta[kept,kept]
mis = mis[kept,kept]
tm = tm[kept,kept]

table3 = read.csv('Table3.csv')
phes = table3[kept,'Phenotype']
types = table3[kept,'Data.type']
phe1 = phe2 = NULL
type1 = type2 = NULL
K = length(kept)
rho12 = se.rho12 = zeta12 = se.zeta12 = obs = timerun = NULL
fracgen12 = se.fracgen12 = fracenv12 = se.fracenv12 = NULL
for (i in 1:K){
  for (j in 1:K){
    if (i>=j) next
    phe1 = append(phe1, phes[i])
    phe2 = append(phe2, phes[j])
    type1 = append(type1, types[i])
    type2 = append(type2, types[j])
    rho12 = append(rho12, rho[j,i])
    se.rho12 = append(se.rho12, se.rho[j,i])
    zeta12 = append(zeta12, zeta[j,i])
    se.zeta12 = append(se.zeta12, se.zeta[j,i])
    fracgen12 = append(fracgen12, fracgen[j,i])
    se.fracgen12 = append(se.fracgen12, se.fracgen[j,i])
    fracenv12 = append(fracenv12, fracenv[j,i])
    se.fracenv12 = append(se.fracenv12, se.fracenv[j,i])
    obs = append(obs, mis[j,i])
    timerun = append(timerun, tm[j,i])
  }
  cat(i,'')
}
rhozeta = round(cbind(rho12,se.rho12,zeta12,se.zeta12),3)
frac = round(cbind(fracgen12,se.fracgen12,fracenv12,se.fracenv12),3)
co = data.frame(cbind(phe1,type1,phe2,type2,obs,
                      rhozeta,frac,round(timerun/60,2)))
colnames(co) = c('Phenotype 1','Data type 1',
                 'Phenotype 2','Data type 2', 'Observations',
                 'Coheritability', 'Std Err',
                 'Environmental correlation', 'Std Err',
                 'Frac gen eff', 'Std Err',
                 'Frac env eff', 'Std Err', 'CPU minutes')
write.csv(co, file='Table4.csv')

df = data.frame(value=na.omit(rho12),
                variable=rep('Coheritability',length(na.omit(rho12))))
df$variable = as.factor(df$variable)
ggplot(df, aes(x=value)) + xlim(-1,1) +
  #geom_histogram(position="identity") +
  geom_histogram(colour="black", fill="lightblue")+
  #geom_density(alpha=.2, color='darkblue', lwd=1) +
  labs(title='(A) Coheritability (288x287/2 pairs)', 
       x='value')
df = data.frame(value=na.omit(zeta12),
                variable=rep('Environmental',length(na.omit(zeta12))))
df$variable = as.factor(df$variable)
ggplot(df, aes(x=value)) + xlim(-1,1) +
  #geom_histogram(position="identity") +
  geom_histogram(binwidth=0.1, colour="black", fill="lightgreen")+
  #geom_density(alpha=.2, color='darkgreen', lwd=1) +
  labs(title='(B) Environmental correlation (288x287/2 pairs)', 
       x='value')
df = data.frame(value=na.omit(fracgen12),
                variable=rep('Fracgen',length(na.omit(fracgen12))))
df$variable = as.factor(df$variable)
ggplot(df, aes(x=value)) + xlim(-1,1) +
  #geom_histogram(position="identity") +
  geom_histogram(colour="black", fill="pink")+
  #geom_density(alpha=.2, color='darkred', lwd=1) +
  labs(title='(B) Fraction of genetic effect (288x287/2 pairs)', 
       x='value')

#QQ plot
pval = pnorm(-abs(rho/se.rho))*2
pvals = -log(sort(na.omit(as.numeric(pval))),10)
Ep = -log(seq(0,1,length=length(pvals)),10)
df = data.frame(Op=pvals[pvals<Inf],Ep=Ep[pvals<Inf])
ggplot(df, aes(x=Ep,y=Op)) + xlim(0,2.4) + ylim(0,180) +
  geom_line(aes(x=Ep,y=Ep), color='lightblue') +
  geom_point(col='darkblue') +
  labs(title='(C) Q-Q plot of coheritability', 
       x='Expected -lg(P value)', y='Observed -lg(P value)')
mean(pval[!is.na(pval)]<0.05,na.rm=TRUE)

h2 = as.numeric(table3[kept,'Heritability.1'])
seh = as.numeric(table3[kept,'Std.Err.1'])
x2 = as.numeric(table3[kept,'Environmental.effect.1'])
flag = !is.na(h2)
K = length(h2)
colnames(rho) = rownames(rho) = rep('',K)#phes
colnames(zeta) = rownames(zeta) = rep('',K)#phes
rho[is.na(rho)] = 0
rho = rho + t(rho)
diag(rho) = h2
se.rho[is.na(se.rho)] = 0
se.rho = se.rho + t(se.rho)
diag(se.rho) = seh
zeta[is.na(zeta)] = 0
zeta = zeta + t(zeta)
diag(zeta) = x2
heatmap.2(rho, scale = "none", col = bluered(100), breaks=seq(-1,1,0.02),
          trace = "none", density.info = "none", #Rowv=FALSE, Colv=FALSE,
          main='Coheritability by MPCH')
heatmap.2(zeta, scale = "none", col = bluered(100), breaks=seq(-1,1,0.02),
          trace = "none", density.info = "none", #Rowv=FALSE, Colv=FALSE,
          main='Environmental correlation')
library(RColorBrewer)
rho2 = rho[flag,flag]^2
heatmap.2(rho2, scale = "none", 
          col=colorRampPalette(brewer.pal(9, "Reds"))(100), 
          breaks=seq(0,0.1,0.1*0.01),
          trace = "none", density.info = "none", #Rowv=FALSE, Colv=FALSE,
          main='(C) Squared coheritability by MPCH')
colnames(rho2) = rownames(rho2) = phes
heatmap.2(rho2, scale = "none", 
          col=colorRampPalette(brewer.pal(9, "Reds"))(100), 
          breaks=seq(0,0.6,0.6*0.01),
          trace = "none", density.info = "none", #Rowv=FALSE, Colv=FALSE,
          main='Squared coheritability by MPCH',margins=c(20,20))
