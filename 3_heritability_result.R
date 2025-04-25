## Single-trait heritability
library(grid)
library(gplots)
library(ggplot2)
library(forestploter)
expit <- function(x) {1/(1+exp(-x))}
logit <- function(x) {log(x/(1-x))}

setwd('/home/ukb')
table2 = read.csv('Table2.csv')
table2 = table2[,-1]
uid = table2[,'Field.UID']
phes = table2[,'Phenotype']
types = table2[,'Data.type']
K1 = 152
K2 = 97
K3 = 27
K4 = 14
K = K1+K2+K3+K4

setwd('/home/ukb/results_alldat/')

fracgen = rep(NA, K)
heritability = rep(NA, K)
fracenv = rep(NA, K)
environmental = rep(NA, K)
se.f = se.rho = se.zeta = rep(NA, K)
alpha = se.alpha = matrix(NA, K, 15)
iterrun = timerun = rep(NA, K)
fit_con = rep(list(NULL),K1)
for (k in 1:K1){
  tr = try(load(paste0('fit_con',k,'.Rdata')))
  if (!'try-error' %in% class(tr)){
    gammak = fit$gamma
    thetak = fit$theta
    sigma2 = fit$sigmau2
    fracgen[k] = gammak/(gammak+thetak^2)
    heritability[k] = gammak/(gammak+thetak^2+sigma2)
    fracenv[k] = thetak^2/(gammak+thetak^2)
    environmental[k] = thetak^2/(gammak+thetak^2+sigma2)
    est = fit$est
    se = fit$se
    IF = fit$phi
    iterrun[k] = fit$iterrun
    timerun[k] = fit$tm
    se.f[k] = fit$se.fracgen
    se.rho[k] = fit$se.rho
    se.zeta[k] = fit$se.zeta
    alpha[k,] = fit$alpha
    se.alpha[k,] = se[4:length(est)-1]
    fit_con[[k]] = fit
  } 
}

fit_bin = rep(list(NULL),K2)
for (k in 1:K2){
  tr = try(load(paste0('fit_bin',k,'.Rdata')))
  if (!'try-error' %in% class(tr)){
    gammak = fit$gamma
    thetak = fit$theta
    fracgen[k+K1] = gammak/(gammak+thetak^2)
    heritability[k+K1] = gammak/(gammak+thetak^2+1)
    fracenv[k+K1] = thetak^2/(gammak+thetak^2)
    environmental[k+K1] = thetak^2/(gammak+thetak^2+1)
    iterrun[k+K1] = fit$iterrun
    timerun[k+K1] = fit$tm
    est = fit$est
    se = fit$se
    IF = fit$phi
    alpha[k+K1,] = c(fit$delta,fit$alpha)
    se.alpha[k+K1,] = se[3:length(est)]
    se.f[k+K1] = fit$se.fracgen
    se.rho[k+K1] = fit$se.rho
    se.zeta[k+K1] = fit$se.zeta
    fit$sigmau2 = 1
    fit_bin[[k]] = fit
  } 
}

fit_ord = rep(list(NULL),K3)
for (k in 1:K3){
  tr = try(load(paste0('fit_ord',k,'.Rdata')))
  p = 14
  if (!'try-error' %in% class(tr)){
    gammak = fit$gamma
    thetak = fit$theta
    fracgen[k+K1+K2] = gammak/(gammak+thetak^2)
    heritability[k+K1+K2] = gammak/(gammak+thetak^2+1)
    fracenv[k+K1+K2] = thetak^2/(gammak+thetak^2)
    environmental[k+K1+K2] = thetak^2/(gammak+thetak^2+1)
    est = fit$est
    se = fit$se
    IF = fit$phi
    iterrun[k+K1+K2] = fit$iterrun
    timerun[k+K1+K2] = fit$tm
    alpha[k+K1+K2,2:15] = fit$alpha
    se.alpha[k+K1+K2,2:15] = se[2+1:p]
    se.f[k+K1+K2] = fit$se.fracgen
    se.rho[k+K1+K2] = fit$se.rho
    se.zeta[k+K1+K2] = fit$se.zeta
    fit$sigmau2 = 1
    fit_ord[[k]] = fit
  } 
}

fit_sur = rep(list(NULL),K4)
for (k in 1:K4){
  tr = try(load(paste0('fit_sur',k,'.Rdata')))
  if (!'try-error' %in% class(tr)){
    gammak = fit$gamma
    thetak = fit$theta
    fracgen[k+K1+K2+K3] = gammak/(gammak+thetak^2)
    heritability[k+K1+K2+K3] = gammak/(gammak+thetak^2+1.644934)
    fracenv[k+K1+K2+K3] = thetak^2/(gammak+thetak^2)
    environmental[k+K1+K2+K3] = thetak^2/(gammak+thetak^2+1.644934)
    iterrun[k+K1+K2+K3] = fit$iterrun
    timerun[k+K1+K2+K3] = fit$tm
    est = fit$est
    se = fit$se
    IF = fit$phi
    alpha[k+K1+K2+K3,2:15] = c(fit$alpha[1],NA,fit$alpha[2:length(fit$alpha)])
    se.alpha[k+K1+K2+K3,2:15] = c(fit$se[3],NA,se[4:length(fit$se)])
    se.f[k+K1+K2+K3] = fit$se.fracgen
    se.rho[k+K1+K2+K3] = fit$se.rho
    se.zeta[k+K1+K2+K3] = fit$se.zeta
    fit$sigmau2 = pi^2/6
    fit_sur[[k]] = fit
  } 
}

names(heritability) = names(se.rho) = table2[,'Phenotype']
setwd('/home/ukb/')
alphase = round(cbind(alpha,se.alpha),4)
table2 = table2[,c('Field.UID','Phenotype','Data.type')]
table3 = cbind(table2, round(heritability,4), round(se.rho,4),
               round(environmental,4), round(se.zeta,4),
               round(fracgen,4), round(fracenv,4), round(se.f,4),
               iterrun, round(timerun/60,3))
colnames(table3) = c('Field UID','Phenotype','Data type',
                     'Heritability', 'Std Err',
                     'Environmental effect', 'Std Err',
                     'Frac gen eff', 'Frac env eff', 'Std Err', 
                     'Iteration','CPU minutes')
rownames(table3) = 1:nrow(table3)
write.csv(table3, file='Table3.csv')

kept = (1:K)[!(1:K)%in%c(193,195)]
heritability = heritability[kept]
environmental = environmental[kept]
fracgen = fracgen[kept]
se.rho = se.rho[kept]

#Histogram
df = data.frame(value=heritability,
                variable=rep('Heritability',length(heritability)))
df$variable = as.factor(df$variable)
ggplot(df, aes(x=value)) + xlim(0,1) +
  #geom_histogram(position="identity") +
  geom_histogram(colour="black", fill="lightblue")+
  #geom_density(alpha=.2, color='darkblue', lwd=1) +
  labs(title='(A) Heritability (K=288 phenotypes)', 
       x='value')
df = data.frame(value=environmental,
                variable=rep('Environmental',length(environmental)))
df$variable = as.factor(df$variable)
ggplot(df, aes(x=value)) + xlim(0,1) +
  #geom_histogram(position="identity") +
  geom_histogram(colour="black", fill="lightgreen")+
  #geom_density(alpha=.2, color='darkgreen', lwd=1) +
  labs(title='(B) Environmental effect (K=288 phenotypes)', 
       x='value')
df = data.frame(value=fracgen,
                variable=rep('Genetic',length(fracgen)))
df$variable = as.factor(df$variable)
ggplot(df, aes(x=value)) + xlim(0,1) +
  #geom_histogram(position="identity") +
  geom_histogram(colour="black", fill="pink")+
  #geom_density(alpha=.2, color='darkred', lwd=1) +
  labs(title='(A) Fraction of genetic effect (K=288 phenotypes)', 
       x='value')
pval = pnorm(-abs(heritability/se.rho))
Op = -log(sort(na.omit(as.numeric(pval))),10)
mean(pval<0.05/288)
Ep = -log(seq(0,1,length=length(Op)),10)
df = data.frame(Op=Op[Op<Inf],Ep=Ep[Op<Inf])
ggplot(df, aes(x=Ep,y=Op)) + xlim(0,0.5) + ylim(0,60) +
  geom_line(aes(x=Ep,y=Ep), color='lightblue') +
  geom_point(col='darkblue') +
  labs(title='(C) Q-Q plot of heritability', 
       x='Expected -lg(P value)', y='Observed -lg(P value)')