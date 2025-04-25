## Generate data

Gi_pc = matrix(c(1,.5,.5,1),2,2)
Gi_fmc = matrix(c(1,0,.5,0,1,.5,.5,.5,1),3,3)
Gi_pcs = matrix(c(1,.5,.5,.5,1,.25,.5,.25,1),3,3)
Gi_fmcs = cbind(rbind(Gi_fmc,c(.5,.5,.25)),c(.5,.5,.25,1))

N_ind = 300
N_pc = 300
N_fmc = 300
N_pcs = 300
N_fmcs = 300
Fam_Num = N_ind + N_pc + N_fmc + N_pcs + N_fmcs
type = c(rep(0,N_ind), rep(1,N_pc), rep(2,N_fmc), rep(3,N_pcs), rep(4,N_fmcs))
familyID = c(1:N_ind, rep(1:N_pc, each=2)+N_ind, rep(1:(N_fmc+N_pcs), each=3)+N_ind+N_pc,
             rep(1:N_fmcs, each=4)+N_ind+N_pc+N_fmc+N_pcs) - 1
N = length(familyID)
Fam_size = c(rep(1,N_pc), rep(2,N_pc), rep(3,N_fmc), rep(3,N_pcs), rep(4,N_fmcs))
KM = matrix(0, nrow=N, ncol=4)
KM[1:N_ind,1] = 1
for (i in 1:N_pc){
  KM[N_ind+(2*i-1):(2*i),1:2] = Gi_pc
}
for (i in 1:N_fmc){
  KM[N_ind+2*N_pc+(3*i-2):(3*i),1:3] = Gi_fmc
}
for (i in 1:N_pcs){
  KM[N_ind+2*N_pc+3*N_fmc+(3*i-2):(3*i),1:3] = Gi_pcs
}
for (i in 1:N_fmcs){
  KM[N_ind+2*N_pc+3*N_fmc+3*N_pcs+(4*i-3):(4*i),1:4] = Gi_fmcs
}

#X0 = rep(1, N)
X1 = runif(N, -1, 1)
X2 = runif(N, -1, 1)
X3 = runif(N, -1, 1)
X = cbind(X1,X2,X3)
e = c(rep(rnorm(N_ind), each=1), rep(rnorm(N_pc), each=2), rep(rnorm(N_fmc+N_pcs), each=3),
      rep(rnorm(N_fmcs), each=4))

Gamma = diag(c(0.64,0.49,0.64,0.49,0.49,0.36))
for (i in 1:6){
  for (j in 1:6){
    if (i!=j) Gamma[i,j] = sqrt(Gamma[i,i]*Gamma[j,j])*0.5
  }
}
eps = c(as.numeric(t(mvrnorm(N_ind,rep(0,6),Gamma))),
        as.numeric(t(mvrnorm(N_pc,rep(0,12),kronecker(Gi_pc,Gamma)))),
        as.numeric(t(mvrnorm(N_fmc,rep(0,18),kronecker(Gi_fmc,Gamma)))),
        as.numeric(t(mvrnorm(N_pcs,rep(0,18),kronecker(Gi_pcs,Gamma)))),
        as.numeric(t(mvrnorm(N_fmcs,rep(0,24),kronecker(Gi_fmcs,Gamma)))))
p = ncol(X)

# Two Normal
alpha1 = c(1,1,1)/2
alpha2 = c(1,2,2)/3
theta1 = 0.8
theta2 = 0.6
gamma1 = Gamma[1,1]
gamma2 = Gamma[2,2]
sig1 = 1
sig2 = 0.8
true1_con = c(theta1,gamma1,alpha1,sig1^2)
true2_con = c(theta2,gamma2,alpha2,sig2^2)
eps1 = eps[6*(1:N)-5]
eps2 = eps[6*(1:N)-4]
Y1_con = as.numeric(X%*%alpha1 + theta1*e + eps1 + rnorm(N,0,sig1))
Y2_con  = as.numeric(X%*%alpha2 + theta2*e + eps2 + rnorm(N,0,sig2))

# Two Ordinal (including one binary)
alpha1 = c(1,1,1)/2
alpha2 = c(1,2,2)/3
theta1 = 0.6
theta2 = 0.5
gamma1 = Gamma[3,3]
gamma2 = Gamma[4,4]
delta1 = 1
delta2 = c(-1,0,1)
true1_ord = c(theta1,gamma1,alpha1,delta1)
true2_ord = c(theta2,gamma2,alpha2,delta2)
eps1 = eps[6*(1:N)-3]
eps2 = eps[6*(1:N)-2]
Z1 = as.numeric(X%*%alpha1) + theta1*e + eps1 + rnorm(N)
Z2 = as.numeric(X%*%alpha2) + theta2*e + eps2 + rnorm(N)
Y1_ord = rep(2, N)
Y1_ord[Z1<=delta1] = 1
Y2_ord = rep(4, N)
Y2_ord[Z2<=delta2[3]] = 3
Y2_ord[Z2<=delta2[2]] = 2
Y2_ord[Z2<=delta2[1]] = 1

# Two Survival
alpha1 = c(1,1,1)/2
alpha2 = c(1,2,2)/3
theta1 = 0.4
theta2 = 0.5
gamma1 = Gamma[5,5]
gamma2 = Gamma[6,6]
true1_sur = c(theta1,gamma1,alpha1)
true2_sur = c(theta2,gamma2,alpha2)
eps1 = eps[6*(1:N)-1]
eps2 = eps[6*(1:N)]
Y1 = exp(-10*log(runif(N))/exp(as.numeric(X%*%alpha1)+theta1*e+eps1))-1
Y1[is.infinite(Y1)] = 99
Y1[is.nan(Y1)] = 99
C1 = runif(N,5,10)
D1_sur = as.numeric(Y1 <= C1)
T1_sur = D1_sur*Y1 + (1-D1_sur)*C1
Y2 = sqrt(-log(runif(N))/exp(as.numeric(X%*%alpha2)+theta2*e+eps2))/0.06
Y2[is.infinite(Y2)] = 99
Y2[is.nan(Y2)] = 99
C2 = runif(N,5,10)
D2_sur = as.numeric(Y2 <= C2)
T2_sur = D2_sur*Y2 + (1-D2_sur)*C2
