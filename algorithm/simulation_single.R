#setwd("/home/yuhaoden/jobs/")


lid = Sys.getenv("SLURM_ARRAY_TASK_ID")
if (lid=="") lid=0

## Load data

set.seed(lid)
source('head.R')
source('generate.R')
source('continuous.R')
source('probit.R')
source('ordinal.R')
source('survival.R')
source('con_con.R')
source('con_ord.R')
source('con_sur.R')
source('ord_ord.R')
source('ord_sur.R')
source('sur_sur.R')

## Estimate diagonal elements, save parameters and influence functions

t0 = Sys.time()
res1_con = cal_con(Y1_con,X,familyID,KM,true1_con)
t1 = Sys.time()
res2_con = cal_con(Y2_con,X,familyID,KM,true2_con)
t2 = Sys.time()
res1_ord = cal_bin(Y1_ord-1,X,familyID,KM,true1_ord)
t3 = Sys.time()
res2_ord = cal_ord(Y2_ord,X,familyID,KM,true2_ord)
t4 = Sys.time()
res1_sur = cal_sur(T1_sur,D1_sur,X,familyID,KM,true1_sur)
t5 = Sys.time()
res2_sur = cal_sur(T2_sur,D2_sur,X,familyID,KM,true2_sur)
t6 = Sys.time()

save(res1_con, res2_con, res1_ord, res2_ord, 
    res1_sur, res2_sur, 
    Y1_con, Y2_con, Y1_ord, Y2_ord, T1_sur, D1_sur, T2_sur, D2_sur,
    X, familyID,
    file=paste0('res',lid,'.Rdata'))


kinship = rep(1,length(familyID))
t7 = Sys.time()
res12 = con_con(Y1_con,Y2_con,X,X,familyID,res1_con,res2_con,kinship)
res13 = con_ord(Y1_con,Y1_ord,X,X,familyID,res1_con,res1_ord,kinship)
res14 = con_ord(Y1_con,Y2_ord,X,X,familyID,res1_con,res2_ord,kinship)
res15 = con_sur(Y1_con,T1_sur,D1_sur,X,X,familyID,res1_con,res1_sur,kinship)
res16 = con_sur(Y1_con,T2_sur,D2_sur,X,X,familyID,res1_con,res2_sur,kinship)
res23 = con_ord(Y2_con,Y1_ord,X,X,familyID,res2_con,res1_ord,kinship)
res24 = con_ord(Y2_con,Y2_ord,X,X,familyID,res2_con,res2_ord,kinship)
res25 = con_sur(Y2_con,T1_sur,D1_sur,X,X,familyID,res2_con,res1_sur,kinship)
res26 = con_sur(Y2_con,T2_sur,D2_sur,X,X,familyID,res2_con,res2_sur,kinship)
res34 = ord_ord(Y1_ord,Y2_ord,X,X,familyID,res1_ord,res2_ord,kinship)
res35 = ord_sur(Y1_ord,T1_sur,D1_sur,X,X,familyID,res1_ord,res1_sur,kinship)
res36 = ord_sur(Y1_ord,T2_sur,D2_sur,X,X,familyID,res1_ord,res2_sur,kinship)
res45 = ord_sur(Y2_ord,T1_sur,D1_sur,X,X,familyID,res2_ord,res1_sur,kinship)
res46 = ord_sur(Y2_ord,T2_sur,D2_sur,X,X,familyID,res2_ord,res2_sur,kinship)
res56 = sur_sur(T1_sur,D1_sur,T2_sur,D2_sur,X,X,familyID,res1_sur,res2_sur,kinship)
t8 = Sys.time()

t.1 = as.numeric(difftime(t1, t0, units = "secs"))
t.2 = as.numeric(difftime(t2, t1, units = "secs"))
t.3 = as.numeric(difftime(t3, t2, units = "secs"))
t.4 = as.numeric(difftime(t4, t3, units = "secs"))
t.5 = as.numeric(difftime(t5, t4, units = "secs"))
t.6 = as.numeric(difftime(t6, t5, units = "secs"))
t.7 = as.numeric(difftime(t8, t7, units = "secs"))
times = c(t.1,t.2,t.3,t.4,t.5,t.6,t.7)

est_rho0 = c(res1_con$rho,res2_con$rho,res1_ord$rho,res2_ord$rho,res1_sur$rho,res2_sur$rho)
se_rho0 = c(res1_con$se.rho,res2_con$se.rho,res1_ord$se.rho,res2_ord$se.rho,
            res1_sur$se.rho,res2_sur$se.rho)
est_zeta0 = c(res1_con$zeta,res2_con$zeta,res1_ord$zeta,res2_ord$zeta,res1_sur$zeta,res2_sur$zeta)
se_zeta0 = c(res1_con$se.zeta,res2_con$se.zeta,res1_ord$se.zeta,res2_ord$se.zeta,
             res1_sur$se.zeta,res2_sur$se.zeta)
est_rho = c(res12$rho12,res13$rho12,res14$rho12,res15$rho12,res16$rho12,
            res23$rho12,res24$rho12,res25$rho12,res26$rho12,
            res34$rho12,res35$rho12,res36$rho12,
            res45$rho12,res46$rho12,res56$rho12)
se_rho = c(res12$se.rho12,res13$se.rho12,res14$se.rho12,res15$se.rho12,res16$se.rho12,
           res23$se.rho12,res24$se.rho12,res25$se.rho12,res26$se.rho12,
           res34$se.rho12,res35$se.rho12,res36$se.rho12,
           res45$se.rho12,res46$se.rho12,res56$se.rho12)
est_zeta = c(res12$zeta12,res13$zeta12,res14$zeta12,res15$zeta12,res16$zeta12,
            res23$zeta12,res24$zeta12,res25$zeta12,res26$zeta12,
            res34$zeta12,res35$zeta12,res36$zeta12,
            res45$zeta12,res46$zeta12,res56$zeta12)
se_zeta = c(res12$se.zeta12,res13$se.zeta12,res14$se.zeta12,res15$se.zeta12,res16$se.zeta12,
            res23$se.zeta12,res24$se.zeta12,res25$se.zeta12,res26$se.zeta12,
            res34$se.zeta12,res35$se.zeta12,res36$se.zeta12,
            res45$se.zeta12,res46$se.zeta12,res56$se.zeta12)
            
est_rho0_p = c(res12$rho1_p,res12$rho2_p,res34$rho1_p,res34$rho2_p,res56$rho1_p,res56$rho2_p)
se_rho0_p = c(res12$se.rho1_p,res12$se.rho2_p,res34$se.rho1_p,res34$se.rho2_p,
            res56$se.rho1_p,res56$se.rho2_p)
est_rho_p = c(res12$rho12_p,res13$rho12_p,res14$rho12_p,res15$rho12_p,res16$rho12_p,
            res23$rho12_p,res24$rho12_p,res25$rho12_p,res26$rho12_p,
            res34$rho12_p,res35$rho12_p,res36$rho12_p,
            res45$rho12_p,res46$rho12_p,res56$rho12_p)
se_rho_p = c(res12$se.rho12_p,res13$se.rho12_p,res14$se.rho12_p,res15$se.rho12_p,res16$se.rho12_p,
           res23$se.rho12_p,res24$se.rho12_p,res25$se.rho12_p,res26$se.rho12_p,
           res34$se.rho12_p,res35$se.rho12_p,res36$se.rho12_p,
           res45$se.rho12_p,res46$se.rho12_p,res56$se.rho12_p)
est_zeta_p = c(res12$zeta12_p,res13$zeta12_p,res14$zeta12_p,res15$zeta12_p,res16$zeta12_p,
            res23$zeta12_p,res24$zeta12_p,res25$zeta12_p,res26$zeta12_p,
            res34$zeta12_p,res35$zeta12_p,res36$zeta12_p,
            res45$zeta12_p,res46$zeta12_p,res56$zeta12_p)
se_zeta_p = c(res12$se.zeta12_p,res13$se.zeta12_p,res14$se.zeta12_p,res15$se.zeta12_p,res16$se.zeta12_p,
            res23$se.zeta12_p,res24$se.zeta12_p,res25$se.zeta12_p,res26$se.zeta12_p,
            res34$se.zeta12_p,res35$se.zeta12_p,res36$se.zeta12_p,
            res45$se.zeta12_p,res46$se.zeta12_p,res56$se.zeta12_p)

save(res1_con, res2_con, res1_ord, res2_ord, res1_sur, res2_sur, times,
    est_rho0, est_zeta0, est_rho, est_zeta, se_rho0, se_zeta0, se_rho, se_zeta,
    est_rho0_p, est_rho_p, est_zeta_p, se_rho0_p, se_rho_p, se_zeta_p,
    file=paste0('res',lid,'.Rdata'))

