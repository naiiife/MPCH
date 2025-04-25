## Extract family
library(data.table)
library(dplyr)
library(ukbtools)
library(gplots)
library(survival)
setwd('/home/ukb')
df = read.table("ukb_rel.dat",head=TRUE)
ukb_gen_rel_count(df, plot=TRUE)
matchid = ukb_gen_related_with_data(df,unique(df[,1]),cutoff=0.0884)
matchid[,5] = matchid[,5]*2
kinpairs = matchid[,c(1,2,5)]
n = nrow(matchid)
family_2 = NULL
allmatch = c(matchid[,1],matchid[,2])
length(unique(allmatch))
for (i in 1:n){
  case = matchid[i,]
  case1 = as.numeric(case[1])
  case2 = as.numeric(case[2])
  ks = as.numeric(case[5])
  if (sum(allmatch==case1)+sum(allmatch==case2)==2){
    family_2 = rbind(family_2, case)
    matchid[i,] = NA
  }
}
matchid = na.omit(matchid)
allmatch1 = matchid[,1]
allmatch2 = matchid[,2]
n = nrow(matchid)
family2 = family3 = family4 = family5 = family6 = family7 = NULL
for (i in 1:n){
  case = matchid[i,]
  case1 = as.numeric(case[1])
  case2 = as.numeric(case[2])
  if (case1==0 | case2==0) next
  rows = i
  rows = append(rows, which(allmatch1==case1))
  rows = append(rows, which(allmatch2==case1))
  rows = append(rows, which(allmatch1==case2))
  rows = append(rows, which(allmatch2==case2))
  rows = unique(rows)
  matched = unique(as.numeric(c(matchid[rows,1],matchid[rows,2])))
  ni = length(matched)
  for (j in 1:ni){
    case3 = matched[j]
    if (case3==case1 | case3==case2) next
    rows_c = which(allmatch1==case3)
    for (k in rows_c){
      if (k%in%rows) next
      rows = append(rows, k)
      case4 = as.numeric(matchid[k,2])
      rows = append(rows, which(allmatch1==case4))
      rows = append(rows, which(allmatch2==case4))
      rows = unique(rows)
    }
    rows_c = which(allmatch2==case3)
    for (k in rows_c){
      if (k%in%rows) next
      rows = append(rows, k)
      case4 = as.numeric(matchid[k,1])
      rows = append(rows, which(allmatch1==case4))
      rows = append(rows, which(allmatch2==case4))
      rows = unique(rows)
    }
  }
  rows = unique(rows)
  if (length(rows)==2) family2 = rbind(family2, matchid[rows,])
  if (length(rows)==3) family3 = rbind(family3, matchid[rows,])
  if (length(rows)==4) family4 = rbind(family4, matchid[rows,])
  if (length(rows)==5) family5 = rbind(family5, matchid[rows,])
  if (length(rows)==6) family6 = rbind(family6, matchid[rows,])
  if (length(rows)==7) family7 = rbind(family7, matchid[rows,])
  matchid[rows,] = 0
}

n2 = nrow(family2)/2
family2 = as.matrix(family2)
n3 = nrow(family3)/3
family3 = as.matrix(family3)
n4 = nrow(family4)/4
family4 = as.matrix(family4)
n5 = nrow(family5)/5
family5 = as.matrix(family5)
n6 = nrow(family6)/6
family6 = as.matrix(family6)
n7 = nrow(family7)/7
family7 = as.matrix(family7)

## Families with 2 members
n1 = nrow(family_2)
family_2 = as.matrix(family_2)
family_2member = NULL
for (i in 1:n1){
  members = family_2[i,1:2]
  c = family_2[i,5]
  members_cov = matrix(c(1,c,c,1),2,2)
  family_2member = rbind(family_2member, cbind(members,members_cov))
}

## Families with 3 members
family_3member = NULL
for (i in 1:n2){
  pairs = family2[(2*i-1):(2*i),]
  members = unique(c(pairs[,1:2]))
  if (length(members)==3){
    members_cov = diag(rep(1,3))
    for (k in 1:2){
      cov1 = which(members==pairs[k,1])
      cov2 = which(members==pairs[k,2])
      members_cov[cov1,cov2] = pairs[k,5]
      members_cov[cov2,cov1] = pairs[k,5]
    }
    family_3member = rbind(family_3member,cbind(members,members_cov))
  }
}
for (i in 1:n3){
  pairs = family3[(3*i-2):(3*i),]
  members = unique(c(pairs[,1:2]))
  if (length(members)==3){
    members_cov = diag(rep(1,3))
    for (k in 1:3){
      cov1 = which(members==pairs[k,1])
      cov2 = which(members==pairs[k,2])
      members_cov[cov1,cov2] = pairs[k,5]
      members_cov[cov2,cov1] = pairs[k,5]
    }
    family_3member = rbind(family_3member,cbind(members,members_cov))
  }
}

## Families with 4 members
family_4member = NULL
for (i in 1:n3){
  pairs = family3[(3*i-2):(3*i),]
  members = unique(c(pairs[,1:2]))
  if (length(members)==4){
    members_cov = diag(rep(1,4))
    for (k in 1:3){
      cov1 = which(members==pairs[k,1])
      cov2 = which(members==pairs[k,2])
      members_cov[cov1,cov2] = pairs[k,5]
      members_cov[cov2,cov1] = pairs[k,5]
    }
    family_4member = rbind(family_4member,cbind(members,members_cov))
  }
}
for (i in 1:n4){
  pairs = family4[(4*i-3):(4*i),]
  members = unique(c(pairs[,1:2]))
  if (length(members)==4){
    members_cov = diag(rep(1,4))
    for (k in 1:4){
      cov1 = which(members==pairs[k,1])
      cov2 = which(members==pairs[k,2])
      members_cov[cov1,cov2] = pairs[k,5]
      members_cov[cov2,cov1] = pairs[k,5]
    }
    family_4member = rbind(family_4member,cbind(members,members_cov))
  }
}
for (i in 1:n5){
  pairs = family5[(5*i-4):(5*i),]
  members = unique(c(pairs[,1:2]))
  if (length(members)==4){
    members_cov = diag(rep(1,4))
    for (k in 1:5){
      cov1 = which(members==pairs[k,1])
      cov2 = which(members==pairs[k,2])
      members_cov[cov1,cov2] = pairs[k,5]
      members_cov[cov2,cov1] = pairs[k,5]
    }
    family_4member = rbind(family_4member,cbind(members,members_cov))
  }
}
for (i in 1:n6){
  pairs = family6[(6*i-5):(6*i),]
  members = unique(c(pairs[,1:2]))
  if (length(members)==4){
    members_cov = diag(rep(1,4))
    for (k in 1:6){
      cov1 = which(members==pairs[k,1])
      cov2 = which(members==pairs[k,2])
      members_cov[cov1,cov2] = pairs[k,5]
      members_cov[cov2,cov1] = pairs[k,5]
    }
    family_4member = rbind(family_4member,cbind(members,members_cov))
  }
}

## Families with 5 members
family_5member = NULL
for (i in 1:n4){
  pairs = family4[(4*i-3):(4*i),]
  members = unique(c(pairs[,1:2]))
  if (length(members)==5){
    members_cov = diag(rep(1,5))
    for (k in 1:4){
      cov1 = which(members==pairs[k,1])
      cov2 = which(members==pairs[k,2])
      members_cov[cov1,cov2] = pairs[k,5]
      members_cov[cov2,cov1] = pairs[k,5]
    }
    family_5member = rbind(family_5member,cbind(members,members_cov))
  }
}
for (i in 1:n5){
  pairs = family5[(5*i-4):(5*i),]
  members = unique(c(pairs[,1:2]))
  if (length(members)==5){
    members_cov = diag(rep(1,5))
    for (k in 1:5){
      cov1 = which(members==pairs[k,1])
      cov2 = which(members==pairs[k,2])
      members_cov[cov1,cov2] = pairs[k,5]
      members_cov[cov2,cov1] = pairs[k,5]
    }
    family_5member = rbind(family_5member,cbind(members,members_cov))
  }
}
for (i in 1:n6){
  pairs = family6[(6*i-5):(6*i),]
  members = unique(c(pairs[,1:2]))
  if (length(members)==5){
    members_cov = diag(rep(1,5))
    for (k in 1:6){
      cov1 = which(members==pairs[k,1])
      cov2 = which(members==pairs[k,2])
      members_cov[cov1,cov2] = pairs[k,5]
      members_cov[cov2,cov1] = pairs[k,5]
    }
    family_5member = rbind(family_5member,cbind(members,members_cov))
  }
}
for (i in 1:n7){
  pairs = family7[(7*i-6):(7*i),]
  members = unique(c(pairs[,1:2]))
  if (length(members)==5){
    members_cov = diag(rep(1,5))
    for (k in 1:7){
      cov1 = which(members==pairs[k,1])
      cov2 = which(members==pairs[k,2])
      members_cov[cov1,cov2] = pairs[k,5]
      members_cov[cov2,cov1] = pairs[k,5]
    }
    family_5member = rbind(family_5member,cbind(members,members_cov))
  }
}
n5 = nrow(family_5member)/5
min.kinship = keep.kinship = rep(0,n5)
for (i in 1:n5){
  ind = 5*(i-1) + 1:5
  mm = family_5member[ind,1]
  mc = family_5member[ind,2:6]
  dc = which.min(rowSums(mc))
  min.kinship[i] = (rowSums(mc)[dc]-1)/4
  keep.kinship[i] = (sum(rowSums(mc)[-dc])-1)/4/4
  family_4member = rbind(family_4member,cbind(mm[-dc],mc[-dc,-dc]))
}

family_all = c(family_2member[,1],family_3member[,1],family_4member[,1])
N = length(family_all)
n2 = nrow(family_2member)
n3 = nrow(family_3member)
n4 = nrow(family_4member)
KM = matrix(0,N,4)
KM[1:n2,1:2] = family_2member[,2:3]
KM[(n2+1):(n2+n3),1:3] = family_3member[,2:4]
KM[(n2+n3+1):(n2+n3+n4),1:4] = family_4member[,2:5]
n2 = n2/2
n3 = n3/3
n4 = n4/4
familyID = c(rep(1:n2, each=2), rep((n2+1):(n2+n3), each=3), 
             rep((n2+n3+1):(n2+n3+n4), each=4)) - 1
n = n2+n3+n4

# n2=5091, n3=632, n4=114
# 4 blocks
familyID1 = c(rep(1:(n2/4), each=2), rep((n2+1):(n2+n3/4), each=3), 
              rep(floor(n2+n3+1):(n2+n3+n4/4), each=4)) - 1
familyID2 = c(rep(floor(n2/4+1):(n2/4*2), each=2), 
              rep(floor(n2+1+n3/4):(n2+n3/4*2), each=3), 
              rep(floor(n2+n3+1+n4/4):(n2+n3+n4/4*2), each=4)) - 1
familyID3 = c(rep(floor(n2/4*2+1):(n2/4*3), each=2), 
              rep(floor(n2+1+n3/4*2):(n2+n3/4*3), each=3), 
              rep(floor(n2+n3+1+n4/4*2):(n2+n3+n4/4*3), each=4)) - 1
familyID4 = c(rep(floor(n2/4*3+1):(n2), each=2), 
              rep(floor(n2+1+n3/4*3):(n2+n3/4*4), each=3), 
              rep(floor(n2+n3+1+n4/4*3):(n2+n3+n4), each=4)) - 1
pairs = kinpairs[kinpairs[,1]%in%family_all&kinpairs[,2]%in%family_all,]

## load UKB data

rm(df)
dat <- read.csv('ukbdata.csv',head=TRUE,sep=',')
subdat = dat[dat$eid%in%family_all,]
inddat = dat[!dat$eid%in%allmatch,]
rm(dat)
n1 = nrow(inddat)
KM_ind = matrix(0,n1,4)
KM_ind[,1] = 1
KM = rbind(KM, KM_ind)
familyID_ind = max(familyID) + 1:n1
familyID_ind1 = max(familyID) + 1:(n1/4)
familyID_ind2 = max(familyID) + floor(n1/4+1):(n1/4*2)
familyID_ind3 = max(familyID) + floor(n1/4*2+1):(n1/4*3)
familyID_ind4 = max(familyID) + floor(n1/4*3+1):(n1)
eid = data.frame(family_all)
colnames(eid) = 'eid'
famdat = inner_join(eid, subdat, by='eid')
famdat$familyID = familyID
inddat$familyID = familyID_ind
familyID = c(familyID, familyID_ind)
familyID1 = c(familyID1, familyID_ind1)
familyID2 = c(familyID2, familyID_ind2)
familyID3 = c(familyID3, familyID_ind3)
familyID4 = c(familyID4, familyID_ind4)
cn = colnames(famdat)
alldat = famdat
#alldat = rbind(famdat, inddat)
#alldat$familyID = familyID
n = n1+n2+n3+n4
size = nrow(alldat)
eid = alldat$eid

rm(famdat)
rm(inddat)
# covariates
# sex, age, education, income
sex = alldat[,'X31.0.0']
education = as.numeric(alldat[,"X845.0.0"]>16)
education[is.na(education)] = 0
income = as.numeric(alldat[,"X738.0.0"]>=4)
income[is.na(income)] = 0
age = as.numeric(alldat[,"X21022.0.0"])
age[is.na(age)] = mean(age,na.rm=TRUE)
gpc = alldat[,c("X22009.0.1","X22009.0.2","X22009.0.3","X22009.0.4","X22009.0.5",
                "X22009.0.6","X22009.0.7","X22009.0.8","X22009.0.9","X22009.0.10")]
gpc[is.na(gpc)] = 0
X = cbind(sex, age, education, income, gpc)
p = ncol(X)

birthyear = alldat[,'X34.0.0']
birthmonth = alldat[,'X52.0.0']
birthday = as.Date(paste0(birthyear,'-',birthmonth,'-',1))


## phenotypes

library(readxl)
library(stringr)
library(RNOmni)
dict1 = read_xlsx('dict.xlsx', sheet='continuous')
dict1.name = (dict1[,2])[[1]]
dict1 = (dict1[,1])[[1]]
dict2 = read_xlsx('dict.xlsx', sheet='binary')
dict2.name = (dict2[,2])[[1]]
dict2 = (dict2[,1])[[1]]
dict3 = read_xlsx('dict.xlsx', sheet='survival')
dict3.name = (dict3[,2])[[1]]
dict3 = (dict3[,1])[[1]]
dict4 = read_xlsx('dict.xlsx', sheet='ordinal')
dict4.name = (dict4[,2])[[1]]
dict4 = (dict4[,1])[[1]]
dict0 = read_xlsx('dict.xlsx', sheet='covariates')
dict0 = (dict0[,1])[[1]]

# continuous
con.name = sapply(dict1, function(l) paste0('X',l,'.0'))
con.dat = sapply(con.name, function(l) 
  rowMeans(as.matrix(alldat[,str_detect(cn,l)]), na.rm=TRUE))
colnames(con.dat) = dict1
prop = sapply(1:ncol(con.dat), function(l) mean(is.na(con.dat[,l])))
con.name2 = sapply(dict1, function(l) paste0('X',l,'.2'))
con.dat2 = sapply(con.name2, function(l) 
  rowMeans(as.matrix(alldat[,str_detect(cn,l)]), na.rm=TRUE))
colnames(con.dat2) = dict1
con.dat[,prop==1] = con.dat2[,prop==1]
con.dat = con.dat[,prop<=0.95]
dict1.name = dict1.name[prop<=0.95]
means.con = colMeans(con.dat, na.rm=TRUE)
means.con.fam = colMeans(con.dat[1:N,], na.rm=TRUE)
sds.con = apply(con.dat, 2, sd, na.rm=TRUE)
sds.con.fam = apply(con.dat[1:N,], 2, sd, na.rm=TRUE)
for (j in 1:ncol(con.dat)){
  index = (!is.na(con.dat[,j]))
  con.dat[index,j] = RankNorm(con.dat[index,j])
}
dim(con.dat)
K1 = ncol(con.dat)

# binary
bin.name = sapply(dict2, function(l) paste0('X',l,'.0'))
bin.dat = sapply(bin.name, function(l) 
  rowMeans(as.matrix(alldat[,str_detect(cn,l)]), na.rm=TRUE))
colnames(bin.dat) = dict2
bin.dat[bin.dat<0] = NaN
fts = sapply(1:ncol(bin.dat), function(l) length(na.omit(unique(bin.dat[,l]))))
bincheck = fts==2
bin.dat = bin.dat[,bincheck]
dict2 = dict2[bincheck]
dict2.name = dict2.name[bincheck]
bin.dat = bin.dat - rep(1,nrow(bin.dat))%*%t(apply(bin.dat,2,min,na.rm=TRUE))
prop = sapply(1:ncol(bin.dat), function(l) mean(is.na(bin.dat[,l])))
prop0 = colMeans(bin.dat,na.rm=TRUE)
bin.dat = bin.dat[,prop<=0.95&prop0>=0.05&prop0<=0.95]
dict2.name = dict2.name[prop<=0.95&prop0>=0.05&prop0<=0.95]
K2 = ncol(bin.dat)
dim(bin.dat)

#binary (ICD)
pheicd = read_xlsx('phecode_icd10.xlsx',sheet='phecode_icd10')
icd10.begin = which(cn=='X41270.0.0')
maxicd = sum(str_detect(cn,'X41270.0'))
icd = alldat[,which(cn=='X41270.0.0')-1+1:maxicd]
pheicd_phe = substr(pheicd[,5][[1]],1,3)
pheicd_icd = pheicd[,1][[1]]
phelist = unique(pheicd_phe)
phelist = na.omit(unique(sapply(phelist, function(l) substr(l,1,3))))
disease = matrix(0, nrow=size, ncol=length(phelist))
rownames(disease) = alldat$eid
colnames(disease) = sapply(phelist, function(l) paste0('PHE',l))
disease.des = sapply(phelist, function(l) (pheicd[,6][[1]])[which(pheicd_phe==l)[1]])
icdtophe = function(x){ pheicd_phe[which(pheicd_icd==x)[1]] }
for (i in 1:size){
  for (j in 1:maxicd){
    icdij = icd[i,j]
    if (icdij=='') break
    x = icdtophe(icdij)
    index = which(phelist==x)
    if (length(na.omit(index))==1) disease[i,index] = 1
  }
}
prev = colMeans(disease)
disease = disease[,prev>=0.05&prev<=0.95]
disease.des = disease.des[prev>=0.05&prev<=0.95]
bin.dat = cbind(bin.dat, disease)
dim(bin.dat)
K2 = ncol(bin.dat)

#####
## Ordinal
ord.name = sapply(dict4, function(l) paste0('X',l,'.0'))
ord.dat = sapply(ord.name, function(l) 
  rowMeans(as.matrix(alldat[,str_detect(cn,l)]), na.rm=TRUE))
colnames(ord.dat) = dict4
ord.dat[ord.dat<0] = NaN
ord.dat = round(ord.dat,0)
for (j in 1:ncol(ord.dat)){
  ord.dat[,j] = ord.dat[,j] - min(c(99,ord.dat[,j]),na.rm=TRUE) + 1
}
ord.dat[,'4803'] = ord.dat[,'4803'] - 10
ord.dat[,'4814'] = ord.dat[,'4814'] - 6
ord.dat[,'4825'] = ord.dat[,'4825'] - 10
ord.dat[,'4836'] = ord.dat[,'4836'] - 10
ord.dat[ord.dat<=0] = 1
ord.dat[ord.dat>90] = 3
prop = sapply(1:ncol(ord.dat), function(l) mean(is.na(ord.dat[,l])))
ord.dat = ord.dat[,prop<=0.95]
dict4.name = dict4.name[prop<=0.95]
dim(ord.dat)
K3 = ncol(ord.dat)

# survival
studyend = round((as.Date('2024-01-31')-birthday)/365.25,1)
sur.name = sapply(dict3, function(l) paste0('X',l,'.0.0'))
datetoage = function(x){
  if (length(x)==0) return(rep(NA,nrow(x)))
  if (is.numeric(x)) return(x)
  if (is.character(x)) return(round((as.Date(x) - birthday)/365.25,1))
  if (is.data.frame(x)) return(as.numeric(x[[1]]))
}
sur.dat = sapply(sur.name, function(l) datetoage(alldat[,str_detect(cn,l)]))
colnames(sur.dat) = dict3
sur.dat[sur.dat<0] = NA
sur.D.dat = !is.na(sur.dat)
prev = colMeans(sur.D.dat)
for (j in 1:ncol(sur.dat)){
  sur.dat[is.na(sur.dat[,j]),j] = studyend[is.na(sur.dat[,j])]
}
sur.dat = sur.dat[,prev>=0.01]
sur.T.dat = sur.dat
sur.T.dat[sur.T.dat==0] = sur.T.dat[sur.T.dat==0]+0.1
sur.D.dat = sur.D.dat[,prev>=0.01]
dict3.name = dict3.name[prev>=0.01]
prev = prev[prev>=0.01]
dim(sur.T.dat)
K4 = ncol(sur.T.dat)
K = K1+K2+K3+K4
phes = c(dict1.name,dict2.name,disease.des,dict4.name,dict3.name)
fields = c(colnames(con.dat),colnames(bin.dat),colnames(ord.dat),colnames(sur.T.dat))

uids = read_xlsx('phes.xlsx')[,1]
uids = uids[[1]]
#sub = c(21001,50,4079,4080,3063,30690,30760,30780,30750,30870,
#        30160,30150,30020,30120,30130,30140,30080,30090,30010,30000,30250)
#sub = sapply(sub, function(k) which(as.numeric(colnames(con.dat))==k))
#con.all = con.dat[,sub]
#save(con.all, familyID, KM, pairs, X, n1,n2,n3,n4, eid,
#     file='conall.Rdata')
meaningful = which(fields%in%uids)
phes = phes[(1:K)%in%meaningful]
con.dat = con.dat[,(1:K1)%in%meaningful]
bin.dat = bin.dat[,(1:K2+K1)%in%meaningful]
ord.dat = ord.dat[,(1:K3+K1+K2)%in%meaningful]
sur.T.dat = sur.T.dat[,(1:K4+K1+K2+K3)%in%meaningful]
sur.D.dat = sur.D.dat[,(1:K4+K1+K2+K3)%in%meaningful]
familyID = familyID[1:size]
KM = KM[1:size,]
X = as.matrix(X[1:size,])
save(con.dat, familyID, KM, pairs, eid, X, n1,n2,n3,n4,
     file='condat.Rdata')
save(bin.dat, familyID, KM, pairs, eid, X, n1,n2,n3,n4,
     file='bindat.Rdata')
save(ord.dat, familyID, KM, pairs, eid, X, n1,n2,n3,n4,
     file='orddat.Rdata')
save(sur.T.dat,sur.D.dat, familyID, KM, pairs, eid, X, n1,n2,n3,n4,
     file='surdat.Rdata')

# save data
save(con.dat,bin.dat,ord.dat,sur.T.dat,sur.D.dat,
     X,KM,n1,n2,n3,n4,n,N,size,
     eid,pairs,familyID,
     file='ukb_alldat.Rdata')

## Summary statistics (Table 2)
dat = cbind(con.dat,bin.dat,ord.dat,sur.D.dat)
fields = colnames(dat)
K1 = ncol(con.dat)
K2 = ncol(bin.dat)
K3 = ncol(ord.dat)
K4 = ncol(sur.T.dat)
K = K1+K2+K3+K4
type = c(rep('Continuous',K1),
         rep('Binary',K2),
         rep('Ordinal',K3),
         rep('Survival',K4))
# non-missing proportion
nmp = round(colMeans(!is.na(dat)),3)
nmp[K1+K2+K3+1:K4] = colMeans(sur.D.dat)
# mean in complete cases
means = colMeans(dat, na.rm=TRUE)
means[1:K1] = means.con[1:length(means.con)%in%meaningful]
means[K1+K2+K3+1:K4] = sapply(1:K4,
                              function(k) mean(sur.T.dat[sur.D.dat[,k]==1,k]))
means = round(means, 3)
# variance in complete cases
TV = apply(dat,2,var,na.rm=TRUE)
TV[K1+K2+K3+1:K4] = sapply(1:K4, 
                           function(k) var(sur.T.dat[sur.D.dat[,k]==1,k]))
TV = round(sqrt(TV), 3)
TV[1:K1] = round(sds.con[1:length(sds.con)%in%meaningful], 3)
table1 = data.frame(cbind(fields,phes,type,nmp,means,TV))

comfam1 = comfam2 = comfam3 = comfam4 = rep(0,K)
dat[,K1+K2+K3+1:K4] = sur.D.dat
for (k in 1:K){
  if (k>K1+K2+K3){
    faml = familyID[!is.nan(dat[,k])&dat[,k]==1]
  } else {
    faml = familyID[!is.nan(dat[,k])]
  }
  n.com = length(faml)
  faml = faml[faml<n2+n3+n4]
  fam = unique(faml)
  comfam2[k] = sum(sapply(fam, function(i) sum(faml==i)==2))
  comfam3[k] = sum(sapply(fam, function(i) sum(faml==i)==3))
  comfam4[k] = sum(sapply(fam, function(i) sum(faml==i)==4))
  comfam1[k] = n.com - comfam2[k] - comfam3[k] - comfam4[k]
  cat(k,'')
}
table1 = data.frame(cbind(table1, comfam1,comfam2,comfam3,comfam4))
colnames(table1) = c('Field UID','Phenotype','Data type',
                     'Complete case proportion','Mean','SD',
                     'Family of 1','Family of 2','Family of 3','Family of 4')
rownames(table1) = 1:nrow(table1)
write.csv(table1, file='Table2.csv')
write.csv(table2, file='Table3.csv')
save(n2,n3,n4,N,familyID,eid,pairs,
     file='ukb_famdat.Rdata')
