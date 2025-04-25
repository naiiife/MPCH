# Ordinal

cppFunction('
  List cal_ord(IntegerVector Y, NumericMatrix X,
              IntegerVector familyID, NumericMatrix KM,
              NumericVector par0){
  Environment glob = Environment::global_env();
  Function ginv("ginv");
  Function mprod("mprod");
  NumericVector w_list = glob.get("w_list");
  NumericVector b_list = glob.get("b_list");
  int ngrid = w_list.size();
  int p = X(0,_).size();
  int n0 = min(familyID);
  int n = max(familyID) + 1;
  int N = familyID.size();
  IntegerVector Fam_size;
  IntegerVector Fam_list;
  for (int i=n0; i<n; i++){
    int ni = sum(ifelse(familyID==i,1,0));
    if(ni>0) {
    Fam_list.push_back(i);
    Fam_size.push_back(ni);
    }
  }
  n = Fam_list.size();
  int L = 1;
  for (int i=0; i<N; i++){
  int Yi = Y[i];
    if(std::isnan(Yi)) {Y[i]=0;}
    if(Yi>L) {L=Yi;}
  }
  Y = Y - 1;
  std::cout << L << std::endl;
  
  NumericVector res (p+1+L);
  NumericVector d1 (N);
  NumericVector d2 (N);
  NumericVector d3 (N);
  NumericVector d4 (N);
  NumericVector d5 (N);
  NumericVector Ee2 (n);
  NumericVector EepsGeps (n);
  NumericMatrix dl2dalpha (p,p);
  NumericVector dldd (L-1);
  NumericVector dl2dd (L-1);
  
  double theta_k = par0[0];
  double gamma_kk = par0[1];
  NumericVector alpha_k (p);
  for (int k=0; k<p; k++){
  alpha_k[k] = par0[2+k];
  }
  NumericVector delta_k (L-1);
  for (int k=0; k<L-1; k++){
  delta_k[k] = par0[2+p+k];
  }
  
  int iterrun = 0;
  double tol = 1;
  while (tol > 0.0001){
  NumericVector res0;
  res0.assign(res.begin(),res.end());
  NumericVector delta (L+1);
  delta[0] = -999.0;
  delta[L] = 999.0;
  for (int l=0; l<L-1; l++){
    delta[l+1] = delta_k[l];
  }
  int id_begin = 0;
  for (int i=0; i<n; i++){
  int ni = Fam_size[i];
  if(ni==0) {continue;}
  IntegerVector Yi (ni);
  NumericMatrix Xi (ni, p);
  NumericVector dli (ni);
  NumericVector dlil (ni);
  NumericVector Xb (ni);
  for (int j=0; j<ni; j++){
    Xb[j] = sum(X(id_begin+j,_)*alpha_k);
    int Yij = Y[id_begin+j];
    Yi[j] = Yij;
    if(Yij<0){
    dli[j] = 999.0;
    dlil[j] = -999.0;
    } else {
    dli[j] = delta[Yij+1];
    dlil[j] = delta[Yij];
    }
  }
  
  NumericVector num1 (ni);
  NumericVector num2 (ni);
  NumericVector num3 (ni);
  NumericVector num4 (ni);
  NumericVector num5 (ni);
  double num6 = 0.0;
  double num7 = 0.0;
  double den = 0.0;
  NumericMatrix KMi (ni,ni);
  NumericMatrix Girt (ni,ni);
  for (int r=0; r<ni; r++){
  for (int s=0; s<ni; s++){
    KMi(r,s) = KM(id_begin+r,s);
  }
  }
  
  if(ni == 1) {
    for (int b1 = 0; b1 < ngrid; b1++) {
      double b_eps1 = b_list[b1];
      double w_eps1 = w_list[b1];
      double eps = b_eps1 * sqrt(2.0*gamma_kk);
        for (int b = 0; b < ngrid; b++) {
          double b_e = b_list[b];
          double w_e = w_list[b];
          double e = b_e*1.4142136*theta_k;
          NumericVector linear = dli - Xb - e - eps;
          NumericVector linearl = dlil - Xb - e - eps;
          NumericVector Plinear (ni);
          NumericVector plinear = exp(-linear*linear*0.5)/2.506628;
          NumericVector plinearl = exp(-linearl*linearl*0.5)/2.506628;
          for (int j=0; j<ni; j++){
          if (Yi[j]<0){
          plinear[j]=0;
          plinearl[j]=0;
          }
          double pnorm1 = erf(linear[j]/1.4142136);
          double pnorm2 = erf(linearl[j]/1.4142136);
          Plinear[j] = (pnorm1-pnorm2)/2;
          }
          NumericVector dlinear = - linear * plinear;
          NumericVector dlinearl = - linearl * plinearl;
          for (int j=0; j<ni; j++){
          if(Plinear[j]==0) { Plinear[j] = 0.0000001; }
          }
          NumericVector loglik = log(Plinear);
          for (int j=0; j<ni; j++){
            if(Yi[j]<0) {loglik[j]=0;}
          }
          double prodlik = exp(sum(loglik));
          double wdens = w_e * w_eps1 * prodlik;
          den = den + wdens;
          num1 = num1 + wdens * plinear / Plinear;
          num2 = num2 + wdens * plinearl / Plinear;
          num3 = num3 + wdens * (dlinear/Plinear - (plinear/Plinear)*(plinear/Plinear));
          num4 = num4 + wdens * (dlinearl/Plinear + (plinearl/Plinear)*(plinearl/Plinear));
          NumericVector ddP = (dlinear-dlinearl)/Plinear;
          NumericVector ppP = (plinear-plinearl)/Plinear;
          num5 = num5 + wdens * (ddP - ppP*ppP);
          num6 = num6 + wdens * e*e;
          num7 = num7 + wdens * eps*eps;
      }
    }
  }
  if(ni == 2) {
    Girt(0,0) = 1.0;
    double g12 = KMi(1,0);
    Girt(1,0) = g12;
    Girt(1,1) = sqrt(1.0-g12*g12);
    for (int b1 = 0; b1 < ngrid; b1++) {
      double b_eps1 = b_list[b1];
      double w_eps1 = w_list[b1];
      for (int b2 = 0; b2 < ngrid; b2++) {
        double b_eps2 = b_list[b2];
        double w_eps2 = w_list[b2];
        NumericVector b_eps = NumericVector::create(b_eps1, b_eps2);
        NumericVector eps (ni);
        for (int j=0; j<ni; j++){
          eps[j] = sum(Girt(j,_)*b_eps);
        }
        eps = eps * sqrt(2.0*gamma_kk);
        for (int b = 0; b < ngrid; b++) {
          double b_e = b_list[b];
          double w_e = w_list[b];
          double e = b_e*1.4142136*theta_k;
          NumericVector linear = dli - Xb - e - eps;
          NumericVector linearl = dlil - Xb - e - eps;
          NumericVector Plinear (ni);
          NumericVector plinear = exp(-linear*linear*0.5)/2.506628;
          NumericVector plinearl = exp(-linearl*linearl*0.5)/2.506628;
          for (int j=0; j<ni; j++){
          if (Yi[j]<0){
          plinear[j]=0;
          plinearl[j]=0;
          }
          double pnorm1 = erf(linear[j]/1.4142136);
          double pnorm2 = erf(linearl[j]/1.4142136);
          Plinear[j] = (pnorm1-pnorm2)/2;
          }
          NumericVector dlinear = - linear * plinear;
          NumericVector dlinearl = - linearl * plinearl;
          for (int j=0; j<ni; j++){
          if(Plinear[j]==0) { Plinear[j] = 0.0000001; }
          }
          NumericVector loglik = log(Plinear);
          for (int j=0; j<ni; j++){
            if(Yi[j]<0) {loglik[j]=0;}
          }
          double prodlik = exp(sum(loglik));
          double wdens = w_e * w_eps1 * w_eps2 * prodlik;
          den = den + wdens;
          num1 = num1 + wdens * plinear / Plinear;
          num2 = num2 + wdens * plinearl / Plinear;
          num3 = num3 + wdens * (dlinear/Plinear - (plinear/Plinear)*(plinear/Plinear));
          num4 = num4 + wdens * (dlinearl/Plinear + (plinearl/Plinear)*(plinearl/Plinear));
          NumericVector ddP = (dlinear-dlinearl)/Plinear;
          NumericVector ppP = (plinear-plinearl)/Plinear;
          num5 = num5 + wdens * (ddP - ppP*ppP);
          num6 = num6 + wdens * e*e;
          num7 = num7 + wdens * 2.0*gamma_kk*sum(b_eps*b_eps);
        }
      }
    }
  }
  if(ni == 3) {
    Girt(0,0) = 1.0;
    double g12 = KMi(1,0);
    double g13 = KMi(2,0);
    Girt(1,0) = g12;
    Girt(2,0) = g13;
    double l22 = sqrt(1.0-g12*g12);
    Girt(1,1) = l22;
    double l23 = (KMi(2,1)-g12*g13)/l22;
    Girt(2,1) = l23;
    Girt(2,2) = sqrt(1.0-g13*g13-l23*l23);
    for (int b1 = 0; b1 < ngrid; b1++) {
      double b_eps1 = b_list[b1];
      double w_eps1 = w_list[b1];
      for (int b2 = 0; b2 < ngrid; b2++) {
        double b_eps2 = b_list[b2];
        double w_eps2 = w_list[b2];
        for (int b3 = 0; b3 < ngrid; b3++) {
        double b_eps3 = b_list[b3];
        double w_eps3 = w_list[b3];
        NumericVector b_eps = NumericVector::create(b_eps1, b_eps2, b_eps3);
        NumericVector eps (ni);
        for (int j=0; j<ni; j++){
          eps[j] = sum(Girt(j,_)*b_eps);
        }
        eps = eps * sqrt(2.0*gamma_kk);
        for (int b = 0; b < ngrid; b++) {
          double b_e = b_list[b];
          double w_e = w_list[b];
          double e = b_e*1.4142136*theta_k;
          NumericVector linear = dli - Xb - e - eps;
          NumericVector linearl = dlil - Xb - e - eps;
          NumericVector Plinear (ni);
          NumericVector plinear = exp(-linear*linear*0.5)/2.506628;
          NumericVector plinearl = exp(-linearl*linearl*0.5)/2.506628;
          for (int j=0; j<ni; j++){
          if (Yi[j]<0){
          plinear[j]=0;
          plinearl[j]=0;
          }
          double pnorm1 = erf(linear[j]/1.4142136);
          double pnorm2 = erf(linearl[j]/1.4142136);
          Plinear[j] = (pnorm1-pnorm2)/2;
          }
          NumericVector dlinear = - linear * plinear;
          NumericVector dlinearl = - linearl * plinearl;
          for (int j=0; j<ni; j++){
          if(Plinear[j]==0) { Plinear[j] = 0.0000001; }
          }
          NumericVector loglik = log(Plinear);
          for (int j=0; j<ni; j++){
            if(Yi[j]<0) {loglik[j]=0;}
          }
          double prodlik = exp(sum(loglik));
          double wdens = w_e * w_eps1 * w_eps2 * w_eps3 * prodlik;
          den = den + wdens;
          num1 = num1 + wdens * plinear / Plinear;
          num2 = num2 + wdens * plinearl / Plinear;
          num3 = num3 + wdens * (dlinear/Plinear - (plinear/Plinear)*(plinear/Plinear));
          num4 = num4 + wdens * (dlinearl/Plinear + (plinearl/Plinear)*(plinearl/Plinear));
          NumericVector ddP = (dlinear-dlinearl)/Plinear;
          NumericVector ppP = (plinear-plinearl)/Plinear;
          num5 = num5 + wdens * (ddP - ppP*ppP);
          num6 = num6 + wdens * e*e;
          num7 = num7 + wdens * 2.0*gamma_kk*sum(b_eps*b_eps);
          }
        }
      }
    }
  }
  if(ni == 4) {
    Girt(0,0) = 1.0;
    double g12 = KMi(1,0);
    double g13 = KMi(2,0);
    Girt(1,0) = g12;
    double l22 = sqrt(1.0-g12*g12);
    Girt(1,1) = l22;
    double l23 = (KMi(2,1)-g12*g13)/l22;
    double l33 = sqrt(1.0-g13*g13-l23*l23);
    Girt(2,0) = g13;
    Girt(2,1) = l23;
    Girt(2,2) = l33;
    double g14 = KMi(3,0);
    double g24 = KMi(3,1);
    double l24 = (g24-g12*g14)/l22;
    double l34 = (KMi(3,2)-g13*g14-l23*l24)/l33;
    Girt(3,0) = g14;
    Girt(3,1) = l24;
    Girt(3,2) = l34;
    Girt(3,3) = sqrt(1.0-g14*g14-l24*l24-l34*l34);
    for (int b1 = 0; b1 < ngrid; b1++) {
      double b_eps1 = b_list[b1];
      double w_eps1 = w_list[b1];
      for (int b2 = 0; b2 < ngrid; b2++) {
        double b_eps2 = b_list[b2];
        double w_eps2 = w_list[b2];
        for (int b3 = 0; b3 < ngrid; b3++) {
        double b_eps3 = b_list[b3];
        double w_eps3 = w_list[b3];
        for (int b4 = 0; b4 < ngrid; b4++) {
        double b_eps4 = b_list[b4];
        double w_eps4 = w_list[b4];
        NumericVector b_eps = NumericVector::create(b_eps1, b_eps2, b_eps3, b_eps4);
        NumericVector eps (ni);
        for (int j=0; j<ni; j++){
          eps[j] = sum(Girt(j,_)*b_eps);
        }
        eps = eps * sqrt(2.0*gamma_kk);
        for (int b = 0; b < ngrid; b++) {
          double b_e = b_list[b];
          double w_e = w_list[b];
          double e = b_e*1.4142136*theta_k;
          NumericVector linear = dli - Xb - e - eps;
          NumericVector linearl = dlil - Xb - e - eps;
          NumericVector Plinear (ni);
          NumericVector plinear = exp(-linear*linear*0.5)/2.506628;
          NumericVector plinearl = exp(-linearl*linearl*0.5)/2.506628;
          for (int j=0; j<ni; j++){
          if (Yi[j]<0){
          plinear[j]=0;
          plinearl[j]=0;
          }
          double pnorm1 = erf(linear[j]/1.4142136);
          double pnorm2 = erf(linearl[j]/1.4142136);
          Plinear[j] = (pnorm1-pnorm2)/2;
          }
          NumericVector dlinear = - linear * plinear;
          NumericVector dlinearl = - linearl * plinearl;
          for (int j=0; j<ni; j++){
          if(Plinear[j]==0) { Plinear[j] = 0.0000001; }
          }
          NumericVector loglik = log(Plinear);
          for (int j=0; j<ni; j++){
            if(Yi[j]<0) {loglik[j]=0;}
          }
          double prodlik = exp(sum(loglik));
          double wdens = w_e * w_eps1 * w_eps2 * w_eps3 * w_eps4 * prodlik;
          den = den + wdens;
          num1 = num1 + wdens * plinear / Plinear;
          num2 = num2 + wdens * plinearl / Plinear;
          num3 = num3 + wdens * (dlinear/Plinear - (plinear/Plinear)*(plinear/Plinear));
          num4 = num4 + wdens * (dlinearl/Plinear + (plinearl/Plinear)*(plinearl/Plinear));
          NumericVector ddP = (dlinear-dlinearl)/Plinear;
          NumericVector ppP = (plinear-plinearl)/Plinear;
          num5 = num5 + wdens * (ddP - ppP*ppP);
          num6 = num6 + wdens * e*e;
          num7 = num7 + wdens * 2.0*gamma_kk*sum(b_eps*b_eps);
          }
        }
      }
    }
  }
  }
  if(den>0){
  for (int j=0; j<ni; j++){
    d1[id_begin+j] = num1[j]/den;
    d2[id_begin+j] = num2[j]/den;
    d3[id_begin+j] = num3[j]/den;
    d4[id_begin+j] = num4[j]/den;
    d5[id_begin+j] = num5[j]/den;
  }
  Ee2[i] = num6/den;
  EepsGeps[i] = num7/den;
  }
  id_begin = id_begin + ni;
  }
  
  for (int l=0; l<L-1; l++){
    NumericVector Yl = ifelse(Y==l, 1.0, 0.0);
    NumericVector Yll = ifelse(Y==l+1, 1.0, 0.0);
    dldd[l] = sum(d1*Yl) - sum(d2*Yll);
    dl2dd[l] = sum(d3*Yl) - sum(d4*Yll);
  }
  NumericVector dldalpha (p);
  for (int k=0; k<p; k++){
    dldalpha[k] = sum(X(_,k)*(d2-d1));
    for (int r=0; r<p; r++){
    dl2dalpha(k,r) = sum(X(_,k)*X(_,r)*d5);
    }
  }
  NumericMatrix d2_1 = ginv(dl2dalpha);
  NumericVector d2d1 (p);
  for (int k=0; k<p; k++){
    d2d1[k] = sum(d2_1(k,_)*dldalpha);
  }
  alpha_k = alpha_k - 0.8 * d2d1;
  for (int r=0; r<p; r++){
    if(std::isnan(alpha_k[r])) {alpha_k[r]=0.0;}
  }
  delta_k = delta_k - 0.8 * dldd / dl2dd;
  for (int r=0; r<L-2; r++){
    if(std::isnan(delta_k[r])) delta_k[r]=-99;
    double difdelta = delta_k[r+1]-delta_k[r];
    if(difdelta<=0) {
    delta_k[r+1] = delta_k[r] + 0.001 + difdelta/3;
    delta_k[r] = delta_k[r] + difdelta/3*2;
    }
  }
  theta_k = sqrt(sum(Ee2) / n);
  gamma_kk = sum(EepsGeps) / N;
  res[0] = theta_k;
  res[1] = gamma_kk;
  for (int k=0; k<p; k++){
    res[k+2] = alpha_k[k];
  }
  for (int k=0; k<L-1; k++){
    res[k+2+p] = delta_k[k];
  }
  tol = max(sqrt((res - res0)*(res - res0)));
  iterrun = iterrun + 1;
  if(iterrun>=5000) {tol = 0.0;}
  std::cout << res[0] << " " << res[1] << " " << res[2] << " " << res[3] << " " << res[4] << std::endl;
  std::cout << delta_k[0] << " " << delta_k[1] << " " << delta_k[2] << " " << tol << std::endl;
  }
  int fail = 0;
  if(iterrun>=5000) {fail=1;}
  if(std::isnan(tol)) {fail=1;}
  
  // Calculate Variance
  
  NumericVector ld_t (n);
  if(theta_k>0){
  ld_t = Ee2/(theta_k*theta_k*theta_k) - 1.0/theta_k;
  }
  NumericVector Famsize = as<NumericVector>(Fam_size);
  NumericVector ld_g = -Famsize/(gamma_kk*2.0);
  ld_g = ld_g + EepsGeps/(gamma_kk*gamma_kk*2.0);
  NumericMatrix Eld (n,p+1+L);
  int id_begin = 0;
  for (int i=0; i<n; i++){
    int ni = Fam_size[i];
    if(ni==0) {continue;}
    Eld(i,0) = ld_t[i];
    Eld(i,1) = ld_g[i];
    for (int j=0; j<ni; j++){
    for (int k=0; k<p; k++){
      Eld(i,k+2) = Eld(i,k+2) + (d2[id_begin+j]-d1[id_begin+j]) * X(id_begin+j,k);
    }
    for (int l=0; l<L-1; l++){
      double d12 = 0.0;
      if (Y[id_begin+j]==l) d12 = d1[id_begin+j];
      if (Y[id_begin+j]==l+1) d12 = -d2[id_begin+j];
      Eld(i,l+p+2) = Eld(i,l+p+2) + d12;
    }
    }
    id_begin = id_begin + ni;
  }
  
  NumericMatrix El2 = mprod(transpose(Eld), Eld);
  double ldd_t = - 2.0 / (theta_k*theta_k) * n;
  if(theta_k==0) {ldd_t=0;}
  double ldd_g = - N / (gamma_kk*gamma_kk*2.0);
  NumericMatrix Ell (p+1+L,p+1+L);
  NumericVector Ddd (N);
  NumericVector delta (L+1);
  delta[0] = -999.0;
  delta[L] = 999.0;
  for (int l=0; l<L-1; l++){
    delta[l+1] = delta_k[l];
  }
  id_begin = 0;
  for (int i=0; i<n; i++){
  int ni = Fam_size[i];
  if(ni==0) {continue;}
  NumericVector Yi (ni);
  NumericMatrix Xi (ni, p);
  NumericVector dli (ni);
  NumericVector dlil (ni);
  NumericVector Xb (ni);
  for (int j=0; j<ni; j++){
    int Yij = Y[id_begin+j];
    Yi[j] = Yij;
    for (int k=0; k<p; k++){
    Xi(j,k) = X(id_begin+j, k);
    }
    Xb[j] = sum(X(id_begin+j,_)*alpha_k);
    if(Yij<0){
    dli[j] = 999.0;
    dlil[j] = -999.0;
    } else {
    dli[j] = delta[Yij+1];
    dlil[j] = delta[Yij];
    }
  }
  NumericMatrix num (p+1+L,p+1+L);
  double den = 0;
  NumericVector num8 (ni);
  NumericMatrix KMi (ni,ni);
  NumericMatrix Girt (ni,ni);
  for (int r=0; r<ni; r++){
  for (int s=0; s<ni; s++){
    KMi(r,s) = KM(id_begin+r,s);
  }
  }
  
  if(ni == 1) {
    for (int b1 = 0; b1 < ngrid; b1++) {
      double b_eps = b_list[b1];
      double w_eps1 = w_list[b1];
      double eps = b_eps * sqrt(2.0*gamma_kk);
        for (int b = 0; b < ngrid; b++) {
          double b_e = b_list[b];
          double w_e = w_list[b];
          double e = b_e*1.4142136*theta_k;
          NumericVector linear = dli - Xb - e - eps;
          NumericVector linearl = dlil - Xb - e - eps;
          NumericVector Plinear (ni);
          NumericVector plinear = exp(-linear*linear*0.5)/2.506628;
          NumericVector plinearl = exp(-linearl*linearl*0.5)/2.506628;
          for (int j=0; j<ni; j++){
          if (Yi[j]<0){
          plinear[j]=0;
          plinearl[j]=0;
          }
          double pnorm1 = erf(linear[j]/1.4142136);
          double pnorm2 = erf(linearl[j]/1.4142136);
          Plinear[j] = (pnorm1-pnorm2)/2;
          }
          NumericVector dlinear = - linear * plinear;
          NumericVector dlinearl = - linearl * plinearl;
          for (int j=0; j<ni; j++){
          if(Plinear[j]==0) { Plinear[j] = 0.0000001; }
          }
          NumericVector loglik = log(Plinear);
          for (int j=0; j<ni; j++){
            if(Yi[j]<0) {loglik[j]=0;}
          }
          double prodlik = exp(sum(loglik));
          double wdens = w_e * w_eps1 * prodlik;
          if(wdens>0){
          NumericVector ld (p+L+1);
          if(theta_k>0){
          ld[0] = e*e/(theta_k*theta_k*theta_k) - 1.0/theta_k;
          }
          ld[1] = b_eps*b_eps/gamma_kk - 1.0/(gamma_kk*2.0);
          for (int k=0; k<p; k++){
          ld[2+k] = sum(Xi(_,k) * (plinearl-plinear)/Plinear);
          }
          for (int l=0; l<L-1; l++){
          NumericVector Yil = ifelse(Yi==l, 1.0, 0.0);
          NumericVector Yill = ifelse(Yi==l+1, 1.0, 0.0);
          ld[2+p+l] = sum(Yil*plinear/Plinear - Yill*plinearl/Plinear);
          }
          for (int mi=0; mi<p+1+L; mi++){
            for (int mj=0; mj<p+1+L; mj++){
              num(mi,mj) = num(mi,mj) + ld[mi] * ld[mj] * wdens;
            }
          }
          num8 = num8 + wdens * (plinear/Plinear)*(plinearl/Plinear);
          den = den + wdens;
          }
      }
    }
  }
  if(ni == 2) {
    Girt(0,0) = 1.0;
    double g12 = KMi(1,0);
    Girt(1,0) = g12;
    double l22 = sqrt(1.0-g12*g12);
    Girt(1,1) = l22;
    for (int b1 = 0; b1 < ngrid; b1++) {
      double b_eps1 = b_list[b1];
      double w_eps1 = w_list[b1];
      for (int b2 = 0; b2 < ngrid; b2++) {
        double b_eps2 = b_list[b2];
        double w_eps2 = w_list[b2];
        NumericVector b_eps = NumericVector::create(b_eps1, b_eps2);
        NumericVector eps (ni);
        for (int j=0; j<ni; j++){
          eps[j] = sum(Girt(j,_)*b_eps);
        }
        eps = eps * sqrt(2.0*gamma_kk);
        for (int b = 0; b < ngrid; b++) {
          double b_e = b_list[b];
          double w_e = w_list[b];
          double e = b_e*1.4142136*theta_k;
          NumericVector linear = dli - Xb - e - eps;
          NumericVector linearl = dlil - Xb - e - eps;
          NumericVector Plinear (ni);
          NumericVector plinear = exp(-linear*linear*0.5)/2.506628;
          NumericVector plinearl = exp(-linearl*linearl*0.5)/2.506628;
          for (int j=0; j<ni; j++){
          if (Yi[j]<0){
          plinear[j]=0;
          plinearl[j]=0;
          }
          double pnorm1 = erf(linear[j]/1.4142136);
          double pnorm2 = erf(linearl[j]/1.4142136);
          Plinear[j] = (pnorm1-pnorm2)/2;
          }
          NumericVector dlinear = - linear * plinear;
          NumericVector dlinearl = - linearl * plinearl;
          for (int j=0; j<ni; j++){
          if(Plinear[j]==0) { Plinear[j] = 0.0000001; }
          }
          NumericVector loglik = log(Plinear);
          for (int j=0; j<ni; j++){
            if(Yi[j]<0) {loglik[j]=0;}
          }
          double prodlik = exp(sum(loglik));
          double wdens = w_e * w_eps1 * w_eps2 * prodlik;
          if(wdens>0){
          NumericVector ld (p+L+1);
          if(theta_k>0){
          ld[0] = e*e/(theta_k*theta_k*theta_k) - 1.0/theta_k;
          }
          ld[1] = sum(b_eps*b_eps)/gamma_kk - 2.0/(gamma_kk*2.0);
          for (int k=0; k<p; k++){
          ld[2+k] = sum(Xi(_,k) * (plinearl-plinear)/Plinear);
          }
          for (int l=0; l<L-1; l++){
          NumericVector Yil = ifelse(Yi==l, 1.0, 0.0);
          NumericVector Yill = ifelse(Yi==l+1, 1.0, 0.0);
          ld[2+p+l] = sum(Yil*plinear/Plinear - Yill*plinearl/Plinear);
          }
          for (int mi=0; mi<p+1+L; mi++){
            for (int mj=0; mj<p+1+L; mj++){
              num(mi,mj) = num(mi,mj) + ld[mi] * ld[mj] * wdens;
            }
          }
          num8 = num8 + wdens * (plinear/Plinear)*(plinearl/Plinear);
          den = den + wdens;
          }
        }
      }
    }
  }
  if(ni == 3) {
    Girt(0,0) = 1.0;
    double g12 = KMi(1,0);
    double g13 = KMi(2,0);
    Girt(1,0) = g12;
    double l22 = sqrt(1.0-g12*g12);
    Girt(1,1) = l22;
    double l23 = (KMi(2,1)-g12*g13)/l22;
    double l33 = sqrt(1.0-g13*g13-l23*l23);
    Girt(2,0) = g13;
    Girt(2,1) = l23;
    Girt(2,2) = l33;
    for (int b1 = 0; b1 < ngrid; b1++) {
      double b_eps1 = b_list[b1];
      double w_eps1 = w_list[b1];
      for (int b2 = 0; b2 < ngrid; b2++) {
        double b_eps2 = b_list[b2];
        double w_eps2 = w_list[b2];
        for (int b3 = 0; b3 < ngrid; b3++) {
        double b_eps3 = b_list[b3];
        double w_eps3 = w_list[b3];
        NumericVector b_eps = NumericVector::create(b_eps1, b_eps2, b_eps3);
        NumericVector eps (ni);
        for (int j=0; j<ni; j++){
          eps[j] = sum(Girt(j,_)*b_eps);
        }
        eps = eps * sqrt(2.0*gamma_kk);
        for (int b = 0; b < ngrid; b++) {
          double b_e = b_list[b];
          double w_e = w_list[b];
          double e = b_e*1.4142136*theta_k;
          NumericVector linear = dli - Xb - e - eps;
          NumericVector linearl = dlil - Xb - e - eps;
          NumericVector Plinear (ni);
          NumericVector plinear = exp(-linear*linear*0.5)/2.506628;
          NumericVector plinearl = exp(-linearl*linearl*0.5)/2.506628;
          for (int j=0; j<ni; j++){
          if (Yi[j]<0){
          plinear[j]=0;
          plinearl[j]=0;
          }
          double pnorm1 = erf(linear[j]/1.4142136);
          double pnorm2 = erf(linearl[j]/1.4142136);
          Plinear[j] = (pnorm1-pnorm2)/2;
          }
          NumericVector dlinear = - linear * plinear;
          NumericVector dlinearl = - linearl * plinearl;
          for (int j=0; j<ni; j++){
          if(Plinear[j]==0) { Plinear[j] = 0.0000001; }
          }
          NumericVector loglik = log(Plinear);
          for (int j=0; j<ni; j++){
            if(Yi[j]<0) {loglik[j]=0;}
          }
          double prodlik = exp(sum(loglik));
          double wdens = w_e * w_eps1 * w_eps2 * w_eps3 * prodlik;
          if(wdens>0){
          NumericVector ld (p+L+1);
          if(theta_k>0){
          ld[0] = e*e/(theta_k*theta_k*theta_k) - 1.0/theta_k;
          }
          ld[1] = sum(b_eps*b_eps)/gamma_kk - 3.0/(gamma_kk*2.0);
          for (int k=0; k<p; k++){
          ld[2+k] = sum(Xi(_,k) * (plinearl-plinear)/Plinear);
          }
          for (int l=0; l<L-1; l++){
          NumericVector Yil = ifelse(Yi==l, 1.0, 0.0);
          NumericVector Yill = ifelse(Yi==l+1, 1.0, 0.0);
          ld[2+p+l] = sum(Yil*plinear/Plinear - Yill*plinearl/Plinear);
          }
          for (int mi=0; mi<p+1+L; mi++){
            for (int mj=0; mj<p+1+L; mj++){
              num(mi,mj) = num(mi,mj) + ld[mi] * ld[mj] * wdens;
            }
          }
          num8 = num8 + wdens * (plinear/Plinear)*(plinearl/Plinear);
          den = den + wdens;
          }
          }
        }
      }
    }
  }
  if(ni == 4) {
    Girt(0,0) = 1.0;
    double g12 = KMi(1,0);
    double g13 = KMi(2,0);
    Girt(1,0) = g12;
    double l22 = sqrt(1.0-g12*g12);
    Girt(1,1) = l22;
    double l23 = (KMi(2,1)-g12*g13)/l22;
    double l33 = sqrt(1.0-g13*g13-l23*l23);
    Girt(2,0) = g13;
    Girt(2,1) = l23;
    Girt(2,2) = l33;
    double g14 = KMi(3,0);
    double g24 = KMi(3,1);
    double l24 = (g24-g12*g14)/l22;
    double l34 = (KMi(3,2)-g13*g14-l23*l24)/l33;
    Girt(3,0) = g14;
    Girt(3,1) = l24;
    Girt(3,2) = l34;
    Girt(3,3) = sqrt(1.0-g14*g14-l24*l24-l34*l34);
    for (int b1 = 0; b1 < ngrid; b1++) {
      double b_eps1 = b_list[b1];
      double w_eps1 = w_list[b1];
      for (int b2 = 0; b2 < ngrid; b2++) {
        double b_eps2 = b_list[b2];
        double w_eps2 = w_list[b2];
        for (int b3 = 0; b3 < ngrid; b3++) {
        double b_eps3 = b_list[b3];
        double w_eps3 = w_list[b3];
        for (int b4 = 0; b4 < ngrid; b4++) {
        double b_eps4 = b_list[b4];
        double w_eps4 = w_list[b4];
        NumericVector b_eps = NumericVector::create(b_eps1, b_eps2, b_eps3, b_eps4);
        NumericVector eps (ni);
        for (int j=0; j<ni; j++){
          eps[j] = sum(Girt(j,_)*b_eps);
        }
        eps = eps * sqrt(2.0*gamma_kk);
        for (int b = 0; b < ngrid; b++) {
          double b_e = b_list[b];
          double w_e = w_list[b];
          double e = b_e*1.4142136*theta_k;
          NumericVector linear = dli - Xb - e - eps;
          NumericVector linearl = dlil - Xb - e - eps;
          NumericVector Plinear (ni);
          NumericVector plinear = exp(-linear*linear*0.5)/2.506628;
          NumericVector plinearl = exp(-linearl*linearl*0.5)/2.506628;
          for (int j=0; j<ni; j++){
          if (Yi[j]<0){
          plinear[j]=0;
          plinearl[j]=0;
          }
          double pnorm1 = erf(linear[j]/1.4142136);
          double pnorm2 = erf(linearl[j]/1.4142136);
          Plinear[j] = (pnorm1-pnorm2)/2;
          }
          NumericVector dlinear = - linear * plinear;
          NumericVector dlinearl = - linearl * plinearl;
          for (int j=0; j<ni; j++){
          if(Plinear[j]==0) { Plinear[j] = 0.0000001; }
          }
          NumericVector loglik = log(Plinear);
          for (int j=0; j<ni; j++){
            if(Yi[j]<0) {loglik[j]=0;}
          }
          double prodlik = exp(sum(loglik));
          double wdens = w_e * w_eps1 * w_eps2 * w_eps3 * w_eps4 * prodlik;
          if(wdens>0){
          NumericVector ld (p+L+1);
          if(theta_k>0){
          ld[0] = e*e/(theta_k*theta_k*theta_k) - 1.0/theta_k;
          }
          ld[1] = sum(b_eps*b_eps)/gamma_kk - 4.0/(gamma_kk*2.0);
          for (int k=0; k<p; k++){
          ld[2+k] = sum(Xi(_,k) * (plinearl-plinear)/Plinear);
          }
          for (int l=0; l<L-1; l++){
          NumericVector Yil = ifelse(Yi==l, 1.0, 0.0);
          NumericVector Yill = ifelse(Yi==l+1, 1.0, 0.0);
          ld[2+p+l] = sum(Yil*plinear/Plinear - Yill*plinearl/Plinear);
          }
          for (int mi=0; mi<p+1+L; mi++){
            for (int mj=0; mj<p+1+L; mj++){
              num(mi,mj) = num(mi,mj) + ld[mi] * ld[mj] * wdens;
            }
          }
          num8 = num8 + wdens * (plinear/Plinear)*(plinearl/Plinear);
          den = den + wdens;
          }
          }
        }
      }
    }
    }
  }
  if(den>0){
  for (int mi=0; mi<1+p+L; mi++){
  for (int mj=0; mj<1+p+L; mj++){
    Ell(mi,mj) = Ell(mi,mj) + num(mi,mj)/den;
  }
  }
  for (int j=0; j<ni; j++){
  Ddd[id_begin+j] = num8[j]/den;
  }
  }
  id_begin = id_begin + ni;
  }
  NumericVector ldd_ad1 = d3 + Ddd;
  NumericVector ldd_ad2 = d4 - Ddd;
  NumericMatrix Icomneg (1+p+L);
  Icomneg(0,0) = ldd_t;
  Icomneg(1,1) = ldd_g;
  for (int mi=0; mi<p; mi++){
  for (int mj=0; mj<p; mj++){
    Icomneg(2+mi,2+mj) = dl2dalpha(mi,mj);
  }
  }
  for (int l=0; l<L-1; l++){
    NumericVector Yl = ifelse(Y==l,1.0,0.0);
    for (int r=0; r<L-1; r++){
      if(r==l) Icomneg(2+p+l,2+p+r) = dl2dd[l];
      if(r==l+1) 
      {
      NumericVector Yr = ifelse(Y==r,1.0,0.0);
      Icomneg(2+p+l,2+p+r) = sum(Yr*Ddd);
      }
      if(r==l-1)
      {
      Icomneg(2+p+l,2+p+r) = sum(Yl*Ddd);
      }
    }
  }
  for (int mi=0; mi<p; mi++){
  for (int mj=0; mj<L-1; mj++){
    NumericVector Yl = ifelse(Y==mj, 1.0, 0.0);
    NumericVector Yll = ifelse(Y==mj+1, 1.0, 0.0);
    Icomneg(2+mi,2+p+mj) = sum(X(_,mi)*(ldd_ad2*Yll-ldd_ad1*Yl));
    Icomneg(2+p+mj,2+mi) = sum(X(_,mi)*(ldd_ad2*Yll-ldd_ad1*Yl));
  }
  }
  NumericMatrix V (p+1+L,p+1+L);
  for (int mi=0; mi<p+1+L; mi++){
    for (int mj=0; mj<p+1+L; mj++){
    V(mi,mj) = - Icomneg(mi,mj) + El2(mi,mj) - Ell(mi,mj);
  }
  }
  V = ginv(V);
  NumericMatrix Vs = mprod(V,El2);
  Vs = mprod(Vs,V);
  NumericVector se (p+1+L);
  for (int k=0; k<p+1+L; k++){
  double sek = V(k,k);
  if(sek > 0)
  {
  se[k] = sqrt(sek);
  } 
  else
  {
  double seks = Vs(k,k);
  se[k] = sqrt(seks);
  }
  }
  NumericMatrix phin = mprod(Eld, V);
  phin = n * phin;
  NumericMatrix phi (n,p+2+L);
  NumericVector meanphi (p+1+L);
  for (int k=0; k<p+1+L; k++){
  meanphi[k] = sum(phin(_,k)) / n;
  }
  for (int i=0; i<n; i++){
  for (int k=0; k<p+1+L; k++){
  phi(i,0) = Fam_list[i];
  phi(i,k+1) = phin(i,k) - meanphi[k];
  }
  }
  double sigmau2 = 1.0;
  double S2 = theta_k*theta_k + gamma_kk + sigmau2;
  double rho = gamma_kk/S2;
  double zeta = theta_k*theta_k/S2;
  double D1 = -2.0*theta_k*gamma_kk/S2/S2;
  double D2 = (theta_k*theta_k+sigmau2)/S2/S2;
  NumericVector IF = phi(_,1)*D1 + phi(_,2)*D2;
  double serho = sqrt(sum(IF*IF)/(n-1)/n);
  D1 = 2.0*theta_k*(gamma_kk+sigmau2)/S2/S2;
  D2 = -theta_k*theta_k/S2/S2;
  IF = phi(_,1)*D1 + phi(_,2)*D2;
  double sezeta = sqrt(sum(IF*IF)/(n-1)/n);
  S2 = theta_k*theta_k + gamma_kk;
  double fracgen = gamma_kk/S2;
  double fracenv = theta_k*theta_k/S2;
  D1 = 2.0*theta_k*gamma_kk/S2/S2;
  D2 = -theta_k*theta_k/S2/S2;
  IF = phi(_,1)*D1 + phi(_,2)*D2;
  double sefrac = sqrt(sum(IF*IF)/(n-1)/n);
  
  // Output
  
  List results;
  results["alpha"] = alpha_k;
  results["theta"] = theta_k;
  results["gamma"] = gamma_kk;
  results["delta"] = delta_k;
  results["se.theta"] = se[0];
  results["se.gamma"] = se[1];
  NumericVector se_alpha (p);
  NumericVector se_delta (L-1);
  for (int k=0; k<p; k++){
  se_alpha[k] = se[2+k];
  }
  for (int k=0; k<L-1; k++){
  se_delta[k] = se[1+p+k];
  }
  results["se.alpha"] = se_alpha;
  results["se.delta"] = se_delta;
  results["est"] = res;
  results["se"] = se;
  results["phi"] = phi;
  results["V"] = V;
  results["L"] = L;
  results["p"] = p;
  results["rho"] = rho;
  results["zeta"] = zeta;
  results["fracgen"] = fracgen;
  results["fracenv"] = fracenv;
  results["se.rho"] = serho;
  results["se.zeta"] = sezeta;
  results["se.fracgen"] = sefrac;
  results["se.fracenv"] = sefrac;
  results["iterrun"] = iterrun;
  results["fail"] = fail;
  return results;
  }
')
