## Survival

cppFunction('
  List cal_sur(NumericVector T, NumericVector D, NumericMatrix X,
              IntegerVector familyID, NumericMatrix KM,
              NumericVector par0){
  Environment glob = Environment::global_env();
  Function ginv("ginv");
  NumericVector w_list = glob.get("w_list");
  NumericVector b_list = glob.get("b_list");
  int ngrid = w_list.size();
  NumericVector X1 = X(0,_);
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
  NumericVector tseq0;
  for (int i=0; i<N; i++){
    if(D[i]>0.5) 
    {
    tseq0.push_back(T[i]);
    }
  }
  std::sort(tseq0.begin(), tseq0.end());
  int L = tseq0.size();
  NumericVector Ts = NumericVector::create(tseq0[0]);
  for (int i=1; i<L; i++){
  if(tseq0[i]>tseq0[i-1]) {Ts.push_back(tseq0[i]);}
  }
  L = Ts.size();
  NumericVector res (p+2+L);
  NumericVector lam (L, 0.001);
  NumericVector Lam (N);
  for (int i=0; i<N; i++){
  double ti = T[i];
  NumericVector TsL = ifelse(Ts<=ti,1.0,0.0);
  Lam[i] = sum(TsL*lam);
  }
  NumericVector Eel (N, 1.0);
  NumericVector Ee2 (n);
  NumericVector EepsGeps (n);
  
  double theta_k = par0[0];
  double gamma_kk = par0[1];
  NumericVector alpha_k (p);
  for (int k=0; k<p; k++){
  alpha_k[k] = par0[2+k];
  }
  
  int iterrun = 0;
  double tol = 1;
  while (tol > 0.0001){
  NumericVector res0;
  NumericVector alpha_k0;
  res0.assign(res.begin(),res.end());
  alpha_k0.assign(alpha_k.begin(),alpha_k.end());
  int id_begin = 0;
  for (int i=0; i<n; i++){
  int ni = Fam_size[i];
  if(ni==0) {continue;}
  NumericVector Di (ni);
  NumericVector Lami (ni);
  NumericVector Xb (ni);
  for (int j=0; j<ni; j++){
    int nij = id_begin+j;
    Di[j] = D[nij];
    Xb[j] = sum(X(nij,_)*alpha_k);
    Lami[j] = Lam[nij];
  }
  NumericVector num1 (ni);
  double num2 = 0;
  double num3 = 0;
  double den = 0;
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
          double e = b_e*sqrt(2)*theta_k;
          NumericVector linear = Xb + e + eps;
          NumericVector El = exp(linear);
          NumericVector loglik = Di*linear - Lami*El;
          double wdens = w_e * w_eps1 * exp(sum(loglik));
          if(std::isnan(wdens)) {wdens=0;}
          den = den + wdens;
          num1 = num1 + wdens * El;
          num2 = num2 + wdens * e*e;
          num3 = num3 + wdens * eps*eps;
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
          double e = b_e*sqrt(2)*theta_k;
          NumericVector linear = Xb + e + eps;
          NumericVector El = exp(linear);
          NumericVector loglik = Di*linear - Lami*El;
          double wdens = w_e * w_eps1 * w_eps2 * exp(sum(loglik));
          if(std::isnan(wdens)) {wdens=0;}
          den = den + wdens;
          num1 = num1 + wdens * El;
          num2 = num2 + wdens * e*e;
          num3 = num3 + wdens * 2.0*gamma_kk*sum(b_eps*b_eps);
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
          double e = b_e*sqrt(2)*theta_k;
          NumericVector linear = Xb + e + eps;
          NumericVector El = exp(linear);
          NumericVector loglik = Di*linear - Lami*El;
          double wdens = w_e * w_eps1 * w_eps2 * w_eps3 * exp(sum(loglik));
          if(std::isnan(wdens)) {wdens=0;}
          den = den + wdens;
          num1 = num1 + wdens * El;
          num2 = num2 + wdens * e*e;
          num3 = num3 + wdens * 2.0*gamma_kk*sum(b_eps*b_eps);
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
          double e = b_e*sqrt(2)*theta_k;
          NumericVector linear = Xb + e + eps;
          NumericVector El = exp(linear);
          NumericVector loglik = Di*linear - Lami*El;
          double wdens = w_e * w_eps1 * w_eps2 * w_eps3 * w_eps4 * exp(sum(loglik));
          if(std::isnan(wdens)) {wdens=0;}
          den = den + wdens;
          num1 = num1 + wdens * El;
          num2 = num2 + wdens * e*e;
          num3 = num3 + wdens * 2.0*gamma_kk*sum(b_eps*b_eps);
          }
          }
        }
      }
    }
  }
  if(den>0){
  for (int j=0; j<ni; j++){
    double eeli = num1[j]/den;
    if(eeli>1000) {eeli = 0;}
    Eel[id_begin+j] = eeli;
  }
  Ee2[i] = num2/den;
  EepsGeps[i] = num3/den;
  }
  id_begin = id_begin + ni;
  }
  
  NumericMatrix E1 (N,p);
  NumericMatrix E2 (p,p);
  for (int i=0; i<N; i++){
  if(D[i]>0.5){
    double ti = T[i];
    NumericVector TL = ifelse(T>=ti,1.0,0.0);
    NumericVector TLE = TL*Eel;
    double SIEel = sum(TLE);
    if(SIEel>0){
    for (int k=0; k<p; k++){
    E1(i,k) = sum(TLE*X(_,k))/SIEel;
    for (int r=0; r<p; r++){
      E2(k,r) = E2(k,r) + sum(TLE*X(_,k)*X(_,r))/SIEel;
    }
    }
    }
  }
  }
  NumericVector dldalpha (p);
  NumericMatrix dl2dalpha (p,p);
  for (int k=0; k<p; k++){
    dldalpha[k] = sum(X(_,k)*D-E1(_,k));
    for (int r=0; r<p; r++){
    dl2dalpha(k,r) = sum(E1(_,k)*E1(_,r)) - E2(k,r);
    }
  }
  NumericMatrix d2_1 (p,p);
  double sdd = sum(dl2dalpha*dl2dalpha);
  if(sdd>0){
  d2_1 = ginv(dl2dalpha);
  }
  NumericVector d2d1 (p);
  for (int k=0; k<p; k++){
    d2d1[k] = sum(d2_1(k,_)*dldalpha);
  }
  for (int r=0; r<p; r++){
    if(std::isnan(d2d1[r])) {
    alpha_k[r] = 0.0;
    d2d1[r] = 0.0;
    }
    if(d2d1[r]>100) {
    alpha_k[r] = 0.0;
    d2d1[r] = 0.0;
    }
  }
  alpha_k = alpha_k - 0.8 * d2d1;
  theta_k = sqrt(sum(Ee2) / n);
  gamma_kk = sum(EepsGeps) / N;
  
  for (int i=0; i<N; i++){
    Eel[i] = Eel[i] * exp(-sum(X(i,_)*0.8*d2d1));
    if(std::isnan(Eel[i])) {Eel[i] = 0.0;}
    if(std::isinf(Eel[i])) {Eel[i] = 0.0;}
  }
  for (int l=0; l<L; l++){
    double Tsl = Ts[l];
    NumericVector TslD = ifelse(T==Tsl,1.0,0.0);
    NumericVector TslL = ifelse(T>=Tsl,1.0,0.0);
    double DTEL = sum(TslL*Eel);
    if(DTEL>0) {lam[l] = sum(D*TslD)/DTEL;}
  }
  for (int i=0; i<N; i++){
    double ti = T[i];
    Lam[i] = sum(lam*ifelse(Ts<=ti,1.0,0.0));
  }
  
  res[0] = theta_k;
  res[1] = gamma_kk;
  for (int k=0; k<p; k++){
    res[k+2] = alpha_k[k];
  }
  for (int l=0; l<L-1; l++){
    res[l+p+2] = lam[l];
  }
  tol = max(sqrt((res - res0)*(res - res0)));
  
  // std::cout << mean(Eel) << " " << mean(alpha_k) << std::endl;
  std::cout << res[0] << " " << res[1] << " " << mean(alpha_k) << " " << tol << std::endl;
  iterrun = iterrun + 1;
  if(iterrun>=5000) {tol = 0.0;}
  }
  int fail = 0;
  if(iterrun>=5000) {fail=1;}
  if(std::isnan(tol)) {fail=1;}
  
  // Calculate Variance
  
  NumericMatrix lam_dh (L, p+2);
  NumericVector h (p+2);
  for (int k=0; k<p+2; k++){
    h[k] = 2/sqrt(n)*res[k];
  }
  NumericVector lamT (N);
  NumericMatrix pfloglik (n,p+2);
  NumericMatrix pfloglik_a (n,p+2);
  for (int i=0; i<N; i++){
    double ti = T[i];
    NumericVector Tsi = ifelse(Ts==ti, 1.0, 0.0);
    double dl = sum(Tsi*lam);
    if (dl>0) {lamT[i] = log(dl);}
  }
  
  for (int q=0; q<2; q++){
  for (int k=0; k<p+2; k++){
  NumericVector lamh;
  NumericVector Lamh;
  NumericVector resh;
  lamh.assign(lam.begin(), lam.end());
  Lamh.assign(Lam.begin(), Lam.end());
  resh.assign(res.begin(), res.end());
  resh[k] = resh[k] + q*h[k];
  theta_k = resh[0];
  gamma_kk = resh[1];
  for (int k=0; k<p; k++){
  alpha_k[k] = resh[k+2];
  }
  
  NumericVector lk (n);
  tol = 1.0;
  while (tol>0.0001) {
  int id_begin = 0;
  NumericVector lam_0;
  lam_0.assign(lamh.begin(), lamh.end());
  for (int i=0; i<n; i++){
  int ni = Fam_size[i];
  if(ni==0) {continue;}
  NumericVector Ti (ni);
  NumericVector Di (ni);
  NumericVector Lami (ni);
  NumericVector lamTi (ni);
  NumericVector Xb (ni);
  for (int j=0; j<ni; j++){
    int nij = id_begin+j;
    Ti[j] = T[nij];
    Di[j] = D[nij];
    Xb[j] = sum(X(nij,_)*alpha_k);
    Lami[j] = Lamh[nij];
    lamTi[j] = lamT[nij];
  }
  NumericVector num1 (ni);
  double den = 0;
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
          double e = b_e*sqrt(2)*theta_k;
          NumericVector linear = Xb + e + eps;
          NumericVector El = exp(linear);
          NumericVector loglik = Di*(linear+lamTi) - Lami*El;
          double wdens = w_e * w_eps1 * exp(sum(loglik));
          if(std::isnan(wdens)) {wdens=0;}
          den = den + wdens;
          num1 = num1 + wdens * El;
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
          double e = b_e*sqrt(2)*theta_k;
          NumericVector linear = Xb + e + eps;
          NumericVector El = exp(linear);
          NumericVector loglik = Di*(linear+lamTi) - Lami*El;
          double wdens = w_e * w_eps1 * w_eps2 * exp(sum(loglik));
          if(std::isnan(wdens)) {wdens=0;}
          den = den + wdens;
          num1 = num1 + wdens * El;
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
          double e = b_e*sqrt(2)*theta_k;
          NumericVector linear = Xb + e + eps;
          NumericVector El = exp(linear);
          NumericVector loglik = Di*(linear+lamTi) - Lami*El;
          double wdens = w_e * w_eps1 * w_eps2 * w_eps3 * exp(sum(loglik));
          if(std::isnan(wdens)) {wdens=0;}
          den = den + wdens;
          num1 = num1 + wdens * El;
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
          double e = b_e*sqrt(2)*theta_k;
          NumericVector linear = Xb + e + eps;
          NumericVector El = exp(linear);
          NumericVector loglik = Di*(linear+lamTi) - Lami*El;
          double wdens = w_e * w_eps1 * w_eps2 * w_eps3 * w_eps4 * exp(sum(loglik));
          if(std::isnan(wdens)) {wdens=0;}
          den = den + wdens;
          num1 = num1 + wdens * El;
          }
          }
        }
      }
    }
  }
  if(den>0){
  for (int j=0; j<ni; j++){
    double eeli = num1[j]/den;
    if(eeli>1000) {eeli = 0;}
    Eel[id_begin+j] = eeli;
  }
  }
  lk[i] = den;
  id_begin = id_begin + ni;
  }
  for (int l=0; l<L; l++){
    double Tsl = Ts[l];
    NumericVector TslD = ifelse(T==Tsl,1.0,0.0);
    NumericVector TslL = ifelse(T>=Tsl,1.0,0.0);
    lamh[l] = sum(D*TslD)/sum(TslL*Eel);
  }
  for (int i=0; i<N; i++){
    double ti = T[i];
    NumericVector TsL = ifelse(Ts<=ti,1.0,0.0);
    NumericVector TsD = ifelse(Ts==ti,1.0,0.0);
    Lamh[i] = sum(lamh*TsL);
    double dl = sum(TsD*lamh);
    if (dl>0) {lamT[i] = log(dl);} 
  }
  tol = max(sqrt((lamh - lam_0)*(lamh - lam_0)));
  }
  std::cout << q << " " << k << std::endl;
  if(q==0) {
  for (int i=0; i<n; i++){
    pfloglik(i,k) = log(lk[i]);
  }
  }
  if(q==1) {
  for (int i=0; i<n; i++){
    pfloglik_a(i,k) = log(lk[i]);
  }
  for (int l=0; l<L; l++){
    lam_dh(l,k) = lamh[l];
  }
  }
  resh[k] = res[k];
  }
  }
  
  NumericMatrix score (n,p+2);
  NumericMatrix Infom (p+2,p+2);
  for (int i=0; i<n; i++){
  for (int k=0; k<p+2; k++){
    if(h[k]>0){
    score(i,k) = (pfloglik_a(i,k) - pfloglik(i,k))/h[k];
    }
  }
  }
  for (int r=0; r<p+2; r++){
  for (int k=0; k<p+2; k++){
    Infom(r,k) = sum(score(_,r)*score(_,k));
  }
  }
  NumericMatrix Iinv = ginv(Infom);
  NumericMatrix phin (n,p+2);
  for (int i=0; i<n; i++){
  for (int k=0; k<p+2; k++){
    phin(i,k) = sum(score(i,_)*Iinv(k,_))*n;
  }
  }
  NumericMatrix phi (n,p+3);
  NumericVector meanphi (p+2);
  for (int k=0; k<p+2; k++){
  meanphi[k] = sum(phin(_,k)) / n;
  }
  for (int i=0; i<n; i++){
  for (int k=0; k<p+2; k++){
  phi(i,0) = Fam_list[i];
  phi(i,k+1) = phin(i,k) - meanphi[k];
  }
  }
  NumericVector se (p+2);
  for (int k=0; k<p+2; k++){
  double sek = Iinv(k,k);
  if(sek>0) {
    se[k] = sqrt(sek);
  } else {
    se[k] = sqrt(sum(phin(_,k)*phin(_,k))/n);
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
  NumericVector se_alpha (p);
  theta_k = res[0];
  gamma_kk = res[1];
  for (int k=0; k<p; k++){
  alpha_k[k] = res[2+k];
  se_alpha[k] = se[2+k];
  }
  NumericVector est (p+2);
  for (int k=0; k<p+2; k++){
    est[k] = res[k];
  }
  results["alpha"] = alpha_k;
  results["theta"] = theta_k;
  results["gamma"] = gamma_kk;
  results["se.alpha"] = se_alpha;
  results["se.theta"] = se[0];
  results["se.gamma"] = se[1];
  results["est"] = est;
  results["se"] = se;
  results["tt"] = Ts;
  results["lam"] = lam;
  results["h"] = h;
  results["lam_dh"] = lam_dh;
  results["phi"] = phi;
  results["V"] = Iinv;
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