# Continuous

cppFunction('
  List cal_con(NumericVector Y, NumericMatrix X,
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
  NumericVector Famsize = as<NumericVector>(Fam_size);
  n = Fam_list.size();
  NumericVector YXb (N);
  NumericVector res (p+3);
  NumericVector Ee (N);
  NumericVector Eeps (N);
  NumericVector Eeps2 (N);
  NumericVector Eeeps (N);
  NumericVector Ee2 (n);
  NumericVector EepsGeps (n);
  
  double theta_k = par0[0];
  double gamma_kk = par0[1];
  NumericVector alpha_k (p);
  for (int k=0; k<p; k++){
  alpha_k[k] = par0[2+k];
  }
  double sigmau2 = par0[p+2];
  NumericMatrix XX (p,p);
  for (int r=0; r<p; r++){
  for (int s=0; s<p; s++){
    XX(r,s) = sum(X(_,r)*X(_,s));
  }
  }
  NumericMatrix XXinv = ginv(XX);
  
  int iterrun = 0;
  double tol = 1;
  while (tol > 0.0001){
  NumericVector res0;
  res0.assign(res.begin(),res.end());
  int id_begin = 0;
  for (int i=0; i<n; i++){
  int ni = Fam_size[i];
  if(ni==0) {continue;}
  NumericVector Yi (ni);
  NumericVector Xb (ni);
  for (int j=0; j<ni; j++){
    Xb[j] = sum(X(id_begin+j,_)*alpha_k);
    double Yij = Y[id_begin+j];
    Yi[j] = Yij;
  }
  NumericMatrix KMi (ni,ni);
  for (int r=0; r<ni; r++){
  for (int s=0; s<ni; s++){
    KMi(r,s) = KM(id_begin+r,s);
  }
  }
  
  NumericMatrix S (ni,ni);
  for (int k=0; k<ni; k++){
  S(k,k) = sigmau2;
  for (int j=0; j<ni; j++){
  S(k,j) = S(k,j) + theta_k*theta_k + gamma_kk*KMi(k,j);
  }
  }
  NumericMatrix Giinv (ni,ni);
  NumericMatrix Si (ni,ni);
  if (ni==1){
  Giinv = 1.0/KMi;
  Si = 1.0/S;
  } else {
  Giinv = ginv(KMi);
  Si = ginv(S);
  }
  NumericVector YXbi = Yi - Xb;
  NumericVector SE (ni);
  NumericVector sumSi (ni);
  NumericMatrix SG (ni,ni);
  for (int k=0; k<ni; k++){
    SE[k] = sum(Si(k,_)*YXbi);
    sumSi[k] = sum(Si(k,_));
    for (int j=0; j<ni; j++){
    SG(k,j) = sum(Si(k,_)*KMi(_,j));
    }
  }
  double Eei = theta_k*theta_k*sum(SE);
  NumericVector Eepsi (ni);
  for (int k=0; k<ni; k++){
    Eepsi[k] = gamma_kk*sum(KMi(k,_)*SE);
  }
  double Vei = theta_k*theta_k - theta_k*theta_k*theta_k*theta_k*sum(sumSi);
  NumericMatrix Vepsi (ni,ni);
  for (int k=0; k<ni; k++){
    for (int j=0; j<ni; j++){
    Vepsi(k,j) = gamma_kk*KMi(k,j) - gamma_kk*gamma_kk*sum(KMi(k,_)*SG(_,j));
    }
  }
  Ee2[i] = Eei*Eei + Vei;
  double DV = 0.0;
  for (int k=0; k<ni; k++){
    DV = DV + Vepsi(k,k);
  }
  Eeps2[i] = sum(Eepsi*Eepsi) + DV;
  NumericVector GE (ni);
  NumericVector GV (ni);
  for (int k=0; k<ni; k++){
  GE[k] = sum(Giinv(k,_)*Eepsi);
  GV[k] = sum(Giinv(k,_)*Vepsi(_,k));
  }
  EepsGeps[i] = sum(Eepsi*GE) + sum(GV);
  for (int j=0; j<ni; j++){
  Ee[id_begin+j] = Eei;
  Eeps[id_begin+j] = Eepsi[j];
  Eeeps[id_begin+j] = Eei*Eepsi[j]-theta_k*theta_k*gamma_kk*sum(KMi(j,_)*sumSi);
  YXb[id_begin+j] = YXbi[j];
  }
  id_begin = id_begin + ni;
  }
  
  double sig1 = sum(YXb*YXb) - 2.0*sum(YXb*(Ee+Eeps));
  double sig2 = sum(Eeps2) + 2.0*sum(Eeeps) + sum(Ee2*Famsize);
  sigmau2 = (sig1 + sig2) / N;
  NumericVector Yee = Y - Ee - Eeps;
  NumericVector XY (p);
  for (int k=0; k<p; k++){
  XY[k] = sum(X(_,k)*Yee);
  }
  for (int k=0; k<p; k++){
  alpha_k[k] = sum(XXinv(k,_)*XY);
  }
  theta_k = sqrt(sum(Ee2)/n);
  gamma_kk = sum(EepsGeps)/N;
  res[0] = theta_k;
  res[1] = gamma_kk;
  for (int k=0; k<p; k++){
    res[k+2] = alpha_k[k];
  }
  res[p+2] = sigmau2;
  tol = max(sqrt((res - res0)*(res - res0)));
  iterrun = iterrun + 1;
  if(iterrun>=5000) {tol = 0.0;}
  std::cout << res[0] << " " << res[1] << " " << mean(alpha_k) << " " << res[p+2] << " " << tol << std::endl;
  }
  int fail = 0;
  if(iterrun>=5000) {fail=1;}
  if(std::isnan(tol)) {fail=1;}
  
  // Calculate Variance
  
  NumericVector ld_t (n);
  if(theta_k>0){
  ld_t = Ee2/(theta_k*theta_k*theta_k) - 1.0/theta_k;
  }
  NumericVector ld_g = EepsGeps/(gamma_kk*gamma_kk*2.0) - Famsize/(gamma_kk*2.0);
  NumericVector Y_ee = YXb - Ee - Eeps;
  NumericMatrix Eld (n,p+3);
  int id_begin = 0;
  for (int i=0; i<n; i++){
    int ni = Fam_size[i];
    if(ni==0) {continue;}
    Eld(i,0) = ld_t[i];
    Eld(i,1) = ld_g[i];
    for (int j=0; j<ni; j++){
    int ij = id_begin + j;
    for (int k=0; k<p; k++){
      Eld(i,k+2) = Eld(i,k+2) + Y_ee[ij]*X(ij,k)/sigmau2;
    }
    double elds1 = (Ee2[i]+Eeps2[ij])/2.0/sigmau2/sigmau2 - 1.0/sigmau2/2.0;
    double elds2 = YXb[ij]*YXb[ij]/2.0 - YXb[ij]*(Ee[ij]+Eeps[ij]) + Eeeps[ij];
    Eld(i,p+2) = Eld(i,p+2) + elds1 + elds2/sigmau2/sigmau2; 
    }
    id_begin = id_begin + ni;
  }
  
  NumericMatrix El2 = mprod(transpose(Eld), Eld);
  NumericMatrix Ell (p+3,p+3);
  id_begin = 0;
  for (int i=0; i<n; i++){
  int ni = Fam_size[i];
  if(ni==0) {continue;}
  NumericMatrix Xi (ni,p);
  NumericVector YXbi (ni);
  for (int j=0; j<ni; j++){
  for (int k=0; k<p; k++){
  Xi(j,k) = X(id_begin+j,k);
  }
  YXbi[j] = Y[id_begin+j] - sum(Xi(j,_)*alpha_k);
  }
  NumericMatrix num (p+3,p+3);
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
      double b_eps = b_list[b1];
      double w_eps1 = w_list[b1];
      double eps = b_eps * sqrt(2.0*gamma_kk);
        for (int b = 0; b < ngrid; b++) {
          double b_e = b_list[b];
          double w_e = w_list[b];
          double e = b_e*1.4142136*theta_k;
          NumericVector linear = YXbi - e - eps;
          NumericVector loglik = -linear*linear/2.0/sigmau2;
          double prodlik = exp(sum(loglik));
          double wdens = w_e * w_eps1 * prodlik;
          if(wdens>0){
          NumericVector ld (p+3);
          if(theta_k>0){
          ld[0] = e*e/(theta_k*theta_k*theta_k) - 1.0/theta_k;
          }
          ld[1] = b_eps*b_eps/gamma_kk - 1.0/(gamma_kk*2.0);
          for (int k=0; k<p; k++){
          ld[2+k] = sum(Xi(_,k)*linear/sigmau2);
          }
          ld[2+p] = -1.0/sigmau2/2.0+sum(linear*linear)/2.0/sigmau2/sigmau2;
          for (int mi=0; mi<p+3; mi++){
            for (int mj=0; mj<p+3; mj++){
              num(mi,mj) = num(mi,mj) + ld[mi] * ld[mj] * wdens;
            }
          }
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
          NumericVector linear = YXbi - e - eps;
          NumericVector loglik = -linear*linear/2.0/sigmau2;
          double prodlik = exp(sum(loglik));
          double wdens = w_e * w_eps1 * w_eps2 * prodlik;
          if(wdens>0){
          NumericVector ld (p+3);
          if(theta_k>0){
          ld[0] = e*e/(theta_k*theta_k*theta_k) - 1.0/theta_k;
          }
          ld[1] = sum(b_eps*b_eps)/gamma_kk - 1.0/gamma_kk;
          for (int k=0; k<p; k++){
          ld[2+k] = sum(Xi(_,k)*linear/sigmau2);
          }
          ld[2+p] = -1.0/sigmau2+sum(linear*linear)/2.0/sigmau2/sigmau2;
          for (int mi=0; mi<p+3; mi++){
            for (int mj=0; mj<p+3; mj++){
              num(mi,mj) = num(mi,mj) + ld[mi] * ld[mj] * wdens;
            }
          }
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
          NumericVector linear = YXbi - e - eps;
          NumericVector loglik = -linear*linear/2.0/sigmau2;
          double prodlik = exp(sum(loglik));
          double wdens = w_e * w_eps1 * w_eps2 * w_eps3 * prodlik;
          if(wdens>0){
          NumericVector ld (p+3);
          if(theta_k>0){
          ld[0] = e*e/(theta_k*theta_k*theta_k) - 1.0/theta_k;
          }
          ld[1] = sum(b_eps*b_eps)/gamma_kk - 3.0/(gamma_kk*2.0);
          for (int k=0; k<p; k++){
          ld[2+k] = sum(Xi(_,k)*linear/sigmau2);
          }
          ld[2+p] = -3.0/sigmau2/2.0+sum(linear*linear)/2.0/sigmau2/sigmau2;
          for (int mi=0; mi<p+3; mi++){
            for (int mj=0; mj<p+3; mj++){
              num(mi,mj) = num(mi,mj) + ld[mi] * ld[mj] * wdens;
            }
          }
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
          NumericVector linear = YXbi - e - eps;
          NumericVector loglik = -linear*linear/2.0/sigmau2;
          double prodlik = exp(sum(loglik));
          double wdens = w_e * w_eps1 * w_eps2 * w_eps3 * w_eps4 * prodlik;
          if(wdens>0){
          NumericVector ld (p+3);
          if(theta_k>0){
          ld[0] = e*e/(theta_k*theta_k*theta_k) - 1.0/theta_k;
          }
          ld[1] = sum(b_eps*b_eps)/gamma_kk - 2.0/gamma_kk;
          for (int k=0; k<p; k++){
          ld[2+k] = sum(Xi(_,k)*linear/sigmau2);
          }
          ld[2+p] = -2.0/sigmau2+sum(linear*linear)/2.0/sigmau2/sigmau2;
          for (int mi=0; mi<p+3; mi++){
            for (int mj=0; mj<p+3; mj++){
              num(mi,mj) = num(mi,mj) + ld[mi] * ld[mj] * wdens;
            }
          }
          den = den + wdens;
          }
        }
        }
      }
    }
    }
  }
  if(den>0){
  for (int mi=0; mi<p+3; mi++){
  for (int mj=0; mj<p+3; mj++){
    Ell(mi,mj) = Ell(mi,mj) + num(mi,mj)/den;
  }
  }
  }
  id_begin = id_begin + ni;
  }
  NumericMatrix Icomneg (p+3);
  double ldd_t = - 2.0 / (theta_k*theta_k) * n;
  if(theta_k==0) {ldd_t=0;}
  Icomneg(0,0) = ldd_t;
  Icomneg(1,1) = - N / (gamma_kk*gamma_kk*2.0);
  for (int mi=0; mi<p; mi++){
  for (int mj=0; mj<p; mj++){
    Icomneg(2+mi,2+mj) = - XX(mi,mj) / sigmau2;
  }
  }
  Icomneg(p+2,p+2) = - N/(sigmau2*sigmau2*2.0);
  NumericMatrix V (p+3,p+3);
  for (int mi=0; mi<p+3; mi++){
    for (int mj=0; mj<p+3; mj++){
    V(mi,mj) = - Icomneg(mi,mj) + El2(mi,mj) - Ell(mi,mj);
  }
  }
  V = ginv(V);
  NumericMatrix Vs = mprod(V,El2);
  Vs = mprod(Vs,V);
  NumericVector se (p+3);
  for (int k=0; k<p+3; k++){
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
  NumericMatrix phi (n,p+4);
  NumericVector meanphi (p+3);
  for (int k=0; k<p+3; k++){
  meanphi[k] = sum(phin(_,k)) / n;
  }
  for (int i=0; i<n; i++){
  for (int k=0; k<p+3; k++){
  phi(i,0) = Fam_list[i];
  phi(i,k+1) = phin(i,k) - meanphi[k];
  }
  }
  double S2 = theta_k*theta_k + gamma_kk + sigmau2;
  double rho = gamma_kk/S2;
  double zeta = theta_k*theta_k/S2;
  double D1 = -2.0*theta_k*gamma_kk/S2/S2;
  double D2 = (theta_k*theta_k+sigmau2)/S2/S2;
  double D3 = -gamma_kk/S2/S2;
  NumericVector IF = phi(_,1)*D1 + phi(_,2)*D2 + phi(_,p+3)*D3;
  double serho = sqrt(sum(IF*IF)/(n-1)/n);
  D1 = 2.0*theta_k*(gamma_kk+sigmau2)/S2/S2;
  D2 = -theta_k*theta_k/S2/S2;
  D3 = D2;
  IF = phi(_,1)*D1 + phi(_,2)*D2 + phi(_,p+3)*D3;
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
  results["sigmau2"] = sigmau2;
  results["se.theta"] = se[0];
  results["se.gamma"] = se[1];
  NumericVector se_alpha (p);
  for (int k=0; k<p; k++){
  se_alpha[k] = se[2+k];
  }
  results["se.alpha"] = se_alpha;
  results["se.sigma"] = se[p+2];
  results["est"] = res;
  results["se"] = se;
  results["phi"] = phi;
  results["V"] = V;
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
