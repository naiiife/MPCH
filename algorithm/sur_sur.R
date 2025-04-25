cppFunction('
  List sur_sur(NumericVector T1, NumericVector D1, NumericVector T2, NumericVector D2,
        NumericMatrix X1, NumericMatrix X2, IntegerVector familyID,
        List res1, List res2, NumericVector kinship, NumericVector par){
  Environment glob = Environment::global_env();
  NumericVector w_list = glob.get("w_list");
  NumericVector b_list = glob.get("b_list");
  int ngrid = w_list.size();
  int N = familyID.size();
  int p1 = X1(0,_).size();
  int p2 = X2(0,_).size();
  double theta_1 = res1["theta"];
  double gamma_1 = res1["gamma"];
  NumericVector alpha_1 = res1["alpha"];
  NumericVector Ts1 = res1["tt"];
  NumericVector lam1 = res1["lam"];
  NumericMatrix lam1_dh = res1["lam_dh"];
  int L1 = Ts1.size();
  NumericVector Lam1 (N);
  for (int i=0; i<N; i++){
  double ti = T1[i];
  NumericVector TsL = ifelse(Ts1<=ti,1.0,0.0);
  Lam1[i] = sum(TsL*lam1);
  }
  NumericMatrix phi1 = res1["phi"];
  double theta_2 = res2["theta"];
  double gamma_2 = res2["gamma"];
  NumericVector alpha_2 = res2["alpha"];
  NumericVector Ts2 = res2["tt"];
  NumericVector lam2 = res2["lam"];
  int L2 = Ts2.size();
  NumericVector Lam2 (N);
  for (int i=0; i<N; i++){
  double ti = T2[i];
  NumericVector TsL = ifelse(Ts2<=ti,1.0,0.0);
  Lam2[i] = sum(TsL*lam2);
  }
  NumericMatrix lam2_dh = res2["lam_dh"];
  NumericMatrix phi2 = res2["phi"];
  int n1 = phi1(_,0).size();
  int n2 = phi2(_,0).size();
  NumericVector nn (n1+n2);
  IntegerVector Fam_list;
  for (int i=0; i<n1; i++){
    nn[i] = phi1(i,0);
  }
  for (int i=0; i<n2; i++){
    nn[i+n1] = phi2(i,0);
  }
  int n = max(nn);
  int n0 = min(nn);
  for (int i=n0; i<n+1; i++){
  int ni = sum(ifelse(phi1(_,0)==i,1,0)) + sum(ifelse(phi2(_,0)==i,1,0));
  if(ni>0) {Fam_list.push_back(i);}
  }
  n = Fam_list.size();
  double maxlk = 0;
  double gamma12_max = par[2];
  double gamma12_min = par[1];
  double gamma12 = par[0];
  double tol = 1.0;
  NumericMatrix Grt (2,2);
  Grt(0,0) = sqrt(gamma_1);
  
  while(tol>0.001){
    double gamma120 = gamma12;
    double gamma12_left = gamma12*0.9+gamma12_min*0.1;
    double gamma12_right = gamma12*0.9+gamma12_max*0.1;
    NumericVector lk;
    double g = gamma12;
    for (int q=0; q<2; q++){
    if(q==0) {
      g = gamma12_left;
    } else{
      g = gamma12_right;
    }
    double loglik = 0.0;
    for (int i=0; i<N; i++){
      double D1i = D1[i];
      double Lam1i = Lam1[i];
      double D2i = D2[i];
      double Lam2i = Lam2[i];
      NumericVector X1i = X1(i,_);
      NumericVector X2i = X2(i,_);
      double ks = kinship[i];
      double den = 0;
      double gamma12k = g * ks;
      Grt(1,0) = gamma12k/sqrt(gamma_1);
      Grt(1,1) = sqrt(gamma_2 - gamma12k*gamma12k/gamma_1);
      for (int b1=0; b1<ngrid; b1++){
        double b_eps1 = 1.4142136 * b_list[b1];
        double w_eps1 = w_list[b1];
        for (int b2=0; b2<ngrid; b2++){
          double b_eps2 = 1.4142136 * b_list[b2];
          double w_eps2 = w_list[b2];
            double eps1 = Grt(0,0) * b_eps1;
            double eps2 = Grt(1,0) * b_eps1 + Grt(1,1) * b_eps2;
          for (int b=0; b<ngrid; b++){
            double e = 1.4142136 * b_list[b];
            double w_e = w_list[b];
            double linear1 = sum(X1i*alpha_1)+theta_1*e+eps1;
            double linear2 = sum(X2i*alpha_2)+theta_2*e+eps2;
            double loglik1 = D1i*linear1 - Lam1i*exp(linear1);
            double loglik2 = D2i*linear2 - Lam2i*exp(linear2);
            double El = exp(loglik1+loglik2);
            double wdens = El * w_eps1 * w_eps2 * w_e;
            den = den + wdens;
          }
        }
      }
      if (den>0){
      loglik = loglik + log(den);
      }
    }
    lk.push_back(loglik);
    }
    if(lk[0]<lk[1]) {
      gamma12_min = gamma12_left;
    } else {
      gamma12_max = gamma12_right;
    }
    gamma12 = (gamma12_min+gamma12_max)/2.0;
    maxlk = max(lk);
    tol = sqrt((gamma12-gamma120)*(gamma12-gamma120));
    std::cout << gamma12 << std::endl;
  }
  double gamma12_alt = gamma12;
  
  theta_2 = -theta_2;
  double maxlk_neg = 0;
  gamma12_max = par[2];
  gamma12_min = par[1];
  gamma12 = par[0];
  tol = 1.0;
  Grt(0,0) = sqrt(gamma_1);
  
  while(tol>0.001){
    double gamma120 = gamma12;
    double gamma12_left = gamma12*0.9+gamma12_min*0.1;
    double gamma12_right = gamma12*0.9+gamma12_max*0.1;
    NumericVector lk;
    double g = gamma12;
    for (int q=0; q<2; q++){
    if(q==0) {
      g = gamma12_left;
    } else{
      g = gamma12_right;
    }
    double loglik = 0.0;
    for (int i=0; i<N; i++){
      double D1i = D1[i];
      double Lam1i = Lam1[i];
      double D2i = D2[i];
      double Lam2i = Lam2[i];
      NumericVector X1i = X1(i,_);
      NumericVector X2i = X2(i,_);
      double ks = kinship[i];
      double den = 0;
      double gamma12k = g * ks;
      Grt(1,0) = gamma12k/sqrt(gamma_1);
      Grt(1,1) = sqrt(gamma_2 - gamma12k*gamma12k/gamma_1);
      for (int b1=0; b1<ngrid; b1++){
        double b_eps1 = 1.4142136 * b_list[b1];
        double w_eps1 = w_list[b1];
        for (int b2=0; b2<ngrid; b2++){
          double b_eps2 = 1.4142136 * b_list[b2];
          double w_eps2 = w_list[b2];
            double eps1 = Grt(0,0) * b_eps1;
            double eps2 = Grt(1,0) * b_eps1 + Grt(1,1) * b_eps2;
          for (int b=0; b<ngrid; b++){
            double e = 1.4142136 * b_list[b];
            double w_e = w_list[b];
            double linear1 = sum(X1i*alpha_1)+theta_1*e+eps1;
            double linear2 = sum(X2i*alpha_2)+theta_2*e+eps2;
            double loglik1 = D1i*linear1 - Lam1i*exp(linear1);
            double loglik2 = D2i*linear2 - Lam2i*exp(linear2);
            double El = exp(loglik1+loglik2);
            double wdens = El * w_eps1 * w_eps2 * w_e;
            den = den + wdens;
          }
        }
      }
      if (den>0){
      loglik = loglik + log(den);
      }
    }
    lk.push_back(loglik);
    }
    if(lk[0]<lk[1]) {
      gamma12_min = gamma12_left;
    } else {
      gamma12_max = gamma12_right;
    }
    gamma12 = (gamma12_min+gamma12_max)/2.0;
    maxlk_neg = max(lk);
    tol = sqrt((gamma12-gamma120)*(gamma12-gamma120));
    std::cout << gamma12 << std::endl;
  }
  
  if(maxlk > maxlk_neg){
    theta_2 = -theta_2;
    gamma12 = gamma12_alt;
  }
  
  // Calculate variance
  
  NumericVector dl_gamma (N);
  NumericVector dl2_gamma (N);
  NumericVector dl_gamma2 (N);
  for (int i=0; i<N; i++){
    double D1i = D1[i];
    double Lam1i = Lam1[i];
    double D2i = D2[i];
    double Lam2i = Lam2[i];
    NumericVector X1i = X1(i,_);
    NumericVector X2i = X2(i,_);
    double ks = kinship[i];
    double gamma12k = gamma12 * ks;
    double dG = gamma_1 * gamma_2 - gamma12k*gamma12k;
    Grt(0,0) = sqrt(gamma_1);
    Grt(1,0) = gamma12k/sqrt(gamma_1);
    Grt(1,1) = sqrt(gamma_2 - gamma12k*gamma12k/gamma_1);
    double den = 0;
    double num1 = 0;
    double num2 = 0;
    double num3 = 0;
    for (int b1=0; b1<ngrid; b1++){
      double b_eps1 = 1.4142136 * b_list[b1];
      double w_eps1 = w_list[b1];
      for (int b2=0; b2<ngrid; b2++){
        double b_eps2 = 1.4142136 * b_list[b2];
        double w_eps2 = w_list[b2];
          double eps1 = Grt(0,0) * b_eps1;
          double eps2 = Grt(1,0) * b_eps1 + Grt(1,1) * b_eps2;
        for (int b=0; b<ngrid; b++){
          double e = 1.4142136 * b_list[b];
          double w_e = w_list[b];
          double linear1 = sum(X1i*alpha_1)+theta_1*e+eps1;
          double linear2 = sum(X2i*alpha_2)+theta_2*e+eps2;
          double loglik1 = D1i*linear1 - Lam1i*exp(linear1);
          double loglik2 = D2i*linear2 - Lam2i*exp(linear2);
          double be2 = b_eps1*b_eps1 + b_eps2*b_eps2;
          double El = exp(loglik1+loglik2);
          double wdens = El * w_eps1 * w_eps2 * w_e;
          if(wdens>0){
          double dldg = gamma12k/dG*(1.0-be2) + eps1*eps2/dG;
          double dl2dg = (1.0-be2)/dG + 2.0*gamma12k*(gamma12k*(1.0-2.0*be2)+2.0*eps1*eps2)/dG/dG;
          double dldg2 = dldg*dldg;
          den = den + wdens;
          num1 = num1 + wdens * dldg;
          num2 = num2 + wdens * dl2dg;
          num3 = num3 + wdens * dldg2;
          }
        }
      }
    }
    if(den>0){
    dl_gamma[i] =  ks*num1/den;
    dl2_gamma[i] = ks*ks*num2/den;
    dl_gamma2[i] = ks*ks*num3/den;
    }
  }
  
  NumericVector res01 (p1+2);
  NumericVector res02 (p2+2);
  res01[0] = theta_1;
  res01[1] = gamma_1;
  for (int k=0; k<p1; k++){
  res01[2+k] = alpha_1[k];
  }
  res02[0] = theta_2;
  res02[1] = gamma_2;
  for (int k=0; k<p2; k++){
  res02[2+k] = alpha_2[k];
  }
  NumericVector h1 = res1["h"];
  NumericVector h2 = res2["h"];
  // NumericVector lam1T (N);
  // NumericVector lam2T (N);
  // for (int i=0; i<N; i++){
  //  double t1i = T1[i];
  //  NumericVector Tsi1 = ifelse(Ts1==t1i, 1.0, 0.0);
  //  double dl1 = sum(Tsi1*lam1);
  //  if(dl1>0) {lam1T[i] = log(dl1);}
  //  double t2i = T2[i];
  //  NumericVector Tsi2 = ifelse(Ts2==t2i, 1.0, 0.0);
  //  double dl2 = sum(Tsi2*lam2);
  //  if(dl2>0) {lam2T[i] = log(dl2);}
  // }
  
  NumericVector dU_eta1 (p1+2);
  NumericVector dU_eta2 (p2+2);
  for (int q=0; q<2; q++){
  for (int k=0; k<p1+2; k++){
  NumericVector lam1h;
  NumericVector Lam1h;
  // NumericVector lam1Th;
  NumericVector resh;
  if(q==0) {
  lam1h.assign(lam1.begin(), lam1.end());
  Lam1h.assign(Lam1.begin(), Lam1.end());
  // lam1Th.assign(lam1T.begin(), lam1T.end());
  } else {
  lam1h.assign(lam1.begin(), lam1.end());
  Lam1h.assign(Lam1.begin(), Lam1.end());
  // lam1Th.assign(lam1T.begin(), lam1T.end());
  for (int l=0; l<L1; l++){
     lam1h[l] = lam1_dh(l,k);
  }
  for (int i=0; i<N; i++){
    double t = T1[i];
    NumericVector TsL = ifelse(Ts1<=t,1.0,0.0);
  //  NumericVector TsD = ifelse(Ts1==t,1.0,0.0);
    Lam1h[i] = sum(TsL*lam1h);
  //  double dl = sum(TsD*lam1h);
  //  if(dl>0) {
  //  lam1Th[i] = log(dl);
  //  } else {
  //  lam1Th[i] = 0;
  //  }
  }
  }
  resh.assign(res01.begin(), res01.end());
  resh[k] = resh[k] + h1[k]*q;
  theta_1 = resh[0];
  gamma_1 = resh[1];
  for (int r=0; r<p1; r++){
  alpha_1[r] = resh[r+2];
  }
  
  NumericVector lk (N);
  for (int i=0; i<N; i++){
  double D1i = D1[i];
  double Lam1i = Lam1h[i];
  // double lam1Ti = lam1Th[i];
  NumericVector X1i = X1(i,_);
  NumericVector X2i = X2(i,_);
  double D2i = D2[i];
  double Lam2i = Lam2[i];
  double ks = kinship[i];
  double den = 0;
  double num7 = 0;
  double gamma12k = gamma12 * ks;
  double dG = gamma_1 * gamma_2 - gamma12k*gamma12k;
  Grt(0,0) = sqrt(gamma_1);
  Grt(1,0) = gamma12k/sqrt(gamma_1);
  Grt(1,1) = sqrt(gamma_2 - gamma12k*gamma12k/gamma_1);
      for (int b1=0; b1<ngrid; b1++){
        double b_eps1 = 1.4142136 * b_list[b1];
        double w_eps1 = w_list[b1];
        for (int b2=0; b2<ngrid; b2++){
          double b_eps2 = 1.4142136 * b_list[b2];
          double w_eps2 = w_list[b2];
            double eps1 = Grt(0,0) * b_eps1;
            double eps2 = Grt(1,0) * b_eps1 + Grt(1,1) * b_eps2;
          for (int b=0; b<ngrid; b++){
            double e = 1.4142136 * b_list[b];
            double w_e = w_list[b];
            double linear1 = sum(X1i*alpha_1)+theta_1*e+eps1;
            double linear2 = sum(X2i*alpha_2)+theta_2*e+eps2;
            double loglik1 = D1i*linear1 - Lam1i*exp(linear1);
            double loglik2 = D2i*linear2 - Lam2i*exp(linear2);
            double be2 = b_eps1*b_eps1 + b_eps2*b_eps2;
            double El = exp(loglik1+loglik2);
            double wdens = El * w_eps1 * w_eps2 * w_e;
            den = den + wdens;
            double dldg = gamma12k/dG*(1.0-be2) + eps1*eps2/dG;
            num7 = num7 + wdens * dldg * ks;
          }
        }
      }
    lk[i] = num7/den;
    }
  
  dU_eta1[k] = dU_eta1[k] + (2.0*q-1.0)*sum(lk);
  resh[k] = res01[k];
  }
  }
  NumericVector B1 = dU_eta1/(h1*n);
  
  alpha_1[p1-1] = res01[p1+1];
  for (int q=0; q<2; q++){
  for (int k=0; k<p2+2; k++){
  NumericVector lam2h;
  NumericVector Lam2h;
  // NumericVector lam2Th;
  NumericVector resh;
  if(q==0) {
  lam2h.assign(lam2.begin(), lam2.end());
  Lam2h.assign(Lam2.begin(), Lam2.end());
  // lam2Th.assign(lam2T.begin(), lam2T.end());
  } else {
  lam2h.assign(lam2.begin(), lam2.end());
  Lam2h.assign(Lam2.begin(), Lam2.end());
  // lam2Th.assign(lam2T.begin(), lam2T.end());
  for (int l=0; l<L2; l++){
    lam2h[l] = lam2_dh(l,k);
  }
  for (int i=0; i<N; i++){
    double t = T2[i];
    NumericVector TsL = ifelse(Ts2<=t,1.0,0.0);
    // NumericVector TsD = ifelse(Ts2==t,1.0,0.0);
    Lam2h[i] = sum(TsL*lam2h);
    // double dl = sum(TsD*lam2h);
    // if(dl>0) {
    // lam2Th[i] = log(dl);
    // } else {
    // lam2Th[i] = 0;
    // }
  }
  }
  resh.assign(res02.begin(), res02.end());
  resh[k] = resh[k] + h2[k]*q;
  theta_2 = resh[0];
  gamma_2 = resh[1];
  for (int r=0; r<p2; r++){
  alpha_2[r] = resh[r+2];
  }
  
  NumericVector lk (N);
  for (int i=0; i<N; i++){
  double D1i = D1[i];
  double Lam1i = Lam1[i];
  NumericVector X1i = X1(i,_);
  NumericVector X2i = X2(i,_);
  double D2i = D2[i];
  double Lam2i = Lam2h[i];
  // double lam2Ti = lam2Th[i];
  double ks = kinship[i];
  double den = 0;
  double num7 = 0;
  double gamma12k = gamma12 * ks;
  double dG = gamma_1 * gamma_2 - gamma12k*gamma12k;
  Grt(1,0) = gamma12k/sqrt(gamma_1);
  Grt(1,1) = sqrt(gamma_2 - gamma12k*gamma12k/gamma_1);
      for (int b1=0; b1<ngrid; b1++){
        double b_eps1 = 1.4142136 * b_list[b1];
        double w_eps1 = w_list[b1];
        for (int b2=0; b2<ngrid; b2++){
          double b_eps2 = 1.4142136 * b_list[b2];
          double w_eps2 = w_list[b2];
            double eps1 = Grt(0,0) * b_eps1;
            double eps2 = Grt(1,0) * b_eps1 + Grt(1,1) * b_eps2;
          for (int b=0; b<ngrid; b++){
            double e = 1.4142136 * b_list[b];
            double w_e = w_list[b];
            double linear1 = sum(X1i*alpha_1)+theta_1*e+eps1;
            double linear2 = sum(X2i*alpha_2)+theta_2*e+eps2;
            double loglik1 = D1i*linear1 - Lam1i*exp(linear1);
            double loglik2 = D2i*linear2 - Lam2i*exp(linear2);
            double be2 = b_eps1*b_eps1 + b_eps2*b_eps2;
            double El = exp(loglik1+loglik2);
            double wdens = El * w_eps1 * w_eps2 * w_e;
            den = den + wdens;
            double dldg = gamma12k/dG*(1.0-be2) + eps1*eps2/dG;
            num7 = num7 + wdens * dldg * ks;
          }
        }
      }
    lk[i] = num7/den;
    }
  
  dU_eta2[k] = dU_eta2[k] + (2.0*q-1.0)*sum(lk);
  resh[k] = res02[k];
  }
  }
  NumericVector B2 = dU_eta2/(h2*n);
  
  NumericVector dl_gammai (n);
  for (int i=0; i<n; i++){
  int fi = Fam_list[i];
  NumericVector fami = ifelse(familyID==fi,1.0,0.0);
  dl_gammai[i] = sum(dl_gamma*fami);
  }
  double ldd_gamma = sum(dl_gamma2+dl2_gamma-dl_gamma*dl_gamma)/n;
  NumericVector M1 (n);
  NumericVector M2 (n);
  B1.push_front(0);
  B2.push_front(0);
  for (int i=0; i<n; i++){
    int fi = Fam_list[i];
    for (int j=0; j<n1; j++){
    if(phi1(j,0)==fi) {M1[i] = sum(B1*phi1(j,_));}
    }
    for (int j=0; j<n2; j++){
    if(phi2(j,0)==fi) {M2[i] = sum(B2*phi2(j,_));}
    }
  }
  NumericVector sc = (dl_gammai - M1 - M2)/ldd_gamma;
  double V = mean(sc*sc);
  double se = sqrt(V/n);

  double sigmau21 = 1.6449341;
  double sigmau22 = 1.6449341;
  double sigma_1 = sqrt(theta_1*theta_1+gamma_1+1.6449341);
  double sigma_2 = sqrt(theta_2*theta_2+gamma_2+1.6449341);
  double rho_1 = gamma_1/(sigma_1*sigma_1);
  double rho_2 = gamma_2/(sigma_2*sigma_2);
  double rho12 = gamma12/(sigma_1*sigma_2);
  double zeta1 = theta_1*theta_1/(sigma_1*sigma_1);
  double zeta2 = theta_2*theta_2/(sigma_2*sigma_2);
  double zeta12 = theta_1*theta_2/(sigma_1*sigma_2);
  NumericVector Drho_1 (p1+3);
  NumericVector Drho_2 (p2+3);
  NumericVector Drho12_1 (p1+3);
  NumericVector Drho12_2 (p2+3);
  NumericVector Dzeta_1 (p1+3);
  NumericVector Dzeta_2 (p2+3);
  NumericVector Dzeta12_1 (p1+3);
  NumericVector Dzeta12_2 (p2+3);
  Drho_1[1] = -2.0*theta_1*gamma_1/(sigma_1*sigma_1*sigma_1*sigma_1);
  Drho_1[2] = (theta_1*theta_1+sigmau21)/(sigma_1*sigma_1*sigma_1*sigma_1);
  Drho_2[1] = -2.0*theta_2*gamma_2/(sigma_2*sigma_2*sigma_2*sigma_2);
  Drho_2[2] = (theta_2*theta_2+sigmau22)/(sigma_2*sigma_2*sigma_2*sigma_2);
  Drho12_1[1] = -rho12*theta_1/(sigma_1*sigma_1);
  Drho12_1[2] = -rho12/(2.0*sigma_1*sigma_1);
  Drho12_2[1] = -rho12*theta_2/(sigma_2*sigma_2);
  Drho12_2[2] = -rho12/(2.0*sigma_2*sigma_2);
  Dzeta_1[1] = 2.0*theta_1*(gamma_1+sigmau21)/(sigma_1*sigma_1*sigma_1*sigma_1);
  Dzeta_1[2] = -theta_1*theta_1/(sigma_1*sigma_1*sigma_1*sigma_1);
  Dzeta_2[1] = 2.0*theta_2*(gamma_2+sigmau22)/(sigma_2*sigma_2*sigma_2*sigma_2);
  Dzeta_2[2] = -theta_2*theta_2/(sigma_2*sigma_2*sigma_2*sigma_2);
  Dzeta12_1[1] = (theta_2*sigma_1-theta_1*theta_1*theta_2/sigma_1)/(sigma_1*sigma_1*sigma_2);
  Dzeta12_1[2] = -zeta12/(2.0*sigma_1*sigma_1);
  Dzeta12_2[1] = (theta_1*sigma_2-theta_2*theta_2*theta_1/sigma_2)/(sigma_1*sigma_2*sigma_2);
  Dzeta12_2[2] = -zeta12/(2.0*sigma_2*sigma_2);
  double Drho12_12 = rho12/gamma12;
  NumericVector IF_rho_1 (n);
  NumericVector IF_rho_2 (n);
  NumericVector IF_rho12 (n);
  NumericVector IF_zeta_1 (n);
  NumericVector IF_zeta_2 (n);
  NumericVector IF_zeta12 (n);
  for (int i=0; i<n; i++){
    int fi = Fam_list[i];
    for (int j=0; j<n1; j++){
      if(phi1(j,0)==fi) {
        IF_rho_1[i] = sum(Drho_1*phi1(j,_));
        IF_rho12[i] = sum(Drho12_1*phi1(j,_));
        IF_zeta_1[i] = sum(Dzeta_1*phi1(j,_));
        IF_zeta12[i] = sum(Dzeta12_1*phi1(j,_));
      }
    }
    for (int j=0; j<n2; j++){
      if(phi2(j,0)==fi) {
        IF_rho_2[i] = sum(Drho_2*phi2(j,_));
        IF_rho12[i] = IF_rho12[i] + sum(Drho12_2*phi2(j,_));
        IF_zeta_2[i] = sum(Dzeta_2*phi2(j,_));
        IF_zeta12[i] = IF_zeta12[i] + sum(Dzeta12_2*phi2(j,_));
      }
    }
    IF_rho12[i] = IF_rho12[i] + Drho12_12*sc[i];
  }
  double se_rho1 = sqrt(sum(IF_rho_1*IF_rho_1)/n1/n1);
  double se_rho2 = sqrt(sum(IF_rho_2*IF_rho_2)/n2/n2);
  double se_rho12 = sqrt(sum(IF_rho12*IF_rho12)/n/n);
  double se_zeta1 = sqrt(sum(IF_zeta_1*IF_zeta_1)/n1/n1);
  double se_zeta2 = sqrt(sum(IF_zeta_2*IF_zeta_2)/n2/n2);
  double se_zeta12 = sqrt(sum(IF_zeta12*IF_zeta12)/n/n);
  
  NumericVector diss1 (p1+3);
  NumericVector diss2 (p2+3);
  diss1[1] = 1;
  diss1[2] = 1;
  diss2[1] = 1;
  diss2[2] = 1;
  Drho_1 = Drho_1 * 0.0;
  Drho_2 = Drho_2 * 0.0;
  Drho12_1 = Drho12_1 * 0.0;
  Drho12_2 = Drho12_2 * 0.0;
  Dzeta_1 = Dzeta_1 * 0.0;
  Dzeta_2 = Dzeta_2 * 0.0;
  Dzeta12_1 = Dzeta12_1 * 0.0;
  Dzeta12_2 = Dzeta12_2 * 0.0;
  IF_rho_1 = IF_rho_1 * 0.0;
  IF_rho_2 = IF_rho_2 * 0.0;
  IF_rho12 = IF_rho12 * 0.0;
  IF_zeta_1 = IF_zeta_1 * 0.0;
  IF_zeta_2 = IF_zeta_2 * 0.0;
  IF_zeta12 = IF_zeta12 * 0.0;
  double sigma_1d = sqrt(theta_1*theta_1+gamma_1);
  double sigma_2d = sqrt(theta_2*theta_2+gamma_2);
  double rho_1_p = gamma_1/(sigma_1d*sigma_1d);
  double rho_2_p = gamma_2/(sigma_2d*sigma_2d);
  double rho12_p = gamma12/(sigma_1d*sigma_2d);
  double zeta1_p = theta_1*theta_1/(sigma_1d*sigma_1d);
  double zeta2_p = theta_2*theta_2/(sigma_2d*sigma_2d);
  double zeta12_p = theta_1*theta_2/(sigma_1d*sigma_2d);
  Drho_1[1] = -2.0*theta_1*gamma_1/(sigma_1d*sigma_1d*sigma_1d*sigma_1d);
  Drho_1[2] = (theta_1*theta_1)/(sigma_1d*sigma_1d*sigma_1d*sigma_1d);
  Drho_2[1] = -2.0*theta_2*gamma_2/(sigma_2d*sigma_2d*sigma_2d*sigma_2d);
  Drho_2[2] = (theta_2*theta_2)/(sigma_2d*sigma_2d*sigma_2d*sigma_2d);
  Drho12_1[1] = -rho12_p*theta_1/(sigma_1d*sigma_1d);
  Drho12_1[2] = -rho12_p/(2.0*sigma_1d*sigma_1d);
  Drho12_2[1] = -rho12_p*theta_2/(sigma_2d*sigma_2d);
  Drho12_2[2] = -rho12_p/(2.0*sigma_2d*sigma_2d);
  Dzeta_1[1] = 2.0*theta_1*gamma_1/(sigma_1d*sigma_1d*sigma_1d*sigma_1d);
  Dzeta_1[2] = -zeta1_p/(sigma_1d*sigma_1d);
  Dzeta_2[1] = 2.0*theta_2*gamma_2/(sigma_2d*sigma_2d*sigma_2d*sigma_2d);
  Dzeta_2[2] = -zeta2/(sigma_2d*sigma_2d);
  Dzeta12_1[1] = (theta_2*sigma_1d-theta_1*theta_1*theta_2/sigma_1d)/(sigma_1d*sigma_1d*sigma_2d);
  Dzeta12_1[2] = -zeta12_p/(2.0*sigma_1d*sigma_1d);
  Dzeta12_2[1] = (theta_1*sigma_2d-theta_2*theta_2*theta_1/sigma_2d)/(sigma_1d*sigma_2d*sigma_2d);
  Dzeta12_2[2] = -zeta12_p/(2.0*sigma_2d*sigma_2d);
  Drho12_12 = rho12_p/gamma12;
  for (int i=0; i<n; i++){
    int fi = Fam_list[i];
    for (int j=0; j<n1; j++){
      if(phi1(j,0)==fi) {
        IF_rho_1[i] = sum(Drho_1*phi1(j,_)*diss1);
        IF_rho12[i] = sum(Drho12_1*phi1(j,_)*diss1);
        IF_zeta_1[i] = sum(Dzeta_1*phi1(j,_)*diss1);
        IF_zeta12[i] = sum(Dzeta12_1*phi1(j,_)*diss1);
      }
    }
    for (int j=0; j<n2; j++){
      if(phi2(j,0)==fi) {
        IF_rho_2[i] = sum(Drho_2*phi2(j,_)*diss2);
        IF_rho12[i] = IF_rho12[i] + sum(Drho12_2*phi2(j,_)*diss2);
        IF_zeta_2[i] = sum(Dzeta_2*phi2(j,_)*diss2);
        IF_zeta12[i] = IF_zeta12[i] + sum(Dzeta12_2*phi2(j,_)*diss2);
      }
    }
    IF_rho12[i] = IF_rho12[i] + Drho12_12*sc[i];
  }
  double se_rho1_p = sqrt(sum(IF_rho_1*IF_rho_1)/n1/n1);
  double se_rho2_p = sqrt(sum(IF_rho_2*IF_rho_2)/n2/n2);
  double se_rho12_p = sqrt(sum(IF_rho12*IF_rho12)/n/n);
  double se_zeta1_p = sqrt(sum(IF_zeta_1*IF_zeta_1)/n1/n1);
  double se_zeta2_p = sqrt(sum(IF_zeta_2*IF_zeta_2)/n2/n2);
  double se_zeta12_p = sqrt(sum(IF_zeta12*IF_zeta12)/n/n);
  
  List res;
  res["gamma12"] = gamma12;
  res["se"] = se;
  res["rho1"] = rho_1;
  res["se.rho1"] = se_rho1;
  res["rho2"] = rho_2;
  res["se.rho2"] = se_rho2;
  res["rho12"] = rho12;
  res["se.rho12"] = se_rho12;
  res["zeta1"] = zeta1;
  res["se.zeta1"] = se_zeta1;
  res["zeta2"] = zeta2;
  res["se.zeta2"] = se_zeta2;
  res["zeta12"] = zeta12;
  res["se.zeta12"] = se_zeta12;
  res["rho1_p"] = rho_1_p;
  res["se.rho1_p"] = se_rho1_p;
  res["rho2_p"] = rho_2_p;
  res["se.rho2_p"] = se_rho2_p;
  res["rho12_p"] = rho12_p;
  res["se.rho12_p"] = se_rho12_p;
  res["zeta1_p"] = zeta1_p;
  res["se.zeta1_p"] = se_zeta1_p;
  res["zeta2_p"] = zeta2_p;
  res["se.zeta2_p"] = se_zeta2_p;
  res["zeta12_p"] = zeta12_p;
  res["se.zeta12_p"] = se_zeta12_p;
  return res;
  }
')