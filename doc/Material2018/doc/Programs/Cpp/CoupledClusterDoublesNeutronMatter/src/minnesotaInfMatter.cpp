#include "minnesotaInfMatter.hpp"
#include <cmath>

// Minnesota Potential for momentum basis, Not antisym
double vint_Minnesota_Momentum(int *qi,int *qj,int *qk,int *ql, double L){
  double vout = 0.;
  int initialMomx, initialMomy, initialMomz;
  int finalMomx, finalMomy, finalMomz;
  // Momentum Conservation Checks.
  initialMomx = qi[0] + qj[0];
  initialMomy = qi[1] + qj[1];
  initialMomz = qi[2] + qj[2];
  finalMomx = qk[0] + ql[0];
  finalMomy = qk[1] + ql[1];
  finalMomz = qk[2] + ql[2];
  
  // maybe add a conservation of spin if here.
  if( initialMomx == finalMomx && initialMomy == finalMomy && initialMomz == finalMomz ){
    
    double V_R, V_T, V_S;
    double V_0R, V_0T, V_0S;
    double kappa_R, kappa_T, kappa_S;
    double *ki = new double[3];  
    double *kj = new double[3];  
    double *kk = new double[3];
    double *kl = new double[3];
    double *relMomBra = new double[3];
    double *relMomKet = new double[3];
    double *relMomTransf = new double[3];
    double qSquared, spinEx, isoSpinEx;
    double IsIt, PsIt, PsPt, IsPt;
    V_0R = 200; //MeV
    V_0T = 178; //MeV
    V_0S = 91.85; //MeV
    kappa_R = 1.487; //fm^-2
    kappa_T = 0.639; //fm^-2
    kappa_S = 0.465; //fm^-2
    
    qSquared = 0.;
    for( int i = 0; i < 3; i++){
      ki[i] = 2*M_PI*qi[i]/L;
      kj[i] = 2*M_PI*qj[i]/L;
      kk[i] = 2*M_PI*qk[i]/L;
      kl[i] = 2*M_PI*ql[i]/L;
      relMomBra[i] = 0.5*(ki[i] - kj[i]);
      relMomKet[i] = 0.5*(kk[i] - kl[i]);
      relMomTransf[i] = relMomBra[i] - relMomKet[i];
      qSquared += relMomTransf[i]*relMomTransf[i];
    }

    V_R = V_0R/(L*L*L)*pow(M_PI/kappa_R,1.5)*exp(-qSquared/(4*kappa_R));
    V_T = -V_0T/(L*L*L)*pow(M_PI/kappa_T,1.5)*exp(-qSquared/(4*kappa_T));
    V_S = -V_0S/(L*L*L)*pow(M_PI/kappa_S,1.5)*exp(-qSquared/(4*kappa_S));
    
    spinEx = spinExchangeMtxEle(qi[3],qj[3],qk[3],ql[3]);
    isoSpinEx = spinExchangeMtxEle(qi[4],qj[4],qk[4],ql[4]);
    
    // 4 terms, IsIt, PsIt, PsPt, IsPt
    // identity spin, identity isospin
    IsIt = kron_del(qi[3],qk[3])*kron_del(qj[3],ql[3])*kron_del(qi[4],qk[4])*kron_del(qj[4],ql[4]);
    // Exchange spin, identity isospin
    PsIt = spinEx*kron_del(qi[4],qk[4])*kron_del(qj[4],ql[4]);
    // Exchange spin, Exchange isospin
    PsPt = spinEx*isoSpinEx;
    // identity spin, Exchange isospin
    IsPt = kron_del(qi[3],qk[3])*kron_del(qj[3],ql[3])*isoSpinEx;
    
    vout = 0.5*(V_R + 0.5*V_T + 0.5*V_S)*IsIt
      + 0.25*(V_T - V_S)*PsIt
      - 0.5*(V_R + 0.5*V_T + 0.5*V_S)*PsPt
      - 0.25*(V_T - V_S)*IsPt;
    
    delete[] ki;  
    delete[] kj;
    delete[] kk;
    delete[] kl;
    delete[] relMomBra;
    delete[] relMomKet;
    delete[] relMomTransf;
  } // end Momentum Conversation Check

  return vout;
 
} // end vint_Minnesota_Momentum


int spinExchangeMtxEle(int i, int j, int k, int l){
  if( i == l && j == k ){
    return 1;
  } else {
    return 0;
  }
} // end spinEx



int kron_del(int i, int j){
  if(i != j){
    return 0;
  }
  return 1;
} // end kron_del


