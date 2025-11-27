#include "electronGasInteraction.hpp"
#include <cmath>

// electron gas interaction for momentum basis, Not antisym
double vint_electronGas_Momentum(int *qi,int *qj,int *qk,int *ql, double L){

  double vout = 0.;

  // spin deltas  
  if( qi[3] == qk[3] && qj[3] == ql[3] ){
  
  int initialMomx, initialMomy, initialMomz;
  int finalMomx, finalMomy, finalMomz;
 
  initialMomx = qi[0] + qj[0];
  initialMomy = qi[1] + qj[1];
  initialMomz = qi[2] + qj[2];
  finalMomx = qk[0] + ql[0];
  finalMomy = qk[1] + ql[1];
  finalMomz = qk[2] + ql[2];
  
   // Momentum Conservation Checks.
  if( initialMomx == finalMomx && initialMomy == finalMomy && initialMomz == finalMomz ){
    
    double mu = 0.0;
    
    double hbarc = 0.1973269788; // eV micron
    hbarc = hbarc*10000/27.21138505; // Hartree Angstrom
    double fine_struc = 1.0/137.035999139;
    double e_sq = fine_struc*hbarc;
    double prefactor = e_sq/(L*L*L); // hard code this later?

    double relMom = 0.0; 
    double term;
    for( int i = 0; i < 3; i++){
      //ki[i] = 2*M_PI*qi[i]/L;
      //kj[i] = 2*M_PI*qj[i]/L;
      // kk[i] = 2*M_PI*qk[i]/L;
      // kl[i] = 2*M_PI*ql[i]/L;
      term = 2*M_PI*(qi[i] - qk[i])/L;
      relMom += term*term;
    }
    
    // get rid of 0 momentum transfer modes
    if(relMom < 1.e-5){
      vout = 0.0;
    } else {
      vout = prefactor*4*M_PI/(relMom + mu*mu);
    }

  } // end Momentum Conversation Check
  } // end spin check
  return vout;
 
} // end vint_electronGas_Momentum

