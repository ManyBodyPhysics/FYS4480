#include "pairingInteraction.hpp"
#include <cmath>


int kronecker_del(int i, int j){
  if(i != j){
    return 0;
  }
  return 1;
} // end kronecker_del


// pairing interaction, already antisym?
double vint_pairing(int *qi,int *qj,int *qk,int *ql, double g){
  double vout = 0.;

  if(kronecker_del(qi[0],qj[0]) * kronecker_del(qk[0],ql[0]) * kronecker_del(qi[1], -qj[1]) * kronecker_del(-qk[1], ql[1]) == 1){
      vout = -0.5*g;
      
      if( qi[1] == -qk[1] ){
	vout = -vout;
      }
  }
  return vout;
  
} // end v_int
